!################################################################################################
!include "Ferdosi.f90"
!include "AlphaVarCallParallelised.f90"
!include "ReadRogerData.f90"

module GlobalPar
	
	use ISO_Fortran_Env
	implicit none

	integer :: nIndSeq                         								! SpecFile - Number of Individuals in the sequence file
	integer :: LenghtSequenceDataFile,nSnp,fistWindow						! SpecFile - Total number of Snps
	
	integer :: InternalEdit                    								! SpecFile - Internal Edit 1==yes or 0==no
	real(kind=8) :: EditingParameter										! SpecFile - 1st Number is the MAF (excluede SNP with MAF=<EditingParameter)
	real(kind=8) :: maxStdForReadsCount,ThresholdMaxReadsCount				! SpecFile - Remove Reads that are above this standard deviation
	real(kind=8) :: ThresholdExcessHetero									! SpecFile - Remove variants with an excess of heterozygotes
	integer 	 :: ThresholdReadsCount										! SpecFile - Remove single/double/n-tones 
	
	integer :: nInd 														! Calculated Internally - Number of Individuals in the Pedigree
	
	real(kind=8) :: GeneProbThresh  										! SpecFile - Threshold to call a genotype from the probabilities First Value
	real(kind=8) :: GeneProbThreshMin										! SpecFile - Threshold to call a genotype from the probabilities Last Value
	real(kind=8) :: ReduceThr 												! SpecFile - Reduce Geno Treshold factor
	integer :: UsePrevGeneProb                      						! SpecFile - Read old results of GeneProb 1==YES, 0==NO
	
	real(kind=8) :: ErrorRate												! SpecFile - Error rates to define the genotypes probabilities
	
	integer :: ChunkLengthA                 								! SpecFile - First value to define Haplotypes length
	integer :: ChunkLengthB                 								! SpecFile - Last value to define Haplotypes length

	character(len=300) :: PedigreeFile      								! SpecFile - Input File Name - Pedigree
	character(len=300) :: ReadsFile             							! SpecFile - Input File Name - Reads Count for the Reference and the Alternative Allele
	character(len=300) :: ReadsType             							! SpecFile - Input File Name - Reads Count for the Reference and the Alternative Allele
	
	character(len=300) :: MapFile             								! SpecFile - Input File Name - Map File - position of the Variants
	character(len=300) :: SnpChipsInformation   							! SpecFile - Input File Name - Snp array to add more information to the Reads
	
	character(len=300) :: GenoFile              							! SpecFile - Control Results File Name - TrueGeno Genotypes to check results 
	character(len=300) :: PhaseFile             							! SpecFile - Control Results File Name - True Phase to check results 

	integer :: IterationNumber                  							! Control Parameter - Define the number of Iterations
	integer(int64) :: CurrentCountFilledPhase          							! Control Parameter - used to finish the program
	integer(int64) :: CurrentCountFilledGenos          							! Control Parameter - used to finish the program
	integer :: SolutionChanged                  							! Control Parameter - used to finish the program 
	integer :: StartSnp,EndSnp
	
	type(PedigreeHolder) :: ped

	integer(int64), dimension(:), allocatable 		:: position
	real(real64), allocatable, dimension(:) 		:: quality

	integer(kind=1),allocatable,dimension(:,:) 		:: TrueGenos			! Control Results - True Genotypes to check results 
	integer(kind=1),allocatable,dimension(:,:,:) 	:: TruePhase			! Control Results - True Phase to check results 

	integer(kind=1),allocatable,dimension(:)		:: MarkersToExclude		! CleanUpTheRawData - 0=use the variant; 1= don't use the variant
	
	character(len=1),allocatable,dimension(:,:) 	:: CheckGenos   		! Control Results - Use character to check True vs Imputed Genotypes
	character(len=1),allocatable,dimension(:,:,:) 	:: CheckPhase 			! Control Results - Use character to check True vs Imputed PhaseFile
	
	integer,allocatable,dimension(:) 				:: GeneProbYesOrNo		! Temporary Array - use gene prob or not
	integer,allocatable,dimension(:,:,:) 			:: FounderAssignment   	! Temporary File - Save the IDs of the grandparents

	! AlphaMPL output
	real(kind=real64),allocatable,dimension(:,:,:) 	:: ReadCounts !< in the format (pedigreeId, snpId, prob) prob is Pr00,Pr01,Pr10,Pr11
    real(kind=real64),allocatable,dimension(:) 		:: Maf
    
	type CountPhase
		integer,allocatable :: old(:)
		integer,allocatable :: diff(:)    
	end type CountPhase

	type(CountPhase) :: CurrentCountID

	integer :: Windows
end module GlobalPar

!################################################################################################

program FamilyPhase
    use ISO_Fortran_Env
    use omp_lib

	use GlobalPar
	use IntelRNGMod
	use AlphaMLPModule
	use CalculateStatisticsForGenoAndPhase
	implicit none

	integer(kind=8) :: OldCount,NewCount
	integer(int32) :: Seed1,Seed2
	real(kind=8) :: InitialGeneProbThresh
	logical:: fileExists
	real(kind=8)::tstart,tend
	
	
	! Use a seed to sample the Haplotypes length of each window and iteration
	! Print out the window/iteratin/haplotype length in a file
	inquire(file="Seed.txt", exist=fileExists)
	if (fileExists) then
		open(99,file="Seed.txt",action="read")
		read(99,*) Seed1
		close(99)
		call IntitialiseIntelRNG(Seed1,"SeedOld.txt",Seed2)

	else
		call IntitialiseIntelRNG(Seedfile="SeedOld.txt",Out=Seed2)
	end if
	
	open(101,file="AlphaFamSeqHaplotypeLengths.txt",status="unknown")
	write(101,'(1a32)') "Window Iter HaplotypesLengthUsed" 

	open(102,file="AlphaFamSeqWindowsInfo.txt",status="unknown")
	write(102,'(1a22)') "Window StartSnp EndSnp" 
	

	call ReadSpecfile ! TODO: the subroutine is not here anymore
	ped = PedigreeHolder(pedigreeFile) ! read pedigree
	call ped%sortPedigreeAndOverwrite() ! sort pedigree


	! TODO check type of genotype info - right now this for sequence 
	! TODO read vcf format
	if (trim(ReadsType)=="AlphaSim") call ped%addSequenceFromFile(readsFile,nSnp) ! read sequence file AlphaSimFormat
	
	! If the nSnp is really big and there are problems with memory allocation
	!  is possible to split the chromosome in multiple windows
	!  this will reduce the amount of memory used, but will print out several
	!  files for each windows. 

	! Multiple windows allows to restart the analysis from the last window if 
	!  the program crash for some reasons.

	Windows=fistWindow
	EndSnp=(Windows*nSnp)
	StartSnp=EndSnp-nSnp+1
	if (EndSnp>LenghtSequenceDataFile) EndSnp=LenghtSequenceDataFile
	InitialGeneProbThresh=GeneProbThresh

	do while(StartSnp.le.LenghtSequenceDataFile)
		!if (StartSnp>ChunkLengthB) StartSnp=StartSnp-ChunkLengthB
		!EndSnp=EndSnp+ChunkLengthB
		if (EndSnp>LenghtSequenceDataFile) EndSnp=LenghtSequenceDataFile
		nSnp=EndSnp-StartSnp+1

		write(102,'(3(1x,i0))') Windows,StartSnp,EndSnp

		CurrentCountFilledPhase=0
		CurrentCountFilledGenos=0
		GeneProbThresh=InitialGeneProbThresh

		!call CleanUpTheRawData

		if (UsePrevGeneProb==0) then
			call runAlphaMLPAlphaImpute(1, nSnp, ped, ReadCounts, Maf)
			!call SaveGeneProbResults
		else if (UsePrevGeneProb==1) then
			tstart = omp_get_wtime()
			call ReadPrevGeneProb
			tend = omp_get_wtime()
  			write(*,*) "Total wall time for Importing Probabilities", tend - tstart
		endif
		print*," "
		write (*,'(1a39)') " Window Iter   ProbThr    %Phase   %Geno"

		IterationNumber=0
		SolutionChanged=1
		do while ((SolutionChanged==1))!.and.(IterationNumber<2))

			OldCount=CurrentCountFilledPhase+CurrentCountFilledGenos
			IterationNumber=IterationNumber+1
			
			if ((IterationNumber>1).and.(GeneProbThresh>GeneProbThreshMin)) then 
				GeneProbThresh=GeneProbThresh-ReduceThr!0.001
				if (GeneProbThresh.lt.GeneProbThreshMin) GeneProbThresh=GeneProbThreshMin
			endif

			call UseGeneProbToSimpleFillInBasedOnOwnReads
			!if (IterationNumber==1) call ReadSamFile
			call SimpleCleanUpFillIn
			!call CurrentCountFilled

			!if (IterationNumber==1) call UseSnpChipInformation 

			call SimpleFillInBasedOnParentsReads
			call SimpleCleanUpFillIn
			!call CurrentCountFilled
			
			call SimpleFillInBasedOnProgenyReads
			call SimpleCleanUpFillIn
			!call CurrentCountFilled

			call CalculateFounderAssignment

			call ChunkDefinition
			
			call BuildConsensus
			call SimpleCleanUpFillIn
			
			call CurrentCountFilled
			
			!if (SuperC==1) then
				!call FerdosiSerap
				!call SimpleCleanUpFillIn
				!call CurrentCountFilled
			!endif

			NewCount=CurrentCountFilledPhase+CurrentCountFilledGenos
			
			if (OldCount/=NewCount) then
				SolutionChanged=1
			else
				SolutionChanged=0
			endif

			if ((maxval(CurrentCountID%diff).gt.(dble(nSnp)*.0001)*2).or.(GeneProbThresh.gt.GeneProbThreshMin)) then
!				print*,maxval(CurrentCountID%diff),(dble(nSnp)*.0001)*2
				SolutionChanged=1
			else
				SolutionChanged=0
			endif

			!write (*,'(2i4,3f10.3)') Windows,IterationNumber,GeneProbThresh,(dble(CurrentCountFilledPhase)/(dble(nInd*nSnp*2))*100),(dble(CurrentCountFilledGenos)/(dble(nInd*nSnp))*100)
			write (*,'(2i4,1f10.3,2i15)') Windows,IterationNumber,GeneProbThresh,CurrentCountFilledPhase,CurrentCountFilledGenos
			!print*,CurrentCountFilledPhase,CurrentCountFilledGenos,nInd,nSnp
		enddo

		call WriteResults

		StartSnp=EndSnp+1
		EndSnp=EndSnp+nSnp
		if (EndSnp>LenghtSequenceDataFile) then 
			EndSnp=LenghtSequenceDataFile
			nSnp=EndSnp-StartSnp+1
		endif
		Windows=Windows+1

		call DeallocateArrays

	enddo
	call UnintitialiseIntelRNG
	close(101)
	close(102)

	! Merge all files and Calculate Statistics

	Windows=Windows-1
	if (Windows>1) then 
		print*," Merge Output Files"
		call MergeResultsFile
	else 
		CALL RENAME("AlphaFamSeqFinalGenos1.txt", "AlphaFamSeqFinalGenos.txt") 
		CALL RENAME("AlphaFamSeqFinalPhase1.txt", "AlphaFamSeqFinalPhase.txt") 
		CALL RENAME("AlphaFamSeqEditingMarkersRemoved1.txt", "AlphaFamSeqEditingMarkersRemoved.txt") 
	endif

	! Compute statistics
	if ((trim(GenoFile)/="None").or.(trim(PhaseFile)/="None")) print*," Calculate Results"
	if (trim(GenoFile)/="None") 	call GetResultsImputation(LenghtSequenceDataFile,"AlphaFamSeqFinalGenos.txt",GenoFile,"AlphaFamSeqEditingMarkersRemoved.txt",1,"Yes","AlphaFamSeq")
	if (trim(PhaseFile)/="None") 	call GetResultsImputation(LenghtSequenceDataFile,"AlphaFamSeqFinalPhase.txt",PhaseFile,"AlphaFamSeqEditingMarkersRemoved.txt",2,"No","AlphaFamSeq")


end program FamilyPhase

!###########################################################################################################################################################


! At the end of each main step of the program, this subroutine fill the gaps
!  i.e. if a genotype is Homozygote but the phase is missing, then it fill the phase
subroutine SimpleCleanUpFillIn

	use GlobalPar
	use omp_lib
	implicit none

	integer i,j,Change 

	! TODO count the number of allele and genotype that are phased. Fill the complement and add the genotypes when no more new info are added

	call ped%phaseComplement()
	call ped%makeGenotype()

end subroutine SimpleCleanUpFillIn

!################################################################################################

! CoreOfTheProgram 5c - Find the consensus between founder-parent-kids haplotype
!  Here we use the chunks previously define, and we try to fill the missing value of these 
!  three individuals if one of the three has information to do that.

subroutine BuildConsensus 
	use GlobalPar

	implicit none

	integer :: i,k,e,j,m,a,ConsensusHaplotype
	integer(2) :: countAllele !0 and 1
	type(individual), pointer :: parent,grandparent
	
	do i=1,nInd
		do k=1,2 ! paternal or maternal gamete
			e=k+1 ! position parents in Pedigree

			if (maxval(FounderAssignment(i,:,k)).ne.0) then
				do j=1,nSnp
					if (FounderAssignment(i,j,k).ne.0) then

						parent => ped%pedigree(i)%getSireDamObjectbyIndex(e)												
						grandparent => parent%getSireDamObjectbyIndex(FounderAssignment(i,j,k))

						ConsensusHaplotype=9
						countAllele=0
						do a=0,1 
							if (grandparent%individualPhase(1)%getPhase+grandparent%individualPhase(2)%getPhase).lt.3) then
								! Avoid to use markers that are not fully phased for the grandparent
								do m=1,2 
									if (grandparent%individualPhase(m)%getPhase.eq.a) countAllele(a+1)=countAllele(a+1)+1
								enddo
							endif
							if (ped%pedigree(i)%individualPhase(k)%getPhase.eq.a) countAllele(a+1)+1
							if (parent%individualPhase((FounderAssignment(i,j,k)-1))%getPhase.eq.a) countAllele(a+1)+1
						enddo

						if ((countAllele(2).gt.1).and.(countAllele(2).gt.countAllele(1))) ConsensusHaplotype=1
						if ((countAllele(1).gt.1).and.(countAllele(1).gt.countAllele(2))) ConsensusHaplotype=0


						if (ConsensusHaplotype.ne.9) then
							call ped%pedigree(i)%individualPhase(e)%setPhase(j,ConsensusHaplotype)
							call parent%individualPhase((FounderAssignment(i,j,k)-1))%setPhase(j,ConsensusHaplotype)
						endif	
						! TODO: fill the complement phase or fill the gentoype
					endif
				enddo
			endif
		enddo
	enddo

end subroutine BuildConsensus 

!################################################################################################

! CoreOfTheProgram 5b - Work Left/Work Rigth to find chunks of haplotypes using the Founder Assignement array.
! 	i.e. From the previous step for the id7 we have the following 
!        seq of founders assignment: 1 0 0 0 1 0 0 0 2 0 0 0 2. In this step we do:
!
!        from Left to Rigth        : 1 1 1 1 1 1 1 1 2 2 2 2 2 
!        from Rigth to Left        : 1 1 1 1 1 2 2 2 2 2 2 2 2 and the consensus of these two vector generates
!        FounderAssignment 		   : 1 1 1 1 1 0 0 0 2 2 2 2 2 
!        where the 0 represent a region where the recombination happens.

! TODO : there can be a lot of small chunks. Find a way to understand if they are a single chunk or if 
!        the recombination happen in the middle of some of them.
!        Be carefull that we can have scenarios with females parents not sequenced.

subroutine ChunkDefinition

	use GlobalPar
	use IntelRNGMod
	implicit none

	integer :: i,j,e,k
	integer :: f1
	integer :: fSnp
	integer :: lSnp
	integer(int32) :: ChunkLength2(1)
	real(real32),allocatable,dimension(:) :: Z
	real(real32) ::	 A,B

	
	integer,allocatable,dimension(:) :: FounderAssignmentF,FounderAssignmentB

	allocate(FounderAssignmentF(1:nSnp))
	allocate(FounderAssignmentB(1:nSnp))

	
	A = ChunkLengthA  
	B = ChunkLengthB
	Z = SampleIntelUniformS(1,0.0,1.0) 
	ChunkLength2(1)=floor((A-1)+(B-(A-1))*Z(1))+1

	write(101,'(3(1x,i0))') Windows,IterationNumber,ChunkLength2(1)
	
	
	!!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (nInd,nSnp,FounderAssignment,ChunkLength2)
	do i=1,nInd
		do k=1,2 ! paternal or maternal gamete
			e=k+1 ! position parents in Pedigree
			
			FounderAssignmentF=0
			FounderAssignmentB=0

			if (maxval(FounderAssignment(i,:,k))/=0) then
				
				fSnp=1
				lSnp=ChunkLength2(1)
				
				
				!print*,"B",i,k,count(FounderAssignment(i,:,k)/=0),fSnp,lSnp
					

				do while (lSnp<=nSnp)

					if ((lSnp+ChunkLength2(1))>nSnp) lSnp=nSnp
					if ((fSnp-ChunkLength2(1))<0) fSnp=1
				
					FounderAssignmentF(fSnp:lSnp)=FounderAssignment(i,fSnp:lSnp,k)
					FounderAssignmentB(fSnp:lSnp)=FounderAssignment(i,fSnp:lSnp,k)

					f1=0
					!!! Forward founder assignment 
					do j=fSnp,lSnp
						if (FounderAssignmentF(j)/=0)then
							f1=FounderAssignmentF(j)
						endif
						if ((FounderAssignmentF(j)==0)) then
							FounderAssignmentF(j)=f1
						endif 
					enddo

					f1=0
					!!! Backward founder assignment 
					do j=lSnp,fSnp,-1
						if (FounderAssignmentB(j)/=0)then
							f1=FounderAssignmentB(j)
						endif
						if ((FounderAssignmentB(j)==0)) then
							FounderAssignmentB(j)=f1
						endif 
					enddo

					fSnp=lSnp+1
					lSnp=lSnp+ChunkLength2(1)
				enddo
				
				do j=1,nSnp
					if (FounderAssignmentF(j)==FounderAssignmentB(j)) FounderAssignment(i,j,k)=FounderAssignmentF(j)
					if (FounderAssignmentF(j)/=FounderAssignmentB(j)) FounderAssignment(i,j,k)=0
				enddo	

			endif
		enddo
	enddo
	!!$OMP END PARALLEL DO

	deallocate(FounderAssignmentF)
	deallocate(FounderAssignmentB)
end subroutine ChunkDefinition

!-------------------------------------------------------------------------------------------------
!> @brief   Assign from which founder the haplotype has been inherited.
!> @detail  This is done using the heterozygous markers of the individual’s parent
!>          i.e. DAD phase = 10, KID Phase=19. The "1" of the kid is coming from the Paternal GrandSire
!>			i.e. DAD phase = 01, KID Phase=19. The "1" of the kid is coming from the Paternal GrandDam
!> @date    August 30, 2017
!--------------------------------------------------------------------------------------------------   

subroutine CalculateFounderAssignment
	use GlobalPar
	use omp_lib
	implicit none

	integer :: i,j,e,k
	integer :: genotype
	integer(2) :: phaseId,phasePar
	type(individual),pointer :: parent

	FounderAssignment(:,:,:)=0
	
   	!$OMP PARALLEL DO ORDERED DEFAULT(PRIVATE) SHARED (ped,FounderAssignment,nSnp,nInd) !collapse(2)	
	do e=2,3 ! Sire and Dam pos in the ped
		do i=1,nInd
			if (.not. ped%pedigree(i)%Founder) then
			parent => ped%pedigree(i)%getSireDamObjectByIndex(e)
			
			do j=1,nSnp
				phase(1) = ped%pedigree(i)%individualPhase(1)%getPhase(j)
				phase(2) = ped%pedigree(i)%individualPhase(2)%getPhase(j)

				genotype = parent%individualGenotype%getGenotype(j)
				phasePar(1) =parent%individualPhase(1)%getPhase(j)
				phasePar(2) =parent%individualPhase(2)%getPhase(j)

				if ((genotype.eq.1).and.(sum(phasePar).lt.3)) then
					if (phasePar(1).eq.PhasePar(e-1)) FounderAssignment(i,j,e-1)=2 !GrandSire 
					if (phasePar(2).eq.PhasePar(e-1)) FounderAssignment(i,j,e-1)=3 !GrandDam
				endif
			enddo
		enddo
	enddo
    !$OMP END PARALLEL DO
end subroutine CalculateFounderAssignment

!-------------------------------------------------------------------------------------------------
!> @brief   Phase individuals' alleles on the basis of his parents/progeny genotyepes
!> @detail  if an individual has one of its two gametes phased and the evidence from 
!>          the progeny indicates that it inherited the other allele, then on the basis
!>			of progeny information the allele can be phased for the relevant gamete of the parent.
!>          i.e., The MGS=00, DAM=0?, KID=11. The genotype of the dam has to be 1 and the phase is 01
!>          !TODO : probably this subroutine is not necessary after the introduction of GeneProb.
!>          Check if it make sense to remove it.
!> @date    August 30, 2017
!--------------------------------------------------------------------------------------------------   

subroutine SimpleFillInBasedOnProgenyReads

	use GlobalPar
	use omp_lib
	implicit none

	integer :: i,j
	integer(2) :: phase,phaseSire,phaseDam
	type(individual),pointer :: sire,dam

	
	!$OMP PARALLEL DO ORDERED DEFAULT(PRIVATE) SHARED (ped,nSnp,nInd) !collapse(2)
	do i=1,nInd

		if (.not. ped%pedigree(i)%Founder) then
			sire => ped%pedigree(i)%sirePointer
			dam => ped%pedigree(i)%damPointer

			do j=1,nSnp
				phase(1) = ped%pedigree(i)%individualPhase(1)%getPhase(j)
				phase(2) = ped%pedigree(i)%individualPhase(2)%getPhase(j)

				phaseSire(1) =sire%individualPhase(1)%getPhase(j)
				phaseSire(2) =sire%individualPhase(2)%getPhase(j)

				phaseDam(1) = dam%individualPhase(1)%getPhase(j)
				phaseDam(2) = dam%individualPhase(2)%getPhase(j)

				if ((sum(phase).eq.0).or.(sum(phase).eq.2)) then
					if (phaseDam(1).ne.9).and.(phaseDam(1).ne.phase(2)) call dam%individualPhase(2)%setPhase(j,phase(2))
					if (phaseDam(2).ne.9).and.(phaseDam(2).ne.phase(2)) call dam%individualPhase(1)%setPhase(j,phase(2))

					if (phaseSire(1).ne.9).and.(phaseSire(1).ne.phase(1)) call sire%individualPhase(2)%setPhase(j,phase(1))
					if (phaseSire(2).ne.9).and.(phaseSire(2).ne.phase(1)) call sire%individualPhase(1)%setPhase(j,phase(2))
				endif

			enddo
		endif
	enddo
	!$OMP END PARALLEL DO

end subroutine SimpleFillInBasedOnProgenyReads

!---------------------------------------------------------------------------
!> @brief   Phasing based on parents informative markers
!> @detail  if a parent is homozygous at a variant, then the phase of the 
!>          progeny for the gamete inherited from that parent can be assigned.
!> @date    August 30, 2017
!---------------------------------------------------------------------------   

subroutine SimpleFillInBasedOnParentsReads

	use GlobalPar
	use omp_lib
	implicit none

	integer :: i,j,e
	type(Individual) , pointer :: parent
	integer :: genotype

    !!$OMP PARALLEL DO ORDERED DEFAULT(PRIVATE) SHARED (ped,nInd,nSnp) !collapse(2)
	do e=1,2
		do i=1,nInd
			parent => ped%pedigree(i)%getSireDamObjectByIndex(e+1)
			if (associated(parent)) then
				do j=1,nSnp
					genotype = parent%individualGenotype%getGenotype(j)
					if (genotype.eq.0) call ped%pedigree(i)%individualPhase(e)%setPhase(j,0)
					if (genotype.eq.2) call ped%pedigree(i)%individualPhase(e)%setPhase(j,1)
				enddo
			endif
		enddo
	enddo
    !!$OMP END PARALLEL DO

end subroutine SimpleFillInBasedOnParentsReads


!---------------------------------------------------------------------------
!> @brief   Use AlphaSLP to fill phase and genotypes
!> @detail  This step can run multiple times, at each iteration the genotype
!>          threshold is relaxed until the minimum value defined by the user
!> @date    August 30, 2017
!---------------------------------------------------------------------------   
 
subroutine UseGeneProbToSimpleFillInBasedOnOwnReads
	
	use GlobalPar
	use omp_lib
	
	implicit none

    integer :: i,j
	integer(2) :: phase
	real(kind=real64) :: p00,p01,p10,p11

	!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (ReadCounts,ped,GeneProbThresh,nSnp,nInd) !collapse(2)
	do i=1,nInd
		do j=1,nSnp

			phase(1) = ped%pedigree(i)%individualPhase(1)%getPhase(j)
			phase(2) = ped%pedigree(i)%individualPhase(2)%getPhase(j)

			p00=ReadCounts(i,j,1)
			p01=ReadCounts(i,j,2)
			p10=ReadCounts(i,j,3)
			p11=ReadCounts(i,j,4)
			
			if ((sum(phase).gt.3).and.(ped%pedigree(i)%individualGenotype%isMissing(j))) then
				if (p00.ge.GeneProbThresh) then
					ped%pedigree(i)%individualGenotype%setGenotype(j,0)
					ped%pedigree(i)%individualPhase(1)%setPhase(j,0)
					ped%pedigree(i)%individualPhase(2)%setPhase(j,0)
				endif

				if (p11.ge.GeneProbThresh) then
					ped%pedigree(i)%individualGenotype%setGenotype(j,2)
					ped%pedigree(i)%individualPhase(1)%setPhase(j,1)
					ped%pedigree(i)%individualPhase(2)%setPhase(j,1)
				endif


				if (p01+p10).ge.GeneProbThresh)) then
					ped%pedigree(i)%individualGenotype%setGenotype(j,1)
					if (p01.ge.GeneProbThresh) then
						ped%pedigree(i)%individualPhase(1)%setPhase(j,0)
						ped%pedigree(i)%individualPhase(2)%setPhase(j,1)
					endif
					if (p10.ge.GeneProbThresh) then
						ped%pedigree(i)%individualPhase(1)%setPhase(j,1)
						ped%pedigree(i)%individualPhase(2)%setPhase(j,0)
					endif
				endif
			endif

			if (((p00+p01).ge.GeneProbThresh).and.(p10.lt.p01).and.(phase(1).eq.9)) ped%pedigree(i)%individualPhase(1)%setPhase(j,0)
			if (((p11+p10).ge.GeneProbThresh).and.(p01.lt.p10).and.(phase(1).eq.9)) ped%pedigree(i)%individualPhase(1)%setPhase(j,1)

			if (((p00+p10).ge.GeneProbThresh).and.(p01.lt.p10).and.(phase(2).eq.9)) ped%pedigree(i)%individualPhase(2)%setPhase(j,0)
			if (((p11+p01).ge.GeneProbThresh).and.(p10.lt.p01).and.(phase(2).eq.9)) ped%pedigree(i)%individualPhase(2)%setPhase(j,1)

		enddo
	enddo
	!$OMP END PARALLEL DO
end subroutine UseGeneProbToSimpleFillInBasedOnOwnReads


!################################################################################################

subroutine AllocateArrays

	use GlobalPar

	implicit none

	! Founders
	allocate(FounderAssignment(nInd,nSnp,2))

	! CountPhase
	allocate(CurrentCountID%old(nInd))
	allocate(CurrentCountID%diff(nInd))

	! Markers to exclude
	allocate(MarkersToExclude(nSnp))

end subroutine AllocateArrays

!################################################################################################

subroutine DeallocateArrays

	use GlobalPar

	implicit none

	deallocate(RawReads)
	IF( ALLOCATED(position)) DEALLOCATE(position) 
	IF( ALLOCATED(quality)) DEALLOCATE(quality) 

	! Founders
	deallocate(FounderAssignment)
	
	! CountPhase
	if (allocated(CurrentCountID%old)) deallocate(CurrentCountID%old)
	if (allocated(CurrentCountID%diff)) deallocate(CurrentCountID%diff)

	! Markers to exclude
	deallocate(MarkersToExclude)

end subroutine DeallocateArrays

!################################################################################################

subroutine CurrentCountFilled

	use GlobalPar

	implicit none

	integer :: i,new

	CurrentCountFilledPhase=count(FilledPhase(:,:,:)/=9)
	CurrentCountFilledGenos=count(FilledGenos(:,:)/=9)

	do i=1,nInd
		new=0
		new=count(FilledPhase(i,:,:)/=9)
		CurrentCountID%diff(i)=new-CurrentCountID%old(i)
		CurrentCountID%old(i)=new
	enddo
end subroutine CurrentCountFilled

!################################################################################################

subroutine ReadData
 	use ISO_Fortran_Env

  	use GlobalPar
	use MaraModule
	use omp_lib

  	implicit none
  
	integer :: i
	integer(int64) :: PosReads
	integer,allocatable,dimension(:) :: TempImput
	integer(int64) :: TmpID
	real(kind=8)::tstart,tend

	!For CurrentCount
	!open (unit=99,file="AlphaFamSeqSummary.log",status="unknown")


	if (trim(ReadsType)=="VcfTools") call readRogerData(ReadsFile, Ids, position, quality, SequenceData,LenghtSequenceDataFile,nSnp,StartSnp,EndSnp,nIndSeq)
	if (trim(ReadsType)=="AlphaSim") call readAlphaSimReads(ReadsFile, Ids,position,SequenceData,LenghtSequenceDataFile,nSnp,StartSnp,EndSnp,nIndSeq)
	
	allocate(RawReads(nInd,nSnp,2))
	RawReads(:,:,:)=0

    tstart = omp_get_wtime()
	!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) PRIVATE(i) SHARED(nIndSeq,Ids,RawReads,SequenceData)
	do i=1, nIndSeq
		PosReads=0
		read(Ids(i),*) TmpID
	    ! call GetID(TmpID, PosReads) 
	    !print*,i,TmpID,PosReads
	    if (PosReads/=0) RawReads(PosReads,:,:)= SequenceData(i,:,:)
	enddo
	!$OMP END PARALLEL DO
	tend = omp_get_wtime()
	write(*,*) "Total wall time for Reads sorting ", tend - tstart
	
	deallocate(SequenceData)
	allocate(TempImput(LenghtSequenceDataFile))
	call AllocateArrays
end subroutine ReadData

!###########################################################################################################################################################

subroutine MergeResultsFile
 	use ISO_Fortran_Env
  	use GlobalPar
	
	use AlphaSortMod

  	implicit none
  
	integer :: i,j,k,l,m,DumI,nSnpWindow,MaxNrMissingSnp
	integer(kind=1),allocatable,dimension(:) :: TempImput,FinalOutput
	integer(int32),allocatable,dimension(:,:) :: MissingSnpTmp,MissingSnp
	integer,allocatable,dimension(:,:) :: WindowsInfo

	character(len=80) :: filout1,filout2,nChar,FmtInt
	!integer :: TmpID
	!real(kind=8)::tstart,tend


	allocate(WindowsInfo(Windows,2))
	open (unit=1,file="AlphaFamSeqWindowsInfo.txt",status="old")

	write(nChar,*) LenghtSequenceDataFile
	FmtInt='(i10,'//trim(adjustl(nChar))//'i2)'
	

	open (unit=5,file="AlphaFamSeqFinalPhase.txt",status="unknown")
	open (unit=6,file="AlphaFamSeqFinalGenos.txt",status="unknown")
	open (unit=7,file="AlphaFamSeqEditingMarkersRemoved.txt",status="unknown")


	read(1,*)

	do i=1,Windows
		read(1,*) DumI,WindowsInfo(i,1) ,WindowsInfo(i,2)
		!print*,WindowsInfo(i,1),WindowsInfo(i,2)
	enddo

	close(1)

	allocate(FinalOutput(LenghtSequenceDataFile))
	
	! Genotypes file
	do i=1,nInd

		FinalOutput=9
		do j=1,Windows
			write (filout1,'("AlphaFamSeqFinalGenos",i0,".txt")') j
			open (unit=2,file=trim(filout1),status="unknown")
			
			if (i>1) then
				do m=1,(i-1)
					read(2,*)
				enddo
			endif


			nSnpWindow=WindowsInfo(j,2)-WindowsInfo(j,1)+1
			
			allocate(TempImput(nSnpWindow))
			TempImput=9
			
			read(2,*), DumI,TempImput(:)

			l=0
			do k=WindowsInfo(j,1),WindowsInfo(j,2)
				l=l+1
				if ((FinalOutput(k)==9).and.(TempImput(l)/=9)) then
					FinalOutput(k)=TempImput(l)
				endif
				if ((FinalOutput(k)/=9).and.(TempImput(l)/=9)) then
					if  (FinalOutput(k)/=TempImput(l)) then
						FinalOutput(k)=9
					endif
				endif
			enddo
			deallocate(TempImput)
			close(2)
		enddo
		
		write(6,FmtInt) DumI,FinalOutput(:)
	enddo
	close(6)

	! Phase File
	do i=1,(nInd*2)
	
		FinalOutput=9
		do j=1,Windows
			write (filout2,'("AlphaFamSeqFinalPhase",i0,".txt")') j
			open (unit=3,file=trim(filout2),status="unknown")
			
			if (i>1) then
				do m=1,(i-1)
					read(3,*)
				enddo
			endif

			nSnpWindow=WindowsInfo(j,2)-WindowsInfo(j,1)+1
			
			allocate(TempImput(nSnpWindow))
			TempImput=9
			
			read(3,*), DumI,TempImput(:)

			l=0
			do k=WindowsInfo(j,1),WindowsInfo(j,2)
				l=l+1
				if ((FinalOutput(k)==9).and.(TempImput(l)/=9)) then
					FinalOutput(k)=TempImput(l)
				endif
				if ((FinalOutput(k)/=9).and.(TempImput(l)/=9)) then
					if  (FinalOutput(k)/=TempImput(l)) then
						FinalOutput(k)=9
					endif
				endif
			enddo
			deallocate(TempImput)
			close(3)	
		enddo
		write(5,FmtInt) DumI,FinalOutput(:)
	enddo
	close(5)

	! Markers with zero reads
	MaxNrMissingSnp=LenghtSequenceDataFile+(Windows*ChunkLengthB*2)
	allocate(MissingSnpTmp(MaxNrMissingSnp,3))
	MissingSnp=0
	i=1
	do j=1,Windows
		write (filout1,'("AlphaFamSeqEditingMarkersRemoved",i0,".txt")') j
		open (unit=4,file=trim(filout1),status="unknown")
		
		do 
			read(4,*,END=904) MissingSnpTmp(i,1),MissingSnpTmp(i,2),MissingSnpTmp(i,3)
			i=i+1
		enddo
		904 continue
		close(4)
	enddo

	i=i-1
	allocate(MissingSnp(1:i,3))
	MissingSnp=MissingSnpTmp(1:i,:)
	deallocate(MissingSnpTmp)
	call HpSortI(i,MissingSnp(:,1))
	write(7,'(1i20,1i10,1i4)') MissingSnp(1,:) ! Problems with array bound if I use a single loop.
	do j=2,i
		if ((j.gt.1).and.(MissingSnp(j,1)/=MissingSnp(j-1,1))) then
			write(7,'(1i20,1i10,1i4)') MissingSnp(j,:)
		endif
	enddo
	close(7)
	
	deallocate(MissingSnp)
	
	do i=1,Windows
		write (filout1,'("AlphaFamSeqFinalGenos",i0,".txt")') i
		open (unit=2,file=trim(filout1),status="unknown")
		close(2,status="delete")
		
		write (filout2,'("AlphaFamSeqFinalPhase",i0,".txt")') i
		open (unit=3,file=trim(filout2),status="unknown")
		close(3,status="delete")
		
		write (filout2,'("AlphaFamSeqEditingMarkersRemoved",i0,".txt")') i
		open (unit=4,file=trim(filout2),status="unknown")
		close(4,status="delete")
	enddo
end subroutine MergeResultsFile

!###########################################################################################################################################################
subroutine ReadPrevGeneProb
	use GlobalPar

	implicit none

	integer(int64) :: i,j,DumI
	character(len=30) :: nChar
	character(len=80) :: FmtReal,filout5
	
	filout5="AlphaFamSeqFinalGeneProb1.bin"
	open (unit=5,file=trim(filout5),status="old",form="unformatted")
	
	do i=1,nInd
	
		read(5) DumI,Pr00(i,:)
		read(5) DumI,Pr01(i,:)
		read(5) DumI,Pr10(i,:)
		read(5) DumI,Pr11(i,:)
	
	enddo

	close (3)
	close (5)
end subroutine ReadPrevGeneProb

!###########################################################################################################################################################

subroutine SaveGeneProbResults
	use GlobalPar

	implicit none

	integer :: i,j,nRow,p,stat,DumI
	character(len=30) :: nChar
	character(len=80) :: FmtReal,filout5,filout6
	
	!integer,allocatable,dimension(:) :: SeqColToKeep
	!real,allocatable,dimension(:) 	 :: ReducedGeneProb

	! WriteOut ReducedFileOf GeneProb
	!open (unit=2,file="SeqColToKeep.txt",status="old")
	
	! nRow = 0
	! do
	!     read(2, *, iostat=stat) DumI
	!     if (stat/=0) exit
	!     nRow = nRow + 1
	! enddo
	! rewind(2)

	!allocate(SeqColToKeep(nRow))
	!allocate(ReducedGeneProb(nRow))

	! SeqColToKeep(:)=0
	! ReducedGeneProb(:)=0

	! do i=1,nRow
	! 	read(2,*) SeqColToKeep(i)
	! enddo
	! close (2)


	! write(nChar,*) nRow
	! FmtReal='(i0,'//trim(adjustl(nChar))//'f7.4)'
	! write (filout6,'("AlphaFamSeqReducedGeneProb",i0,".txt")') Windows
	! open (unit=6,file=trim(filout6),status="unknown")

	! do i=1,nInd
	! 	p=1
	! 	do j=1,nRow
	! 		ReducedGeneProb(p)=Pr00(i,SeqColToKeep(j))
	! 		p=p+1
	! 	enddo
	! 	write (6,FmtReal) Ped(i,1),ReducedGeneProb(:)
	! 	p=1
	! 	do j=1,nRow
	! 		ReducedGeneProb(p)=Pr01(i,SeqColToKeep(j))
	! 		p=p+1
	! 	enddo
	! 	write (6,FmtReal) Ped(i,1),ReducedGeneProb(:)
	! 	p=1
	! 	do j=1,nRow
	! 		ReducedGeneProb(p)=Pr10(i,SeqColToKeep(j))
	! 		p=p+1
	! 	enddo
	! 	write (6,FmtReal) Ped(i,1),ReducedGeneProb(:)
	! 	p=1
	! 	do j=1,nRow
	! 		ReducedGeneProb(p)=Pr11(i,SeqColToKeep(j))
	! 		p=p+1
	! 	enddo
	! 	write (6,FmtReal) Ped(i,1),ReducedGeneProb(:)
	! enddo
	! close(6)

	! deallocate(SeqColToKeep)
	! deallocate(ReducedGeneProb)

	! Write Out Full file of GeneProb

	write(nChar,*) nSnp
	FmtReal='(i0,'//trim(adjustl(nChar))//'f7.4)'
	!write (filout5,'("AlphaFamSeqFinalGeneProb",i0,".txt")') Windows
	! open (unit=5,file=trim(filout5),status="unknown")

	! do i=1,nInd
	! 	write (5,FmtReal) Ped(i,1),Pr00(i,:)
	! 	write (5,FmtReal) Ped(i,1),Pr01(i,:)
	! 	write (5,FmtReal) Ped(i,1),Pr10(i,:)
	! 	write (5,FmtReal) Ped(i,1),Pr11(i,:)
	! enddo

	write (filout5,'("AlphaFamSeqFinalGeneProb",i0,".bin")') Windows
	open (unit=5,file=trim(filout5),status="unknown",form="unformatted")

	do i=1,nInd
		write (5) (Ped(i,1)),Pr00(i,:)
		write (5) (Ped(i,1)),Pr01(i,:)
		write (5) (Ped(i,1)),Pr10(i,:)
		write (5) (Ped(i,1)),Pr11(i,:)
	enddo


	close (5)

end subroutine SaveGeneProbResults

!###########################################################################################################################################################

subroutine WriteResults

	use GlobalPar

	implicit none

	integer :: i,j,nRow,stat,DumI,p
	character(len=30) :: nChar
	character(len=80) :: FmtInt,FmtInt2,FmtCha,FmtReal,filout1,filout2,filout3,filout4,filout6,filout7
		
	! integer,allocatable,dimension(:) 			:: SeqColToKeep
	! integer(kind=1),allocatable,dimension(:) 	:: ReducedGenos

	! ! Save Reduced Genotypes
	! open (unit=6,file="SeqColToKeep.txt",status="old")
	
	! nRow = 0
	! do
	!     read(6, *, iostat=stat) DumI
	!     if (stat/=0) exit
	!     nRow = nRow + 1
	! enddo
	! rewind(6)

	! allocate(SeqColToKeep(nRow))
	! allocate(ReducedGenos(nRow))
	
	! SeqColToKeep(:)=0
	! ReducedGenos(:)=9

	! do i=1,nRow
	! 	read(6,*) SeqColToKeep(i)
	! enddo
	! close (6)


	! write(nChar,*) nRow
	! FmtInt='(i0,'//trim(adjustl(nChar))//'i2)'
	! write (filout7,'("AlphaFamSeqReducedGenos",i0,".txt")') Windows
	! open (unit=7,file=trim(filout7),status="unknown")

	! do i=1,nInd
	! 	p=1
	! 	do j=1,nRow
	! 		ReducedGenos(p)=FilledGenos(i,SeqColToKeep(j))
	! 		p=p+1
	! 	enddo
	! 	write (7,FmtInt) Ped(i,1),ReducedGenos(:)
	! enddo
	! close(7)

	! deallocate(SeqColToKeep)
	! deallocate(ReducedGenos)

	! WriteOut Full Output

	write(nChar,*) nSnp
	FmtInt='(i0,'//trim(adjustl(nChar))//'i2)'
	FmtInt2='(i0,'//trim(adjustl(nChar))//'(1x,i0))'
	FmtCha='(i0,a2,'//trim(adjustl(nChar))//'a2)'
	FmtReal='(i0,'//trim(adjustl(nChar))//'f7.4)'
	
	write (filout1,'("AlphaFamSeqFinalPhase",i0,".txt")') Windows
	write (filout2,'("AlphaFamSeqFinalGenos",i0,".txt")') Windows
	write (filout4,'("AlphaFamSeqFounderAssignment",i0,".txt")') Windows
	
	open (unit=1,file=trim(filout1),status="unknown")
	open (unit=2,file=trim(filout2),status="unknown")
	open (unit=4,file=trim(filout4),status="unknown")
	

	do i=1,nInd
		write (1,FmtInt) Ped(i,1),FilledPhase(i,:,1)
		write (1,FmtInt) Ped(i,1),FilledPhase(i,:,2)
		write (2,FmtInt) Ped(i,1),FilledGenos(i,:)

		if (maxval(FounderAssignment(i,:,:))/=0) then 
			write (4,FmtInt2) Ped(i,1),FounderAssignment(i,:,1)
			write (4,FmtInt2) Ped(i,1),FounderAssignment(i,:,2)
		endif
	enddo

	close (1)
	close (2)
	close (4)

end subroutine WriteResults

!################################################################################################

subroutine ReadsLikelihood(nRef,nAlt,ErrorRate,Pr0,Pr1,Pr2)
	implicit none
	
	real,intent(in) :: nRef,nAlt
	real(kind=8),intent(in) :: ErrorRate
	real,intent(inout) :: Pr0,Pr1,Pr2
	real :: lPr0,lPr1,lPr2

      if ((nRef+nAlt)==0) then
        lPr0=log(1.)
        lPr1=log(1.)
        lPr2=log(1.)
      else
        lPr0=(nRef*log(1-ErrorRate))+(nAlt*log(ErrorRate))
        if (lPr0.lt.log(.000000001)) lPr0=-9999
        lPr1=(nRef*log(0.5))+(nAlt*log(0.5))
        if (lPr1.lt.log(.000000001)) lPr1=-9999
        lPr2=(nAlt*log(1-ErrorRate))+(nRef*log(ErrorRate))
        if (lPr2.lt.log(.000000001)) lPr2=-9999
      endif

      Pr0=exp(lPr0)/(exp(lPr0)+exp(lPr1)+exp(lPr2))
      Pr1=exp(lPr1)/(exp(lPr0)+exp(lPr1)+exp(lPr2))
      Pr2=exp(lPr2)/(exp(lPr0)+exp(lPr1)+exp(lPr2))
end subroutine ReadsLikelihood

!###########################################################################################
subroutine InternalEdititing
	! 	use GlobalPar
	! 	use omp_lib
	! 	implicit none

	! 	integer :: i,j
	! 	real(kind=4),allocatable,dimension(:) :: RefAllele
	! 	real(kind=4) :: maf

	! 	allocate(RefAllele(nSnp))
	! 	allocate(GeneProbYesOrNo(nSnp))
	! 	GeneProbYesOrNo=1
	! 	i=0

	! 	open (unit=1,file="AlphaFamSeqEditedSnp.txt",status="unknown")
	! 	print*,"Internal Editing"
	! 	write(*,'(1a10,1f10.3)'),"(1) MAF",EditingParameter
		
	! 	!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nSnp,RawReads,RefAllele,EditingParameter,GeneProbYesOrNo,i)
	! 	do j=1,nSnp
	! 		RefAllele(j)=real(sum(RawReads(:,j,1)))
	! 		maf=RefAllele(j)/sum(RawReads(:,j,:))
	! 		if ((maf.le.EditingParameter).or.(maf.ge.(1-EditingParameter))) then
	! 			GeneProbYesOrNo(j)=0
	! 			i=i+1
	! 			write(1,'(1i10,1f10.3,2i10)'),j,maf,sum(RawReads(:,j,1)),sum(RawReads(:,j,2))
	! 		endif
	! 	end do
	! 	!$OMP END PARALLEL DO
		
	! 	close(1)

	! 	if (i==0) print*,"All SNPs passed the editing"
	! 	if (i>0)  write(*,'(1i12,1a3,1i12,1a30)'), i,"on",nSnp,"didn't pass the editing"

	! 	deallocate(RefAllele)
end subroutine InternalEdititing

!###########################################################################################################################################################

! This step perform 3 editing (at the moment):
! 1) Markers with no reads (it happens with simulated data or with a subset of real data) 
! 2) Remove variants with less than nUserDefined reads for the referece or alternative allele (i.e., single- and double-tones)
! 3) Remove ReadsCount/invidual that are above a user defined threshold (i.e., excess of reads compare to the average)
! 4) Remove variants with an excess of heterozygotes using the Exact Tests of HWE described in Wigginton et al., Am.J.Hum.Genet. 76:887-883,2005

! Moreover, the coverage/individual and coverage/marker are calculated in this step
subroutine CleanUpTheRawData

	! 	use GlobalPar
	! 	use omp_lib
	! 	implicit none

	! 	integer :: i,j,nReadsRemoved,nTmpInd,e,pos
	! 	real 	:: cov(nInd)
	! 	real 	:: std,covSnp
	! 	character(len=50) :: filout1,filout2,filout3,filout4,filout5
	! 	character(len=30) :: nChar
	! 	character(len=80) :: FmtInt

	! 	integer(kind=2),allocatable,dimension(:) :: ReadsRemoved
	! 	integer(kind=2),allocatable,dimension(:,:) :: tmpReads

	! 	real(kind=8)					 				:: pHetExcess
	! 	integer											:: ObsGenos(3), EstGenos(3) ! observed genotypes


	! 	write (filout1,'("AlphaFamSeqEditingIndividualCoverage",i0,".txt")') Windows
	! 	open (unit=1,file=trim(filout1),status="unknown")

	! 	write (filout2,'("AlphaFamSeqEditingMarkerCoverage",i0,".txt")') Windows
	! 	open (unit=2,file=trim(filout2),status="unknown")

	! 	write (filout3,'("AlphaFamSeqEditingMarkersRemoved",i0,".txt")') Windows
	! 	open (unit=3,file=trim(filout3),status="unknown")

	! 	write (filout4,'("AlphaFamSeqEditingIndividualReadsRemoved",i0,".txt")') Windows
	! 	open (unit=4,file=trim(filout4),status="unknown")

	! 	write (filout5,'("AlphaFamSeqEditingPvalueExcessOfHeterozygotes",i0,".txt")') Windows
	! 	open (unit=5,file=trim(filout5),status="unknown")

	! 	!write (filout3,'("AlphaFamSeqReads",i0,".txt")') Windows
	! 	!open (unit=3,file=trim(filout3),status="unknown")

	! 	allocate(ReadsRemoved(nSnp))
	! 	ReadsRemoved(:)=0
	! 	MarkersToExclude=0 ! 0=Use The marker ; 1=Don't use the marker

	! 	do j=1,nSnp
	! 		covSnp=0
	! 		!e=0
	! 		covSnp=sum(RawReads(:,j,:))/dble(nTmpInd)
	! 		!if (cov.lt.EditingParameter) then
	! 		!	RawReads(:,j,:)=0
	! 		!	e=1
	! 		!endif
	! 		write (2,'(1i0,1x,1f7.3,1x,1i0)') j,covSnp,e
	! 	enddo
	! 	close(2)


	! 	nTmpInd=0
	! 	do i=1,nInd
	! 		cov(i)=0
	! 		do j=1,nSnp
	! 			cov(i)=cov(i)+sum(RawReads(i,j,:))
	! 		enddo
	! 		if (cov(i).gt.0) nTmpInd=nTmpInd+1
	! 		cov(i)=cov(i)/dble(nSnp)

	! 		std=0
	! 		do j=1,nSnp
	! 			std=std+((sum(RawReads(i,j,:))-cov(i))**2)
	! 		enddo
	! 		std=sqrt(std/(dble(nSnp)-1))

	! 		nReadsRemoved=0
	! 		if (cov(i).gt.0.0) then
	! 			do j=1,nSnp
	! 				if (((sum(RawReads(i,j,:))-cov(i))/std).gt.maxStdForReadsCount) then
	! 					write(4,'(2i20,1f10.6,2i4)'),Ped(i,1),j,cov(i),RawReads(i,j,:)
	! 					RawReads(i,j,:)=0
	! 					nReadsRemoved=nReadsRemoved+1
	! 					ReadsRemoved(j)=ReadsRemoved(j)+1
	! 				endif
	! 			enddo
	! 		endif

	! 		write (1,'(1i0,1f7.3,1f10.6)') Ped(i,1),cov(i),dble(nReadsRemoved)/dble(nSnp)*100

	! 		!write (3,FmtInt) Ped(i,1), RawReads(i,:,1)
	! 		!write (3,FmtInt) Ped(i,1), RawReads(i,:,2)
	! 	enddo

	! 	close(1)
	! 	close(4)


	! 	allocate(tmpReads(nTmpInd,2))
		
	! 	!!$OMP PARALLEL DO ORDERED DEFAULT(PRIVATE) SHARED (RawReads,nSnp)
	! 	do j=1,nSnp
			
	! 		! 1) Markers with zero reads	
	! 		if ((maxval(RawReads(:,j,:))==0).or.(ReadsRemoved(j).gt.(dble(nInd)*ThresholdMaxReadsCount))) then
	! 			write (3,'(1i10,1i20,1i4)') j,position(j),0
	! 			MarkersToExclude(j)=1
	! 			RawReads(:,j,:)=0
	! 		endif

	! 		! 2) Remove single/double-tones and excess heterozygotes
			
	! 		tmpReads=9
	! 		ObsGenos=0
	! 		EstGenos=0
	! 		pHetExcess=0.0

	! 		pos=1
	! 		do i=1,nInd
	! 			if (cov(i).gt.0.0) then
	! 				tmpReads(pos,:)=RawReads(i,j,:)
	! 				pos=pos+1
	! 			endif
	! 		enddo

	! 		call ExcessHeterozygotes(tmpReads(:,:),nTmpInd,ObsGenos,EstGenos,pHetExcess)

	! 		if (((ObsGenos(2)+dble(2*ObsGenos(1))).le.ThresholdReadsCount).or.((ObsGenos(2)+dble(2*ObsGenos(3))).le.ThresholdReadsCount)) then
	! 			if (MarkersToExclude(j).ne.1) then
	! 				write (3,'(1i10,1i20,1i4)') j,position(j),2
	! 				MarkersToExclude(j)=1
	! 			endif
	! 		endif

	! 		if (MarkersToExclude(j).ne.1) then
	! 			if (pHetExcess.le.ThresholdExcessHetero) then
	! 				write (3,'(1i10,1i20,1i4)') j,position(j),1
	! 				write(5,'(1i20,3i4,5x,3i4,1f10.5)') position(j),ObsGenos(:),EstGenos(:),pHetExcess

	! 				MarkersToExclude(j)=1
			
	! 			endif
	! 		endif

	! 	enddo
	! 	!!$OMP END PARALLEL DO

	! 	write(*,'(1a10,1x,1i0,1a4,1x,1i0,1a10)'),"Exclude",count(MarkersToExclude(:)==1),"on",nSnp,"Markers"

	! 	deallocate(tmpReads)
		
	! 	close (3)
	! 	close (5)
end subroutine CleanUpTheRawData


!###########################################################################################################################################################
subroutine ExcessHeterozygotes(ReadsCount,n,ObsGenos,EstGenosInit,pHetExcess)

	implicit none

	integer,intent(in) 						:: n ! Nr of individuals
	integer(kind=2),intent(in)				:: ReadsCount(n,2) ! Reads for 1 variant and all individuals
	real(kind=8),intent(inout) 				:: pHetExcess
	integer,intent(inout)					:: ObsGenos(3),EstGenosInit(3) ! observed  & estimated genotypes
	
	integer									:: EstGenos(3) ! estimated genotypes
	
	integer									:: i,nGeno,RareHomo,CommomHomo,RareCopies,mid,posEst,posObs,FirstInd
	real(kind=8),allocatable,dimension(:) 	:: probs,plow,phigh
	real(kind=8)							:: SumProbs

	ObsGenos=0
	EstGenos=0
	nGeno=0

	do i=1,n
		if (sum(ReadsCount(i,:)).gt.0) then
			if ((ReadsCount(i,1).gt.0).and.(ReadsCount(i,2).eq.0)) ObsGenos(1)=ObsGenos(1)+1
			if ((ReadsCount(i,1).gt.0).and.(ReadsCount(i,2).gt.0)) ObsGenos(2)=ObsGenos(2)+1
			if ((ReadsCount(i,1).eq.0).and.(ReadsCount(i,2).gt.0)) ObsGenos(3)=ObsGenos(3)+1
			nGeno=nGeno+1
		endif
	enddo

	RareHomo=3
	CommomHomo=1
	if (ObsGenos(1).lt.ObsGenos(2)) then
		RareHomo=1 ! if the number of HomoReferenceAllele < HomoAlternativeAllele switch the common and rare allele
		CommomHomo=3
	endif

	RareCopies=nint(2.0*dble(ObsGenos(RareHomo))+ObsGenos(2))
	
	if (RareCopies.lt.2) then
		pHetExcess=1
	else if (RareCopies.ge.2) then

		mid=nint(dble(RareCopies)*dble(2*nGeno-RareCopies)/dble(2*nGeno))
		
		if (mid.lt.0) mid=0

		if ((mod(RareCopies,2).eq.0).and.(mod(mid,2).ne.0)) mid=mid+1
		if ((mod(RareCopies,2).ne.0).and.(mod(mid,2).eq.0)) mid=mid+1

		EstGenos(2)=mid
		EstGenos(RareHomo)=floor(dble(RareCopies-mid)/2.0)
		if (EstGenos(RareHomo).lt.0) EstGenos(RareHomo)=0
		EstGenos(CommomHomo)=nGeno-EstGenos(RareHomo)-EstGenos(2)
		if (EstGenos(CommomHomo).lt.0) EstGenos(CommomHomo)=0

		EstGenosInit(:)=EstGenos(:)
		
		! # we observed 21 copies of the minor allele (RareCopies) --> the observed nr of hetero  (ObsGenos(2)) will vary seq(1,21,2)

		! The number of possible heterozygotes is odd if RareCopies is odd, and is even if RareCopies is even
		if (mod(RareCopies,2).eq.0) then
			allocate(probs(0:(RareCopies/2)))
			allocate(plow(0:(RareCopies/2)))
			
			posEst=(mid)/2
			posObs=ObsGenos(2)/2
			FirstInd=0
		else if(mod(RareCopies,2).eq.1) then
			allocate(probs(1:((RareCopies+1)/2)))
			allocate(plow(1:((RareCopies+1)/2)))
			
			posEst=(mid+1)/2
			posObs=(ObsGenos(2)+1)/2
			FirstInd=1
		endif

		probs(posEst)=1.0
		SumProbs=probs(posEst)

		! Start to calculate the probabilities using the equations 2 of Am.J.Hum.Genet.76:887-883,2005
		!
		!!!  P(NAB=nAB-2|N,nA) = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
		!!!  P(NAB=nAB+2|N,nA) = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc	/((curr_hets + 2.0) * (curr_hets + 1.0))
		do i=posEst,(FirstInd+1),-1
			probs(i-1)=probs(i)*dble(EstGenos(2))*dble(EstGenos(2)-1.0)/(4.0*dble(EstGenos(RareHomo)+1.0)*dble(EstGenos(CommomHomo)+1.0))
			EstGenos(2)=EstGenos(2)-2
			EstGenos(1)=EstGenos(1)+1
			EstGenos(3)=EstGenos(3)+1
			SumProbs=SumProbs+probs(i-1)
		enddo

		EstGenos(2)=mid
		EstGenos(RareHomo)=floor(dble(RareCopies-mid)/2.0)
		if (EstGenos(RareHomo).lt.0) EstGenos(RareHomo)=0
		EstGenos(CommomHomo)=nGeno-EstGenos(RareHomo)-EstGenos(2)
		if (EstGenos(CommomHomo).lt.0) EstGenos(CommomHomo)=0

		if (EstGenos(2).lt.RareCopies) then
			do i=posEst,(size(probs)-(2-FirstInd))
				probs(i+1)=probs(i)*4.0*dble(EstGenos(CommomHomo))*dble(EstGenos(RareHomo))/(dble(EstGenos(2)+2.0)*dble(EstGenos(2)+1.0))
				EstGenos(2)=EstGenos(2)+2
				EstGenos(1)=EstGenos(1)-1
				EstGenos(3)=EstGenos(3)-1
				SumProbs=SumProbs+probs(i+1)
			enddo
		endif
		
	!		if ((ObsGenos(1)==611).and.(ObsGenos(2)==0).and.(ObsGenos(3)==1)) print*,probs(:)


		probs(:)=probs(:)/SumProbs
		plow=1
		plow(FirstInd)=probs(FirstInd)

	!		if ((ObsGenos(1)==611).and.(ObsGenos(2)==0).and.(ObsGenos(3)==1)) print*,probs(:)

		do i=(FirstInd+1),(size(probs)-1)
			plow(i)=plow(i-1)+probs(i)
		enddo



	!		if ((ObsGenos(1)==611).and.(ObsGenos(2)==0).and.(ObsGenos(3)==1)) print*,plow(:),posObs
		if (FirstInd.eq.0) then
			if (posObs.gt.0) pHetExcess=1-plow(posObs-1)
			if (posObs.eq.0) pHetExcess=1
		endif

		if (FirstInd.eq.1) then
			if (posObs.gt.1) pHetExcess=1-plow(posObs-1)
			if (posObs.eq.1) pHetExcess=1
		endif

		deallocate(probs)
		deallocate(plow)
	endif
end subroutine ExcessHeterozygotes

!###########################################################################################################################################################
 
subroutine ReadSamFile
	
	! 	use GlobalPar
	! 	use omp_lib
		
	! 	implicit none

	! 	integer :: i,j,k,HapLength,tmpId,wp,wm,cp,cm
	! 	integer,allocatable,dimension(:) :: SamHap,HapPos
	! 	real 	:: LikeP,LikeM
	! 	real(kind=8) :: p0,m0
	! 	character(len=100) :: filout12
	! 	logical :: exist12

	! 	allocate(HapPos(5000))
	! 	allocate(SamHap(5000))


	! 	do i=nInd,1,-1
	! 		tmpId=0
	! 		HapLength=0
	! 		HapPos=0
	! 		SamHap=9

	! 		cp=0
	! 		cm=0
	! 		wp=0
	! 		wm=0
	! 		write (filout12,'(i0,".FilledPhase.txt")') Ped(i,1)
	! 		inquire  (file=trim(filout12),exist=exist12)
	! 		if (exist12) then
	! 			open(unit=12, file=trim(filout12), status="unknown")

	! 			do 
	! 				read(12,*,end=12),tmpId,HapLength,(HapPos(k),k=1,HapLength)
	! 		    	read(12,*),tmpId,HapLength,(SamHap(k),k=1,HapLength)
	! 		    	!print*,tmpId,HapLength,HapPos(1:HapLength)

	! 		    	LikeP=log(.5)
	! 		    	LikeM=log(.5)

	! 			    do j=1,HapLength
	! 			        p0=Pr00(i,HapPos(j))+Pr01(i,HapPos(j))
	! 			        m0=Pr00(i,HapPos(j))+Pr10(i,HapPos(j))
				        
	! 			        if (p0.lt.0.0000001) p0=0.0000001
	! 			        if (m0.lt.0.0000001) m0=0.0000001

	! 			        if (p0.gt.0.9999999) p0=0.9999999
	! 			        if (m0.gt.0.9999999) m0=0.9999999

	! 			      if (SamHap(j)==0) then
	! 			        LikeP=LikeP+log(p0)
	! 			        LikeM=LikeM+log(m0)
	! 			      else if (SamHap(j)==1) then
	! 			        LikeP=LikeP+log(abs(p0-1))
	! 			        LikeM=LikeM+log(abs(m0-1))
	! 			      endif
	! 			    enddo

	! 			    !if (i==1) print*,LikeP,LikeM
	! 			    if (LikeP.gt.LikeM) then ! The haplotype is paternal
	! 			      do j=1,HapLength
	! 			        if (FilledPhase(i,HapPos(j),1)==9) then
	! 			        	cp=cp+1
	! 			        	FilledPhase(i,HapPos(j),1)=SamHap(j)
	! 			        	FilledGenos(i,HapPos(j))=1
	! 			        endif
	! 			        if (FilledPhase(i,HapPos(j),2)==9) then
	! 			        	cm=cm+1
	! 			        	FilledPhase(i,HapPos(j),2)=abs(SamHap(j)-1)
	! 			        	FilledGenos(i,HapPos(j))=1
	! 			        endif

	! 			      	if (FilledPhase(i,HapPos(j),1)/=SamHap(j)) then
	! 			      		wp=wp+1
	! 			      		FilledPhase(i,HapPos(j),1)=9
	! 			      		FilledGenos(i,HapPos(j))=1
	! 			      	endif
	! 			      	if (FilledPhase(i,HapPos(j),2)/=abs(SamHap(j)-1)) then
	! 			      		wm=wm+1
	! 			      		FilledPhase(i,HapPos(j),1)=9
	! 			      		FilledGenos(i,HapPos(j))=1
	! 			      	endif
	! 			      enddo
	! 			    endif
				    
	! 			    if (LikeM.gt.LikeP) then ! The haplotype is maternal
	! 			      do j=1,HapLength
	! 			        if (FilledPhase(i,HapPos(j),2)==9) then
	! 			        	cm=cm+1
	! 			        	FilledPhase(i,HapPos(j),2)=SamHap(j)
	! 			        	FilledGenos(i,HapPos(j))=1
	! 			        endif
	! 			        if (FilledPhase(i,HapPos(j),1)==9) then
	! 			        	cp=cp+1
	! 			        	FilledPhase(i,HapPos(j),1)=abs(SamHap(j)-1)
	! 			        	FilledGenos(i,HapPos(j))=1
	! 			        endif

	! 			      	if (FilledPhase(i,HapPos(j),2)/=SamHap(j)) then
	! 			      		wm=wm+1
	! 			      		FilledPhase(i,HapPos(j),1)=9
	! 			      		FilledGenos(i,HapPos(j))=1
	! 			      	endif
	! 			      	if (FilledPhase(i,HapPos(j),1)/=abs(SamHap(j)-1)) then
	! 			      		wp=wp+1
	! 			      		FilledPhase(i,HapPos(j),1)=9
	! 			      		FilledGenos(i,HapPos(j))=1
	! 			      	endif
	! 			      enddo
	! 			    endif

	! 			   	if (LikeM.eq.LikeP) then
	! 			      do j=1,HapLength
	! 			      	if (FilledGenos(i,HapPos(j))/=1) then
	! 			      		FilledGenos(i,HapPos(j))=1
	! 			      		FilledPhase(i,HapPos(j),:)=9
	! 			      	endif
	! 			      enddo
	! 			    endif
	! 			enddo
	! 			12 close(12)
	! 			if (tmpId.ne.0) print*,tmpId,cp,cm,wp,wm
	! 		endif
	! 	enddo
end subroutine ReadSamFile
