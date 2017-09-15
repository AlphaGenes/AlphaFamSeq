! #ifdef OS_UNIX

! #DEFINE DASH "/"
! #DEFINE MD "mkdir -p"
! #else
! #DEFINE DASH "\"
! #DEFINE MD "md"
! #endif


!################################################################################################

module GlobalPar
	
	use ISO_Fortran_Env
	use pedigreemodule
	implicit none

	! SpecFile Parametes --------------------------------------------------------------------------
	character(len=300) :: PedigreeFile      								! SpecFile - Input File Name - Pedigree
	character(len=300) :: SequenceFile             							! SpecFile - Input File Name - Sequence Data
	character(len=300) :: SequenceFileFormat             					! SpecFile - Input SequenceFile Option - AlphaSim or VcfTools format
	character(len=300) :: SequenceDataType									! SpecFile - Input SequenceFile Option - RC = read counts, GL= genotype likelihood, GT= genotype 

	character(len=300) :: SnpChipFile   									! SpecFile - Input File Name - Snp chip array to add more information to the Reads
	character(len=300) :: MapSnpChipFile   									! SpecFile - Input File Name - Map file for Snp chip array
	
	integer 		   :: nSnp												! SpecFile - Input - Total number of Snps

	character(len=300) :: chr 			  									! SpecFile - Input SequenceFile Option - chromosome ID
	integer 		   :: StartPos,EndPos									! SpecFile - Input SequenceFile Option - first and last position

	real 			:: maxStdForReadsCount,ThresholdMaxReadsCount				! SpecFile - Editing Parametes - Remove Reads that are above this standard deviation
	integer 	 	:: ThresholdReadsCount										! SpecFile - Editing Parametes - Remove single/double/n-tones 
	real			:: ThresholdExcessHetero									! SpecFile - Editing Parametes - Remove variants with an excess of heterozygotes
	
	integer(kind=1)	:: UsePrevGeneProb                      					! SpecFile - Input SingleLocusPeeler - Read old results of GeneProb 1==YES, 0==NO
	real 		 	:: GeneProbThresh  											! SpecFile - Input SingleLocusPeeler - Threshold to call a genotype from the probabilities First Value
	real 			:: GeneProbThreshMin										! SpecFile - Input SingleLocusPeeler - Threshold to call a genotype from the probabilities Last Value
	real 			:: ReduceThr 												! SpecFile - Input SingleLocusPeeler - Reduce Geno Treshold factor
	
	integer 	 :: minWindowSizeHapDefinition      						! SpecFile - Input Build Consensu Haplotype - First value to define Haplotypes length
	integer 	 :: maxWindowSizeHapDefinition      						! SpecFile - Input Build Consensu Haplotype - Last value to define Haplotypes length

	character(len=300) :: GenoFile              							! SpecFile - Control Results File Name - TrueGeno Genotypes to check results 
	character(len=300) :: PhaseFile             							! SpecFile - Control Results File Name - True Phase to check results 


	! Global Parametes --------------------------------------------------------------------------
	type(PedigreeHolder) :: ped
	integer :: nInd 														! Calculated Internally - Number of Individuals in the Pedigree

	integer :: IterationNumber                  							! Control Parameter - Define the number of Iterations
	integer :: StartSnp,EndSnp

	! Sequence Data VcfTools Info
	integer(int64), dimension(:), allocatable 		:: position
	real(real64), allocatable, dimension(:) 		:: quality

	! AlphaMPL output
	real(kind=real64),allocatable,dimension(:,:,:) 	:: ReadCounts !< in the format (pedigreeId, snpId, prob) prob is Pr00,Pr01,Pr10,Pr11
    
    real(kind=real64),allocatable,dimension(:) 		:: Maf


	integer(kind=1),allocatable,dimension(:,:) 		:: TrueGenos			! Control Results - True Genotypes to check results 
	integer(kind=1),allocatable,dimension(:,:,:) 	:: TruePhase			! Control Results - True Phase to check results 

	integer(kind=1),allocatable,dimension(:)		:: MarkersToExclude		! CleanUpTheRawData - 0=use the variant; 1= don't use the variant
	
	integer,allocatable,dimension(:) 				:: GeneProbYesOrNo		! Temporary Array - use gene prob or not
	integer(kind=1),allocatable,dimension(:,:,:) 	:: FounderAssignment   	! Temporary File - Save the IDs of the grandparents

	integer :: Windows
end module GlobalPar

!################################################################################################

program FamilyPhase
    
    use ISO_Fortran_Env
    use specFileModule
    use omp_lib
	use GlobalPar
	use IntelRNGMod
	use AlphaMLPModule
	use CalculateStatisticsForGenoAndPhase
	
	implicit none

	!integer 		:: i,j
	integer(int64) 	:: OldCount,NewCount
	integer(int64)	:: CurrentCountMissingPhase,CurrentCountMissingGenos
	integer(int32) 	:: Seed1,Seed2
	real(kind=8) 	:: InitialGeneProbThresh
	logical			:: fileExists
	real(kind=8)	:: tstart,tend

	! Seed Definition --------------------------------------------------------------------------------------------------
	! Use a seed to sample the Haplotypes length of each window and iteration
	inquire(file="Seed.txt", exist=fileExists)
	if (fileExists) then
		open(99,file="Seed.txt",action="read")
		read(99,*) Seed1
		close(99)
		call IntitialiseIntelRNG(Seed1,"SeedOld.txt",Seed2)
	else
		call IntitialiseIntelRNG(Seedfile="SeedOld.txt",Out=Seed2)
	end if

	!call MakeDirectories
	
	! Print out the window/iteratin/haplotype length in a file
	open(101,file="AlphaFamSeqHaplotypeLengths.txt",status="unknown")
	write(101,'(1a32)') "Window Iter Core StartSnp EndSnp" 

	! User-defined parameters ------------------------------------------------------------------------------------------
	print*,"ReadSpecfile"
	call ReadSpecfile(PedigreeFile, &
                      SequenceFile,SequenceFileFormat,SequenceDataType, &
                      SnpChipFile,MapSnpChipFile, &
                      nSnp,chr, StartPos,EndPos, & 
                      maxStdForReadsCount,ThresholdMaxReadsCount, &
                      ThresholdReadsCount,ThresholdExcessHetero, &
                      UsePrevGeneProb,GeneProbThresh,GeneProbThreshMin,ReduceThr, &
                      minWindowSizeHapDefinition,maxWindowSizeHapDefinition, &
                      GenoFile,PhaseFile)
 
	! Read Pedigree ----------------------------------------------------------------------------------------------------
	print*,"Read Pedigree"
	ped = PedigreeHolder(pedigreeFile,nsnps=nSnp) 
	!call ped%sortPedigreeAndOverwrite()
	nInd=ped%pedigreeSize

	! Read Sequence Data -----------------------------------------------------------------------------------------------
	! TODO split windows to avoid huge memory allocation
	! TODO check type of genotype info - right now this for sequence 
	if (UsePrevGeneProb==0) then
		print*,"Read SequenceData"
		if (trim(SequenceDataType)=="RC") then
			if (trim(SequenceFileFormat)=="AlphaSim") call ped%addSequenceFromFile(SequenceFile,nSnp,startSnp=StartPos,endSnp=EndPos) ! read sequence file AlphaSimFormat
			if (trim(SequenceFileFormat)=="VcfTools") call ped%addSequenceFromVCFFile(seqFile=SequenceFile,nSnpsIn=nSnp,chr=chr,StartPos=StartPos,EndPos=EndPos)
		endif
	endif	
	! Edit The Row Data ------------------------------------------------------------------------------------------------
	! TODO : Check Mendelian Inconsistencies
	! TODO : Check Excess of Reads
	! TODO : Check Single- and Double-tones

	! Run GeneProb -----------------------------------------------------------------------------------------------------
	InitialGeneProbThresh=GeneProbThresh
	GeneProbThresh=InitialGeneProbThresh

	allocate(ReadCounts(4,nSnp,nInd))
	if (UsePrevGeneProb==0) then
		print*,"Run GeneProb"
		allocate(Maf(nSnp))
		tstart = omp_get_wtime()
		call runAlphaMLPAlphaImpute(1,nSnp,ped,ReadCounts,Maf)
		tend = omp_get_wtime()
		write(*,*) "Total wall time for Running SingleLP", tend - tstart
		call SaveGeneProbResults
	else if (UsePrevGeneProb==1) then
		print*,"Read old results of SingleLP"
		tstart = omp_get_wtime()
		call ReadPrevGeneProb
		tend = omp_get_wtime()
		write(*,*) "Total wall time for Importing Probabilities", tend - tstart
	endif

	if (maxWindowSizeHapDefinition.gt.1) then
		allocate(FounderAssignment(ped%pedigreeSize-ped%nDummys,nSnp,2))
	endif
	
	! Iterate on the next steps ----------------------------------------------------------------------------------------
	 print*," "
	 write (*,'(1a39)') " Window Iter   ProbThr    %Phase   %Geno"

	 CurrentCountMissingGenos=0
	 CurrentCountMissingPhase=0
	 NewCount=-1
	 OldCount=0

	 IterationNumber=0
	 do while (NewCount.ne.OldCount)
	 	OldCount=NewCount
	 	IterationNumber=IterationNumber+1
			
	 	if ((IterationNumber>1).and.(GeneProbThresh>GeneProbThreshMin)) then 
	 		GeneProbThresh=GeneProbThresh-ReduceThr
	 		if (GeneProbThresh.lt.GeneProbThreshMin) GeneProbThresh=GeneProbThreshMin
	 	endif
	 	call UseGeneProbToSimpleFillInBasedOnOwnReads
	 	call SimpleCleanUpFillIn ! count Geno is not working

	! 	!if (IterationNumber==1) call ReadSamFile
	! 	call SimpleCleanUpFillIn
	! 	!if (IterationNumber==1) call UseSnpChipInformation 
	 	call SimpleFillInBasedOnParentsReads
	 	call SimpleCleanUpFillIn

	 	call SimpleFillInBasedOnProgenyReads
	 	call SimpleCleanUpFillIn
	 		
	 	if (maxWindowSizeHapDefinition.gt.1) then
!			call CalculateFounderAssignment

		 	call ChunkDefinition
!	 		call BuildConsensus
!	 		call SimpleCleanUpFillIn
	 	endif
	! 	! Count Missing
	 	CurrentCountMissingGenos=ped%CountMissingGenotypesNoDummys()
	 	CurrentCountMissingPhase=ped%CountMissingPhaseNoDummys()

	 	NewCount=CurrentCountMissingPhase + CurrentCountMissingGenos
		write (*,'(2i4,1f10.3,2i15)') Windows,IterationNumber,GeneProbThresh,CurrentCountMissingPhase,CurrentCountMissingGenos

	 enddo

	! ! WriteOut Results -------------------------------------------------------------------------------------------------

	 call WriteResults
	! call DeallocateArrays

	! call UnintitialiseIntelRNG
	! close(101)

	! ! ! Compute statistics
	! ! if ((trim(GenoFile)/="None").or.(trim(PhaseFile)/="None")) print*," Calculate Results"
	! ! if (trim(GenoFile)/="None") 	call GetResultsImputation(LenghtSequenceDataFile,"AlphaFamSeqFinalGenos.txt",GenoFile,"AlphaFamSeqEditingMarkersRemoved.txt",1,"Yes","AlphaFamSeq")
	! ! if (trim(PhaseFile)/="None") 	call GetResultsImputation(LenghtSequenceDataFile,"AlphaFamSeqFinalPhase.txt",PhaseFile,"AlphaFamSeqEditingMarkersRemoved.txt",2,"No","AlphaFamSeq")


end program FamilyPhase


! subroutine MakeDirectories


! 	use GlobalPar
! 	use ifport
! 	implicit none

! 	integer :: CSTAT,ESTAT
! 	character(len=100) :: CMSG
! 	logical :: dirExists

! 	#ifdef OS_UNIX
! 		inquire(directory="Output", exist=dirExists)
! 		if (.not. dirExists) then
! 			CALL EXECUTE_COMMAND_LINE ("mkdir -p Output",EXITSTAT=ESTAT, CMDSTAT=CSTAT, CMDMSG=CMSG)
! 			if (CSTAT.gt.0) then
! 				print*,"ERROR : Command execution failed",trim(CMSG)
! 			else if (CSTAT.lt.0) then
! 				print*,"ERROR : Command execution not supported"
! 			else
! 				print*,"Command completed with status",ESTAT
! 			endif
! 		endif

! 		inquire(directory="Temporary", exist=dirExists)
! 		if (.not. dirExists) then
! 			CALL EXECUTE_COMMAND_LINE ("mkdir -p Output",EXITSTAT=ESTAT, CMDSTAT=CSTAT, CMDMSG=CMSG)
! 			if (CSTAT.gt.0) then
! 				print*,"ERROR : Command execution failed",trim(CMSG)
! 			else if (CSTAT.lt.0) then
! 				print*,"ERROR : Command execution not supported"
! 			else
! 				print*,"Command completed with status",ESTAT
! 			endif
! 		endif

! 		inquire(directory="Stats", exist=dirExists)
! 		if (.not. dirExists) then
! 			CALL EXECUTE_COMMAND_LINE ("mkdir -p Output",EXITSTAT=ESTAT, CMDSTAT=CSTAT, CMDMSG=CMSG)
! 			if (CSTAT.gt.0) then
! 				print*,"ERROR : Command execution failed",trim(CMSG)
! 			else if (CSTAT.lt.0) then
! 				print*,"ERROR : Command execution not supported"
! 			else
! 				print*,"Command completed with status",ESTAT
! 			endif
! 		endif

! 	#endif
	
! end subroutine MakeDirectories

!-------------------------------------------------------------------------------------------------
!> @brief   Get phase complement and make the genotype
!> @detail  Do while thre are no more information to add
!> @date    August 31, 2017
!--------------------------------------------------------------------------------------------------   

subroutine SimpleCleanUpFillIn

	use GlobalPar
	use omp_lib
	implicit none

	integer :: oldMissingGeno,newMissingGeno,oldMissingPhase,newMissingPhase

	oldMissingGeno=ped%CountMissingGenotypesNoDummys()
	oldMissingPhase=ped%CountMissingPhaseNoDummys()
	
	newMissingGeno=0
	newMissingPhase=0

	!do while ((newMissingGeno.ne.oldMissingGeno).or.(newMissingPhase.ne.oldMissingPhase))
	do while ((newMissingPhase.ne.oldMissingPhase))
		oldMissingGeno=newMissingGeno
		oldMissingPhase=newMissingPhase

		call ped%phaseComplement()
		call ped%makeGenotype()
		
		newMissingGeno=ped%CountMissingGenotypesNoDummys()
		newMissingPhase=ped%CountMissingPhaseNoDummys()
	enddo
end subroutine SimpleCleanUpFillIn

!-------------------------------------------------------------------------------------------------
!> @brief   Build the consensus haplotype for grandparent-parent-probands 
!> @detail  Here we use the chunks previously define, and we try to fill the missing value of these 
!>          three (or more, if there are sibs) individuals if two of these individuals have information.
!> @date    September 08, 2017
!
!  REVISION HISTORY:
!  2017.08.30  mbattagin - Initial Version: Build the consensus haplotype for grandparent-parent-proband  
!  2017.09.08  mbattagin - Use all the progenies at the same time to build the consensus
!--------------------------------------------------------------------------------------------------   

! subroutine BuildConsensus 
! 	use GlobalPar

! 	implicit none

! 	integer :: i,k,e,j,m,a,nOffs,o,nFounders
! 	integer(kind=1) :: ConsensusHaplotype
! 	integer,dimension(2) :: countAllele !0 and 1
! 	type(individual), pointer :: grandparent
! 	integer,allocatable :: posOffs(:),founderOffspring(:)

! 	do i=1,ped%pedigreeSize-ped%nDummys ! Parent
! 		nOffs=ped%pedigree(i)%nOffs ! store nr of Offsprings
! 		nFounders=0 ! Store number of offsprings with informative positions to build the consensus
! 		k=ped%pedigree(i)%gender ! store the gender of the parent
				
! 		if ((nOffs.gt.0).and.(.not.ped%pedigree(i)%isDummy).and.(.not.ped%pedigree(i)%Founder)) then ! Build the consensus haplotypes using all the progeny informations
! 			allocate(posOffs(nOffs))
! 			allocate(founderOffspring(nOffs))

! 			do o=1,nOffs
! 				posOffs(o)=ped%pedigree(i)%offsprings(o)%p%id
! 				if (maxval(FounderAssignment(posOffs(o),:,k)).ne.0) nFounders=nFounders+1 !Count how many individuals have informative snps
! 			enddo
! 			if (nFounders.gt.0) then ! Start to build the consensus
! 				founderOffspring=0
! 				do j=1,nSnp
! 					founderOffspring=FounderAssignment(posOffs,j,k)
! 					if (maxval(founderOffspring(:)).gt.0) then ! The snps have a founder, build it's consensus
! 						do e=2,3
! 							grandparent => ped%pedigree(i)%getSireDamObjectbyIndex(e)
						
! 							ConsensusHaplotype=9
! 							countAllele=0
! 							do a=0,1 
! 								if (((grandparent%individualPhase(1)%getPhase(j)+grandparent%individualPhase(2)%getPhase(j)).lt.3).and.(.not. grandparent%isDummy)) then
! 									! Avoid to use markers that are not fully phased for the grandparent
! 									do m=1,2 
! 										if (grandparent%individualPhase(m)%getPhase(j).eq.a) countAllele(a+1)=countAllele(a+1)+1
! 									enddo
! 								endif
	
! 								if (ped%pedigree(i)%individualPhase(e-1)%getPhase(j).eq.a) countAllele(a+1)=countAllele(a+1)+1

! 								do o=1,nOffs
! 									if (founderOffspring(o).eq.e) then
! 										if (ped%pedigree(posOffs(o))%individualPhase(k)%getPhase(j).eq.a) countAllele(a+1)=countAllele(a+1)+1
! 									endif
! 								enddo
! 							enddo
! 							! Fill the phases
! 	 						if ((countAllele(2).gt.1).and.(countAllele(2).gt.countAllele(1))) ConsensusHaplotype=1
! 	 						if ((countAllele(1).gt.1).and.(countAllele(1).gt.countAllele(2))) ConsensusHaplotype=0
	 						
! 	 						if (ConsensusHaplotype.ne.9) then
! 	 							!if ((minval(countAllele).gt.0).and.(nFounders.gt.1)) write(*,'(1a10,7(1x,i0),10x,100i1)'),ped%pedigree(i)%originalID,i,j,e,nOffs,nFounders,countAllele, founderOffspring
! 								call ped%pedigree(i)%individualPhase(e-1)%setPhase(j,ConsensusHaplotype)
! 								do o=1,nOffs
! 	 								if (founderOffspring(o).eq.e) then
! 			 							call ped%pedigree(posOffs(o))%individualPhase(k)%setPhase(j,ConsensusHaplotype)
! 									endif
! 	 							enddo
! 	 						endif	
	
! 						enddo
! 					endif
! 				enddo
! 				call ped%pedigree(i)%makeIndividualPhaseCompliment()
! 				call ped%pedigree(i)%makeIndividualGenotypeFromPhase()
! 				do o=1,nOffs
! 					call ped%pedigree(posOffs(o))%makeIndividualPhaseCompliment()
! 					call ped%pedigree(posOffs(o))%makeIndividualGenotypeFromPhase()
! 				enddo

! 			endif	
! 			deallocate(posOffs)
! 			deallocate(founderOffspring)

! 		endif
! 	enddo

! 	call ped%cleangenotypesbasedonhaplotypes()

! end subroutine BuildConsensus 

!-------------------------------------------------------------------------------------------------
!> @brief   Work Left/Work Rigth to find chunks of haplotypes using the Founder Assignement array
!> @detail  i.e. From the previous step for the id7 we have the following 
!>				seq of founders assignment: 1 0 0 0 1 0 0 0 2 0 0 0 2. In this step we do:
!>		        from Left to Rigth        : 1 1 1 1 1 1 1 1 2 2 2 2 2 
!>  		    from Rigth to Left        : 1 1 1 1 1 2 2 2 2 2 2 2 2 and the consensus of these two vector generates
!>
!>        		FounderAssignment   	  : 1 1 1 1 1 0 0 0 2 2 2 2 2 
!>        where the 0 represent a region where the recombination happens.

!> TODO : there can be a lot of small chunks. Find a way to understand if they are a single chunk or if 
!>        the recombination happen in the middle of some of them.
!>        Be carefull that we can have scenarios with females parents not sequenced.
!> @date    August 31, 2017

!  REVISION HISTORY:
!  2017.08.31  mbattagin - Initial Version: Work Left/Work Rigth to find chunks of haplotypes using the Founder Assignement array  
!  2017.09.12  mbattagin - Use core index to have windows with uniform numbers of snps (i.e., avoid few snps in the tail)
! 
!--------------------------------------------------------------------------------------------------   

subroutine ChunkDefinition

	use GlobalPar
	use IntelRNGMod
	use coreutils
	implicit none

	integer 			:: i,j,e,k,c,nCores
	integer 			:: f1,p1,p2,founder
	integer 			:: fSnp
	integer 			:: lSnp
	integer 			:: ChunkLength
	real,dimension(1)   :: Z
	integer,allocatable,dimension (:,:) :: CoreIndex
	
	integer,allocatable,dimension(:) :: FounderAssignmentF,FounderAssignmentB


	real,dimension(2) :: f

	! Check window size
	if (minWindowSizeHapDefinition.gt.nSnp) then
		print*,"WARNING: min window size > nSnp"
		print*,"         use nSnp as min window size"
		minWindowSizeHapDefinition=nSnp
	endif

	if (maxWindowSizeHapDefinition.gt.nSnp) then
		print*,"WARNING: max window size > nSnp"
		print*,"         use nSnp as max window size"
		maxWindowSizeHapDefinition=nSnp
	endif

	if (minWindowSizeHapDefinition.gt.maxWindowSizeHapDefinition) then
		print*,"WARNING: min window size > max window size"
		print*,"         use max window size as min window size"
		minWindowSizeHapDefinition=maxWindowSizeHapDefinition
	endif

	Z = SampleIntelUniformS(n=1,a=0.0,b=1.0) 
	ChunkLength=floor((dble(minWindowSizeHapDefinition)-1)+(dble(maxWindowSizeHapDefinition)-(dble(minWindowSizeHapDefinition)-1))*Z(1))+1
	nCores=nSnp/ChunkLength
	allocate(CoreIndex(nCores,2))
	
	CoreIndex=calculatecores(nSnp,ChunkLength,.false.)
	nCores=size(CoreIndex)/2

	do i=1,nCores
		write(101,'(6(1x,i0))') Windows,IterationNumber,i,CoreIndex(i,:)
	enddo


	do c=1,nCores
		fSnp=CoreIndex(c,1)
		lSnp=CoreIndex(c,2)

		call CalculateFounderAssignment(fSnp,lSnp)
	enddo

! 	!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE(c,i,k,e,j,f1,ChunkLength,fSnp,lSnp,FounderAssignmentF,FounderAssignmentB)
! 	do c=1,nCores
! 		ChunkLength=CoreIndex(c,2)-CoreIndex(c,1)+1
! 		do i=1,ped%pedigreeSize-ped%nDummys
! 			do k=1,2 ! paternal or maternal gamete
! 				e=k+1 ! position parents in Pedigree
				
				
! 				fSnp=CoreIndex(c,1)
! 				lSnp=CoreIndex(c,2)
				
! 				!if (count(FounderAssignment(i,fSnp:lSnp,k).ne.0).gt.1) then
! 				if (dble(count(FounderAssignment(i,fSnp:lSnp,k).ne.0))/dble(ChunkLength).gt.0.05) then

! 					allocate(FounderAssignmentF(ChunkLength))
! 					allocate(FounderAssignmentB(ChunkLength))

! 					FounderAssignmentF=0
! 					FounderAssignmentB=0

! 					FounderAssignmentF(:)=FounderAssignment(i,fSnp:lSnp,k)
! 					FounderAssignmentB(:)=FounderAssignment(i,fSnp:lSnp,k)

					
! 					f1=0
! 					!!! Forward founder assignment 
! 					do j=1,ChunkLength
! 						if (FounderAssignmentF(j)/=0)then
! 							f1=FounderAssignmentF(j)
! 						endif
! 						if ((FounderAssignmentF(j)==0)) then
! 							FounderAssignmentF(j)=f1
! 						endif 
! 					enddo

! 					f1=0
! 					!!! Backward founder assignment 
! 					do j=ChunkLength,1,-1
! 						if (FounderAssignmentB(j)/=0)then
! 							f1=FounderAssignmentB(j)
! 						endif
! 						if ((FounderAssignmentB(j)==0)) then
! 							FounderAssignmentB(j)=f1
! 						endif 
! 					enddo
					
! 					do j=1,ChunkLength
! 						if (FounderAssignmentF(j)==FounderAssignmentB(j)) FounderAssignment(i,(fSnp+j-1),k)=FounderAssignmentF(j)
! 						if (FounderAssignmentF(j)/=FounderAssignmentB(j)) FounderAssignment(i,(fSnp+j-1),k)=0
! 					enddo	
! 					deallocate(FounderAssignmentF)
! 					deallocate(FounderAssignmentB)
! 				endif
! 			enddo
! 		enddo
! 	enddo
! 	!$OMP END PARALLEL DO

 end subroutine ChunkDefinition

!-------------------------------------------------------------------------------------------------
!> @brief   Assign from which founder the haplotype has been inherited.
!> @detail  This is done using the heterozygous markers of the individual’s parent
!>          i.e. DAD phase = 10, KID Phase=19. The "1" of the kid is coming from the Paternal GrandSire
!>			i.e. DAD phase = 01, KID Phase=19. The "1" of the kid is coming from the Paternal GrandDam
!> @date    August 30, 2017
!--------------------------------------------------------------------------------------------------   

subroutine CalculateFounderAssignment(fSnp,lSnp)
	use GlobalPar
	use omp_lib
	implicit none

	integer,intent(in) :: fSnp,lSnp
	integer :: i,j,e,k,phaseId,geno
	integer,dimension(2) :: phasePar
	type(individual),pointer :: parent

	real :: ParentProb,ProgenyProb
	real,allocatable,dimension(:) :: HapProb
	character(len=80) :: filout4
	
	filout4 ="AlphaFamSeqHapProb.txt"

	open (unit=4,file=trim(filout4),status="unknown")


	allocate(HapProb(2))
	
	FounderAssignment(:,:,:)=0

    	do i=1,ped%pedigreeSize-ped%nDummys	
    		if (.not. ped%pedigree(i)%Founder) then
			    do e=2,3
		    		HapProb(:)=.5

			    	k=e-1
					parent => ped%pedigree(i)%getSireDamObjectByIndex(e)
	    			if (associated(parent)) then
		    			do j=fSnp,lSnp
		    				ParentProb=0
		    				ProgenyProb=0

		    				! The parent is heterozygous
		    				if ((ReadCounts(2,j,parent%id).eq.maxval(ReadCounts(:,j,parent%id))).or.(ReadCounts(3,j,parent%id).eq.maxval(ReadCounts(:,j,parent%id)))) then
		    					ParentProb=1-ReadCounts(2,j,parent%id)
		    					ProgenyProb=ReadCounts(1,j,i)+ReadCounts(2,j,i)
		    					HapProb(1)=HapProb(1)+(1-ProgenyProb)*ParentProb
		    					HapProb(2)=HapProb(2)+(ProgenyProb*ParentProb)
		    				endif
		    			enddo
		    			if (HapProb(1)/HapProb(2).ge.1.5) then
		    				write(4,'(1a20,1i2,2i10,1i2)') ped%pedigree(i)%originalID,k,fSnp,lSnp,2
		    			else if (HapProb(1)/HapProb(2).le.0.5) then
		    				write(4,'(1a20,1i2,2i10,1i2)') ped%pedigree(i)%originalID,k,fSnp,lSnp,3
		    			else
		    				write(4,'(1a20,1i2,2i10,1i2)') ped%pedigree(i)%originalID,k,fSnp,lSnp,0
	    				endif
	    			endif
	    		enddo
	    	endif
    	enddo
    close(4)







	! OldFounder Assignement

 	!   	!$OMP PARALLEL DO ORDERED DEFAULT(SHARED)  PRIVATE(e,i,j,parent,phaseId,geno,phasePar) !collapse(2)	
	! do e=2,3 ! Sire and Dam pos in the ped
	! 	do i=1,ped%pedigreeSize-ped%nDummys
	! 		if (.not. ped%pedigree(i)%Founder) then
	! 			parent => ped%pedigree(i)%getSireDamObjectByIndex(e)
				
	! 			if (associated(parent)) then
	! 	 			do j=1,nSnp
	! 					phaseId = ped%pedigree(i)%individualPhase(e-1)%getPhase(j)
						
	! 					geno = parent%individualGenotype%getGenotype(j)
	! 					phasePar(1) =parent%individualPhase(1)%getPhase(j)
	! 	 				phasePar(2) =parent%individualPhase(2)%getPhase(j)

	! 					if ((geno.eq.1).and.(sum(phasePar).lt.3)) then
	! 						if (phasePar(1).eq.phaseId) FounderAssignment(i,j,(e-1))=2 !GrandSire 
	! 						if (phasePar(2).eq.phaseId) FounderAssignment(i,j,(e-1))=3 !GrandDam
	! 					endif

	! 				enddo
	! 			endif
	! 		endif
	! 	enddo
	! enddo
 !    !$OMP END PARALLEL DO
end subroutine CalculateFounderAssignment

!-------------------------------------------------------------------------------------------------
!> @brief   Phase individuals' alleles on the basis of his parents/progeny genotyepes
!> @detail  If an individual has one of its two gametes phased and the evidence from 
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
	integer(kind=1),dimension(2) :: phase,phaseSire,phaseDam
	type(individual),pointer :: sire,dam

	
	!$OMP PARALLEL DO ORDERED DEFAULT(PRIVATE) SHARED (ped,nSnp,nInd) !collapse(2)
	do i=1,ped%pedigreeSize-ped%nDummys

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
					if ((phaseDam(1).ne.9).and.(phaseDam(1).ne.phase(2))) call dam%individualPhase(2)%setPhase(j,phase(2))
					if ((phaseDam(2).ne.9).and.(phaseDam(2).ne.phase(2))) call dam%individualPhase(1)%setPhase(j,phase(2))

					if ((phaseSire(1).ne.9).and.(phaseSire(1).ne.phase(1))) call sire%individualPhase(2)%setPhase(j,phase(1))
					if ((phaseSire(2).ne.9).and.(phaseSire(2).ne.phase(1))) call sire%individualPhase(1)%setPhase(j,phase(2))
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
	integer(kind=1) :: geno

    !$OMP PARALLEL DO ORDERED DEFAULT(PRIVATE) SHARED (ped,nInd,nSnp) !collapse(2)
	do e=1,2
		do i=1,ped%pedigreeSize-ped%nDummys
			parent => ped%pedigree(i)%getSireDamObjectByIndex(e+1)
			if (associated(parent)) then
				do j=1,nSnp
					geno = parent%individualGenotype%getGenotype(j)
					if (geno.eq.0) call ped%pedigree(i)%individualPhase(e)%setPhase(j,0)
					if (geno.eq.2) call ped%pedigree(i)%individualPhase(e)%setPhase(j,1)
				enddo
			endif
		enddo
	enddo
    !$OMP END PARALLEL DO
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
	integer(kind=1),dimension(2) :: phase
	real(kind=real64) :: p00,p01,p10,p11

	!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (ReadCounts,ped,GeneProbThresh,nSnp,nInd) !collapse(2)
	do i=1,ped%pedigreeSize-ped%nDummys
		do j=1,nSnp
			phase(1) = ped%pedigree(i)%individualPhase(1)%getPhase(j)
			phase(2) = ped%pedigree(i)%individualPhase(2)%getPhase(j)
			
			p00=ReadCounts(1,j,i)
			p01=ReadCounts(2,j,i)
			p10=ReadCounts(3,j,i)
			p11=ReadCounts(4,j,i)
			
			if ((sum(phase).gt.3).and.(ped%pedigree(i)%individualGenotype%isMissing(j))) then
				if (p00.ge.GeneProbThresh) then
					call ped%pedigree(i)%individualGenotype%setGenotype(j,0)
					call ped%pedigree(i)%individualPhase(1)%setPhase(j,0)
					call ped%pedigree(i)%individualPhase(2)%setPhase(j,0)
				endif

				if (p11.ge.GeneProbThresh) then
					call ped%pedigree(i)%individualGenotype%setGenotype(j,2)
					call ped%pedigree(i)%individualPhase(1)%setPhase(j,1)
					call ped%pedigree(i)%individualPhase(2)%setPhase(j,1)
				endif


				if ((p01+p10).ge.GeneProbThresh) then
					call ped%pedigree(i)%individualGenotype%setGenotype(j,1)
					if (p01.ge.GeneProbThresh) then
						call ped%pedigree(i)%individualPhase(1)%setPhase(j,0)
						call ped%pedigree(i)%individualPhase(2)%setPhase(j,1)
					endif
					if (p10.ge.GeneProbThresh) then
						call ped%pedigree(i)%individualPhase(1)%setPhase(j,1)
						call ped%pedigree(i)%individualPhase(2)%setPhase(j,0)
					endif
				endif
			endif

			if (((p00+p01).ge.GeneProbThresh).and.(p10.lt.p01).and.(phase(1).eq.9)) call ped%pedigree(i)%individualPhase(1)%setPhase(j,0)
			if (((p11+p10).ge.GeneProbThresh).and.(p01.lt.p10).and.(phase(1).eq.9)) call ped%pedigree(i)%individualPhase(1)%setPhase(j,1)

			if (((p00+p10).ge.GeneProbThresh).and.(p01.lt.p10).and.(phase(2).eq.9)) call ped%pedigree(i)%individualPhase(2)%setPhase(j,0)
			if (((p11+p01).ge.GeneProbThresh).and.(p10.lt.p01).and.(phase(2).eq.9)) call ped%pedigree(i)%individualPhase(2)%setPhase(j,1)

		enddo
	enddo
	!$OMP END PARALLEL DO
end subroutine UseGeneProbToSimpleFillInBasedOnOwnReads

!################################################################################################


!---------------------------------------------------------------------------
!> @brief   Read the results of GeneProb if we already run it
!> @detail  Avoid to run GeneProb again if we don't change the input data
!>          TODO: print out a warning in the err file:
!>  			   WARNING: Reading in old output of GeneProb. If the input
!>							data (pedigree and sequence reads) changed then 
!>							GeneProb must be run again.
!> @date    August 31, 2017
!---------------------------------------------------------------------------   

subroutine ReadPrevGeneProb
	use GlobalPar
	use constantModule, only : IDLENGTH
	use alphahousemod, only : int2char
	implicit none

	integer									:: i,p,id
	logical 								:: exist
	character(len=IDLENGTH)					:: DumC
	character(len=:), allocatable			:: filout5
	real(real64),dimension(:), allocatable 	:: temp
	integer 								:: unit, count
	integer(4) recl_at_start, recl_at_end
	
	filout5 ="AlphaFamSeqFinalSLP_Chr" // trim(adjustl(chr)) // "_StartSnp" // int2char(StartPos) // "_EndSnp" // int2char(EndPos) // ".bin"
	
	inquire(file=trim(filout5),exist=exist)
	
	count = 0
	if (exist) then
		allocate(temp(nSnp))
		open (newunit=unit,file=trim(filout5),status="old", FORM='binary')
		do i=1,nInd

			if (ped%pedigree(i)%isDummy) cycle
			do p=1,4 !4 probabilities
				DUMC = ""
				read(unit) recl_at_start
				read(unit) DumC(1:recl_at_start)
				read(unit) recl_at_end
				read(unit) recl_at_start
				read(unit) temp
				read(unit) recl_at_end
				id = ped%dictionary%getValue(DumC)
				if (id /= DICT_NULL) then
					ReadCounts(p,:,id) = temp
				endif
			enddo
		enddo
		close (unit)
		deallocate(temp)
	else if (.not. exist) then
		print*,"ERROR : ",trim(filout5)," doesn't exist"
		print*,"        use RunSingleLocusPeeler = 0 to run Single Locus Peeler"
		stop
	endif

end subroutine ReadPrevGeneProb

!---------------------------------------------------------------------------
!> @brief   Save GeneProb results in a binary format
!> @detail  TODO: save reduced file in a text format
!> @date    August 31, 2017
!---------------------------------------------------------------------------   

subroutine SaveGeneProbResults
	use GlobalPar
	use alphahousemod, only : int2char
	implicit none

	integer :: i,p,tmpId, unit, count
	character(len=30) :: nChar
	character(len=80) :: FmtInt2
	character(len=:), allocatable:: filout5,filout6
	real(real64),dimension(:), allocatable 	:: temp
	
	! Write Out Full file of GeneProb
	

	filout5 ="AlphaFamSeqFinalSLP_Chr" // trim(adjustl(chr)) // "_StartSnp" // int2char(StartPos) // "_EndSnp" // int2char(EndPos) // ".bin"
	filout6 ="AlphaFamSeqFinalSLP_Chr" // trim(adjustl(chr)) // "_StartSnp" // int2char(StartPos) // "_EndSnp" // int2char(EndPos) // ".txt"
	
	open (newunit=unit,file=trim(filout5),status="unknown", FORM='UNFORMATTED')

	count =  0
	print *, nInd
	do i=1,nInd

		if (ped%pedigree(i)%isDummy) cycle
		tmpId = ped%inputMap(i)
		do p=1,4 ! 4 probabilities
			temp = ReadCounts(p,:,tmpId)
			count = count +1
			write(unit) ped%pedigree(tmpId)%originalID
			write(unit) temp
		enddo
	enddo
	close (unit)
	
	! Write Out Full file of GeneProb
	write(nChar,*) nSnp
	FmtInt2='(1a20,'//trim(adjustl(nChar))//'f7.4)'

	open (newunit=unit,file=trim(filout6),status="unknown")

	do i=1,nInd

		if (ped%pedigree(i)%isDummy) cycle
		tmpId = ped%inputMap(i)
		do p=1,4 ! 4 probabilities
			write(unit,FmtInt2) ped%pedigree(tmpId)%originalID,ReadCounts(p,:,tmpId)
		enddo
	enddo
	close (unit)


end subroutine SaveGeneProbResults

!---------------------------------------------------------------------------
!> @brief   Write Genotypes, phase and founder assignment
!> @detail  TODO: in the final version don't write out FounderAssignment
!> @date    August 31, 2017
!---------------------------------------------------------------------------   

subroutine WriteResults

	use GlobalPar
	use alphahousemod, only : int2char
	implicit none

	integer :: i,tmpId 
	character(len=30) :: nChar
	character(len=80) :: FmtInt2
	character(len=:), allocatable:: filout1,filout2,filout4

	! WriteOut Full Output

	write(nChar,*) nSnp
	FmtInt2='(a20,1x,'//trim(adjustl(nChar))//'(i0))'
	
	filout1 ="AlphaFamSeqFinalPhase_Chr" // trim(adjustl(chr)) // "_StartSnp" // int2char(StartPos) // "_EndSnp" // int2char(EndPos) // ".txt"
	filout2 ="AlphaFamSeqFinalGenos_Chr" // trim(adjustl(chr)) // "_StartSnp" // int2char(StartPos) // "_EndSnp" // int2char(EndPos) // ".txt"
	
	call ped%writeoutphase(trim(filout1))
	call ped%writeoutgenotypes(trim(filout2))

	close (1)
	close (2)

	if (maxWindowSizeHapDefinition.gt.1) then
		filout4 ="AlphaFamSeqFinalFounderAssignment_Chr" // trim(adjustl(chr)) // "_StartSnp" // int2char(StartPos) // "_EndSnp" // int2char(EndPos) // ".txt"

		open (unit=4,file=trim(filout4),status="unknown")
		do i=1,ped%pedigreeSize-ped%nDummys
			!tmpId = ped%inputMap(i)
			!if (maxval(FounderAssignment(i,:,:))/=0) then 
				!print*,i,tmpId," ",ped%pedigree(i)," ",ped%pedigree(tmpId)%originalID
				write (4,FmtInt2) ped%pedigree(i)%originalID,FounderAssignment(i,:,1)
				write (4,FmtInt2) ped%pedigree(i)%originalID,FounderAssignment(i,:,2)
			!endif
		enddo
		close (4)
	endif
	
end subroutine WriteResults

!################################################################################################

subroutine ReadsLikelihood!(nRef,nAlt,ErrorRate,Pr0,Pr1,Pr2)
	! implicit none
	
	! real,intent(in) :: nRef,nAlt
	! real(kind=8),intent(in) :: ErrorRate
	! real,intent(inout) :: Pr0,Pr1,Pr2
	! real :: lPr0,lPr1,lPr2

 !      if ((nRef+nAlt)==0) then
 !        lPr0=log(1.)
 !        lPr1=log(1.)
 !        lPr2=log(1.)
 !      else
 !        lPr0=(nRef*log(1-ErrorRate))+(nAlt*log(ErrorRate))
 !        if (lPr0.lt.log(.000000001)) lPr0=-9999
 !        lPr1=(nRef*log(0.5))+(nAlt*log(0.5))
 !        if (lPr1.lt.log(.000000001)) lPr1=-9999
 !        lPr2=(nAlt*log(1-ErrorRate))+(nRef*log(ErrorRate))
 !        if (lPr2.lt.log(.000000001)) lPr2=-9999
 !      endif

 !      Pr0=exp(lPr0)/(exp(lPr0)+exp(lPr1)+exp(lPr2))
 !      Pr1=exp(lPr1)/(exp(lPr0)+exp(lPr1)+exp(lPr2))
 !      Pr2=exp(lPr2)/(exp(lPr0)+exp(lPr1)+exp(lPr2))
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
subroutine ExcessHeterozygotes!(ReadsCount,n,ObsGenos,EstGenosInit,pHetExcess)

	! implicit none

	! integer,intent(in) 						:: n ! Nr of individuals
	! integer(kind=2),intent(in)				:: ReadsCount(n,2) ! Reads for 1 variant and all individuals
	! real(kind=8),intent(inout) 				:: pHetExcess
	! integer,intent(inout)					:: ObsGenos(3),EstGenosInit(3) ! observed  & estimated genotypes
	
	! integer									:: EstGenos(3) ! estimated genotypes
	
	! integer									:: i,nGeno,RareHomo,CommomHomo,RareCopies,mid,posEst,posObs,FirstInd
	! real(kind=8),allocatable,dimension(:) 	:: probs,plow,phigh
	! real(kind=8)							:: SumProbs

	! ObsGenos=0
	! EstGenos=0
	! nGeno=0

	! do i=1,n
	! 	if (sum(ReadsCount(i,:)).gt.0) then
	! 		if ((ReadsCount(i,1).gt.0).and.(ReadsCount(i,2).eq.0)) ObsGenos(1)=ObsGenos(1)+1
	! 		if ((ReadsCount(i,1).gt.0).and.(ReadsCount(i,2).gt.0)) ObsGenos(2)=ObsGenos(2)+1
	! 		if ((ReadsCount(i,1).eq.0).and.(ReadsCount(i,2).gt.0)) ObsGenos(3)=ObsGenos(3)+1
	! 		nGeno=nGeno+1
	! 	endif
	! enddo

	! RareHomo=3
	! CommomHomo=1
	! if (ObsGenos(1).lt.ObsGenos(2)) then
	! 	RareHomo=1 ! if the number of HomoReferenceAllele < HomoAlternativeAllele switch the common and rare allele
	! 	CommomHomo=3
	! endif

	! RareCopies=nint(2.0*dble(ObsGenos(RareHomo))+ObsGenos(2))
	
	! if (RareCopies.lt.2) then
	! 	pHetExcess=1
	! else if (RareCopies.ge.2) then

	! 	mid=nint(dble(RareCopies)*dble(2*nGeno-RareCopies)/dble(2*nGeno))
		
	! 	if (mid.lt.0) mid=0

	! 	if ((mod(RareCopies,2).eq.0).and.(mod(mid,2).ne.0)) mid=mid+1
	! 	if ((mod(RareCopies,2).ne.0).and.(mod(mid,2).eq.0)) mid=mid+1

	! 	EstGenos(2)=mid
	! 	EstGenos(RareHomo)=floor(dble(RareCopies-mid)/2.0)
	! 	if (EstGenos(RareHomo).lt.0) EstGenos(RareHomo)=0
	! 	EstGenos(CommomHomo)=nGeno-EstGenos(RareHomo)-EstGenos(2)
	! 	if (EstGenos(CommomHomo).lt.0) EstGenos(CommomHomo)=0

	! 	EstGenosInit(:)=EstGenos(:)
		
	! 	! # we observed 21 copies of the minor allele (RareCopies) --> the observed nr of hetero  (ObsGenos(2)) will vary seq(1,21,2)

	! 	! The number of possible heterozygotes is odd if RareCopies is odd, and is even if RareCopies is even
	! 	if (mod(RareCopies,2).eq.0) then
	! 		allocate(probs(0:(RareCopies/2)))
	! 		allocate(plow(0:(RareCopies/2)))
			
	! 		posEst=(mid)/2
	! 		posObs=ObsGenos(2)/2
	! 		FirstInd=0
	! 	else if(mod(RareCopies,2).eq.1) then
	! 		allocate(probs(1:((RareCopies+1)/2)))
	! 		allocate(plow(1:((RareCopies+1)/2)))
			
	! 		posEst=(mid+1)/2
	! 		posObs=(ObsGenos(2)+1)/2
	! 		FirstInd=1
	! 	endif

	! 	probs(posEst)=1.0
	! 	SumProbs=probs(posEst)

	! 	! Start to calculate the probabilities using the equations 2 of Am.J.Hum.Genet.76:887-883,2005
	! 	!
	! 	!!!  P(NAB=nAB-2|N,nA) = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
	! 	!!!  P(NAB=nAB+2|N,nA) = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc	/((curr_hets + 2.0) * (curr_hets + 1.0))
	! 	do i=posEst,(FirstInd+1),-1
	! 		probs(i-1)=probs(i)*dble(EstGenos(2))*dble(EstGenos(2)-1.0)/(4.0*dble(EstGenos(RareHomo)+1.0)*dble(EstGenos(CommomHomo)+1.0))
	! 		EstGenos(2)=EstGenos(2)-2
	! 		EstGenos(1)=EstGenos(1)+1
	! 		EstGenos(3)=EstGenos(3)+1
	! 		SumProbs=SumProbs+probs(i-1)
	! 	enddo

	! 	EstGenos(2)=mid
	! 	EstGenos(RareHomo)=floor(dble(RareCopies-mid)/2.0)
	! 	if (EstGenos(RareHomo).lt.0) EstGenos(RareHomo)=0
	! 	EstGenos(CommomHomo)=nGeno-EstGenos(RareHomo)-EstGenos(2)
	! 	if (EstGenos(CommomHomo).lt.0) EstGenos(CommomHomo)=0

	! 	if (EstGenos(2).lt.RareCopies) then
	! 		do i=posEst,(size(probs)-(2-FirstInd))
	! 			probs(i+1)=probs(i)*4.0*dble(EstGenos(CommomHomo))*dble(EstGenos(RareHomo))/(dble(EstGenos(2)+2.0)*dble(EstGenos(2)+1.0))
	! 			EstGenos(2)=EstGenos(2)+2
	! 			EstGenos(1)=EstGenos(1)-1
	! 			EstGenos(3)=EstGenos(3)-1
	! 			SumProbs=SumProbs+probs(i+1)
	! 		enddo
	! 	endif
		
	! !		if ((ObsGenos(1)==611).and.(ObsGenos(2)==0).and.(ObsGenos(3)==1)) print*,probs(:)


	! 	probs(:)=probs(:)/SumProbs
	! 	plow=1
	! 	plow(FirstInd)=probs(FirstInd)

	! !		if ((ObsGenos(1)==611).and.(ObsGenos(2)==0).and.(ObsGenos(3)==1)) print*,probs(:)

	! 	do i=(FirstInd+1),(size(probs)-1)
	! 		plow(i)=plow(i-1)+probs(i)
	! 	enddo



	! !		if ((ObsGenos(1)==611).and.(ObsGenos(2)==0).and.(ObsGenos(3)==1)) print*,plow(:),posObs
	! 	if (FirstInd.eq.0) then
	! 		if (posObs.gt.0) pHetExcess=1-plow(posObs-1)
	! 		if (posObs.eq.0) pHetExcess=1
	! 	endif

	! 	if (FirstInd.eq.1) then
	! 		if (posObs.gt.1) pHetExcess=1-plow(posObs-1)
	! 		if (posObs.eq.1) pHetExcess=1
	! 	endif

	! 	deallocate(probs)
	! 	deallocate(plow)
	! endif
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
