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

	integer :: SuperC                           							! SpecFile - Parameter to use/not use Ferdosi module

	character(len=300) :: PedigreeFile      								! SpecFile - Input File Name - Pedigree
	character(len=300) :: ReadsFile             							! SpecFile - Input File Name - Reads Count for the Reference and the Alternative Allele
	character(len=300) :: ReadsType             							! SpecFile - Input File Name - Reads Count for the Reference and the Alternative Allele
	
	character(len=300) :: MapFile             								! SpecFile - Input File Name - Map File - position of the Variants
	character(len=300) :: SnpChipsInformation   							! SpecFile - Input File Name - Snp array to add more information to the Reads
	
	character(len=300) :: GenoFile              							! SpecFile - Control Results File Name - TrueGeno Genotypes to check results 
	character(len=300) :: PhaseFile             							! SpecFile - Control Results File Name - True Phase to check results 

	integer :: IterationNumber                  							! Control Parameter - Define the number of Iterations
	integer(kind=8) :: CurrentCountFilledPhase          							! Control Parameter - used to finish the program
	integer(kind=8) :: CurrentCountFilledGenos          							! Control Parameter - used to finish the program
	integer :: SolutionChanged                  							! Control Parameter - used to finish the program 
	integer :: StartSnp,EndSnp
	
	integer(int64),allocatable,dimension(:,:) 		:: Ped         			! Input File - Pedigree
	integer(int64),allocatable,dimension(:,:) 		:: RecPed				! Temporary File - Pedigree Recoded
	integer(int64),allocatable,dimension(:) 		:: Id           		! Read Data - used to read unsorted data

	integer(kind=2),allocatable,dimension(:,:,:) 	:: SequenceData			! Input File - Snp array to add more information to the Reads
	integer(kind=2),allocatable,dimension(:,:,:) 	:: RawReads				! Input File - Snp array to add more information to the Reads
	
	character(len=100), allocatable, dimension(:) 	:: Ids
	integer(int64), dimension(:), allocatable 		:: position
	real(real64), allocatable, dimension(:) 		:: quality

	integer(kind=1),allocatable,dimension(:,:) 		:: TrueGenos			! Control Results - True Genotypes to check results 
	integer(kind=1),allocatable,dimension(:,:,:) 	:: TruePhase			! Control Results - True Phase to check results 

	integer(kind=1),allocatable,dimension(:)		:: MarkersToExclude		! CleanUpTheRawData - 0=use the variant; 1= don't use the variant
	
	character(len=1),allocatable,dimension(:,:) 	:: CheckGenos   		! Control Results - Use character to check True vs Imputed Genotypes
	character(len=1),allocatable,dimension(:,:,:) 	:: CheckPhase 			! Control Results - Use character to check True vs Imputed PhaseFile
	
	integer,allocatable,dimension(:) 				:: GeneProbYesOrNo		! Temporary Array - use gene prob or not
	integer,allocatable,dimension(:,:,:) 			:: FounderAssignment   	! Temporary File - Save the IDs of the grandparents

	real(kind=4),allocatable,dimension(:,:) 		:: Pr00   				! Output GeneProb - Probabilities for Ind(i) and Spn(j) to be Homozygote for Reference Allele 
    real(kind=4),allocatable,dimension(:,:) 		:: Pr01	  				! Output GeneProb - Probabilities for Ind(i) and Spn(j) to be Heterozygote (0 from dad, 1 from mum)
    real(kind=4),allocatable,dimension(:,:) 		:: Pr10					! Output GeneProb - Probabilities for Ind(i) and Spn(j) to be Heterozygote (1 from dad, 0 from mum)
    real(kind=4),allocatable,dimension(:,:) 		:: Pr11					! Output GeneProb - Probabilities for Ind(i) and Spn(j) to be Homozygote for Alternative Allele 
    
	integer(kind=1),allocatable,dimension(:,:) 		:: FilledGenos 			! Output - Imputed Genotypes
	integer(kind=1),allocatable,dimension(:,:,:) 	:: FilledPhase  		! Output - Imputed Phase

	!												:: StatBySnp            ! Save coverage, nr 
	!												:: StatByInd            ! 

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
	

	! Read SpecFile and Pedigree. Those files are in common for all the windows
	!  (if there are multiple Windows)
	call ReadSpecfile
	call ReadPedigree
	
	
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

		call ReadData
		call InitialiseArrays
		call CleanUpTheRawData

		if (UsePrevGeneProb==0) then
			call RunGeneProb
			call SaveGeneProbResults
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
			if (IterationNumber==1) call ReadSamFile
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
!			call CountFounder
			
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

			write (*,'(2i4,3f10.3)') Windows,IterationNumber,GeneProbThresh,(dble(CurrentCountFilledPhase)/(dble(nInd*nSnp*2))*100),(dble(CurrentCountFilledGenos)/(dble(nInd*nSnp))*100)
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
! STOLEN FROM ALPHASIM, WRITTEN BY DAVID WILSON 
! Function returns a character (Of 512 Bytes for compatibility) that is a completely lowercase copy of input str
function TLC(str)
    
    character(*), intent(in) :: str
    character(len=512) :: TLC
    integer :: i
    TLC = trim(str)
    do i = 1, len(TLC)
        select case(TLC(i:i))
            case("A":"Z")
                TLC(i:i) = achar(iachar(TLC(i:i))+32)
        end select
    enddo
    return
end function TLC

!###########################################################################################################################################################

function probscore(x1)
	implicit none
	integer(kind=2)         :: probscore
	real(kind=8),intent(in) :: x1
	real(kind=8) :: x2

	x2=x1
	if (x2.le.0.0001) x2=0.0001
	if (x2.ge.0.9999) x2=0.9999
	probscore=nint(-10*log10(x2)*100)

	return
end function probscore

!###########################################################################################################################################################

function score2prob(x1)
	implicit none
	real(kind=8)         		:: score2prob
	integer(kind=2),intent(in)  :: x1

	score2prob=10**(-x1/10*100) 
	return
end function score2prob

!###########################################################################################################################################################

! Reads in and initialises specfile parameters 
subroutine ReadSpecfile

    use GlobalPar

    implicit none

    integer :: FileLength, stat, i
    !character(len=256) :: Var
    character(len=30) :: SpecParam
   	character (len=512) :: TLC

    open(unit=1, file="AlphaFamSeqSpec.txt", status="old")

    FileLength = 0

    do
        read(1, *, iostat=stat) SpecParam
        if (stat/=0) exit
        FileLength = FileLength + 1
    enddo

    rewind(1)

    do i=1, FileLength
        read(1,'(a30,A)', advance='NO', iostat=stat) SpecParam 
        
        select case(trim(TLC(SpecParam)))

        	case('numberofindividuals')
                read(1, *, iostat=stat) nIndSeq
                if (stat /= 0) then
                    print *, "NumberOfIndividuals not set properly in spec file"
                    stop 2
                endif 
               
			case('numberofsnps')
                read(1, *, iostat=stat) LenghtSequenceDataFile,nSnp,fistWindow
                if (stat /= 0) then
                    print *, "NumberOfSnps not set properly in spec file"
                    print *, LenghtSequenceDataFile,nSnp
                    stop 2
                endif 

			case('internaledit')
                read(1, *, iostat=stat) InternalEdit
                if (stat /= 0) then
                    print *, "InternalEdit not set properly in spec file"
                    print *, InternalEdit
                    stop 2
                endif 

			case('editingparameter')
                read(1, *, iostat=stat) EditingParameter
                if (stat /= 0) then
                    print *, "EditingParameter not set properly in spec file"
                    print *, EditingParameter
                    stop 2
                endif 

			case('removeoutliersreadscount')
                read(1, *, iostat=stat) maxStdForReadsCount
                if (stat /= 0) then
                    print *, "RemoveOutliersReadsCount not set properly in spec file"
                    print *, maxStdForReadsCount,ThresholdMaxReadsCount
                    stop 2
                endif 

			case('removemarkerslownrreads')
                read(1, *, iostat=stat) ThresholdReadsCount
                if (stat /= 0) then
                    print *, "RemoveMarkersLowNrReads not set properly in spec file"
                    print *, ThresholdReadsCount
                    stop 2
                endif 

			case('removeexcessheteropvalue')
                read(1, *, iostat=stat) ThresholdExcessHetero
                if (stat /= 0) then
                    print *, "RemoveExcessHeteroPvalue not set properly in spec file"
                    print *, ThresholdExcessHetero
                    stop 2
                endif 



            case('genotypeprobability')
                read(1, *, iostat=stat) GeneProbThresh,GeneProbThreshMin,ReduceThr,UsePrevGeneProb !nIter 
                if (stat /= 0) then
                    print *, "GenotypeProbability not set properly in spec file"
                    stop 8
                endif   

            case('errorrate')
                read(1, *, iostat=stat) ErrorRate
                if (stat /= 0) then
                    print *, "ErrorRate not set properly in spec file"
                    stop 8
                endif   

            case('rangechunklenght')
                read(1, *, iostat=stat) ChunkLengthA,ChunkLengthB
                if (stat /= 0) then
                    print *, "RangeChunkLenght not set properly in spec file"
                    stop 8
                endif   

			case('superconsensus')
                read(1, *, iostat=stat) SuperC
                if (stat /= 0) then
                    print *, "SuperConsensus not set properly in spec file"
                    stop 8
                endif   
                
			case('pedigreefile')
                read(1, *, iostat=stat) PedigreeFile
                if (stat /= 0) then
                    print *, "PedigreeFile not set properly in spec file"
                    stop 8
                endif   

			case('readsfile')
                read(1, *, iostat=stat) ReadsFile,ReadsType
                if (stat /= 0) then
                    print *, "ReadsFile not set properly in spec file"
                    stop 10
                endif   

			case('mapfile')
                read(1, *, iostat=stat) MapFile
                if (stat /= 0) then
                    print *, "MapFile not set properly in spec file"
                    stop 10
                endif   
                
			case('genofile')
                read(1, *, iostat=stat) GenoFile
                if (stat /= 0) then
                    print *, "GenoFile not set properly in spec file"
                    stop 10 
                endif   
            
            case('snpchipsinformation')
            	read(1, *, iostat=stat) SnpChipsInformation
            	if (stat /= 0) then
                    print *, "SnpChipsInformation not set properly in spec file"
                    stop 10
                endif
                
			case('phasefile')
                read(1, *, iostat=stat) PhaseFile
                if (stat /= 0) then
                    print *, "PhaseFile not set properly in spec file"
                    stop 10
                endif
                
			case default
               print *, "Error in specfile, please check", SpecParam
                stop 16
    	end select
    enddo

    close(1)
end subroutine ReadSpecfile

!###########################################################################################

! Not used at the moment. It allows to use the Snp Chip info, if available
! TODO: istead filling the genotyepes and the phase, merge these info with the one coming 
!       from the Reads.

subroutine UseSnpChipInformation

	use GlobalPar

	implicit none

	integer:: FileLength,stat,DumI,i,j
	integer(int64) :: PosGeno
	integer,allocatable,dimension(:) :: TempImput
	integer,allocatable,dimension(:,:) :: IdsSnpChip
	
	if (trim(SnpChipsInformation)/="None") then

		allocate(TempImput(nSnp))
		
		open (unit=11,file=trim(SnpChipsInformation),status="old")

		FileLength = 0
	    do
	        read(11, *, iostat=stat) DumI
	        if (stat/=0) exit
	        FileLength = FileLength + 1
	    enddo
	    
	    rewind(11)

	    allocate(IdsSnpChip(nInd,nSnp))
	    
	    IdsSnpChip=9

	    do i=1, FileLength
			PosGeno=0
	    	read(11,*) Id(i), TempImput(:) 
	    	call GetID(Id(i), PosGeno)
	    	IdsSnpChip(PosGeno,:) = TempImput(:)
	    enddo

		close(11)
	    
	    
	    do i=1,nInd
			do j=1,nSnp
				if (IdsSnpChip(i,j)==0) then 
					FilledPhase(i,j,:)=0
					FilledGenos(i,j)=0
				endif
				if (IdsSnpChip(i,j)==2) then
					FilledPhase(i,j,:)=1
					FilledGenos(i,j)=2
				endif
				if (IdsSnpChip(i,j)==1) then
					FilledGenos(i,j)=1
				endif
			enddo
		enddo

		call CurrentCountFilled

		deallocate(TempImput)
		deallocate(IdsSnpChip)
	
	endif
end subroutine UseSnpChipInformation

!###########################################################################################

! This step perform 3 editing (at the moment):
! 1) Markers with no reads (it happens with simulated data or with a subset of real data) 
! 2) Remove variants with less than nUserDefined reads for the referece or alternative allele (i.e., single- and double-tones)
! 3) Remove ReadsCount/invidual that are above a user defined threshold (i.e., excess of reads compare to the average)
! 4) Remove variants with an excess of heterozygotes using the Exact Tests of HWE described in Wigginton et al., Am.J.Hum.Genet. 76:887-883,2005

! Moreover, the coverage/individual and coverage/marker are calculated in this step

subroutine CleanUpTheRawData

	use GlobalPar
	use omp_lib
	implicit none

	integer :: i,j,nReadsRemoved,nTmpInd,e,pos
	real 	:: cov(nInd)
	real 	:: std,covSnp
	character(len=50) :: filout1,filout2,filout3,filout4,filout5
	character(len=30) :: nChar
	character(len=80) :: FmtInt

	integer(kind=2),allocatable,dimension(:) :: ReadsRemoved
	integer(kind=2),allocatable,dimension(:,:) :: tmpReads

	real(kind=8)					 				:: pHetExcess
	integer											:: ObsGenos(3), EstGenos(3) ! observed genotypes


	write (filout1,'("AlphaFamSeqEditingIndividualCoverage",i0,".txt")') Windows
	open (unit=1,file=trim(filout1),status="unknown")

	write (filout2,'("AlphaFamSeqEditingMarkerCoverage",i0,".txt")') Windows
	open (unit=2,file=trim(filout2),status="unknown")

	write (filout3,'("AlphaFamSeqEditingMarkersRemoved",i0,".txt")') Windows
	open (unit=3,file=trim(filout3),status="unknown")

	write (filout4,'("AlphaFamSeqEditingIndividualReadsRemoved",i0,".txt")') Windows
	open (unit=4,file=trim(filout4),status="unknown")

	write (filout5,'("AlphaFamSeqEditingPvalueExcessOfHeterozygotes",i0,".txt")') Windows
	open (unit=5,file=trim(filout5),status="unknown")

	!write (filout3,'("AlphaFamSeqReads",i0,".txt")') Windows
	!open (unit=3,file=trim(filout3),status="unknown")

	allocate(ReadsRemoved(nSnp))
	ReadsRemoved(:)=0
	MarkersToExclude=0 ! 0=Use The marker ; 1=Don't use the marker

	do j=1,nSnp
		covSnp=0
		!e=0
		covSnp=sum(RawReads(:,j,:))/dble(nTmpInd)
		!if (cov.lt.EditingParameter) then
		!	RawReads(:,j,:)=0
		!	e=1
		!endif
		write (2,'(1i0,1x,1f7.3,1x,1i0)') j,covSnp,e
	enddo
	close(2)


	nTmpInd=0
	do i=1,nInd
		cov(i)=0
		do j=1,nSnp
			cov(i)=cov(i)+sum(RawReads(i,j,:))
		enddo
		if (cov(i).gt.0) nTmpInd=nTmpInd+1
		cov(i)=cov(i)/dble(nSnp)

		std=0
		do j=1,nSnp
			std=std+((sum(RawReads(i,j,:))-cov(i))**2)
		enddo
		std=sqrt(std/(dble(nSnp)-1))

		nReadsRemoved=0
		if (cov(i).gt.0.0) then
			do j=1,nSnp
				if (((sum(RawReads(i,j,:))-cov(i))/std).gt.maxStdForReadsCount) then
					write(4,'(2i20,1f10.6,2i4)'),Ped(i,1),j,cov(i),RawReads(i,j,:)
					RawReads(i,j,:)=0
					nReadsRemoved=nReadsRemoved+1
					ReadsRemoved(j)=ReadsRemoved(j)+1
				endif
			enddo
		endif

		write (1,'(1i0,1f7.3,1f10.6)') Ped(i,1),cov(i),dble(nReadsRemoved)/dble(nSnp)*100

		!write (3,FmtInt) Ped(i,1), RawReads(i,:,1)
		!write (3,FmtInt) Ped(i,1), RawReads(i,:,2)
	enddo

	close(1)
	close(4)


	allocate(tmpReads(nTmpInd,2))
	
	!!$OMP PARALLEL DO ORDERED DEFAULT(PRIVATE) SHARED (RawReads,nSnp)
	do j=1,nSnp
		
		! 1) Markers with zero reads	
		if ((maxval(RawReads(:,j,:))==0).or.(ReadsRemoved(j).gt.(dble(nInd)*ThresholdMaxReadsCount))) then
			write (3,'(1i10,1i20,1i4)') j,position(j),0
			MarkersToExclude(j)=1
			RawReads(:,j,:)=0
		endif

		! 2) Remove single/double-tones and excess heterozygotes
		
		tmpReads=9
		ObsGenos=0
		EstGenos=0
		pHetExcess=0.0

		pos=1
		do i=1,nInd
			if (cov(i).gt.0.0) then
				tmpReads(pos,:)=RawReads(i,j,:)
				pos=pos+1
			endif
		enddo

		call ExcessHeterozygotes(tmpReads(:,:),nTmpInd,ObsGenos,EstGenos,pHetExcess)

		if (((ObsGenos(2)+dble(2*ObsGenos(1))).le.ThresholdReadsCount).or.((ObsGenos(2)+dble(2*ObsGenos(3))).le.ThresholdReadsCount)) then
			if (MarkersToExclude(j).ne.1) then
				write (3,'(1i10,1i20,1i4)') j,position(j),2
				MarkersToExclude(j)=1
			endif
		endif

		if (MarkersToExclude(j).ne.1) then
			if (pHetExcess.le.ThresholdExcessHetero) then
				write (3,'(1i10,1i20,1i4)') j,position(j),1
				write(5,'(1i20,3i4,5x,3i4,1f10.5)') position(j),ObsGenos(:),EstGenos(:),pHetExcess

				MarkersToExclude(j)=1
		
			endif
		endif

	enddo
	!!$OMP END PARALLEL DO

	write(*,'(1a10,1x,1i0,1a4,1x,1i0,1a10)'),"Exclude",count(MarkersToExclude(:)==1),"on",nSnp,"Markers"

	deallocate(tmpReads)
	
	close (3)
	close (5)
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

! At the end of each main step of the program, this subroutine fill the gaps
!  i.e. if a genotype is Homozygote but the phase is missing, then it fill the phase
subroutine SimpleCleanUpFillIn

	use GlobalPar
	use omp_lib
	implicit none

	integer i,j,Change !,e,k,g1,g2

	!Change=1
	!do while (Change==1)
		Change=0
		!$OMP PARALLEL DO ORDERED DEFAULT(FIRSTPRIVATE) SHARED (FilledGenos,FilledPhase,nInd,Change, nSnp) !collapse(2)
		do i=1,nInd
			do j=1,nSnp
				if (FilledGenos(i,j)==9) then
					if ((FilledPhase(i,j,1)/=9).and.(FilledPhase(i,j,2)/=9)) then
						FilledGenos(i,j)=sum(FilledPhase(i,j,:))
						Change=1
					endif
				endif
				if ((FilledPhase(i,j,1)/=9).and.(FilledPhase(i,j,2)==9).and.(FilledGenos(i,j)==1)) then
					if (FilledPhase(i,j,1)==0) FilledPhase(i,j,2)=1
					if (FilledPhase(i,j,1)==1) FilledPhase(i,j,2)=0
					Change=1
				endif
				if ((FilledPhase(i,j,2)/=9).and.(FilledPhase(i,j,1)==9).and.(FilledGenos(i,j)==1)) then
					if (FilledPhase(i,j,2)==0) FilledPhase(i,j,1)=1
					if (FilledPhase(i,j,2)==1) FilledPhase(i,j,1)=0
					Change=1
				endif

				if (FilledGenos(i,j)/=9) then
					if ((FilledPhase(i,j,1)/=9).and.(FilledPhase(i,j,2)/=9)) then
						if (FilledGenos(i,j)/=sum(FilledPhase(i,j,:))) then
							FilledGenos(i,j)=9
							FilledPhase(i,j,:)=9
							Change=1
							!print*,"CheckError",i,j
						endif
					endif
				endif

			enddo
			
		enddo
		!$OMP END PARALLEL DO
	!enddo
end subroutine SimpleCleanUpFillIn

!###########################################################################################

 !subroutine FerdosiSerap

	! 	use GlobalPar
	! 	use Ferdosi
	! 	implicit none

	! 	integer :: nHalfSib,i,k,p,hs,ParentGender,j,ChunckStart,ChunckStop,count1,count2,j2,nHalfSibOwnReads
	! 	integer, allocatable, dimension(:) :: HalfSibID,ParentGeno,HasOwnReads
	! 	integer, allocatable, dimension(:,:) :: ParentPhase, HalfSibGeno
	! 	integer, allocatable, dimension(:,:,:) :: HalfSibPhase
		

	! 	nHalfSib=0
	! 	do i=nInd,1,-1 ! count how many time is a parent
	! 		do k=2,3 ! sire or dam
	! 			ParentGender=k-1
	! 			nHalfSib=count(RecPed(:,k)==i)
		
	! 			allocate(HalfSibID(nHalfSib))
	! 			allocate(HasOwnReads(nHalfSib))

	! 			HasOwnReads=0
	! 			nHalfSibOwnReads=0

	! 			p=1
	! 			do hs=1,nInd
	! 				if (RecPed(hs,k)==i) then
	! 					HalfSibID(p)=hs
	! 					if ((maxval(RawReads(hs,:,:))>0).and.(count(FilledGenos(hs,:)==9)<20000)) then
	! 						HasOwnReads(p)=1
	! 						nHalfSibOwnReads=nHalfSibOwnReads+1
	! 					endif
	! 					p=p+1
	! 				endif
	! 			enddo

	! 			if (nHalfSibOwnReads > 19 ) then

	! 				allocate(HalfSibGeno(nHalfSibOwnReads, nSnp))
	! 				allocate(HalfSibPhase(nHalfSibOwnReads, nSnp, 2))

	! 				allocate(ParentPhase(nSnp,2))
	! 				allocate(ParentGeno(nSnp))

	! 				p=1
	! 				do hs=1,nHalfSib
	! 					if (HasOwnReads(hs)==1) then
	! 						HalfSibGeno(p,:)=FilledGenos(HalfSibID(hs),:)
	! 						p=p+1
	! 					endif
	! 				enddo

	! 				ParentGeno = 9
	! 				ParentPhase = 9
	! 				HalfSibPhase = 9

	! 				call FerdosiRunner(nHalfSibOwnReads, nSnp, ParentGender, ParentPhase, HalfSibGeno, HalfSibPhase)

	! 				do j=1,nSnp ! Update Parents genotypes
	! 					if (sum(ParentPhase(j,:))<3) then 

	! 						ParentGeno=sum(ParentPhase(j,:))

	! 						if (FilledGenos(i,j)==9) then
	! 							FilledGenos(i,j)=ParentGeno(j)
	! 						endif

	! 						if (FilledGenos(i,j)==0) then
	! 							if ((FilledPhase(i,j,1)==9).and.(FilledPhase(i,j,2)==9)) then
	! 								FilledPhase(i,j,:)=0
	! 							endif
	! 						endif
							
	! 						if (FilledGenos(i,j)==2) then
	! 							if ((FilledPhase(i,j,1)==9).and.(FilledPhase(i,j,2)==9)) then
	! 								FilledPhase(i,j,:)=1
	! 							endif
	! 						endif

	! 	!						if (FilledGenos(i,j)==1) then
	! 	!							if ((FilledPhase(i,j,1)/=9).or.(FilledPhase(i,j,2)/=9)) then
	! 	!								if (FilledPhase(i,j,1)==0) FilledPhase(i,j,2)=1
	! 	!								if (FilledPhase(i,j,1)==1) FilledPhase(i,j,2)=0
	! 	!								if (FilledPhase(i,j,2)==0) FilledPhase(i,j,1)=1
	! 	!								if (FilledPhase(i,j,2)==1) FilledPhase(i,j,1)=0
	! 	!							endif
	! 	!						endif
						
	! 					! check
	! 					if (i==100) then
	! 						if ((ParentGeno(j)/=TrueGenos(i,j)).and.(ParentGeno(j)==FilledGenos(i,j))) then
	! 							write(*,'(6i6)'),j,RawReads(i,j,:),TrueGenos(i,j),FilledGenos(i,j),ParentGeno(j)
	! 						endif
	! 					endif


	! 					endif
	! 				enddo


	! 	!				p=1
	! 	!				do hs=1,nInd
	! 	!					if (RecPed(hs,k)==i) then
	! 	!						if (sum(HalfSibPhase(p,j,:))<3) then
								
	! 	!							if (FilledGenos(hs,j)==9) then
	! 	!								FilledGenos(hs,j)=HalfSibGeno(p,j)
	! 	!							endif

	! 							!if (HalfSibGeno(p,j) /= FilledGenos(hs,j)) then
	! 							!	FilledGenos(hs,j)=9
	! 							!endif

	! 	!							if (sum(FilledPhase(hs,j,:))==9) then
	! 	!								FilledPhase(hs,j,:)=HalfSibPhase(p,j,:)
	! 	!							endif

	! 							!if (HalfSibPhase(p,j,1) /= FilledPhase(hs,j,1)) then
	! 							!	FilledPhase(hs,j,:)=9
	! 							!endif

	! 							!if (HalfSibPhase(p,j,2) /= FilledPhase(hs,j,2)) then
	! 							!	FilledPhase(hs,j,:)=9
	! 							!endif
	! 	!							p=p+1
	! 	!						endif
	! 	!					endif
	! 	!				enddo


	! 				deallocate(ParentPhase)
	! 				deallocate(ParentGeno)
	! 				deallocate(HalfSibGeno)
	! 				deallocate(HalfSibPhase)

	! 			endif

	! 			deallocate(HasOwnReads)
	! 			deallocate(HalfSibID)


	! 		enddo
	! 	enddo
 !end subroutine FerdosiSerap
!################################################################################################

! CoreOfTheProgram 5c - Find the consensus between founder-parent-kids haplotype
!  Here we use the chunks previously define, and we try to fill the missing value of these 
!  three individuals if one of the three has information to do that.

subroutine BuildConsensus 
	use GlobalPar

	implicit none

	integer :: i,k,e,j,m,Count1,Count0

	integer,allocatable,dimension(:) :: ConsensusHaplotype
	integer(int64),allocatable,dimension(:,:,:) :: ConsensusIds

	allocate(ConsensusIds(nInd,4,2))
	allocate(ConsensusHaplotype(nSnp))

	do i=1,nInd
		do k=1,2 ! paternal or maternal gamete
			e=k+1 ! position parents in Pedigree

			if (maxval(FounderAssignment(i,:,k))/=0) then
				do j=1,nSnp
					if (MarkersToExclude(j).eq.0) then
						
						ConsensusIds(i,:,:)=0
						if (FounderAssignment(i,j,k)/=0) then

							! Save Ids of the trio (founder-parent-proband)

							ConsensusIds(i,4,1)=RecPed(i,1) 													! Id proband
							ConsensusIds(i,3,1)=RecPed(i,e) 													! Id parent
							
							ConsensusIds(i,1,1)=RecPed(ConsensusIds(i,3,1),FounderAssignment(i,j,k))			! Id Founder
							ConsensusIds(i,2,1)=RecPed(ConsensusIds(i,3,1),FounderAssignment(i,j,k))			! Id Founder

							! Save gamete that the trio shared

							ConsensusIds(i,1,2)=1 	! gamete1 Founder
							ConsensusIds(i,2,2)=2 	! gamete2 Founder

							if (FounderAssignment(i,j,k)==2) ConsensusIds(i,3,2)=1	! From male 
							if (FounderAssignment(i,j,k)==3) ConsensusIds(i,3,2)=2  ! From female
							
							ConsensusIds(i,4,2)=k 	
								
							ConsensusHaplotype(j)=9

							Count1=0
							Count0=0

							if (sum(FilledPhase(ConsensusIds(i,1,1),j,:))<3) then
								do m=1,2 ! members and gamets ! Avoid to use markers that are not fully phased
									if (FilledPhase(ConsensusIds(i,m,1),j,ConsensusIds(i,m,2))==1) Count1=Count1+1
									if (FilledPhase(ConsensusIds(i,m,1),j,ConsensusIds(i,m,2))==0) Count0=Count0+1
								enddo
							endif

							do m=3,4 ! members and gamets
								if (FilledPhase(ConsensusIds(i,m,1),j,ConsensusIds(i,m,2))==1) Count1=Count1+1
								if (FilledPhase(ConsensusIds(i,m,1),j,ConsensusIds(i,m,2))==0) Count0=Count0+1
							enddo

							if (Count1>1 .and. Count1.gt.Count0) ConsensusHaplotype(j)=1 !Count0==0
							if (Count0>1 .and. Count0.gt.Count1) ConsensusHaplotype(j)=0 !Count1==0


							if (ConsensusHaplotype(j)/=9) then
								do m=3,4
									FilledPhase(ConsensusIds(i,m,1),j,ConsensusIds(i,m,2))=ConsensusHaplotype(j)
									if ((FilledGenos(ConsensusIds(i,m,1),j)/=9).and.(sum(FilledPhase(ConsensusIds(i,m,1),j,:))<3)) then
										if (FilledGenos(ConsensusIds(i,m,1),j)/=sum(FilledPhase(ConsensusIds(i,m,1),j,:))) then
											FilledPhase(ConsensusIds(i,m,1),j,ConsensusIds(i,m,2))=9
										endif
									endif
								enddo
							endif


							!AlternativeGamete=Abs((ConsensusIds(i,4,2)-1)-1)+1 ! gamete for the other founder (i.e. if i=male --> AlternativeGamete=female)
							!if (ConsensusHaplotype(j)==0) then
							!	if ((FilledPhase(ConsensusIds(i,4,1),j,AlternativeGamete)==9).and.(RawReads(ConsensusIds(i,4,1),j,2)>=ReadCountMinThresh)) then
							!		FilledPhase(ConsensusIds(i,4,1),j,AlternativeGamete)=1
							!	endif
							!endif

							!if (ConsensusHaplotype(j)==1) then
							!	if ((FilledPhase(ConsensusIds(i,4,1),j,AlternativeGamete)==9).and.(RawReads(ConsensusIds(i,4,1),j,1)>=ReadCountMinThresh)) then
							!		FilledPhase(ConsensusIds(i,4,1),j,AlternativeGamete)=0
							!	endif
							!endif

						
						endif
					endif
				enddo
			endif
		enddo
	enddo

	deallocate(ConsensusIds)
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

!#####################################################################################################################

subroutine CountFounder
	use GlobalPar
	implicit none

	integer :: i,j,k,f1,f2,c1,c2

	do i=1,nInd
		do k=1,2
			if (maxval(FounderAssignment(i,:,k))/=0) then
				f1=0
				f2=0
				c1=0
				c2=0
				do j=1,nSnp
					if (FounderAssignment(i,j,k)/=0) then
						if (FounderAssignment(i,j,k)==f1) c1=c1+1
						if (FounderAssignment(i,j,k)==f2) c2=c2+1
						if (f1==0) then
							f1=FounderAssignment(i,j,k)
							c1=c1+1
						endif
						if ((f1/=0).and.(FounderAssignment(i,j,k)/=f1).and.(f2==0)) then
							f2=FounderAssignment(i,j,k)
							c2=c2+1
						endif
					endif
				enddo
				write(*,'(7i10)'),i,Ped(i,1),k,f1,c1,f2,c2
			endif
		enddo
	enddo
end subroutine CountFounder

!#####################################################################################################################

! CoreOfTheProgram 5a - Assign from which founder the haplotype has been inherited.
!  This is done using the heterozygous markers of the individual’s parent
!  i.e. DAD phase = 10, KID Phase=19. The "1" of the kid is coming from the Paternal GrandSire
!  i.e. DAD phase = 01, KID Phase=19. The "1" of the kid is coming from the Paternal GrandDam

subroutine CalculateFounderAssignment
	use GlobalPar
	use omp_lib
	implicit none

	integer :: i,j,e,k
	
	FounderAssignment(:,:,:)=0
	
   	!$OMP PARALLEL DO ORDERED DEFAULT(PRIVATE) SHARED (FilledGenos,FilledPhase,RecPed,FounderAssignment,nSnp,nInd) !collapse(2)	
	do e=2,3 ! Sire and Dam pos in the ped
		do i=1,nInd
			do j=1,nSnp
				if (MarkersToExclude(j).eq.0) then
					!k=e-1 ! Sire and Dam gamete
					!if (RecPed(i,e-1)==0) exit
					if (FilledGenos(RecPed(i,e),j)==1) then
						if (sum(FilledPhase(RecPed(i,e),j,:))<3) then
							!if (FilledPhase(RecPed(i,e),j,1)==FilledPhase(i,j,e-1)) FounderAssignment(i,j,e-1)=RecPed(RecPed(i,e),2) !Sire of Parent !FounderId(k,1)
							!if (FilledPhase(RecPed(i,e),j,2)==FilledPhase(i,j,e-1)) FounderAssignment(i,j,e-1)=RecPed(RecPed(i,e),3) !Dam  of Parent !FounderId(k,2)
							if (FilledPhase(RecPed(i,e),j,1)==FilledPhase(i,j,e-1)) FounderAssignment(i,j,e-1)=2 !Sire of Parent !FounderId(k,1)
							if (FilledPhase(RecPed(i,e),j,2)==FilledPhase(i,j,e-1)) FounderAssignment(i,j,e-1)=3 !Dam  of Parent !FounderId(k,2)

						endif
					endif
				endif
			enddo
		enddo
	enddo
    !$OMP END PARALLEL DO
end subroutine CalculateFounderAssignment

!#####################################################################################################################

! CoreOfTheProgram 4 - Phasing of heterozygous genotypes for individuals 
!  that have their parent and their progeny with opposing homozygotes
!  i.e. The MGS=0, DAM=?, KID=2. The genotype of the dam has to be 1 and the phase is 01
!  from that parent at that marker can be solved.
!  TODO : probably this subroutine is not necessary after the introduction of GeneProb.
!         Check if it make sense to remove it.

subroutine SimpleFillInBasedOnProgenyReads
	! Fill missing gamete according to the one of progeny and the alternative gamete 

	use GlobalPar
	use omp_lib
	implicit none

	integer :: i,j
	integer(int64) :: IdSir,IdDam

	!$OMP PARALLEL DO ORDERED DEFAULT(PRIVATE) SHARED (FilledPhase,RecPed,nSnp, nInd) !collapse(2)
	do i=1,nInd
		do j=1,nSnp
			if ((MarkersToExclude(j).eq.0).and.(sum(FilledPhase(i,j,:))==0).or.(sum(FilledPhase(i,j,:))==2)) then
				IdSir=RecPed(i,2)
				IdDam=RecPed(i,3)
				!if (IdDam/=0) then
					if ((FilledPhase(IdDam,j,1)/=9).and.(FilledPhase(IdDam,j,1)/=FilledPhase(i,j,2))) FilledPhase(IdDam,j,2)=FilledPhase(i,j,2)
					if ((FilledPhase(IdDam,j,2)/=9).and.(FilledPhase(IdDam,j,2)/=FilledPhase(i,j,2))) FilledPhase(IdDam,j,1)=FilledPhase(i,j,2)
				!endif

				!if (IdSir/=0) then				
					if ((FilledPhase(IdSir,j,1)/=9).and.(FilledPhase(IdSir,j,1)/=FilledPhase(i,j,1))) FilledPhase(IdSir,j,2)=FilledPhase(i,j,1)
					if ((FilledPhase(IdSir,j,2)/=9).and.(FilledPhase(IdSir,j,2)/=FilledPhase(i,j,1))) FilledPhase(IdSir,j,1)=FilledPhase(i,j,1)
				!endif
			endif
		enddo
	enddo
	!$OMP END PARALLEL DO

end subroutine SimpleFillInBasedOnProgenyReads

!################################################################################################

! CoreOfTheProgram 3 - Phasing based on parents informative markers
!  If a parent has the homozygote markers known (filled), then the phase of the progeny for the gamete inherited 
!  from that parent at that marker can be solved.

subroutine SimpleFillInBasedOnParentsReads

	use GlobalPar
	use omp_lib
	implicit none

	integer :: i,j,e,k

    !$OMP PARALLEL DO ORDERED DEFAULT(PRIVATE) SHARED (FilledGenos,FilledPhase,RecPed,nInd,nSnp) !collapse(2)
	do e=1,2
		do i=1,nInd
			do j=1,nSnp
				!if (RecPed(i,e+1)==0) exit
				if ((MarkersToExclude(j).eq.0).and.(maxval(FilledPhase(i,j,:))==9)) then
					k=Abs((e-1)-1)+1
					if (FilledGenos(RecPed(i,e+1),j)==0) then
						FilledPhase(i,j,e)=0
					endif

					if (FilledGenos(RecPed(i,e+1),j)==2) then
						FilledPhase(i,j,e)=1
					endif
				endif
			enddo
		enddo
	enddo
    !$OMP END PARALLEL DO

end subroutine SimpleFillInBasedOnParentsReads

!################################################################################################

subroutine ReadSamFile
	
	use GlobalPar
	use omp_lib
	
	implicit none

	integer :: i,j,k,HapLength,tmpId,wp,wm,cp,cm
	integer,allocatable,dimension(:) :: SamHap,HapPos
	real 	:: LikeP,LikeM
	real(kind=8) :: p0,m0
	character(len=100) :: filout12

	allocate(HapPos(5000))
	allocate(SamHap(5000))

	HapPos=0
	SamHap=9

	do i=nInd,1,-1
		cp=0
		cm=0
		wp=0
		wm=0
		write (filout12,'(i0,".FilledPhase.txt")') Ped(i,1)
		open(unit=12, file=trim(filout12), status="unknown")

		do 
			read(12,*,end=12),tmpId,HapLength,(HapPos(k),k=1,HapLength)
	    	read(12,*),tmpId,HapLength,(SamHap(k),k=1,HapLength)
	    	!print*,tmpId,HapLength,HapPos(1:HapLength)

	    	LikeP=log(.5)
	    	LikeM=log(.5)

		    do j=1,HapLength
		        p0=Pr00(i,HapPos(j))+Pr01(i,HapPos(j))
		        m0=Pr00(i,HapPos(j))+Pr10(i,HapPos(j))
		        
		        if (p0.lt.0.0000001) p0=0.0000001
		        if (m0.lt.0.0000001) m0=0.0000001

		        if (p0.gt.0.9999999) p0=0.9999999
		        if (m0.gt.0.9999999) m0=0.9999999

		      if (SamHap(j)==0) then
		        LikeP=LikeP+log(p0)
		        LikeM=LikeM+log(m0)
		      else if (SamHap(j)==1) then
		        LikeP=LikeP+log(abs(p0-1))
		        LikeM=LikeM+log(abs(m0-1))
		      endif
		    enddo

		    !if (i==1) print*,LikeP,LikeM
		    if (LikeP.gt.LikeM) then ! The haplotype is paternal
		      do j=1,HapLength
		        if (FilledPhase(i,HapPos(j),1)==9) then
		        	cp=cp+1
		        	FilledPhase(i,HapPos(j),1)=SamHap(j)
		        	FilledGenos(i,HapPos(j))=1
		        endif
		        if (FilledPhase(i,HapPos(j),2)==9) then
		        	cm=cm+1
		        	FilledPhase(i,HapPos(j),2)=abs(SamHap(j)-1)
		        	FilledGenos(i,HapPos(j))=1
		        endif

		      	if (FilledPhase(i,HapPos(j),1)/=SamHap(j)) then
		      		wp=wp+1
		      		FilledPhase(i,HapPos(j),1)=9
		      		FilledGenos(i,HapPos(j))=1
		      	endif
		      	if (FilledPhase(i,HapPos(j),2)/=abs(SamHap(j)-1)) then
		      		wm=wm+1
		      		FilledPhase(i,HapPos(j),1)=9
		      		FilledGenos(i,HapPos(j))=1
		      	endif
		      enddo
		    endif
		    
		    if (LikeM.gt.LikeP) then
		      do j=1,HapLength
		        if (FilledPhase(i,HapPos(j),2)==9) then
		        	cm=cm+1
		        	FilledPhase(i,HapPos(j),2)=SamHap(j)
		        	FilledGenos(i,HapPos(j))=1
		        endif
		        if (FilledPhase(i,HapPos(j),1)==9) then
		        	cp=cp+1
		        	FilledPhase(i,HapPos(j),1)=abs(SamHap(j)-1)
		        	FilledGenos(i,HapPos(j))=1
		        endif

		      	if (FilledPhase(i,HapPos(j),2)/=SamHap(j)) then
		      		wm=wm+1
		      		FilledPhase(i,HapPos(j),1)=9
		      		FilledGenos(i,HapPos(j))=1
		      	endif
		      	if (FilledPhase(i,HapPos(j),1)/=abs(SamHap(j)-1)) then
		      		wp=wp+1
		      		FilledPhase(i,HapPos(j),1)=9
		      		FilledGenos(i,HapPos(j))=1
		      	endif
		      enddo
		    endif

		   	if (LikeM.eq.LikeP) then
		      do j=1,HapLength
		      	if (FilledGenos(i,HapPos(j))/=1) then
		      		FilledGenos(i,HapPos(j))=1
		      		FilledPhase(i,HapPos(j),:)=9
		      	endif
		      enddo
		    endif


		enddo
		12 close(12)
		print*,tmpId,cp,cm,wp,wm
	enddo

end subroutine ReadSamFile

!################################################################################################

! CoreOfTheProgram 2 - Use the output of geneprob to fill the genotypes and the phase of the individuals
!  This step run multiple time. At each time the genotype Threshold (i.e., threshold that measures the 
!  ... of calling a genotype) is decreased until the value specified by the user.
!  Input : Probabilities from the RunGeneProb
!  Output : FilledPhase and FilledGenos

subroutine UseGeneProbToSimpleFillInBasedOnOwnReads
	
	use GlobalPar
	use omp_lib
	
	implicit none

    integer :: i,j
    !integer(kind=2) :: probscore

	!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (Pr00,Pr11,Pr01,Pr10,FilledGenos,FilledPhase,GeneProbThresh,nSnp,nInd) !collapse(2)
	do i=1,nInd
		do j=1,nSnp

			if ((Pr00(i,j).ge.GeneProbThresh).and.(sum(FilledPhase(i,j,:))>3)) then
				FilledGenos(i,j)=0
				FilledPhase(i,j,:)=0
			endif

			if (((Pr01(i,j)+Pr10(i,j)).ge.GeneProbThresh).and.(sum(FilledPhase(i,j,:))>3)) then
				FilledGenos(i,j)=1
				if (Pr01(i,j).ge.GeneProbThresh) then
					FilledPhase(i,j,1)=0
					FilledPhase(i,j,2)=1
				endif
				if (Pr10(i,j).ge.GeneProbThresh) then
					FilledPhase(i,j,1)=1
					FilledPhase(i,j,2)=0
				endif
			endif

			if ((Pr11(i,j).ge.GeneProbThresh).and.(sum(FilledPhase(i,j,:))>3)) then
				FilledGenos(i,j)=2
				FilledPhase(i,j,:)=1
			endif

			if (((Pr00(i,j)+Pr01(i,j)).ge.GeneProbThresh).and.(Pr10(i,j).lt.Pr01(i,j)).and.(FilledPhase(i,j,1)==9)) FilledPhase(i,j,1)=0
			if (((Pr11(i,j)+Pr10(i,j)).ge.GeneProbThresh).and.(Pr01(i,j).lt.Pr10(i,j)).and.(FilledPhase(i,j,1)==9)) FilledPhase(i,j,1)=1

			if (((Pr00(i,j)+Pr10(i,j)).ge.GeneProbThresh).and.(Pr01(i,j).lt.Pr10(i,j)).and.(FilledPhase(i,j,2)==9)) FilledPhase(i,j,2)=0
			if (((Pr11(i,j)+Pr01(i,j)).ge.GeneProbThresh).and.(Pr10(i,j).lt.Pr01(i,j)).and.(FilledPhase(i,j,2)==9)) FilledPhase(i,j,2)=1

		enddo

	enddo
	!$OMP END PARALLEL DO
end subroutine UseGeneProbToSimpleFillInBasedOnOwnReads

!################################################################################################

! CoreOfTheProgram 1 - Use geneprob to call the variants.
!  This step run only one time. 
!  Input : Pedigree and the RawReads to call the variants.
!  Output : 4 array of dimension(nInd,nSnp) with probabilities for the genotype 00, 01, 10 and 11
subroutine RunGeneProb
	
	use GlobalPar
	use AlphaVarCallFuture
	!use omp_lib
	
	implicit none

	integer(int64),allocatable,dimension(:) :: SeqSire,SeqDam !SeqId
	
	!integer :: i,j
	!real(kind=8),allocatable,dimension(:,:,:) :: ReadCounts
    
	
	!allocate(SeqId(nInd))
    allocate(SeqSire(nInd))
    allocate(SeqDam(nInd))

    !SeqId=RecPed(1:nInd,1)
	SeqSire=RecPed(1:nInd,2)
	SeqDam=RecPed(1:nInd,3)
	
	
	call AlphaVarCall(nInd,nSnp,1,nSnp,ErrorRate,0,SeqSire,SeqDam,RawReads,FilledGenos(1:nInd,nSnp),Pr00,Pr01,Pr10,Pr11)

	!deallocate(SeqId)
	deallocate(SeqSire)
	deallocate(Seqdam)
end subroutine RunGeneProb

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

!################################################################################################

subroutine AllocateArrays

	use GlobalPar

	implicit none

	! New phase/geno
	allocate(FilledGenos(0:nInd,nSnp))
	allocate(FilledPhase(0:nInd,nSnp,2))

	! GeneProb
	allocate(Pr00(nInd,nSnp))
	allocate(Pr01(nInd,nSnp))
	allocate(Pr10(nInd,nSnp))
	allocate(Pr11(nInd,nSnp))

	! Founders
	allocate(FounderAssignment(nInd,nSnp,2))
	!allocate(ConsensusFounders(nInd,2))

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

	! Input Files
	!deallocate(Id)
 	!if (trim(GenoFile)/="None")  deallocate(TrueGenos)
	!if (trim(PhaseFile)/="None")  deallocate(TruePhase)

	deallocate(RawReads)
	IF( ALLOCATED(position)) DEALLOCATE(position) 
	IF( ALLOCATED(quality)) DEALLOCATE(quality) 

	! New phase/geno
	deallocate(FilledGenos)
	deallocate(FilledPhase)

	! GeneProb
	deallocate(Pr00)
	deallocate(Pr01)
	deallocate(Pr10)
	deallocate(Pr11)

	! Founders
	deallocate(FounderAssignment)
	
	! CountPhase
	if (allocated(CurrentCountID%old)) deallocate(CurrentCountID%old)
	if (allocated(CurrentCountID%diff)) deallocate(CurrentCountID%diff)

	! Markers to exclude
	deallocate(MarkersToExclude)

	! Checks
	!deallocate(CheckPhase)
	!deallocate(CheckGenos)
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

subroutine InitialiseArrays

	use GlobalPar

	implicit none

	!integer(kind=2) :: probscore


	FilledPhase=9
	FilledGenos=9

	Pr00=0.0
	Pr01=0.0
	Pr10=0.0
	Pr11=0.0

	!ConsensusFounders=0
end subroutine InitialiseArrays

!################################################################################################

subroutine InternalEdititing
	use GlobalPar
	use omp_lib
	implicit none

	integer :: i,j
	real(kind=4),allocatable,dimension(:) :: RefAllele
	real(kind=4) :: maf

	allocate(RefAllele(nSnp))
	allocate(GeneProbYesOrNo(nSnp))
	GeneProbYesOrNo=1
	i=0

	open (unit=1,file="AlphaFamSeqEditedSnp.txt",status="unknown")
	print*,"Internal Editing"
	write(*,'(1a10,1f10.3)'),"(1) MAF",EditingParameter
	
	!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nSnp,RawReads,RefAllele,EditingParameter,GeneProbYesOrNo,i)
	do j=1,nSnp
		RefAllele(j)=real(sum(RawReads(:,j,1)))
		maf=RefAllele(j)/sum(RawReads(:,j,:))
		if ((maf.le.EditingParameter).or.(maf.ge.(1-EditingParameter))) then
			GeneProbYesOrNo(j)=0
			i=i+1
			write(1,'(1i10,1f10.3,2i10)'),j,maf,sum(RawReads(:,j,1)),sum(RawReads(:,j,2))
		endif
	end do
	!$OMP END PARALLEL DO
	
	close(1)

	if (i==0) print*,"All SNPs passed the editing"
	if (i>0)  write(*,'(1i12,1a3,1i12,1a30)'), i,"on",nSnp,"didn't pass the editing"

	deallocate(RefAllele)
end subroutine InternalEdititing

!################################################################################################

subroutine ReadPedigree

	use GlobalPar
	use MaraModule

	implicit none

	integer(int64) :: i,DumI,j,stat
	!integer,allocatable,dimension(:) :: TempImput
	!integer :: TmpID

	open (unit=2,file=trim(PedigreeFile),status="old")
	
	nInd = 0
	do
	    read(2, *, iostat=stat) DumI
	    if (stat/=0) exit
	    nInd = nInd + 1
	enddo
	print*,"PedigreeLenght",nInd	 
	rewind(2)

	allocate(Ped(0:nInd,3))
	allocate(RecPed(0:nInd,3))
	RecPed=0
	Ped(0,:)=0

	do i=1,nInd
		read(2,*) Ped(i,:)
	enddo

	close (2)


	do i=1,nInd ! Recode the pedigree in a new array
		do j=1,3
			call GetID(Ped(i,j),RecPed(i,j))
		enddo
	enddo
end subroutine ReadPedigree

!###########################################################################################################################################################

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
	    call GetID(TmpID, PosReads)
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

subroutine GetID(InputId, PosId)

    use GlobalPar
    implicit none

    integer(int64), intent(in) :: InputId
    integer(int64), intent(out) :: PosId

    integer :: i,check

    PosId = 0
    check = 0

    do i=1, nInd 
        if (Ped(i,1) == InputId) then !Ped is in global module
            PosId = i
        endif
    enddo

    !if (PosId == 0) then ! just a check, you don't need this really but may be useful for testing!
    !    print*, "individual not present, please check file"
    !   stop
    !check=check+1
    !endif
end subroutine GetID

!###########################################################################################################################################################
subroutine ReadPrevGeneProb
	use GlobalPar

	implicit none

	integer :: i,j,DumI
	character(len=30) :: nChar
	character(len=80) :: FmtReal,filout5
	
	filout5="AlphaFamSeqFinalGeneProb1.txt"
	open (unit=5,file=trim(filout5),status="old")
	
	do i=1,nInd
	
		read(5,*)DumI,Pr00(i,:)
		read(5,*)DumI,Pr01(i,:)
		read(5,*)DumI,Pr10(i,:)
		read(5,*)DumI,Pr11(i,:)
	
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
	
	integer,allocatable,dimension(:) :: SeqColToKeep
	real,allocatable,dimension(:) 	 :: ReducedGeneProb

	! WriteOut ReducedFileOf GeneProb
	open (unit=2,file="SeqColToKeep.txt",status="old")
	
	nRow = 0
	do
	    read(2, *, iostat=stat) DumI
	    if (stat/=0) exit
	    nRow = nRow + 1
	enddo
	rewind(2)

	allocate(SeqColToKeep(nRow))
	allocate(ReducedGeneProb(nRow))

	SeqColToKeep(:)=0
	ReducedGeneProb(:)=0

	do i=1,nRow
		read(2,*) SeqColToKeep(i)
	enddo
	close (2)


	write(nChar,*) nRow
	FmtReal='(i0,'//trim(adjustl(nChar))//'f7.4)'
	write (filout6,'("AlphaFamSeqReducedGeneProb",i0,".txt")') Windows
	open (unit=6,file=trim(filout6),status="unknown")

	do i=1,nInd
		p=1
		do j=1,nRow
			ReducedGeneProb(p)=Pr00(i,SeqColToKeep(j))
			p=p+1
		enddo
		write (6,FmtReal) Ped(i,1),ReducedGeneProb(:)
		p=1
		do j=1,nRow
			ReducedGeneProb(p)=Pr01(i,SeqColToKeep(j))
			p=p+1
		enddo
		write (6,FmtReal) Ped(i,1),ReducedGeneProb(:)
		p=1
		do j=1,nRow
			ReducedGeneProb(p)=Pr10(i,SeqColToKeep(j))
			p=p+1
		enddo
		write (6,FmtReal) Ped(i,1),ReducedGeneProb(:)
		p=1
		do j=1,nRow
			ReducedGeneProb(p)=Pr11(i,SeqColToKeep(j))
			p=p+1
		enddo
		write (6,FmtReal) Ped(i,1),ReducedGeneProb(:)
	enddo
	close(6)

	deallocate(SeqColToKeep)
	deallocate(ReducedGeneProb)

	! Write Out Full file of GeneProb

	write(nChar,*) nSnp
	FmtReal='(i0,'//trim(adjustl(nChar))//'f7.4)'
	write (filout5,'("AlphaFamSeqFinalGeneProb",i0,".txt")') Windows
	open (unit=5,file=trim(filout5),status="unknown")

	do i=1,nInd
		write (5,FmtReal) Ped(i,1),Pr00(i,:)
		write (5,FmtReal) Ped(i,1),Pr01(i,:)
		write (5,FmtReal) Ped(i,1),Pr10(i,:)
		write (5,FmtReal) Ped(i,1),Pr11(i,:)
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
		
	integer,allocatable,dimension(:) 			:: SeqColToKeep
	integer(kind=1),allocatable,dimension(:) 	:: ReducedGenos

	! Save Reduced Genotypes
	open (unit=6,file="SeqColToKeep.txt",status="old")
	
	nRow = 0
	do
	    read(6, *, iostat=stat) DumI
	    if (stat/=0) exit
	    nRow = nRow + 1
	enddo
	rewind(6)

	allocate(SeqColToKeep(nRow))
	allocate(ReducedGenos(nRow))
	
	SeqColToKeep(:)=0
	ReducedGenos(:)=9

	do i=1,nRow
		read(6,*) SeqColToKeep(i)
	enddo
	close (6)


	write(nChar,*) nRow
	FmtInt='(i0,'//trim(adjustl(nChar))//'i2)'
	write (filout7,'("AlphaFamSeqReducedGenos",i0,".txt")') Windows
	open (unit=7,file=trim(filout7),status="unknown")

	do i=1,nInd
		p=1
		do j=1,nRow
			ReducedGenos(p)=FilledGenos(i,SeqColToKeep(j))
			p=p+1
		enddo
		write (7,FmtInt) Ped(i,1),ReducedGenos(:)
	enddo
	close(7)

	deallocate(SeqColToKeep)
	deallocate(ReducedGenos)







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

