!################################################################################################
!include "Ferdosi.f90"
!include "AlphaVarCallParallelised.f90"
!include "ReadRogerData.f90"

module GlobalPar
	
	use ISO_Fortran_Env
	implicit none

	integer :: nIndSeq                         					! SpecFile - Number of Individuals in the sequence file
	integer :: LenghtSequenceDataFile,nSnp,fistWindow			! SpecFile - Total number of Snps
	integer :: nInd 											! number of Individuals in the Pedigree
	integer :: InternalEdit                    					! SpecFile - Internal Edit 1==yes or 0==no
	real(kind=8) :: EditingParameter							! SpecFile - 1st Number is the MAF (excluede SNP with MAF=<EditingParameter)


	real(kind=8) :: GeneProbThresh  							! SpecFile - Threshold to call a genotype from the probabilities First Value
	real(kind=8) :: GeneProbThreshMin							! SpecFile - Threshold to call a genotype from the probabilities Last Value
	real(kind=8) :: ReduceThr 									! SpecFile - Reduce Geno Treshold factor
	!integer :: nIter                        					! SpecFile - Every nIter the ReadContMaxThresh decrease of 1 value until GeneProbThresh
	
	real(kind=8) :: ErrorRate									! SpecFile - Error rates to define the genotypes probabilities
	
	integer :: ChunkLengthA                 					! SpecFile - First value to define Haplotypes length
	integer :: ChunkLengthB                 					! SpecFile - Last value to define Haplotypes length

	integer :: SuperC                           				! SpecFile - Parameter to use/not use Ferdosi module

	character(len=300) :: PedigreeFile      					! Input File Name - Pedigree
	character(len=300) :: ReadsFile             				! Input File Name - Reads Count for the Reference and the Alternative Allele
	character(len=300) :: ReadsType             				! Input File Name - Reads Count for the Reference and the Alternative Allele
	
	character(len=300) :: MapFile             					! Input File Name - Map File - position of the Variants
	character(len=300) :: SnpChipsInformation   				! Input File Name - Snp array to add more information to the Reads
	
	character(len=300) :: GenoFile              				! Control Results File Name - TrueGeno Genotypes to check results 
	character(len=300) :: PhaseFile             				! Control Results File Name - True Phase to check results 

	integer :: IterationNumber                  				! Control Parameter - Define the number of Iterations
	integer :: CurrentCountFilledPhase          				! Control Parameter - used to finish the program
	integer :: CurrentCountFilledGenos          				! Control Parameter - used to finish the program
	integer :: SolutionChanged                  				! Control Parameter - used to finish the program 
	integer :: StartSnp,EndSnp
	
	integer,allocatable,dimension(:,:) :: Ped               	! Input File - Pedigree
	integer,allocatable,dimension(:,:) :: RecPed				! Temporary File - Pedigree Recoded
	integer,allocatable,dimension(:) :: Id                      ! Read Data - used to read unsorted data

	integer(int32),allocatable,dimension(:,:,:) :: SequenceData	! Input File - Snp array to add more information to the Reads
	real(kind=8),allocatable,dimension(:,:,:) :: RawReads		! Input File - Snp array to add more information to the Reads
	character(len=100), allocatable, dimension(:) :: Ids
	integer(int32), dimension(:), allocatable :: position
	real(real64), allocatable, dimension(:) :: quality

	integer,allocatable,dimension(:,:) :: TrueGenos				! Control Results - True Genotypes to check results 
	integer,allocatable,dimension(:,:,:) :: TruePhase				! Control Results - True Phase to check results 

	
	character(len=1),allocatable,dimension(:,:) :: CheckGenos   ! Control Results - Use character to check True vs Imputed Genotypes
	character(len=1),allocatable,dimension(:,:,:) :: CheckPhase ! Control Results - Use character to check True vs Imputed Phase
	character(len=1),allocatable,dimension(:,:) :: CharPhase    ! Control Results - Use character to check True vs Imputed Genotypes
	
	integer,allocatable,dimension(:) :: GeneProbYesOrNo			! Temporary Array - use gene prob or not
	integer,allocatable,dimension(:,:,:) :: FounderAssignment   ! Temporary File - Save the IDs of the grandparents

	real(kind=4),allocatable,dimension(:,:) :: Pr00   			! Output GeneProb - Probabilities for Ind(i) and Spn(j) to be Homozygote for Reference Allele 
    real(kind=4),allocatable,dimension(:,:) :: Pr01	  			! Output GeneProb - Probabilities for Ind(i) and Spn(j) to be Heterozygote (0 from dad, 1 from mum)
    real(kind=4),allocatable,dimension(:,:) :: Pr10				! Output GeneProb - Probabilities for Ind(i) and Spn(j) to be Heterozygote (1 from dad, 0 from mum)
    real(kind=4),allocatable,dimension(:,:) :: Pr11				! Output GeneProb - Probabilities for Ind(i) and Spn(j) to be Homozygote for Alternative Allele 


	integer,allocatable,dimension(:,:) :: FilledGenos 			! Output - Imputed Genotypes
	integer,allocatable,dimension(:,:,:) :: FilledPhase  		! Output - Imputed Phase



	integer :: Windows
end module GlobalPar

!################################################################################################

program FamilyPhase

	use GlobalPar

	implicit none

	integer :: OldCount,NewCount
	real(kind=8) :: InitialGeneProbThresh

	call ReadSpecfile
	call ReadPedigree
	
	Windows=fistWindow
	EndSnp=(Windows*nSnp)
	StartSnp=EndSnp-nSnp+1
	if (EndSnp>LenghtSequenceDataFile) EndSnp=LenghtSequenceDataFile
	
	InitialGeneProbThresh=GeneProbThresh

	do while(StartSnp.le.LenghtSequenceDataFile)
		write(*,'(1a10,2i10,1a1,1i10)')"Windows",Windows,StartSnp,"-",EndSnp

		CurrentCountFilledPhase=0
		CurrentCountFilledGenos=0
		GeneProbThresh=InitialGeneProbThresh

		!call AllocateArrays
		call ReadData
		call InitialiseArrays
		call CheckMissingData
		call RunGeneProb
		write (*,'(1a39)') "      Iter   ProbThr    %Phase   %Geno"


		IterationNumber=0
		SolutionChanged=1
		do while ((SolutionChanged==1))!.and.(IterationNumber<2))

			OldCount=CurrentCountFilledPhase+CurrentCountFilledGenos
			IterationNumber=IterationNumber+1
			
			if ((IterationNumber>1).and.(GeneProbThresh>GeneProbThreshMin)) then 
				GeneProbThresh=GeneProbThresh-ReduceThr!0.001
			endif

			!print*,IterationNumber,GeneProbThresh,GeneProbThreshMin
			!if (InternalEdit==1) then
			!	call SimpleFillInBasedOnOwnReads   ! Old Subroutine that don't use AlphaVarCall
			!	call SimpleCleanUpFillIn
			!	call CurrentCountFilled
			!endif

			!if ((GeneProbThresh>GeneProbThreshMin).or.(IterationNumber==1)) then
			!	if (GeneProbThresh<GeneProbThreshMin) GeneProbThresh=GeneProbThreshMin
				call UseGeneProbToSimpleFillInBasedOnOwnReads
			!endif
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
				!print*,"Run Gene Prob Filled Genos"
				!call RunGeneProbWithFilledGenos
				!call UseGeneProbToSimpleFillInBasedOnOwnReads
				!call CurrentCountFilled
				!NewCount=CurrentCountFilledPhase+CurrentCountFilledGenos
				!if (OldCount/=NewCount) SolutionChanged=1
			endif

			write (*,'(1i10,3f10.5)') IterationNumber,GeneProbThresh,(dble(CurrentCountFilledPhase)/(dble(nInd*nSnp*2))*100),(dble(CurrentCountFilledGenos)/(dble(nInd*nSnp))*100)
		enddo

		if ((trim(PhaseFile)/="None").or.(trim(GenoFile)/="None")) then 
			call Checker
		endif
		call WriteResults
		call WriteStatistics

		StartSnp=EndSnp+1
		EndSnp=EndSnp+nSnp
		if (EndSnp>LenghtSequenceDataFile) then 
			EndSnp=LenghtSequenceDataFile
			nSnp=EndSnp-StartSnp+1
		endif
		Windows=Windows+1

		call DeallocateArrays

	enddo
end program FamilyPhase

!###########################################################################################################################################################

subroutine ProcessLine(Var)

	use GlobalPar

	implicit none

	integer :: nChrom
    integer, allocatable, dimension(:) :: nSnpPerChrom
    
	character(len=256), intent(inout) :: Var
    character (len=512) :: TLC
    character(len=8) :: temp
    character(len=1) :: comma
    character(len=512) :: snps

    integer :: i

    allocate(nSnpPerChrom(nChrom))

	read(Var,'(A,a8)') comma, temp

	if (trim(TLC(temp))=='constant') then
		read(Var, '(A,a8,A,i)') comma, temp, comma, nSnp
		do i=2, nChrom
			nSnpPerChrom(i) = nSnpPerChrom(1)
		enddo
	else
		read(Var, '(A,a8,A,a512)') comma, temp, comma, snps
		read(snps, *) nSnpPerChrom(1:nChrom)
	endif
end subroutine ProcessLine

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

! reads in and initialises specfile parameters 
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

                
            case('genotypeprobability')
                read(1, *, iostat=stat) GeneProbThresh,GeneProbThreshMin,ReduceThr !nIter 
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

subroutine UseSnpChipInformation

	use GlobalPar

	implicit none

	integer:: FileLength,stat,DumI,i,PosGeno,j
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

subroutine CheckMissingData

	use GlobalPar

	implicit none

	integer :: j
	character(len=50) :: filout1

	write (filout1,'("AlphaFamSeqMarkersWithZeroReads",i0,".txt")') Windows
	open (unit=1,file=trim(filout1),status="unknown")
		
	!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (RawReads,nSnp)
		do j=1,nSnp
			if (maxval(RawReads(:,j,:))==0) then
				write (1,'(i10)') j
			endif
		enddo
	!$OMP END PARALLEL DO

	close (1)

end subroutine CheckMissingData

!###########################################################################################

subroutine SimpleCleanUpFillIn

	use GlobalPar

	implicit none

	integer i,j,Change !,e,k,g1,g2

	Change=1
	do while (Change==1)
		Change=0
		!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (FilledGenos,FilledPhase,nInd,nSnp)
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

			enddo
		enddo
		!$OMP END PARALLEL DO
	enddo
end subroutine SimpleCleanUpFillIn

!###########################################################################################

! subroutine FerdosiSerap

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
! end subroutine FerdosiSerap

!###########################################################################################

subroutine BuildConsensus 
	use GlobalPar

	implicit none

	integer :: i,k,e,j,m,Count1,Count0

	integer,allocatable,dimension(:) :: ConsensusHaplotype
	integer,allocatable,dimension(:,:,:) :: ConsensusIds

	allocate(ConsensusIds(nInd,4,2))
	allocate(ConsensusHaplotype(nSnp))

	do i=1,nInd
		do k=1,2 ! paternal or maternal gamete
			e=k+1 ! position parents in Pedigree

			if (maxval(FounderAssignment(i,:,k))/=0) then
				do j=1,nSnp
						
					ConsensusIds(i,:,:)=0
					if (FounderAssignment(i,j,k)/=0) then

						! Save Ids of the trio (founder-parent-proband)

						ConsensusIds(i,1,1)=FounderAssignment(i,j,k)   			! Id Founder
						ConsensusIds(i,2,1)=FounderAssignment(i,j,k)			! Id Founder
						ConsensusIds(i,3,1)=RecPed(i,e) 						! Id parent
						ConsensusIds(i,4,1)=RecPed(i,1) 						! Id proband

						! Save gamete that the trio shared

						ConsensusIds(i,1,2)=1 	! gamete1 Founder
						ConsensusIds(i,2,2)=2 	! gamete2 Founder

						if(RecPed(ConsensusIds(i,3,1),2)==FounderAssignment(i,j,k)) then
							ConsensusIds(i,3,2)=1	! From male 
						endif

						if(RecPed(ConsensusIds(i,3,1),3)==FounderAssignment(i,j,k)) then
							ConsensusIds(i,3,2)=2  ! From female
						endif

						ConsensusIds(i,4,2)=k 	
							
						ConsensusHaplotype(j)=9

						Count1=0
						Count0=0

						do m=1,2 ! members and gamets ! Avoid to use markers that are not fully phased
							if (sum(FilledPhase(ConsensusIds(i,m,1),j,:))<3) then 
								if (FilledPhase(ConsensusIds(i,m,1),j,ConsensusIds(i,m,2))==1) Count1=Count1+1
								if (FilledPhase(ConsensusIds(i,m,1),j,ConsensusIds(i,m,2))==0) Count0=Count0+1
							endif
						enddo

						do m=3,4 ! members and gamets
							if (FilledPhase(ConsensusIds(i,m,1),j,ConsensusIds(i,m,2))==1) Count1=Count1+1
							if (FilledPhase(ConsensusIds(i,m,1),j,ConsensusIds(i,m,2))==0) Count0=Count0+1
						enddo

						if (Count1>0 .and. Count0==0) ConsensusHaplotype(j)=1
						if (Count0>0 .and. Count1==0) ConsensusHaplotype(j)=0


						if (ConsensusHaplotype(j)/=9) then
							do m=3,4								
								FilledPhase(ConsensusIds(i,m,1),j,ConsensusIds(i,m,2))=ConsensusHaplotype(j)
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

				enddo
			endif
		enddo
	enddo

	deallocate(ConsensusIds)
end subroutine BuildConsensus 

!################################################################################################

subroutine ChunkDefinition

	use GlobalPar
	implicit none

	integer :: i,j,e,k
	integer :: f1
	integer :: fSnp,ChunkLength2(1)
	integer :: lSnp
	real 	 Y, Z(1),A,B
	
	integer,allocatable,dimension(:) :: FounderAssignmentF,FounderAssignmentB

	allocate(FounderAssignmentF(1:nSnp))
	allocate(FounderAssignmentB(1:nSnp))

	CALL RANDOM_SEED()
	CALL RANDOM_NUMBER (HARVEST = Y)
	CALL RANDOM_NUMBER (Z)

	A = ChunkLengthA  
	B = ChunkLengthB 
	ChunkLength2=floor((A-1)+(B-(A-1))*Z)+1

	!print*,"CL",ChunkLength2(1)
	
	
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

subroutine CalculateFounderAssignment
	use GlobalPar

	implicit none

	integer :: i,j,e,k
	
	FounderAssignment(:,:,:)=0 
	
	!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (FilledGenos,FilledPhase,RecPed,FounderAssignment,nInd,nSnp)
	do i=1,nInd
		do e=2,3 ! Sire and Dam pos in the ped
			k=e-1 ! Sire and Dam gamete
			do j=1,nSnp
				if (FilledGenos(RecPed(i,e),j)==1) then
					if (sum(FilledPhase(RecPed(i,e),j,:))<3) then
						if (FilledPhase(RecPed(i,e),j,1)==FilledPhase(i,j,k)) FounderAssignment(i,j,k)=RecPed(RecPed(i,e),2) !Sire of Parent !FounderId(k,1)
						if (FilledPhase(RecPed(i,e),j,2)==FilledPhase(i,j,k)) FounderAssignment(i,j,k)=RecPed(RecPed(i,e),3) !Dam  of Parent !FounderId(k,2)
					endif
				endif
			enddo
		enddo
	enddo
	!$OMP END PARALLEL DO
end subroutine CalculateFounderAssignment

!#####################################################################################################################

subroutine SimpleFillInBasedOnProgenyReads
	! Fill missing gamete according to the one of progeny and the alternative gamete 

	use GlobalPar

	implicit none

	integer :: i,j,IdSire,IdDam

	do i=1,nInd

		!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (RawReads,FilledGenos,FilledPhase,RecPed,i,nSnp)
		do j=1,nSnp

			if ((FilledGenos(i,j)==0).or.(FilledGenos(i,j)==2)) then
				IdSire=RecPed(i,2)
				IdDam=RecPed(i,3)


				if ((FilledPhase(i,j,2)/=9).and.(FilledPhase(IdDam,j,1)/=9)) then
					if ((FilledPhase(IdDam,j,2)==9).and.(FilledPhase(i,j,2)/=FilledPhase(IdDam,j,1))) then
						FilledPhase(IdDam,j,2)=FilledPhase(i,j,2)
					endif
				endif

				if ((FilledPhase(i,j,2)/=9).and.(FilledPhase(IdDam,j,2)/=9)) then
					if ((FilledPhase(IdDam,j,1)==9).and.(FilledPhase(i,j,2)/=FilledPhase(IdDam,j,2))) then
						FilledPhase(IdDam,j,1)=FilledPhase(i,j,2)
					endif
				endif


				if ((FilledPhase(i,j,1)/=9).and.(FilledPhase(IdSire,j,1)/=9)) then
					if ((FilledPhase(IdSire,j,2)==9).and.(FilledPhase(i,j,1)/=FilledPhase(IdSire,j,1))) then
						FilledPhase(IdSire,j,2)=FilledPhase(i,j,1)
					endif
				endif

				if ((FilledPhase(i,j,1)/=9).and.(FilledPhase(IdSire,j,2)/=9)) then
					if ((FilledPhase(IdSire,j,1)==9).and.(FilledPhase(i,j,1)/=FilledPhase(IdSire,j,2))) then
						FilledPhase(IdSire,j,1)=FilledPhase(i,j,1)
					endif
				endif

				if ((FilledGenos(IdSire,j)==9).and.(sum(FilledPhase(IdSire,j,:))<3)) FilledGenos(IdSire,j)= sum(FilledPhase(IdSire,j,:))
				if ((FilledGenos(IdDam,j)==9).and.(sum(FilledPhase(IdDam,j,:))<3)) FilledGenos(IdDam,j)= sum(FilledPhase(IdDam,j,:))
			endif
		enddo
		!$OMP END PARALLEL DO
		
	enddo
end subroutine SimpleFillInBasedOnProgenyReads

!################################################################################################

subroutine SimpleFillInBasedOnParentsReads

	use GlobalPar

	implicit none

	integer :: i,j,e,k
	real :: nRef,nAlt,Pr0,Pr1,Pr2

	do i=1,nInd
	!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (RawReads,FilledGenos,FilledPhase,RecPed,ErrorRate,GeneProbThresh,i,nSnp)
		do j=1,nSnp
			if (maxval(FilledPhase(i,j,:))==9) then
				nRef=RawReads(i,j,1)
				nAlt=RawReads(i,j,2)
				
				call ReadsLikelihood(nRef,nAlt,ErrorRate,Pr0,Pr1,Pr2)

				do e=1,2
					k=Abs((e-1)-1)+1
					if (FilledGenos(RecPed(i,e+1),j)==0) then
						FilledPhase(i,j,e)=0
						if (Pr1.ge.GeneProbThresh) then
							FilledPhase(i,j,k)=1
							FilledGenos(i,j)=1
						endif
					endif

					if (FilledGenos(RecPed(i,e+1),j)==2) then
						FilledPhase(i,j,e)=1
						if (Pr1.ge.GeneProbThresh) then
							FilledPhase(i,j,k)=0
							FilledGenos(i,j)=1
						endif
					endif
				enddo
			endif
		enddo
	!$OMP END PARALLEL DO
	enddo
end subroutine SimpleFillInBasedOnParentsReads

!################################################################################################

subroutine UseGeneProbToSimpleFillInBasedOnOwnReads
	
	use GlobalPar
	use omp_lib
	
	implicit none

    integer :: i,j

    !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (Pr00,Pr11,Pr01,Pr10,FilledGenos,FilledPhase,nInd,nSnp,GeneProbThresh)
	do i=1,nInd
		do j=1,nSnp
			!if (sum(RawReads(i,j,:)).ge.10) write(*,'(3i10,6f10.4)'),i,j,Ped(i,1),RawReads(i,j,:),Pr00(i,j),Pr01(i,j),Pr10(i,j),Pr11(i,j)
			!if (i==7) write(*,'(3i10,6f10.4)'),i,j,Ped(i,1),RawReads(i,j,:),Pr00(i,j),Pr01(i,j),Pr10(i,j),Pr11(i,j)
			if ((Pr00(i,j)>=GeneProbThresh).and.(sum(FilledPhase(i,j,:))>3)) then
				FilledGenos(i,j)=0
				FilledPhase(i,j,:)=0
			endif

			if (((Pr01(i,j)+Pr10(i,j))>=GeneProbThresh).and.(sum(FilledPhase(i,j,:))>3)) then
				FilledGenos(i,j)=1
				if (Pr01(i,j)>=GeneProbThresh) then
					FilledPhase(i,j,1)=0
					FilledPhase(i,j,2)=1
				endif
				if (Pr10(i,j)>=GeneProbThresh) then
					FilledPhase(i,j,1)=1
					FilledPhase(i,j,2)=0
				endif
			endif

			if ((Pr11(i,j)>=GeneProbThresh).and.(sum(FilledPhase(i,j,:))>3)) then
				FilledGenos(i,j)=2
				FilledPhase(i,j,:)=1
			endif

		enddo
	enddo
	!$OMP END PARALLEL DO
end subroutine UseGeneProbToSimpleFillInBasedOnOwnReads


!################################################################################################


subroutine RunGeneProb
	
	use GlobalPar
	use AlphaVarCallFuture
	use omp_lib
	
	implicit none

	integer,allocatable,dimension(:) :: SeqSire,SeqDam !SeqId
	
	!integer :: i,j
	!real(kind=8),allocatable,dimension(:,:,:) :: ReadCounts
    
	
	!allocate(SeqId(nInd))
    allocate(SeqSire(nInd))
    allocate(SeqDam(nInd))

    !SeqId=RecPed(1:nInd,1)
	SeqSire=RecPed(1:nInd,2)
	SeqDam=RecPed(1:nInd,3)
	
	! do i=1,nInd
	! 	do j=1,nSnp
	! 		if (FilledGenos(i,j)==1) then
	! 			RawReads(i,j,:)=15.0
	! 		endif
			
	! 		if (FilledGenos(i,j)==0) then 
	! 			RawReads(i,j,1)=30.0
	! 			RawReads(i,j,2)=0.0
	! 		endif
			
	! 		if (FilledGenos(i,j)==2) then 
	! 			RawReads(i,j,1)=0.0
	! 			RawReads(i,j,2)=30.0
	! 		endif

	! 	enddo
	! enddo

	
	! Pr00=0.0
	! Pr01=0.0
	! Pr10=0.0
	! Pr11=0.0
	
	call AlphaVarCall(nInd,nSnp,1,nSnp,ErrorRate,0,SeqSire,SeqDam,RawReads,FilledGenos(1:nInd,nSnp),Pr00,Pr01,Pr10,Pr11)

	!deallocate(SeqId)
	deallocate(SeqSire)
	deallocate(Seqdam)

end subroutine RunGeneProb

!################################################################################################

subroutine RunGeneProbWithFilledGenos
	
	use GlobalPar
	use AlphaVarCallFuture
	use omp_lib
	
	implicit none

	integer,allocatable,dimension(:) :: SeqSire,SeqDam !SeqId,
	
	!allocate(SeqId(nInd))
    allocate(SeqSire(nInd))
    allocate(SeqDam(nInd))

    !SeqId=RecPed(1:nInd,1)
	SeqSire=RecPed(1:nInd,2)
	SeqDam=RecPed(1:nInd,3)
	
	! Pr00=0.0
	! Pr01=0.0
	! Pr10=0.0
	! Pr11=0.0
	 call AlphaVarCall(nInd,nSnp,1,nSnp,ErrorRate,1,SeqSire,SeqDam,RawReads,FilledGenos,Pr00,Pr01,Pr10,Pr11)

	!deallocate(SeqId)
	deallocate(SeqSire)
	deallocate(Seqdam)
end subroutine RunGeneProbWithFilledGenos

!################################################################################################

subroutine UseGeneProbToSimpleFillInBasedOnFilledGenos
	
	use GlobalPar
	use omp_lib
	
	implicit none

    integer :: i,j

    !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (Pr00,Pr11,Pr01,Pr10,FilledGenos,FilledPhase,nInd,nSnp,GeneProbThresh)
	do i=1,nInd
		do j=1,nSnp
			!if (sum(RawReads(i,j,:)).ge.10) write(*,'(3i10,6f10.4)'),i,j,Ped(i,1),RawReads(i,j,:),Pr00(i,j),Pr01(i,j),Pr10(i,j),Pr11(i,j)
			!if (i==7) write(*,'(3i10,6f10.4)'),i,j,Ped(i,1),RawReads(i,j,:),Pr00(i,j),Pr01(i,j),Pr10(i,j),Pr11(i,j)
			if ((Pr00(i,j)>=GeneProbThresh).and.(sum(FilledPhase(i,j,:))>3)) then
				FilledGenos(i,j)=0
				FilledPhase(i,j,:)=0
			endif

			if (((Pr01(i,j)+Pr10(i,j))>=GeneProbThresh).and.(sum(FilledPhase(i,j,:))>3)) then
				FilledGenos(i,j)=1
				if (Pr01(i,j)>=GeneProbThresh) then
					FilledPhase(i,j,1)=0
					FilledPhase(i,j,2)=1
				endif
				if (Pr10(i,j)>=GeneProbThresh) then
					FilledPhase(i,j,1)=1
					FilledPhase(i,j,2)=0
				endif
			endif

			if ((Pr11(i,j)>=GeneProbThresh).and.(sum(FilledPhase(i,j,:))>3)) then
				FilledGenos(i,j)=2
				FilledPhase(i,j,:)=1
			endif

		enddo
	enddo
	!$OMP END PARALLEL DO
end subroutine UseGeneProbToSimpleFillInBasedOnFilledGenos


!################################################################################################


subroutine SimpleFillInBasedOnOwnReads

	use GlobalPar

	implicit none

	integer :: i,j
	real :: nRef,nAlt,Pr0,Pr1,Pr2

	nRef=0.00
	nAlt=0.00
	!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (RawReads,FilledGenos,FilledPhase,nInd,nSnp,ErrorRate,GeneProbThresh,GeneProbYesOrNo)
	do i=1,nInd
		do j=1,nSnp
			
			!if (GeneProbYesOrNo(j)==0) then ! Only Snps that are Fixed
				if ((maxval(RawReads(i,j,:))/=0).and.(FilledGenos(i,j)/=9)) then
					nRef=RawReads(i,j,1)
					nAlt=RawReads(i,j,2)
					
					call ReadsLikelihood(nRef,nAlt,ErrorRate,Pr0,Pr1,Pr2)
					
					if (FilledGenos(i,j)==9) then 

						if (Pr1.ge.GeneProbThresh) then
							FilledGenos(i,j)=1
						endif

						if (Pr0.ge.GeneProbThresh) then  
							FilledPhase(i,j,:)=0
							FilledGenos(i,j)=0
						endif

						if (Pr2.ge.GeneProbThresh) then
							FilledPhase(i,j,:)=1
							FilledGenos(i,j)=2
						endif
						
					endif
				endif
			!endif
		enddo
	enddo
	!$OMP END PARALLEL DO
end subroutine SimpleFillInBasedOnOwnReads

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

	! Input Files
	allocate(Id(0:nInd))
 	if (trim(GenoFile)/="None")  allocate(TrueGenos(nInd,nSnp))
	if (trim(PhaseFile)/="None")  allocate(TruePhase(nInd,nSnp,2))

	! New phase/geno
	allocate(FilledGenos(0:nInd,nSnp))
	allocate(FilledPhase(0:nInd,nSnp,2))

	! GeneProb
	allocate(Pr00(nInd,nSnp))
	allocate(Pr01(nInd,nSnp))
	allocate(Pr10(nInd,nSnp))
	allocate(Pr11(nInd,nSnp))

	Pr00=0.0
	Pr01=0.0
	Pr10=0.0
	Pr11=0.0
	

	! Founders
	allocate(FounderAssignment(nInd,nSnp,2))
	!allocate(ConsensusFounders(nInd,2))

	! Checks
	allocate(CheckPhase(nInd,nSnp,2))
	allocate(CheckGenos(nInd,nSnp))
	allocate(CharPhase(nSnp,2))
end subroutine AllocateArrays

!################################################################################################

subroutine DeallocateArrays

	use GlobalPar

	implicit none

	! Input Files
	deallocate(Id)
 	if (trim(GenoFile)/="None")  deallocate(TrueGenos)
	if (trim(PhaseFile)/="None")  deallocate(TruePhase)

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
	
	! Checks
	deallocate(CheckPhase)
	deallocate(CheckGenos)
	deallocate(CharPhase)
end subroutine DeallocateArrays

!################################################################################################

subroutine CurrentCountFilled

	use GlobalPar

	implicit none

	!integer :: i,j,e,p
	character(len=30) :: FmtCha,nChar
	
	write(nChar,*) nSnp
	FmtCha='(i10,a2,'//trim(adjustl(nChar))//'a2)'

	CurrentCountFilledPhase=count(FilledPhase(:,:,:)/=9)
	CurrentCountFilledGenos=count(FilledGenos(:,:)/=9)


	!write (*,'(a20,5i10)') "Current Counts ", IterationNumber,CurrentCountFilledPhase,CurrentCountFilledGenos

	! ! Paste here checker
	! open (unit=99,file="AlphaFamSeqSummary.log",status="unknown")
	! !open (unit=98,file="AlphaFamSeqPhaseCorrectnessIters.txt",status="unknown")

	! CheckPhase='_'
	! CheckGenos='_'

	! do i=1,nInd
	! 	do j=1,nSnp
	! 		do e=1,2
	! 			if (FilledPhase(i,j,e)/=9) then
	! 				! write (98,'(10i6)') i,j,e,FilledPhase(i,j,e),TruePhase(i,j,e),RawReads(i,j,:),FilledGenos(i,j),TrueGenos(i,j)
	! 				if (trim(PhaseFile)/="None") then
	! 					if (FilledPhase(i,j,e)==TruePhase(i,j,e)) CheckPhase(i,j,e)='*'
	! 					if (FilledPhase(i,j,e)/=TruePhase(i,j,e)) then
	! 						CheckPhase(i,j,e)='/'
	! 					endif
	! 				endif

	! 				if (trim(PhaseFile)=="None") then	
	! 					CheckPhase(i,j,e)='*'
	! 				endif
					
	! 			endif
	! 		enddo
				
	! 		if (FilledGenos(i,j)/=9) then
	! 			if (trim(GenoFile)/="None") then
	! 				if (FilledGenos(i,j)==TrueGenos(i,j)) CheckGenos(i,j)='*'
	! 				if (FilledGenos(i,j)/=TrueGenos(i,j)) then
	! 					CheckGenos(i,j)='/'
	! 				endif
	! 			endif

	! 			if (trim(GenoFile)=="None") then
	! 				CheckGenos(i,j)='*'
	! 			endif
	! 		endif
				
	! 	enddo
	! enddo


	! do i=1,nInd
	! 	if ((trim(PhaseFile)/="None").and.(trim(GenoFile)/="None")) then
	! 		write (99,'(i20,4i10)') Ped(i,1),count(CheckPhase(i,:,1)=='*'),count(CheckPhase(i,:,1)=='/'),count(CheckPhase(i,:,1)=='_') 
	! 		write (99,'(i20,4i10)') Ped(i,1),count(CheckPhase(i,:,2)=='*'),count(CheckPhase(i,:,2)=='/'),count(CheckPhase(i,:,2)=='_')
	! 		write (99,'(i20,4i10)') Ped(i,1),count(CheckGenos(i,:)=='*'),count(CheckGenos(i,:)=='/'),count(CheckGenos(i,:)=='_')
	! 	endif
	
	! 	!write (98,FmtCha) Ped(i,1),CheckPhase(i,:,1)
	! 	!write (98,FmtCha) Ped(i,1),CheckPhase(i,:,2)

	! enddo

	! close(99)
	! close(98)
end subroutine CurrentCountFilled

!################################################################################################

subroutine InitialiseArrays

	use GlobalPar

	implicit none

	FilledPhase=9
	FilledGenos=9
	!ConsensusFounders=0
end subroutine InitialiseArrays

!################################################################################################

subroutine InternalEdititing
	use GlobalPar

	implicit none

	integer :: i,j
	real(kind=8),allocatable,dimension(:) :: RefAllele
	real(kind=8) :: maf

	allocate(RefAllele(nSnp))
	allocate(GeneProbYesOrNo(nSnp))
	GeneProbYesOrNo=1
	i=0

	open (unit=1,file="AlphaFamSeqEditedSnp.txt",status="unknown")
	print*,"Internal Editing"
	write(*,'(1a10,1f10.3)'),"(1) MAF",EditingParameter
	
	!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(nSnp,RawReads,RefAllele,EditingParameter,GeneProbYesOrNo,i)
	do j=1,nSnp
		RefAllele(j)=sum(RawReads(:,j,1))
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

	integer :: i,DumI,j,stat
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

  	implicit none
  
	integer :: i,PosReads,PosGeno,PosPhase
	integer,allocatable,dimension(:) :: TempImput
	integer :: TmpID

	if (trim(GenoFile)/="None") open (unit=3,file=trim(GenoFile),status="old")
	if (trim(PhaseFile)/="None") open (unit=4,file=trim(PhaseFile),status="old")
	!For CurrentCount
	!open (unit=99,file="AlphaFamSeqSummary.log",status="unknown")


	if (trim(ReadsType)=="VcfTools") call readRogerData(ReadsFile, Ids, position, quality, SequenceData,LenghtSequenceDataFile,nSnp,StartSnp,EndSnp,nIndSeq)
	if (trim(ReadsType)=="AlphaSim") call readAlphaSimReads(ReadsFile, Ids, SequenceData,LenghtSequenceDataFile,nSnp,StartSnp,EndSnp,nIndSeq)
	
	allocate(RawReads(nInd,nSnp,2))
	RawReads(:,:,:)=0.0
	do i=1, nIndSeq
		PosReads=0
		read(Ids(i),*) TmpID
	    call GetID(TmpID, PosReads)
	    !print*,i,TmpID,PosReads
	    if (PosReads/=0) RawReads(PosReads,:,:)= real(SequenceData(i,:,:))
	enddo

	deallocate(SequenceData)
	allocate(TempImput(LenghtSequenceDataFile))
	call AllocateArrays

	do i=1,nInd
		
	    if (trim(GenoFile)/="None") then
		    PosGeno=0
		    read(3,*) Id(i), TempImput(:) 
		    call GetID(Id(i), PosGeno)
		    TrueGenos(PosGeno,:) = TempImput(StartSnp:EndSnp)
		endif

		if (trim(PhaseFile)/="None") then
		 	PosPhase=0
		    read(4,*) Id(i), TempImput(:) 
		    call GetID(Id(i), PosPhase)
		    TruePhase(PosPhase,:,1)= TempImput(StartSnp:EndSnp)

		    PosPhase=0
			read(4,*) Id(i), TempImput(:) 
		    call GetID(Id(i), PosPhase)
		    TruePhase(PosPhase,:,2)= TempImput(StartSnp:EndSnp)

		endif

	enddo

	deallocate(TempImput)
	
	if (trim(GenoFile)/="None") close (3)
	if (trim(PhaseFile)/="None") close (4)
end subroutine ReadData

!###########################################################################################################################################################

subroutine GetID(InputId, PosId)

    use GlobalPar
    implicit none

    integer, intent(in) :: InputId
    integer, intent(out) :: PosId

    integer :: i,check

    PosId = 0
    check = 0

    do i=1, nInd 
        if (Ped(i,1) == InputId) then !Ped is in global module
            PosId = i
        endif
    enddo

    if (PosId == 0) then ! just a check, you don't need this really but may be useful for testing!
    !    print*, "individual not present, please check file"
    !   stop
    check=check+1
    endif
end subroutine GetID

!###########################################################################################################################################################

subroutine Checker

	use GlobalPar

	implicit none

	integer :: i,j,e
	character(len=50) :: filout1

	write (filout1,'("AlphaFamSeqMistakeIdentifiers",i0,".txt")') Windows
	open (unit=1,file=trim(filout1),status="unknown")
	write (1,'(a95)') "i,j,e,FilledPhase(i,j,e),TruePhase(i,j,e),RawReads(i,j,:),FilledGenos(i,j),TrueGenos(i,j)"

	CheckPhase='_'
	CheckGenos='_'

	do i=1,nInd
		do j=1,nSnp
			if (trim(PhaseFile)/="None") then
				do e=1,2
					if (FilledPhase(i,j,e)/=9) then
						if (FilledPhase(i,j,e)==TruePhase(i,j,e)) CheckPhase(i,j,e)='*'
						if (FilledPhase(i,j,e)/=TruePhase(i,j,e)) then
							CheckPhase(i,j,e)='/'
							write (1,'(5i10,2f6.0,2i6)') i,j,e,FilledPhase(i,j,e),TruePhase(i,j,e),RawReads(i,j,:),FilledGenos(i,j),TrueGenos(i,j)
						endif
					endif
				enddo
			endif
		
			if (trim(GenoFile)/="None") then
				if (FilledGenos(i,j)/=9) then
					if (FilledGenos(i,j)==TrueGenos(i,j)) CheckGenos(i,j)='*'
					if (FilledGenos(i,j)/=TrueGenos(i,j)) then
						CheckGenos(i,j)='/'
					endif
				endif
			endif

		enddo
	enddo

	close (1)
end subroutine Checker

!################################################################################################

subroutine WriteResults

	use GlobalPar

	implicit none

	integer :: i,j
	real(kind=4),allocatable,dimension(:) :: AlleleDosage
	character(len=30) :: nChar
	character(len=80) :: FmtInt,FmtCha,FmtReal,FmtIntF,filout1,filout2,filout3,filout4,filout5

	
	write(nChar,*) nSnp
	FmtInt='(i10,'//trim(adjustl(nChar))//'i2)'
	FmtCha='(i10,a2,'//trim(adjustl(nChar))//'a2)'
	FmtReal='(i10,'//trim(adjustl(nChar))//'f7.4)'
	FmtIntF='(i10,'//trim(adjustl(nChar))//'i6)'

	write (filout1,'("AlphaFamSeqFinalPhase",i0,".txt")') Windows
	write (filout2,'("AlphaFamSeqFinalGenos",i0,".txt")') Windows
	write (filout3,'("AlphaFamSeqFinalAlleleDosage",i0,".txt")') Windows
	write (filout5,'("AlphaFamSeqFinalGeneProb",i0,".txt")') Windows
	write (filout4,'("AlphaFamSeqFounderAssignment",i0,".txt")') Windows
	
	open (unit=1,file=trim(filout1),status="unknown")
	open (unit=2,file=trim(filout2),status="unknown")
	open (unit=3,file=trim(filout3),status="unknown")
	open (unit=5,file=trim(filout5),status="unknown")
	open (unit=4,file=trim(filout4),status="unknown")
	
	allocate(AlleleDosage(nSnp))

	do i=1,nInd
		write (1,FmtInt) Ped(i,1),FilledPhase(i,:,1)
		write (1,FmtInt) Ped(i,1),FilledPhase(i,:,2)
		write (2,FmtInt) Ped(i,1),FilledGenos(i,:)

		do j=1,nSnp
			AlleleDosage(j)=(0.0*Pr00(i,j)+(Pr01(i,j)+Pr10(i,j))+2.0*Pr11(i,j))
		enddo

		write (3,FmtReal) Ped(i,1),AlleleDosage(:)
		write (5,FmtReal) Ped(i,1),Pr00(i,:)
		write (5,FmtReal) Ped(i,1),Pr01(i,:)
		write (5,FmtReal) Ped(i,1),Pr10(i,:)
		write (5,FmtReal) Ped(i,1),Pr11(i,:)

		write (4,FmtIntF) Ped(i,1),FounderAssignment(i,:,1)
		write (4,FmtIntF) Ped(i,1),FounderAssignment(i,:,2)

	enddo

	deallocate(AlleleDosage)


	close (1)
	close (2)
	close (3)
	close (4)
end subroutine WriteResults

!################################################################################################

subroutine WriteStatistics
	use GlobalPar
	use AlphaStatMod
	
	implicit none

	REAL(8), PARAMETER :: D_QNAN = &
	TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8) ! NaN value

	integer :: i,j,p,gYield,pYield1,pYield2,nZeroReads,nOneReads,meanCoverage
	integer,allocatable,dimension(:) :: trueVector,impVector
	real :: gCorrect,pCorrect1,pCorrect2,AFtrue
	!real,allocatable,dimension(:) :: meanTG,sdTG,meanFG,sdFG
	character(len=30) :: nChar
	character(len=80) :: FmtInt,FmtCha,FmtIntF,filout1,filout2
	type(CorrelationReal32) :: gCor,pCor1,pCor2
	!type(DescStatReal32) :: stat

	write(nChar,*) nSnp
	FmtInt='(i10,'//trim(adjustl(nChar))//'i2)'
	FmtCha='(i10,a2,'//trim(adjustl(nChar))//'a2)'
	FmtIntF='(i10,'//trim(adjustl(nChar))//'i6)'

	write (filout1,'("AlphaFamSeqStatByIndividual",i0,".txt")') Windows
	write (filout2,'("AlphaFamSeqStatByMarker",i0,".txt")') Windows
	!write (filout3,'("AlphaFamSeqStatSummary",i0,".txt")') Windows

	open (unit=1,file=trim(filout1),status="unknown")
	open (unit=2,file=trim(filout2),status="unknown")
	!open (unit=3,file=trim(filout3),status="unknown")


	write (1,'(a110)') "Id x %Markers0Reads %Markers<2Reads Geno%Yield Geno%Correct GenoCor Phase%Yield Phase%Correct PhaseCor"
	write (2,'(a110)') "Snp MAFtrue x %Ids0Reads %Ids<2Reads Geno%Yield Geno%Correct GenoCor Phase%Yield Phase%Correct PhaseCor"
	!write (3,'(a110)') "What x Geno%Yield Geno%Correct GenoCor Phase%Yield Phase%Correct PhaseCor"


	! allocate(meanTG(nSnp))
	! allocate(sdTG(nSnp))
	! allocate(meanFG(nSnp))
	! allocate(sdFG(nSnp))

	! do j=1,nSnp
	! 	stat=DescStat(TrueGenos(:,j))
	! 	meanTG(j)=stat%Mean
	! 	sdTG(j)=stat%SD

	! 	meanFG=D_QNAN
	! 	sdFG=D_QNAN
	! 	if (count(FilledGenos(:,j)==9)<(nInd-1)) then
	! 		stat=DescStat(RemoveMissing(FilledGenos(:,j), 9))
	! 		meanFG(j)=stat%Mean
	! 		sdFG(j)=stat%SD
	! 	endif
		
	! enddo

	do i=1,nInd

		if (trim(GenoFile)/="None") then ! Check Genotypes
			gYield=nSnp-count(CheckGenos(i,:)=='_')
			gCorrect=dble(count(CheckGenos(i,:)=='*'))/dble(gYield)*100

			allocate(trueVector(gYield))
			allocate(impVector(gYield))

			if (gYield==nSnp) then
				gCor = Cor(TrueGenos(i,:),FilledGenos(i,:))
			else if (gYield<nSnp) then
				p=1
				do j=1,nSnp
					if (FilledGenos(i,j)/=9) then
						trueVector(p)=TrueGenos(i,j)
						impVector(p)=FilledGenos(i,j)
						p=p+1
					endif
				enddo
				
				gCor=Cor(trueVector,impVector)
			endif

			IF (ALLOCATED(trueVector)) deallocate(trueVector)
			IF (ALLOCATED(impVector)) deallocate(impVector)


		endif




		if (trim(PhaseFile)/="None") then ! Check Phase
			pYield1=(nSnp)-(count(CheckPhase(i,:,1)=='_'))
			pYield2=(nSnp)-(count(CheckPhase(i,:,2)=='_'))
			pCorrect1=dble(count(CheckPhase(i,:,1)=='*'))/dble(pYield1)*100
			pCorrect2=dble(count(CheckPhase(i,:,2)=='*'))/dble(pYield2)*100
			
			if (pYield1==nSnp) then
				pCor1 = Cor(TruePhase(i,:,1),FilledPhase(i,:,1))
				
			else if (pYield1<nSnp) then
				
				allocate(trueVector(pYield1))
				allocate(impVector(pYield1))
				p=1
				do j=1,nSnp
					if (FilledPhase(i,j,1)/=9) then
						trueVector(p)=TruePhase(i,j,1)
						impVector(p)=FilledPhase(i,j,1)
						p=p+1
					endif
				enddo
				pCor1=Cor(trueVector,impVector)
				IF (ALLOCATED(trueVector)) deallocate(trueVector)
				IF (ALLOCATED(impVector)) deallocate(impVector)
			endif

			if (pYield2==nSnp) then
				pCor2 = Cor(TruePhase(i,:,2),FilledPhase(i,:,2))
				
			else if (pYield2<nSnp) then
				
				allocate(trueVector(pYield2))
				allocate(impVector(pYield2))
				p=1
				do j=1,nSnp
					if (FilledPhase(i,j,2)/=9) then
						trueVector(p)=TruePhase(i,j,2)
						impVector(p)=FilledPhase(i,j,2)
						p=p+1
					endif
				enddo
				pCor2=Cor(trueVector,impVector)
				IF (ALLOCATED(trueVector)) deallocate(trueVector)
				IF (ALLOCATED(impVector)) deallocate(impVector)
			endif

		endif

		IF (ALLOCATED(trueVector)) deallocate(trueVector)
		IF (ALLOCATED(impVector)) deallocate(impVector)


		nZeroReads=0
		nOneReads=0
		meanCoverage=0
		do j=1,nSnp
			if (sum(RawReads(i,j,:))==0) nZeroReads=nZeroReads+1
			if (sum(RawReads(i,j,:))<2) nOneReads=nOneReads+1 
			meanCoverage=meanCoverage+sum(RawReads(i,j,:))
		enddo
		
		write(1,'(1i10,9f10.5)') Ped(i,1),(dble(meanCoverage)/dble(nSnp)),(dble(nZeroReads)/dble(nSnp)*100),(dble(nOneReads)/dble(nSnp)*100),(dble(gYield)/(dble(nSnp))*100),gCorrect,gCor%Cor,(dble(pYield1+pYield2)/dble(2*nSnp)*100),(pCorrect1+pCorrect2)/2,((pCor1%Cor+pCor2%Cor)/2)
	enddo

	IF (ALLOCATED(trueVector)) deallocate(trueVector)
	IF (ALLOCATED(impVector)) deallocate(impVector)


	do j=1,nSnp

		if (trim(GenoFile)/="None") then ! Check Genotypes
			gYield=nInd-count(CheckGenos(:,j)=='_')
			gCorrect=dble(count(CheckGenos(:,j)=='*'))/dble(gYield)*100
			
			AFtrue=dble(sum(TrueGenos(:,j)))/dble(nInd*2)
			if ((AFtrue)>0.5) AFtrue=1-AFtrue
			gCor%Cor=D_QNAN !IEEE_VALUE(IEEE_QUIET_NAN) 

			if (gYield>1) then
				if (gYield==nInd) then
					gCor = Cor(TrueGenos(:,j),FilledGenos(:,j))
				else if (gYield<nInd) then
					allocate(trueVector(gYield))
					allocate(impVector(gYield))

					p=1
					do i=1,nInd
						if (FilledGenos(i,j)/=9) then
							trueVector(p)=TrueGenos(i,j)
							impVector(p)=FilledGenos(i,j)
							p=p+1
						endif
					enddo

					gCor=Cor(trueVector,impVector)
					IF (ALLOCATED(trueVector)) deallocate(trueVector)
					IF (ALLOCATED(impVector)) deallocate(impVector)
				endif
			endif
		endif


		if (trim(PhaseFile)/="None") then ! Check Phase
			pYield1=(nInd)-(count(CheckPhase(:,j,1)=='_'))
			pYield2=(nInd)-(count(CheckPhase(:,j,2)=='_'))
			pCorrect1=dble(count(CheckPhase(:,j,1)=='*'))/dble(pYield1)*100
			pCorrect2=dble(count(CheckPhase(:,j,2)=='*'))/dble(pYield2)*100
			
			pCor1%Cor=D_QNAN 
			if (pYield1==nInd) then
				pCor1 = Cor(TruePhase(:,j,1),FilledPhase(:,j,1))
				
			else if ((pYield1<nInd).and.(pYield1>1)) then
				
				allocate(trueVector(pYield1))
				allocate(impVector(pYield1))
				p=1
				do i=1,nInd
					if (FilledPhase(i,j,1)/=9) then
						trueVector(p)=TruePhase(i,j,1)
						impVector(p)=FilledPhase(i,j,1)
						p=p+1
					endif
				enddo
				
				pCor1=Cor(trueVector,impVector)
				IF (ALLOCATED(trueVector)) deallocate(trueVector)
				IF (ALLOCATED(impVector)) deallocate(impVector)
			endif

			pCor2%Cor=D_QNAN 
			if (pYield2==nInd) then
				pCor2 = Cor(TruePhase(:,j,2),FilledPhase(:,j,2))
				
			else if ((pYield2<nInd).and.(pYield2>1)) then
				
				allocate(trueVector(pYield2))
				allocate(impVector(pYield2))
				p=1
				do i=1,nInd
					if (FilledPhase(i,j,2)/=9) then
						trueVector(p)=TruePhase(i,j,2)
						impVector(p)=FilledPhase(i,j,2)
						p=p+1
					endif
				enddo
				
				pCor2=Cor(trueVector,impVector)
				IF (ALLOCATED(trueVector)) deallocate(trueVector)
				IF (ALLOCATED(impVector)) deallocate(impVector)
			endif

		endif

		IF (ALLOCATED(trueVector)) deallocate(trueVector)
		IF (ALLOCATED(impVector)) deallocate(impVector)


		nZeroReads=0
		nOneReads=0
		meanCoverage=0
		do i=1,nInd
			if (sum(RawReads(i,j,:))==0) nZeroReads=nZeroReads+1
			if (sum(RawReads(i,j,:))<2) nOneReads=nOneReads+1
			meanCoverage=meanCoverage+sum(RawReads(i,j,:))
		enddo
		
		write(2,'(1i10,10f10.5)') j,AFtrue,(dble(meanCoverage)/dble(nInd)),(dble(nZeroReads)/dble(nInd)*100),(dble(nOneReads)/dble(nInd)*100),(dble(gYield)/(dble(nInd))*100),gCorrect,gCor%Cor,(dble(pYield1+pYield2)/dble(2*nInd)*100),(pCorrect1+pCorrect2)/2,((pCor1%Cor+pCor2%Cor)/2)
	enddo


	close(1)
	close(2)
	!close(3)

	! write (filout5,'("AlphaFamSeqPhaseCorrectness",i0,".txt")') Windows
	! write (filout6,'("AlphaFamSeqSummary",i0,".txt")') Windows
	! write (filout7,'("AlphaFamSeqFinalPhaseCharFormat",i0,".txt")') Windows
	! write (filout8,'("AlphaFamSeqAlleleFreqAndGenoCorrectnessByMarker",i0,".txt")') Windows

	! open (unit=5,file=trim(filout5),status="unknown")
	! open (unit=6,file=trim(filout6),status="unknown")
	! open (unit=7,file=trim(filout7),status="unknown")
	! open (unit=7,file=trim(filout7),status="unknown")
	! open (unit=8,file=trim(filout8),status="unknown")
	
	! write (6,'(a54)') "Id PhaseOrGeno Correct Wrong Missing"
	! write (8,'(a66)') "Snp AlleleFreq Correct Wrong Missing TrueHet cHet wHet nMhomo mHet"





! 	write (6,'(1i20,1a5,3i10)') Ped(i,1),"G ",count(CheckGenos(i,:)=='*'),count(CheckGenos(i,:)=='/'),count(CheckGenos(i,:)=='_')
		




! 		if (trim(PhaseFile)/="None") then
! 			!write (5,FmtCha) Ped(i,1),CheckPhase(i,:,1)
! 			!write (5,FmtCha) Ped(i,1),CheckPhase(i,:,2)
! 			write (6,'(1i20,1a5,3i10)') Ped(i,1),"P1",count(CheckPhase(i,:,1)=='*'),count(CheckPhase(i,:,1)=='/'),count(CheckPhase(i,:,1)=='_')
! 			write (6,'(1i20,1a5,3i10)') Ped(i,1),"P2",count(CheckPhase(i,:,2)=='*'),count(CheckPhase(i,:,2)=='/'),count(CheckPhase(i,:,2)=='_')

! 			CharPhase='_'
! 			do j=1,nSnp
! 				do e=1,2
! 					if (FilledPhase(i,j,e)==1) CharPhase(j,e)='1'
! 					if (FilledPhase(i,j,e)==0) CharPhase(j,e)='0'
! 				enddo
! 			enddo
! 			write (7,FmtCha) Ped(i,1),CharPhase(:,1)
! 			write (7,FmtCha) Ped(i,1),CharPhase(:,2)
! 		endif


! 			if (trim(GenoFile)/="None") then
! 		do j=1,nSnp
! 			Af=0
! 			nC=0
! 			nW=0
! 			nM=0
! 			nThet=0
! 			nChet=0
! 			nWhet=0
! 			nMhet=0
! 			nMhomo=0
! 			Af=real(sum(TrueGenos(:,j)))/real(nInd*2)
! 			nC=count(CheckGenos(:,j)=='*')
! 			nW=count(CheckGenos(:,j)=='/')
! 			nM=count(CheckGenos(:,j)=='_')
! 			nThet=count(TrueGenos(:,j)==1)
! 			do i=1,nInd
! 				if (TrueGenos(i,j)==1) then
! 					if (FilledGenos(i,j)==1)  nChet=nChet+1
! 					if ((FilledGenos(i,j)/=1).and.(FilledGenos(i,j)/=9))  nWhet=nWhet+1
! 					!if (FilledGenos(i,j)==9)  nMhet=nMhet+1
! 				endif

! 				if (FilledGenos(i,j)==9) then
! 					if (TrueGenos(i,j)==1)  nMhet=nMhet+1
! 					if (TrueGenos(i,j)/=1)  nMhomo=nMhomo+1
! 				endif

! 			enddo
! 			write (8,'(1i20,1f10.5,8i10)') j,Af,nC,nW,nM,nThet,nChet,nWhet,nMhomo,nMhet
! 		enddo
! 	endif

! enddo

! 	close (5)
! 	close (6)
! 	close (7)

end subroutine WriteStatistics
