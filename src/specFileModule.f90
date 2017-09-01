module specFileModule
    contains
    
    ! Reads in and initialises specfile parameters 
    subroutine ReadSpecfile(LenghtSequenceDataFile,nSnp,fistWindow,maxStdForReadsCount, &
                            ThresholdMaxReadsCount,ThresholdReadsCount,ThresholdExcessHetero, &
                            GeneProbThresh,GeneProbThreshMin,ReduceThr,UsePrevGeneProb, &
                            minWindowSizeHapDefinition,maxWindowSizeHapDefinition, &
                            PedigreeFile,ReadsFile,ReadsType,GenoFile,SnpChipsInformation,PhaseFile)

    use AlphaHouseMod, only : countLines

    implicit none

    integer,intent(inout) :: LenghtSequenceDataFile,nSnp,fistWindow                       ! SpecFile - Total number of Snps
    
    real(kind=8),intent(inout) :: maxStdForReadsCount,ThresholdMaxReadsCount              ! SpecFile - Remove Reads that are above this standard deviation
    real(kind=8),intent(inout) :: ThresholdExcessHetero                                   ! SpecFile - Remove variants with an excess of heterozygotes
    integer,intent(inout)      :: ThresholdReadsCount                                     ! SpecFile - Remove single/double/n-tones 
    
    real(kind=8),intent(inout) :: GeneProbThresh                                          ! SpecFile - Threshold to call a genotype from the probabilities First Value
    real(kind=8),intent(inout) :: GeneProbThreshMin                                       ! SpecFile - Threshold to call a genotype from the probabilities Last Value
    real(kind=8),intent(inout) :: ReduceThr                                               ! SpecFile - Reduce Geno Treshold factor
    integer,intent(inout)      :: UsePrevGeneProb                                         ! SpecFile - Read old results of GeneProb 1==YES, 0==NO
    
    integer,intent(inout)      :: minWindowSizeHapDefinition                              ! SpecFile - First value to define Haplotypes length
    integer,intent(inout)      :: maxWindowSizeHapDefinition                              ! SpecFile - Last value to define Haplotypes length

    character(len=300),intent(inout) :: PedigreeFile                                      ! SpecFile - Input File Name - Pedigree
    character(len=300),intent(inout) :: ReadsFile                                         ! SpecFile - Input File Name - Reads Count for the Reference and the Alternative Allele
    character(len=300),intent(inout) :: ReadsType                                         ! SpecFile - Input File Name - Reads Count for the Reference and the Alternative Allele
    
    character(len=300),intent(inout) :: SnpChipsInformation                               ! SpecFile - Input File Name - Snp array to add more information to the Reads
    
    character(len=300),intent(inout) :: GenoFile                                          ! SpecFile - Control Results File Name - TrueGeno Genotypes to check results 
    character(len=300),intent(inout) :: PhaseFile                                         ! SpecFile - Control Results File Name - True Phase to check results 



    integer :: stat, i
    character(len=30) :: SpecParam
   	character (len=512) :: TLC

    open(unit=1, file="AlphaFamSeqSpec.txt", status="old")


    do i=1, countLines("AlphaFamSeqSpec.txt")
        read(1,'(a30,A)', advance='NO', iostat=stat) SpecParam 
        
        select case(trim(TLC(SpecParam)))

			case('numberofsnps')
                read(1, *, iostat=stat) LenghtSequenceDataFile,nSnp,fistWindow
                if (stat /= 0) then
                    print *, "NumberOfSnps not set properly in spec file"
                    print *, LenghtSequenceDataFile,nSnp
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

            case('rangechunklenght')
                read(1, *, iostat=stat) minWindowSizeHapDefinition,maxWindowSizeHapDefinition
                if (stat /= 0) then
                    print *, "RangeChunkLenght not set properly in spec file"
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



end module specFileModule

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
