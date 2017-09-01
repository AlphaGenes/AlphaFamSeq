module specFileModule
    contains
    
    ! Reads in and initialises specfile parameters 
    subroutine ReadSpecfile(PedigreeFile, &
                            SequenceFile,SequenceFileFormat,SequenceDataType, &
                            SnpChipFile,MapSnpChipFile, &
                            nSnp,chr, StartPos,EndPos, & 
                            maxStdForReadsCount,ThresholdMaxReadsCount, &
                            ThresholdReadsCount,ThresholdExcessHetero, &
                            UsePrevGeneProb,GeneProbThresh,GeneProbThreshMin,ReduceThr, &
                            minWindowSizeHapDefinition,maxWindowSizeHapDefinition, &
                            GenoFile,PhaseFile)
        
    use AlphaHouseMod, only : countLines

    implicit none

    character(len=300),intent(inout) :: PedigreeFile                                      ! SpecFile - Input File Name - Pedigree
    character(len=300),intent(inout) :: SequenceFile                                      ! SpecFile - Input File Name - Sequence Data
    character(len=300),intent(inout) :: SequenceFileFormat                                ! SpecFile - Input SequenceFile Option - AlphaSim or VcfTools format
    character(len=300),intent(inout) :: SequenceDataType                                  ! SpecFile - Input SequenceFile Option - RC = read counts, GL= genotype likelihood, GT= genotype 

    character(len=300),intent(inout) :: SnpChipFile                                       ! SpecFile - Input File Name - Snp chip array to add more information to the Reads
    character(len=300),intent(inout) :: MapSnpChipFile                                    ! SpecFile - Input File Name - Map file for Snp chip array
    
    integer,intent(inout)            :: nSnp                                              ! SpecFile - Input - Total number of Snps

    character(len=300),intent(inout) :: chr                                               ! SpecFile - Input SequenceFile Option - chromosome ID
    character(len=300),intent(inout) :: StartPos,EndPos                                   ! SpecFile - Input SequenceFile Option - first and last position

    real(kind=8),intent(inout) :: maxStdForReadsCount,ThresholdMaxReadsCount              ! SpecFile - Editing Parametes - Remove Reads that are above this standard deviation
    real(kind=8),intent(inout) :: ThresholdExcessHetero                                   ! SpecFile - Editing Parametes - Remove variants with an excess of heterozygotes
    integer,intent(inout)      :: ThresholdReadsCount                                     ! SpecFile - Editing Parametes - Remove single/double/n-tones 
    
    integer,intent(inout)      :: UsePrevGeneProb                                         ! SpecFile - Input SingleLocusPeeler - Read old results of GeneProb 1==YES, 0==NO
    real(kind=8),intent(inout) :: GeneProbThresh                                          ! SpecFile - Input SingleLocusPeeler - Threshold to call a genotype from the probabilities First Value
    real(kind=8),intent(inout) :: GeneProbThreshMin                                       ! SpecFile - Input SingleLocusPeeler - Threshold to call a genotype from the probabilities Last Value
    real(kind=8),intent(inout) :: ReduceThr                                               ! SpecFile - Input SingleLocusPeeler - Reduce Geno Treshold factor
    
    integer,intent(inout)      :: minWindowSizeHapDefinition                              ! SpecFile - Input Build Consensu Haplotype - First value to define Haplotypes length
    integer,intent(inout)      :: maxWindowSizeHapDefinition                              ! SpecFile - Input Build Consensu Haplotype - Last value to define Haplotypes length

    character(len=300),intent(inout) :: GenoFile                                          ! SpecFile - Control Results File Name - TrueGeno Genotypes to check results 
    character(len=300),intent(inout) :: PhaseFile                                         ! SpecFile - Control Results File Name - True Phase to check results 



    integer :: stat, i
    character(len=30) :: SpecParam
   	character (len=512) :: TLC

    open(unit=1, file="AlphaFamSeqSpec.txt", status="old")


    do i=1, countLines("AlphaFamSeqSpec.txt")
        read(1,'(a30,A)', advance='NO', iostat=stat) SpecParam 
        
        if (SpecParam(1:1)=="=" .or. len(trim(SpecParam))==0) then
            cycle
        else
        
            select case(trim(TLC(SpecParam)))

                case('pedigreefile')
                    read(1, *, iostat=stat) PedigreeFile
                    if (stat /= 0) then
                        print *, "PedigreeFile not set properly in spec file"
                        stop 8
                    endif   

                case('sequencedatafile')
                    read(1, *, iostat=stat) SequenceFile,SequenceFileFormat,SequenceDataType
                    if (stat /= 0) then
                        print *, "SequenceDataFile not set properly in spec file"
                        stop 10
                    endif   

                case('genotypedatafileandmap')
                    read(1, *, iostat=stat) SnpChipFile,MapSnpChipFile
                    if (stat /= 0) then
                        print *, "GenotypeDataFileAndMap not set properly in spec file"
                        stop 10
                    endif

    			case('numberofsnps')
                    read(1, *, iostat=stat) nSnp
                    if (stat /= 0) then
                        print *, "NumberOfSnps not set properly in spec file"
                        stop 2
                    endif 

                case('chromosomeandinterval')
                    read(1, *, iostat=stat) chr, StartPos,EndPos
                    if (stat /= 0) then
                        print *, "ChromosomeAndInterval not set properly in spec file"
                        stop 2
                    endif 
    			
                case('removeoutliersreadscount')
                    read(1, *, iostat=stat) maxStdForReadsCount,ThresholdMaxReadsCount
                    if (stat /= 0) then
                        print *, "RemoveOutliersReadsCount not set properly in spec file"
                        stop 2
                    endif 

    			case('removemarkerslownrreads')
                    read(1, *, iostat=stat) ThresholdReadsCount
                    if (stat /= 0) then
                        print *, "RemoveMarkersLowNrReads not set properly in spec file"
                        stop 2
                    endif 

    			case('removeexcessheteropvalue')
                    read(1, *, iostat=stat) ThresholdExcessHetero
                    if (stat /= 0) then
                        print *, "RemoveExcessHeteroPvalue not set properly in spec file"
                        stop 2
                    endif 

                case('runsinglelocuspeeler')
                    read(1, *, iostat=stat) UsePrevGeneProb !nIter 
                    if (stat /= 0) then
                        print *, "RunSingleLocusPeeler not set properly in spec file"
                        stop 8
                    endif   

                case('allelethresholdprobability')
                    read(1, *, iostat=stat) GeneProbThresh,GeneProbThreshMin,ReduceThr
                    if (stat /= 0) then
                        print *, "AlleleThresholdProbability not set properly in spec file"
                        stop 8
                    endif   

                case('WindowSizeHaplotype')
                    read(1, *, iostat=stat) minWindowSizeHapDefinition,maxWindowSizeHapDefinition
                    if (stat /= 0) then
                        print *, "windowsizehaplotype not set properly in spec file"
                        stop 8
                    endif   

    			case('truegenofile')
                    read(1, *, iostat=stat) GenoFile
                    if (stat /= 0) then
                        print *, "TrueGenoFile not set properly in spec file"
                        stop 10 
                    endif   
                    
    			case('truephasefile')
                    read(1, *, iostat=stat) PhaseFile
                    if (stat /= 0) then
                        print *, "TruePhaseFile not set properly in spec file"
                        stop 10
                    endif
                    
    			case default
                   print *, "Error in specfile, please check", SpecParam
                    stop 16
        	end select
        endif
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
