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
    integer,intent(inout) :: StartPos,EndPos                                              ! SpecFile - Input SequenceFile Option - first and last position

    real,intent(inout)          :: maxStdForReadsCount,ThresholdMaxReadsCount              ! SpecFile - Editing Parametes - Remove Reads that are above this standard deviation
    integer,intent(inout)       :: ThresholdReadsCount                                     ! SpecFile - Editing Parametes - Remove single/double/n-tones 
    real,intent(inout)          :: ThresholdExcessHetero                                   ! SpecFile - Editing Parametes - Remove variants with an excess of heterozygotes
    
    integer(kind=1),intent(inout)   :: UsePrevGeneProb                                         ! SpecFile - Input SingleLocusPeeler - Read old results of GeneProb 1==YES, 0==NO
    real,intent(inout)              :: GeneProbThresh                                          ! SpecFile - Input SingleLocusPeeler - Threshold to call a genotype from the probabilities First Value
    real,intent(inout)              :: GeneProbThreshMin                                       ! SpecFile - Input SingleLocusPeeler - Threshold to call a genotype from the probabilities Last Value
    real,intent(inout)              :: ReduceThr                                               ! SpecFile - Input SingleLocusPeeler - Reduce Geno Treshold factor
    
    integer,intent(inout)      :: minWindowSizeHapDefinition                              ! SpecFile - Input Build Consensu Haplotype - First value to define Haplotypes length
    integer,intent(inout)      :: maxWindowSizeHapDefinition                              ! SpecFile - Input Build Consensu Haplotype - Last value to define Haplotypes length

    character(len=300),intent(inout) :: GenoFile                                          ! SpecFile - Control Results File Name - TrueGeno Genotypes to check results 
    character(len=300),intent(inout) :: PhaseFile                                         ! SpecFile - Control Results File Name - True Phase to check results 



    integer :: stat, i
    character(len=90) :: SpecParam
   	character (len=512) :: TLC


    ! PedigreeFile="Pedigree.txt"
    ! SequenceFile="SequeceData.txt"
    ! SequenceFileFormat="ASim"
    ! SequenceDataType="RC"
    ! SnpChipFile="None"
    ! MapSnpChipFile="None"
    ! UsePrevGeneProb=0
    ! GeneProbThresh=1.0
    ! GeneProbThreshMin=0.98
    ! ReduceThr=0.1
    ! minWindowSizeHapDefinition=1
    ! maxWindowSizeHapDefinition=1
    ! GenoFile="None"
    ! PhaseFile="None"

    open(unit=1, file="AlphaFamSeqSpec.txt", status="old")

    do i=1, countLines("AlphaFamSeqSpec.txt")
        read(1,'(a30,A)', advance='NO', iostat=stat) SpecParam 
        
        ! if (SpecParam(1:1)=="=" .or. len(trim(SpecParam))==0) then
        !     print*,SpecParam(1:9)
        !     cycle
        ! else
        
            select case(trim(TLC(SpecParam)))

                case('pedigreefile')
                    read(1, *, iostat=stat) PedigreeFile
                    if (stat /= 0) then
                        print*, "ERROR - PedigreeFile not set properly in spec file"
                        print*,"default",PedigreeFile
                        stop 8
                    endif
 
                case('sequencedatafile')
                    read(1, *, iostat=stat) SequenceFile,SequenceFileFormat,SequenceDataType
                    if (stat /= 0) then
                        print*, "ERROR - SequenceDataFile not set properly in spec file"
                        stop 10
                    endif

                    if ((trim(SequenceFileFormat)/="ASim").and.(trim(SequenceFileFormat)/="GATK").and.(trim(SequenceFileFormat)/="None")) then
                        print*,"ERROR - Sequence file format not recognized"
                        print*,"        Valid formats are: ASim or GATK"
                        print*,"        Please check the parameters file"
                        print*,"        ",trim(SequenceFileFormat)
                        stop
                    endif

                    if ((trim(SequenceDataType)/="RC").and.(trim(SequenceDataType)/="None")) then
                        print*,"ERROR - Sequence data type not recognized or not supported"
                        print*,"        Supported type is: RC"
                        print*,"        Please check the parameters file"
                        print*,"        ",trim(SequenceDataType)
                        stop
                    endif

                case('genotypedatafileandmap')
                    read(1, *, iostat=stat) SnpChipFile,MapSnpChipFile
                    if (stat /= 0) then
                        print*, "ERROR - GenotypeDataFileAndMap not set properly in spec file"
                        stop 10
                    endif

    			case('numberofsnps')
                    read(1, *, iostat=stat) nSnp
                    if (stat /= 0) then
                        print*, "ERROR - NumberOfSnps not set properly in spec file"
                        stop 2
                    endif 

                case('chromosomeandinterval')
                    read(1, *, iostat=stat) chr, StartPos,EndPos
                    if (stat /= 0) then
                        print*, "ERROR - ChromosomeAndInterval not set properly in spec file"
                        stop 2
                    endif 

                    if (StartPos.gt.EndPos) then
                        print*,"ERROR - Start POS or SNP greater than Last POS or SNP"
                        print*,"        Please check the parameters file"
                        print*,"        StartPos:",StartPos
                        print*,"        EndPos  :",EndPos
                        stop
                    endif
    			
                case('removeoutliersreadscount')
                    read(1, *, iostat=stat) maxStdForReadsCount,ThresholdMaxReadsCount
                    if (stat /= 0) then
                        print*, "ERROR - RemoveOutliersReadsCount not set properly in spec file"
                        stop 2
                    endif 

    			case('removemarkerslownrreads')
                    read(1, *, iostat=stat) ThresholdReadsCount
                    if (stat /= 0) then
                        print*, "ERROR - RemoveMarkersLowNrReads not set properly in spec file"
                        stop 2
                    endif 

    			case('removeexcessheteropvalue')
                    read(1, *, iostat=stat) ThresholdExcessHetero
                    if (stat /= 0) then
                        print*, "ERROR - RemoveExcessHeteroPvalue not set properly in spec file"
                        stop 2
                    endif 

                case('runsinglelocuspeeler')
                    read(1, *, iostat=stat) UsePrevGeneProb !nIter 
                    if (stat /= 0) then
                        print*, "ERROR - RunSingleLocusPeeler not set properly in spec file"
                        stop 8
                    endif

                    if ((UsePrevGeneProb.ne.0).and.(UsePrevGeneProb.ne.1)) then
                        print*,"ERROR - RunSingleLocusPeeler not recognized"
                        print*,"        Use 0 to run the Single Locus Peeler or 1 to use previous results"
                        print*,"        Please check the parameters file:"
                        print*,"        RunSingleLocusPeeler:",UsePrevGeneProb
                        stop
                    endif
     

                case('allelethresholdprobability')
                    read(1, *, iostat=stat) GeneProbThresh,GeneProbThreshMin,ReduceThr
                    if (stat /= 0) then
                        print*, "ERROR - AlleleThresholdProbability not set properly in spec file"
                        stop 8
                    endif   

                   if (GeneProbThresh.gt.1) then
                        print*,"WARNING - First allele threshold greater than 1.00"
                        print*,"          Using default value of 1.00 instead of ",GeneProbThresh
                        GeneProbThresh=1.00
                    endif
 
                   if (GeneProbThreshMin.le.0.5) then
                        print*,"WARNING - Second allele threshold lower than 0.5"
                        print*,"          Using default value of 0.51 instead of ",GeneProbThreshMin
                        GeneProbThreshMin=0.51
                    endif

                   if (GeneProbThresh-GeneProbThreshMin.lt.ReduceThr) then
                        print*,"WARNING - Decrement of the allele probability threshold is greater than the difference of the two thresholds"
                        print*,"          Using value of ",GeneProbThresh-GeneProbThreshMin," instead of ",GeneProbThreshMin
                        ReduceThr=GeneProbThresh-GeneProbThreshMin
                    endif


                case('windowsizehaplotype')
                    read(1, *, iostat=stat) minWindowSizeHapDefinition,maxWindowSizeHapDefinition
                    if (stat /= 0) then
                        print*, "ERROR - WindowSizeHaplotype not set properly in spec file"
                        stop 8
                    endif

                    if (minWindowSizeHapDefinition.gt.nSnp) then
                        print*,"WARNING: min window size > nSnp"
                        print*,"         use 1 as min window size"
                        minWindowSizeHapDefinition=1
                    endif

                    if (maxWindowSizeHapDefinition.gt.nSnp) then
                        print*,"WARNING: max window size > nSnp"
                        print*,"         use nSnp as max window size"
                        maxWindowSizeHapDefinition=nSnp
                    endif

                    if (minWindowSizeHapDefinition.gt.maxWindowSizeHapDefinition) then
                        print*,"WARNING: min window size > max window size"
                        print*,"         use 1 as min window size"
                        minWindowSizeHapDefinition=1
                    endif
   

    			case('truegenofile')
                    read(1, *, iostat=stat) GenoFile
                    if (stat /= 0) then
                        print*, "ERROR - TrueGenoFile not set properly in spec file"
                        stop 10 
                    endif   
                    
    			case('truephasefile')
                    read(1, *, iostat=stat) PhaseFile
                    if (stat /= 0) then
                        print*, "ERROR - TruePhaseFile not set properly in spec file"
                        stop 10
                    endif
                    
    			case default
                    print*, "ERROR - Error in specfile, please check", SpecParam
                    stop 16
        	end select
        !endif
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
