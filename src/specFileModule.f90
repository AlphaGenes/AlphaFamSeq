module specFileModule
    contains
    
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



end module specFileModule