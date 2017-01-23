!-----------------------------------------------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-----------------------------------------------------------------------------------------------------------------------
!
! MODULE: PhaseRounds
! 
!> @file        HaplotypeBits.f90
!
! DESCRIPTION:
!> @brief       Module to bit-wise operations
!>
!> @details     This MODULE contains a single routine to read in the data that Roger produces for Mara's program.   Input is the
!filename, the total number of Snps, and the total number of Individuals.    Output is the Ids (from the first line), the position
!of the snp, the snp Quality and the sequence reads.   
!filename: input.   Name of file to be read. character, dimension(1)
!Ids:: output.   Id of the Snps. character, dimension(nSnp, 100) (i.e. character length 100)
!position: output. Position of individuals? Row 2. integer, dimension (nIndiv)
!quality: output. The quality of the reads. Row 5. real(real64), dimension(nIndiv)
!rawReads: The sequence reads. integer(int32). dimension(nIndiv, nSnp, 2)
!nSnpIn: input. number of Snps
!nIndivIn: input, number of Individuals
!
!> @author      Diarmaid de Burca, diarmaid.deburca@.ed.ac.uk
!
!> @date        September 23, 2016
!
!> @version     0.0.1 (alpha)
!
! REVISION HISTORY:
! 2016.09.23  Diarmaid de Burca - Initial Version
!
!---------------------------------------------------------------------------------------------------------------------
module MaraModule
  use ISO_Fortran_Env
  implicit none
  contains
subroutine readRogerData(filename, Ids, position, quality, SequenceData,nSnpIn,SnpUsed,StartSnp,EndSnp,nIndivIn)
  use omp_lib
  implicit none
  !filename has the name of the file to be read
  character(len=*), intent(in)::filename
  character(len=100), allocatable, dimension(:), intent(out):: Ids
  !info holds the quality and the poisition (

  character(len=100), dimension(:), allocatable::dumE, dumC
  real(real64), allocatable, dimension(:), intent(out):: quality
  integer(int32), dimension(:), allocatable, intent(out):: position
  integer(int32), dimension(:,:,:), allocatable, intent(out):: SequenceData

  integer(int32), intent(in):: SnpUsed,nSnpIn,StartSnp,EndSnp,nIndivIn
  integer(int32)::nSnp,pos,fileUnit, nIndiv
  integer(int32)::i,j,k


  ! Open the file
  open(newunit=fileUnit, file=filename, action="read")

  nSnp = nSnpIn
  nIndiv = nIndivIn
  
  allocate(SequenceData(nIndiv, SnpUsed, 2))
  allocate(position(SnpUsed))
  allocate(quality(SnpUsed))
  allocate(Ids(nIndiv))
  allocate(dumE(5+2*nIndiv))
  allocate(dumC(5+nIndiv))

  read(fileUnit, *) dumC 
  write(*,"(5A)") "STUFF", trim(dumC(1)), trim(dumC(2)), trim(dumC(3)), "ENDSTUFF"
  do i =1, nIndiv
    write(Ids(i), *) dumC(i+5)
  end do

  pos=1
  do j = 1, nSnp
    read(fileUnit, *) dumE
    if ((j.ge.StartSnp).and.(j.le.EndSnp)) then
      read(dumE(2), *) position(pos)
      read(dumE(5), *) quality(pos)
      !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED (nIndiv,j,pos,dumE,SequenceData,position)
      do i = 1, nIndiv
        k = i*2+4
        read(dumE(k), *) SequenceData(i,pos, 1)
        read(dumE(k+1), *) SequenceData(i, pos, 2)
      end do
      !$OMP END PARALLEL DO
      pos=pos+1
    endif
  end do

  end subroutine readRogerData

  subroutine readAlphaSimReads(filename, Ids,SequenceData,nSnpIn,SnpUsed,StartSnp,EndSnp,nIndivIn)
  use ISO_Fortran_Env
  implicit none
  
  character(len=*), intent(in)::filename
  character(len=100), allocatable, dimension(:), intent(out):: Ids
  
  character(len=100), dimension(:), allocatable::dumE
  integer(int32), dimension(:,:,:), allocatable::SequenceData

  integer(int32), intent(in)::SnpUsed,nSnpIn,StartSnp,EndSnp,nIndivIn
  integer(int32)::nSnp,pos,fileUnit, nIndiv
  integer(int32)::i,j, k,dumI

  
  ! Open the file
  open(newunit=fileUnit, file=filename, action="read")

  nSnp = nSnpIn
  nIndiv = nIndivIn
 
  allocate(dumE(nSnp))
  allocate(SequenceData(nIndiv, SnpUsed, 2))
  allocate(Ids(nIndiv))
  
  do i = 1, nIndiv
    read(fileUnit, *) Ids(i),dumE
    read(dumE(StartSnp:EndSnp), *) SequenceData(i,:, 1)
    read(fileUnit, *) dumI,dumE
    read(dumE(StartSnp:EndSnp), *) SequenceData(i,:, 2)
  end do

end subroutine readAlphaSimReads


end module MaraModule



  !So we want to read in and inver


