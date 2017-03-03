!-----------------------------------------------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-----------------------------------------------------------------------------------------------------------------------
!
! MODULE: PhaseRounds
! 
!> @file        ComputeStat.f90
!
! DESCRIPTION:
!> @brief       Module to calculate accuracy of imputation for imputed Genotypes/Phase 
!>
!> @details     This module contains subtoutines to calculate:
!               Yield(%) = Markers imputed on the total markers 
!				CorrectRate(%) = Markers well imputed on the total markers imputed
!				ErrorRate(%) = Markers wrongly imputed on the total markers imputed
!				Correlation by individual = correlation between true and imputed markers per each individual
!				Correlation by Markers = correlation between true and imputed markers per each marker
! 					
!
!> @author      Mara Battagin, mara.battagin@roslin.ed.ac.uk
!
!> @date        February 09, 2017
!
!> @version     0.0.1 (alpha)
!
! REVISION HISTORY:
! 2017.02.13  Mara Battagin - Calculate Yield removing markers with missing information (in the raw data)
! 2017.02.09  Mara Battagin - Initial Version
!
!---------------------------------------------------------------------------------------------------------------------

module CalculateStatisticsForGenoAndPhase
  use ISO_Fortran_Env
  implicit none
  contains

subroutine GetResultsImputation(nSnp,ImpFile,TrueFile,ExclueSnpFile,Geno1orPhase2,MistakeIdentifier,prefix)
	use AlphaStatMod
  	
	implicit none
	
	character(len=*), intent(in):: ImpFile
	character(len=*), intent(in):: TrueFile
	character(len=*), intent(in):: ExclueSnpFile
	character(len=3), intent(in) :: MistakeIdentifier
    
	character(len=*), intent(in):: prefix
	
	integer(int32),   intent(in):: nSnp
	integer(kind=1),  intent(in):: Geno1orPhase2

	
	integer :: nIndImp,nIndTrue,nInd
	integer :: nSnpUsed
	integer :: gam
	
	integer(int64), allocatable,dimension(:) :: Id
	integer(int32), allocatable,dimension(:) :: MarkerToExclude
	integer, allocatable,dimension(:,:) :: Yield,Correct
	integer(int32), allocatable,dimension(:,:,:) :: ImpSnp, TrueSnp
	
	real(real32), allocatable,dimension(:,:) :: MAF
	real(real32), allocatable,dimension(:,:) :: FinalCor

	character(len=5) :: fileKind
	character(len=300) :: filout1,filout2,filout3

	call GetLengthOfInputFiles(ImpFile,nIndImp,TrueFile,nIndTrue,nInd)
	call SetSomeParametes(gam,Geno1orPhase2,fileKind,prefix,nInd,nIndTrue,nIndImp,filout1,filout2,filout3)
	call AllocateInputDataArrays(nInd,nSnp,gam,Id,ImpSnp,TrueSnp)
	if (nIndImp.le.nIndTrue) call ReadDataIn(gam,nSnp,ImpFile,TrueFile,nIndImp,nIndTrue,Id,ImpSnp,TrueSnp)
	if (nIndImp.gt.nIndTrue) call ReadDataIn(gam,nSnp,TrueFile,ImpFile,nIndTrue,nIndImp,Id,TrueSnp,ImpSnp)

	call ReadMarkersToExclude(MarkerToExclude,ExclueSnpFile,nSnpUsed,nSnp,TrueSnp)
	! CalculareResultsBySnp
	call AllocateResultsArrays(nSnp,gam,Yield,Correct,FinalCor,MAF)
	call CalculateResultsBySnp	(nSnp,gam,nInd,ImpSnp,TrueSnp,MAF,Yield,Correct,FinalCor)
	call WriteResultsBySnp(nSnp,nInd,gam,MarkerToExclude,MAF,Yield,Correct,FinalCor,filout1)
	call DeallocateResultsArrays(Yield,MAF,Correct,FinalCor)
	! CalculateResultsByIndividual
	call AllocateResultsArrays(nInd,gam,Yield,Correct,FinalCor,MAF)
	call CalculateResultsByIndividual(nSnp,gam,nInd,ImpSnp,TrueSnp,MarkerToExclude,Yield,Correct,FinalCor)
	call WriteResultsByIndividual(nInd,gam,nSnpUsed,Id,Yield,Correct,FinalCor,filout2)	
	call DeallocateResultsArrays(Yield,MAF,Correct,FinalCor)
	! Print Mistakes
	print*,trim(MistakeIdentifier)
	if (trim(MistakeIdentifier)=="Yes") then
		call PrintMistakeIdentifier(nInd,nSnp,gam,MarkerToExclude,ImpSnp,TrueSnp,filout3)
	endif

	call DeallocateInputDataArrays(Id,ImpSnp,TrueSnp)

end subroutine GetResultsImputation

!###########################################################################################################################################################

subroutine PrintMistakeIdentifier(nInd,nSnp,gam,MarkerToExclude,ImpSnp,TrueSnp,filout3)

	implicit none

	integer, 		intent(in) :: nInd,gam
	integer(int32), intent(in) :: nSnp

	integer(int32),intent(in),allocatable,dimension(:,:,:) :: ImpSnp, TrueSnp
	integer(int32),	intent(in), allocatable,dimension(:) :: MarkerToExclude

	character(len=300),	intent(in) :: filout3

	integer i,g,j

	open(103, file=trim(filout3), status="unknown")
	write(103,'(1a42)') "Id gam Snp True Imputed"


	do i=1,nInd
		do g=1,gam
			do j=1,nSnp
				if ((ImpSnp(i,j,g)/=9).and.(TrueSnp(i,j,g)/=ImpSnp(i,j,g))) then
					write(103,'(5(1x,i0))') i,g,j,TrueSnp(i,j,g),ImpSnp(i,j,g)
				endif
			enddo
		enddo
	enddo

	close(103)
end subroutine PrintMistakeIdentifier

!###########################################################################################################################################################

subroutine WriteResultsByIndividual(nInd,gam,nSnpUsed,Id,Yield,Correct,FinalCor,filout2)

	implicit none

	integer,			intent(in) :: nInd,gam,nSnpUsed
	integer(int64),		intent(in), allocatable,dimension(:) :: Id
	integer,			intent(in),allocatable,dimension(:,:) :: Yield,Correct
	real(real32),		intent(in),allocatable,dimension(:,:) :: FinalCor

	character(len=300),	intent(in) :: filout2

	integer :: i,g

	open(102, file=trim(filout2), status="unknown")
	write(102,'(1a42)') "Id Yield CorrectRate ErrorRate Correlation"


	do i=1,nInd
		do g=1,gam
			write(102,'(1i0,1x,5f10.4)') Id(i), 100*(dble(Yield(i,g))/dble(nSnpUsed)), 100*(dble(Correct(i,g))/dble(Yield(i,g))), 100*(dble(Yield(i,g)-Correct(i,g))/dble(Yield(i,g))), FinalCor(i,g)		
		enddo
	enddo
	
	close(102)
end subroutine WriteResultsByIndividual

!###########################################################################################################################################################

subroutine CalculateResultsByIndividual(nSnp,gam,nInd,ImpSnp,TrueSnp,MarkerToExclude,Yield,Correct,FinalCor)

	use AlphaStatMod
	implicit none


	integer(int32),intent(in) :: nSnp
	integer,intent(in) :: gam,nInd
	integer(int32),intent(in),allocatable,dimension(:,:,:) :: ImpSnp, TrueSnp
	integer(int32),intent(in), allocatable,dimension(:) :: MarkerToExclude
	
	integer,intent(inout),allocatable,dimension(:,:) :: Yield,Correct
	real(real32),intent(inout),allocatable,dimension(:,:) :: FinalCor


	integer :: i,g,j
	type(CorrelationReal32) :: CorTrueImp
	REAL(8), PARAMETER :: D_QNAN = &
	TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8) ! NaN value


	do i=1,nInd
		do g=1,gam
			Yield(i,g)=0
			Correct(i,g)=0
			FinalCor(i,g)=D_QNAN

			do j=1,nSnp
				if (MarkerToExclude(j)==0) then
					if ((ImpSnp(i,j,g)/=9).and.(TrueSnp(i,j,g)/=9)) then
						Yield(i,g)=Yield(i,g)+1
						if (ImpSnp(i,j,g)==TrueSnp(i,j,g)) Correct(i,g)=Correct(i,g)+1
					endif
				endif
			enddo

			if (Yield(i,g).gt.1) then
				call CalculateCorrelation(Yield(i,g),nSnp,TrueSnp(i,:,g),ImpSnp(i,:,g),CorTrueImp)
				FinalCor(i,g)=CorTrueImp%Cor
			endif
			
		enddo
	enddo
end subroutine CalculateResultsByIndividual

!###########################################################################################################################################################

subroutine WriteResultsBySnp(nSnp,nInd,gam,MarkerToExclude,MAF,Yield,Correct,FinalCor,filout1)

	implicit none

	integer,			intent(in) :: nSnp,nInd,gam
	integer(int32),		intent(in), allocatable,dimension(:) :: MarkerToExclude
	
	integer,			intent(in),allocatable,dimension(:,:) :: Yield,Correct
	real(real32),		intent(in),allocatable,dimension(:,:) :: MAF
	real(real32), allocatable,dimension(:,:) :: FinalCor

	
	character(len=300),	intent(in) :: filout1

	integer :: j,g

	open(101, file=trim(filout1), status="unknown")
	write(101,'(1a47)') "Snp MAF Yield CorrectRate ErrorRate Correlation"


	do j=1,nSnp
		do g=1,gam
			if (MarkerToExclude(j)==0) then
				write(101,'(1i0,1x,5f10.4)') j, MAF(j,g), 100*(dble(Yield(j,g))/dble(nInd)), 100*(dble(Correct(j,g))/dble(Yield(j,g))), 100*(dble(Yield(j,g)-Correct(j,g))/dble(Yield(j,g))), FinalCor(j,g)
			endif
		enddo
	enddo

	close(101)
end subroutine WriteResultsBySnp

!###########################################################################################################################################################

subroutine CalculateResultsBySnp(nSnp,gam,nInd,ImpSnp,TrueSnp,MAF,Yield,Correct,FinalCor)

	use AlphaStatMod
	implicit none

	integer,intent(in) :: nSnp,gam,nInd
	integer(int32),intent(in),allocatable,dimension(:,:,:) :: ImpSnp, TrueSnp
	
	integer,intent(inout),allocatable,dimension(:,:) :: Yield,Correct
	real(real32),intent(inout),allocatable,dimension(:,:) :: MAF
	real(real32),intent(inout),allocatable,dimension(:,:) :: FinalCor


	integer :: j,g,i
	type(CorrelationReal32) :: CorTrueImp
	REAL(8), PARAMETER :: D_QNAN = &
	TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8) ! NaN value

	do j=1,nSnp
		do g=1,gam
			Yield(j,g)=0
			MAF(j,g)=0
			Correct(j,g)=0
			FinalCor(j,g)=D_QNAN

			do i=1,nInd
				if ((ImpSnp(i,j,g)/=9).and.(TrueSnp(i,j,g)/=9)) then
					Yield(j,g)=Yield(j,g)+1
					MAF(j,g)=MAF(j,g)+dble(TrueSnp(i,j,g))
					if (ImpSnp(i,j,g)==TrueSnp(i,j,g)) Correct(j,g)=Correct(j,g)+1
				endif
			enddo

			MAF(j,g)=MAF(j,g)/dble(nInd*2)

			if (Yield(j,g).gt.1) then
				call CalculateCorrelation(Yield(j,g),nInd,TrueSnp(:,j,g),ImpSnp(:,j,g),CorTrueImp)
				FinalCor(j,g)=CorTrueImp%Cor
			endif

		enddo
	enddo
end subroutine CalculateResultsBySnp

!###########################################################################################################################################################

subroutine ReadMarkersToExclude(MarkerToExclude,ExclueSnpFile,nSnpUsed,nSnp,TrueSnp)

	implicit none

	integer(int32),		intent(inout), allocatable,dimension(:) :: MarkerToExclude
	integer,			intent(inout) :: nSnpUsed
	integer(int32),   	intent(in):: nSnp
	character(len=*), 	intent(in):: ExclueSnpFile
	integer(int32),		intent(in),allocatable,dimension(:,:,:) :: TrueSnp

	integer :: DumI,i,a

	allocate(MarkerToExclude(nSnp))
	MarkerToExclude=0
	
	open(12,file=trim(ExclueSnpFile),action="read")
	nSnpUsed=0
	do
		read(12,*,end=912) DumI
		MarkerToExclude(DumI)=1
	enddo
	912 continue
	close(12)

	a=0
	do i=1,nSnp
		if (minval(TrueSnp(:,i,:))==9) then
			MarkerToExclude(i)=1
			a=a+1
		endif
	enddo
	print*,a

	nSnpUsed=count(MarkerToExclude(:)==0)
	print*,nSnpUsed
end subroutine ReadMarkersToExclude

!###########################################################################################################################################################

subroutine ReadDataIn(gam,nSnp,ShortFile,LongFile,nIndShort,nIndLong,Id,ShortSnp,LongSnp)

	implicit none
	
	character(len=*), 	intent(in):: ShortFile
	character(len=*), 	intent(in):: LongFile
	integer(int32),		intent(in) :: nSnp
	integer,  			intent(in):: gam,nIndShort,nIndLong
	
	integer(int64),		intent(inout), allocatable,dimension(:) :: Id
	integer(int32),	intent(inout),allocatable,dimension(:,:,:) :: ShortSnp, LongSnp
	
	integer :: i,g
	integer(int64) :: DumI,PosId
	integer,allocatable,dimension (:) :: TmpInput
					
	open(10, file=trim(ShortFile), action="read")
	open(11, file=trim(LongFile), action="read")

	do i=1,nIndShort
		do g=1,gam
			read(10,*) Id(i), ShortSnp(i,:,g)
		enddo
	enddo

	allocate(TmpInput(nSnp))

	do i=1,nIndLong
		PosId=0
		do g=1,gam
			read(11,*) DumI,TmpInput
			call GetID(Id,DumI,PosId)
			if (PosId/=0) LongSnp(PosId,:,g)=TmpInput(:)
		enddo
	enddo
	
	deallocate(TmpInput)

	close(10)
	close(11)
end subroutine ReadDataIn

!###########################################################################################################################################################

subroutine SetSomeParametes(gam,Geno1orPhase2,fileKind,prefix,nInd,nIndTrue,nIndImp,filout1,filout2,filout3)
	implicit none

	integer(kind=1),intent(in) :: Geno1orPhase2
	integer,intent(inout) :: gam
	integer,intent(inout) :: nIndImp,nIndTrue,nInd

	character(len=5),intent(inout) :: fileKind
	character(len=*), intent(in):: prefix
	
	character(len=300),intent(inout) :: filout1,filout2,filout3

	
	gam=Geno1orPhase2

	if (Geno1orPhase2==1) then
		fileKind="Genos"
	else
		fileKind="Phase"
		nInd=nInd/2
		nIndTrue=nIndTrue/2
		nIndImp=nIndImp/2
	endif

	filout1=trim(adjustl(prefix))//"Stat"//trim(adjustl(fileKind))//"ByMarker.txt"
	filout2=trim(adjustl(prefix))//"Stat"//trim(adjustl(fileKind))//"ByIndividual.txt"
	filout3=trim(adjustl(prefix))//"MistakeIdentifiers"//trim(adjustl(fileKind))//".txt"
end subroutine SetSomeParametes

!###########################################################################################################################################################

subroutine GetLengthOfInputFiles(ImpFile,nIndImp,TrueFile,nIndTrue,nInd)
	implicit none

	character(len=*), intent(in):: ImpFile
	character(len=*), intent(in):: TrueFile
	
	integer,intent(inout) :: nIndImp,nIndTrue,nInd

	! Count number of individuals Imputed file
	open(10, file=trim(ImpFile), action="read")
	nIndImp=0
	do 
		read(10,*,end=910)
		nIndImp=nIndImp+1
	enddo
	910 continue
	close(10)

	! Count number of individuals True file
	open(11, file=trim(TrueFile), action="read")
	nIndTrue=0
	do 
		read(11,*,end=911)
		nIndTrue=nIndTrue+1
	enddo
	911 continue
	close(11)

	! start to do stuff
	if (nIndImp.le.nIndTrue) nInd=nIndImp
	if (nIndImp.gt.nIndTrue) nInd=nIndTrue
end subroutine GetLengthOfInputFiles

!###########################################################################################################################################################

subroutine CalculateCorrelation(Yield,n,TrueSnp,ImpSnp,CorTrueImp)
	use AlphaStatMod
  	
	implicit none
	
	
	integer,intent(in) :: n,Yield
	
	integer(int32),intent(in),dimension(n) :: ImpSnp, TrueSnp
	type(CorrelationReal32),intent(inout) :: CorTrueImp
	
	integer(int32), allocatable,dimension(:) :: TrueTmp,ImpTmp
	integer :: i,p
	
	if (Yield==n) then
		CorTrueImp = Cor(TrueSnp,ImpSnp)
	endif
	
	if (Yield.lt.n) then
		allocate(TrueTmp(Yield))
		allocate(ImpTmp(Yield))
		p=1
		do i=1,n
			if ((ImpSnp(i)/=9).and.(TrueSnp(i)/=9)) then
				TrueTmp(p)=TrueSnp(i)
				ImpTmp(p)=ImpSnp(i)
				p=p+1
			endif
		enddo
		
		CorTrueImp = Cor(TrueTmp,ImpTmp)
		if (allocated(TrueTmp)) deallocate(TrueTmp)
		if (allocated(ImpTmp)) deallocate(ImpTmp)
	endif
end subroutine CalculateCorrelation

!###########################################################################################################################################################

subroutine GetID(Id,InputId,PosId)

    implicit none

    integer(int64), intent(in),dimension(:) :: Id
    integer(int64), intent(in) :: InputId
    integer(int64), intent(out) :: PosId

    integer :: nInd,i

    PosId = 0
    nInd=size(Id)

    do i=1, nInd 
        if (Id(i) == InputId) then 
            PosId = i
        endif
    enddo
end subroutine GetID

!###########################################################################################################################################################

subroutine DeallocateResultsArrays(Yield,MAF,Correct,FinalCor)
	
	implicit none
	
	integer,		intent(inout),allocatable,dimension(:,:) :: Yield,Correct
	real(real32),	intent(inout),allocatable,dimension(:,:) :: MAF,FinalCor
	

	if(allocated(Yield)) deallocate(Yield)
	if(allocated(MAF)) deallocate(MAF)
	if(allocated(Correct)) deallocate(Correct)
	if(allocated(FinalCor)) deallocate(FinalCor)
end subroutine DeallocateResultsArrays

!###########################################################################################################################################################

subroutine AllocateResultsArrays(nRow,gam,Yield,Correct,FinalCor,MAF)
	implicit none
	
	integer,intent(in) :: nRow,gam
	integer,intent(inout),allocatable,dimension(:,:) :: Yield,Correct
	
	real(real32),intent(inout),allocatable,dimension(:,:) :: FinalCor
	real(real32),intent(inout),optional,allocatable,dimension(:,:) :: MAF

	allocate(Yield(nRow,gam))
	allocate(MAF(nRow,gam)) ! GoOut
	allocate(Correct(nRow,gam))
	allocate(FinalCor(nRow,gam))
end subroutine AllocateResultsArrays

!###########################################################################################################################################################

subroutine DeallocateInputDataArrays(Id,ImpSnp,TrueSnp)
	implicit none
	
	integer(int64),intent(inout), allocatable,dimension(:) :: Id
	integer(int32),intent(inout),allocatable,dimension(:,:,:) :: ImpSnp, TrueSnp

	deallocate(Id)
	deallocate(ImpSnp)
	deallocate(TrueSnp)
end subroutine DeallocateInputDataArrays

!###########################################################################################################################################################

subroutine AllocateInputDataArrays(nRow,nCol,gam,Id,ImpSnp,TrueSnp)
	implicit none
	
	integer,intent(in) :: nRow,nCol,gam
	integer(int64),intent(inout), allocatable,dimension(:) :: Id
	integer(int32),intent(inout),allocatable,dimension(:,:,:) :: ImpSnp, TrueSnp

	allocate(Id(nRow))
	allocate(ImpSnp(nRow,nCol,gam))
	allocate(TrueSnp(nRow,nCol,gam))
end subroutine AllocateInputDataArrays

!###########################################################################################################################################################

end module CalculateStatisticsForGenoAndPhase

!###########################################################################################################################################################










