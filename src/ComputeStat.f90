!-----------------------------------------------------------------------------------------------------------------------
! The Roslin Institute, The University of Edinburgh - AlphaGenes Group
!-----------------------------------------------------------------------------------------------------------------------
!
! MODULE: PhaseRounds
! 
!> @file        HaplotypeBits.f90
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
! 2017.02.09  Mara Battagin - Initial Version
!
!---------------------------------------------------------------------------------------------------------------------

module CalculateStatisticsForGenoAndPhase
  use ISO_Fortran_Env
  implicit none
  contains

subroutine GetResultsImputation(nTotSnp,ImpFile,TrueFile,Geno1orPhase2,prefix)
	use AlphaStatMod
  	use omp_lib
	
	implicit none
	
	character(len=*), intent(in):: ImpFile
	character(len=*), intent(in):: TrueFile
	character(len=*), intent(in):: prefix
	
	integer(int32),   intent(in):: nTotSnp
	integer(kind=1),  intent(in):: Geno1orPhase2

	
	integer :: nIndImp,nIndTrue,nInd,nSnp,nChunk,EndSnp,StartSnp
	integer :: i,j,g,k,t,p,gam
	integer(int32), allocatable,dimension(:) :: Id
	integer, allocatable,dimension(:,:) :: Yield,Correct
	integer(kind=1), allocatable,dimension(:,:,:) :: ImpSnp, TrueSnp
	
	real(real64), allocatable,dimension(:,:) :: MAF,FinalCor
	
	character(len=5) :: fileKind
	character(len=300) :: filout1,filout2
	type(CorrelationReal32) :: CorTrueImp

	REAL(8), PARAMETER :: D_QNAN = &
	TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8) ! NaN value

	
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

	nSnp=100000
	gam=Geno1orPhase2

	if (Geno1orPhase2==1) then
		fileKind="Genos"
	else
		fileKind="Phase"
		nInd=nInd/2
		nIndTrue=nIndTrue/2
		nIndImp=nIndImp/2

	endif
	
	allocate(Id(nInd))

	filout1=trim(adjustl(prefix))//"Stat"//trim(adjustl(fileKind))//"ByMarker.txt"
	filout2=trim(adjustl(prefix))//"Stat"//trim(adjustl(fileKind))//"ByIndividual.txt"

	open(101, file=trim(filout1), status="unknown")
	write(101,'(1a47)') "Snp MAF Yield CorrectRate ErrorRate Correlation"

	! Calculate Stat By SnpUsed - Genotypes
	i=0
	do while (EndSnp.lt.nTotSnp)
		i=i+1

		EndSnp=(i*nSnp)
		if (EndSnp>nTotSnp) EndSnp=nTotSnp
		StartSnp=EndSnp-nSnp+1
		
		call AllocateArrays(nInd,nSnp,nSnp,gam,ImpSnp,TrueSnp,Yield,MAF,Correct,FinalCor)

		if (nIndImp.le.nIndTrue) call ReadDataSnp(StartSnp,nSnp,gam,ImpFile,TrueFile,nIndImp,nIndTrue,Id,ImpSnp,TrueSnp)
		if (nIndImp.gt.nIndTrue) call ReadDataSnp(StartSnp,nSnp,gam,TrueFile,ImpFile,nIndTrue,nIndImp,Id,TrueSnp,ImpSnp)
		
		do j=1,nSnp
			do g=1,gam
				Yield(j,g)=0
				MAF(j,g)=0
				Correct(j,g)=0
				FinalCor(j,g)=D_QNAN

				do k=1,nInd
					if ((ImpSnp(k,j,g).ge.0).and.(ImpSnp(k,j,g).le.2)) Yield(j,g)=Yield(j,g)+1
					MAF(j,g)=MAF(j,g)+dble(TrueSnp(k,j,g))
					if (ImpSnp(k,j,g)==TrueSnp(k,j,g)) Correct(j,g)=Correct(j,g)+1
				enddo

				MAF(j,g)=MAF(j,g)/dble(nInd*2)

				if ((Correct(j,g).lt.Yield(j,g)).and.(Yield(j,g).gt.0)) then
					call CalculateCorrelation(Yield(j,g),nInd,TrueSnp(:,j,g),ImpSnp(:,j,g),CorTrueImp)
					FinalCor(j,g)=CorTrueImp%Cor
				endif

				write(101,'(1i0,1x,5f10.4)') StartSnp+j-1,MAF(j,g),100*(dble(Yield(j,g))/dble(nInd)),100*(dble(Correct(j,g))/dble(Yield(j,g))),100*(dble(Yield(j,g)-Correct(j,g))/dble(Yield(j,g))),FinalCor(j,g)
			enddo
		enddo
		
		call Fottiti(ImpSnp,TrueSnp,Yield,MAF,Correct,FinalCor)
	enddo
	close(101)

	
	! Calculate Stat By individual - Genotypes
	open(102, file=trim(filout2), status="unknown")
	write(102,'(1a42)') "Id Yield CorrectRate ErrorRate Correlation"

	call AllocateArrays(1,nTotSnp,nInd,gam,ImpSnp,TrueSnp,Yield,MAF,Correct,FinalCor)
	

	do i=1,nInd
		
		if (nIndImp.le.nIndTrue) call ReadDataInd(i,gam,ImpFile,TrueFile,nIndImp,nIndTrue,Id,ImpSnp,TrueSnp)
		if (nIndImp.gt.nIndTrue) call ReadDataInd(i,gam,TrueFile,ImpFile,nIndTrue,nIndImp,Id,TrueSnp,ImpSnp)


		do g=1,gam
			Yield(i,g)=0
			Correct(i,g)=0
			FinalCor(i,g)=D_QNAN

			do k=1,nTotSnp
				if ((ImpSnp(1,k,g).ge.0).and.(ImpSnp(1,k,g).le.2)) Yield(i,g)=Yield(i,g)+1
				if (ImpSnp(1,k,g)==TrueSnp(1,k,g)) Correct(i,g)=Correct(i,g)+1
			enddo

			if ((Correct(i,g).lt.Yield(i,g)).and.(Yield(i,g).gt.0)) then
				call CalculateCorrelation(Yield(i,g),nTotSnp,TrueSnp(1,:,g),ImpSnp(1,:,g),CorTrueImp)
				FinalCor(i,g)=CorTrueImp%Cor
			endif
				
			write(102,'(1i0,1x,5f10.4)') Id(i),100*(dble(Yield(i,g))/dble(nTotSnp)),100*(dble(Correct(i,g))/dble(Yield(i,g))),100*(dble(Yield(i,g)-Correct(i,g))/dble(Yield(i,g))),FinalCor(i,g)
		enddo
	enddo
	close(102)	

	call Fottiti(ImpSnp,TrueSnp,Yield,MAF,Correct,FinalCor)



end subroutine GetResultsImputation

!###########################################################################################################################################################

subroutine CalculateCorrelation(Yield,n,TrueSnp,ImpSnp,CorTrueImp)
	use AlphaStatMod
  	
	implicit none
	
	
	integer,intent(in) :: n,Yield
	
	integer(kind=1),intent(in),dimension(n) :: ImpSnp, TrueSnp
	type(CorrelationReal32),intent(inout) :: CorTrueImp
	
	integer(int32), allocatable,dimension(:) :: TrueTmp,ImpTmp
	integer :: i,p
	
	if (Yield.gt.1) then
		if (Yield==n) then
			CorTrueImp = Cor(TrueSnp,ImpSnp)
		else if (Yield<n) then
			allocate(TrueTmp(Yield))
			allocate(ImpTmp(Yield))

			p=1
			do i=1,n
				if (ImpSnp(i)/=9) then
					TrueTmp(p)=TrueSnp(i)
					ImpTmp(p)=ImpSnp(i)
					p=p+1
				endif
			enddo

			CorTrueImp=Cor(TrueTmp,ImpTmp)
			if (allocated(TrueTmp)) deallocate(TrueTmp)
			if (allocated(ImpTmp)) deallocate(ImpTmp)
		endif
	endif
end subroutine CalculateCorrelation

!###########################################################################################################################################################

subroutine ReadDataInd(skip,gam,ShortFile,LongFile,nIndShort,nIndLong,Id,ShortSnp,LongSnp)

	implicit none
	
	character(len=*), intent(in):: ShortFile
	character(len=*), intent(in):: LongFile
	
	integer,  intent(in):: skip,gam,nIndShort,nIndLong
	!integer(kind=1),  intent(in):: Geno1orPhase2
	integer(int32),intent(inout), allocatable,dimension(:) :: Id
	
	integer(kind=1),intent(inout),allocatable,dimension(:,:,:) :: ShortSnp, LongSnp
	
	integer :: i,g,j,k,DumI,Pos
	
					
	open(10, file=trim(ShortFile), action="read")
	open(11, file=trim(LongFile), action="read")

	if (skip==1) then
		! do nothing
	else
		do i=1,(skip-1)
			do g=1,gam
				read(10,*)
			enddo
		enddo
	endif

	do g=1,gam
		read(10,*) Id(skip), ShortSnp(1,:,g)
	enddo

	do k=1,nIndLong
		Pos=0
		do g=1,gam
			read(11,*) DumI,LongSnp(1,:,g)
		enddo
		if (DumI==Id(skip)) exit
	enddo
	
	close(10)
	close(11)
end subroutine ReadDataInd

!###########################################################################################################################################################

subroutine ReadDataSnp(StartSnp,nSnp,gam,ShortFile,LongFile,nIndShort,nIndLong,Id,ShortSnp,LongSnp)

	implicit none
	
	character(len=*), intent(in):: ShortFile
	character(len=*), intent(in):: LongFile
	
	integer,  intent(in):: StartSnp,nSnp,nIndShort,nIndLong,gam
	!integer(kind=1),  intent(in):: Geno1orPhase2
	integer(int32),intent(inout), allocatable,dimension(:) :: Id
	
	integer(kind=1),intent(inout),allocatable,dimension(:,:,:) :: ShortSnp, LongSnp
	integer(kind=1), allocatable,dimension(:) :: DumLeft,tmpSnp


	
	integer :: j,g,t,DumI,Pos
	
					
	open(10, file=trim(ShortFile), action="read")
	open(11, file=trim(LongFile), action="read")

	allocate(dumLeft(StartSnp-1))
	allocate(tmpSnp(1:nSnp))
	
	do j=1,nIndShort
		do g=1,gam
			if (StartSnp.eq.1) read(10,*) Id(j), 						  		ShortSnp(j,1:nSnp,g)
			if (StartSnp.gt.1) read(10,*) Id(j), (DumLeft(t),t=1,(StartSnp-1)), ShortSnp(j,1:nSnp,g)
		enddo		
	enddo
	close(10)

	do j=1,nIndLong
		Pos=0
		do g=1,gam
			if (StartSnp.eq.1) read(11,*) DumI,tmpSnp(1:nSnp)
			if (StartSnp.gt.1) read(11,*) DumI, (DumLeft(t),t=1,(StartSnp-1)), tmpSnp(1:nSnp)
			if (g==1) call GetID(Id,nIndShort,DumI,Pos)
			if (Pos/=0) LongSnp(Pos,1:nSnp,g) = tmpSnp(1:nSnp)
		enddo
		
	enddo
	close(11)

	if(allocated(DumLeft)) 	deallocate(dumLeft)
	if(allocated(tmpSnp)) 	deallocate(tmpSnp)
end subroutine ReadDataSnp

!###########################################################################################################################################################

subroutine GetID(Id,nInd,InputId, PosId)

    implicit none

    integer(int32), intent(in),dimension(:) :: Id
    integer, intent(in) :: InputId,nInd
    integer, intent(out) :: PosId

    integer :: i,check

    PosId = 0
    check = 0

    do i=1, nInd 
        if (Id(i) == InputId) then !Ped is in global module
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

subroutine Fottiti(ImpSnp,TrueSnp,Yield,MAF,Correct,FinalCor)
	
	implicit none
	
	integer,intent(inout),allocatable,dimension(:,:) :: Yield,Correct
	integer(kind=1),intent(inout),allocatable,dimension(:,:,:) :: ImpSnp, TrueSnp
	
	real(real64),intent(inout),allocatable,dimension(:,:) :: MAF,FinalCor
	

	if(allocated(ImpSnp)) deallocate(ImpSnp)
	if(allocated(TrueSnp)) deallocate(TrueSnp)
	if(allocated(Yield)) deallocate(Yield)
	if(allocated(MAF)) deallocate(MAF)
	if(allocated(Correct)) deallocate(Correct)
	if(allocated(FinalCor)) deallocate(FinalCor)
	
end subroutine Fottiti

!###########################################################################################################################################################

subroutine AllocateArrays(nRow,nCol,LengthVector,gam,ImpSnp,TrueSnp,Yield,MAF,Correct,FinalCor)
	implicit none
	
	integer,intent(in) :: nRow,nCol,LengthVector,gam
	integer,intent(inout),allocatable,dimension(:,:) :: Yield,Correct
	integer(kind=1),intent(inout),allocatable,dimension(:,:,:) :: ImpSnp, TrueSnp
	
	real(real64),intent(inout),allocatable,dimension(:,:) :: MAF,FinalCor
	

	allocate(ImpSnp(nRow,nCol,gam))
	allocate(TrueSnp(nRow,nCol,gam))
	allocate(Yield(LengthVector,gam))
	allocate(MAF(LengthVector,gam)) ! GoOut
	allocate(Correct(LengthVector,gam))
	allocate(FinalCor(LengthVector,gam))

	
end subroutine AllocateArrays


end module CalculateStatisticsForGenoAndPhase

!###########################################################################################################################################################

! program test

! 	use CalculateStatisticsForGenoAndPhase
! 	implicit none
	
! 	character(len=100) :: A
! 	character(len=100) :: B
! 	character(len=11) :: prefix
! 	integer(int32)  :: nSnp
! 	integer(kind=1)  :: Geno1orPhase2
	
	
! 	A="AlphaFamSeqFinalGenos.txt"
! 	B="Geno.txt"
! 	prefix="AlphaFamSeq"
! 	nSnp=700000
! 	Geno1orPhase2=1

! 	call GetResultsImputation(nSnp,A,B,Geno1orPhase2,prefix)

! 	A="AlphaFamSeqFinalPhase.txt"
! 	B="Phase.txt"
! 	Geno1orPhase2=2


! 	call GetResultsImputation(nSnp,A,B,Geno1orPhase2,prefix)
	
! end program test

!###########################################################################################################################################################









