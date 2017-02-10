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

subroutine GetResultsImputation(nTotSnp,ImpFile,TrueFile,Geno1orPhase2)
	use AlphaStatMod
  	use omp_lib
	
	implicit none
	
	character(len=*), intent(in):: ImpFile
	character(len=*), intent(in):: TrueFile
	
	integer(int32),   intent(in):: nTotSnp
	integer(kind=1),  intent(in):: Geno1orPhase2

	
	integer :: nIndImp,nIndTrue,nInd,nSnp,nChunk,EndSnp,StartSnp
	integer :: i,j,k,t,p,DumI,Pos
	integer(int32), allocatable,dimension(:) :: Id,v1,v2
	integer(int32), allocatable,dimension(:,:) :: TrueTmp,ImpTmp
	integer, allocatable,dimension(:) :: Yield,Correct
	integer(kind=1), allocatable,dimension(:,:) :: ImpSnp, TrueSnp
	integer(kind=1), allocatable,dimension(:) :: dumLeft,tmpSnp

	real(real64), allocatable,dimension(:) :: MAF,FinalCor

	type(CorrelationReal32) :: CorTrueImp
	REAL(8), PARAMETER :: D_QNAN = &
	TRANSFER((/ Z'00000000', Z'7FF80000' /),1.0_8) ! NaN value

	
	! Count number of individuals Imputed file
	open(10, file=ImpFile, action="read")
	nIndImp=0
	do 
		read(10,*,end=910)
		nIndImp=nIndImp+1
	enddo
	910 continue
	rewind(10)

	! Count number of individuals True file
	open(11, file=TrueFile, action="read")
	nIndTrue=0
	do 
		read(11,*,end=911)
		nIndTrue=nIndTrue+1
	enddo
	911 continue
	rewind(11)


	! start to do stuff
	if (nIndImp.le.nIndTrue) nInd=nIndImp
	if (nIndImp.gt.nIndTrue) nInd=nIndTrue

	nSnp=10000
	
	allocate(Id(nInd))
	
	i=0
	do while (EndSnp.lt.nTotSnp)
		i=i+1

		EndSnp=(i*nSnp)
		if (EndSnp>nTotSnp) EndSnp=nTotSnp
		StartSnp=EndSnp-nSnp+1
		
		allocate(ImpSnp(nInd,EndSnp-StartSnp+1))
		allocate(TrueSnp(nInd,EndSnp-StartSnp+1))
		allocate(tmpSnp(EndSnp-StartSnp+1))
		allocate(Yield(EndSnp-StartSnp+1))
		allocate(MAF(EndSnp-StartSnp+1))
		allocate(Correct(EndSnp-StartSnp+1))
		allocate(FinalCor(EndSnp-StartSnp+1))


		!print*,i,StartSnp,EndSnp
		
		! Calculate Stat By SnpUsed - Genotypes

		open(101, file="AlphaFamSeqStatGenoByMarker.txt", status="unknown")
		write(101,'(1a47)') "Snp MAF Yield CorrectRate ErrorRate Correlation"

		!StartSnp+j-1,MAF(j),100*(dble(Yield(j))/dble(nInd)),100*(dble(Correct(j))/dble(Yield(j))),100*(dble(Yield(j)-Correct(j))/dble(Yield(j))),FinalCor(j)
		allocate(dumLeft(StartSnp-1))
			
		if (nIndImp.le.nIndTrue) then
			
			do j=1,nIndImp
				if (StartSnp.eq.1) read(10,*) Id(j), 						  		 ImpSnp(j,1:nSnp)
				if (StartSnp.gt.1) read(10,*) Id(j), (DumLeft(t),t=1,(StartSnp-1)), ImpSnp(j,1:nSnp)
				
			enddo
			rewind(10)

			do j=1,nIndTrue
				Pos=0
				if (StartSnp.eq.1) read(11,*) DumI,tmpSnp(1:nSnp)
	 		    if (StartSnp.gt.1) read(11,*) DumI, (DumLeft(t),t=1,(StartSnp-1)), tmpSnp(1:nSnp)
	 		    call GetID(Id,nIndImp,DumI,Pos)
	 		    if (Pos/=0) TrueSnp(Pos,1:nSnp) = tmpSnp(1:nSnp)
				
			enddo
			rewind(11)
		endif
		
		if(allocated(dumLeft)) 	deallocate(dumLeft)


		do j=1,nSnp

			Yield(j)=0
			MAF(j)=0
			Correct(j)=0
			FinalCor(j)=D_QNAN

			do k=1,nInd
				if ((ImpSnp(k,j).ge.0).and.(ImpSnp(k,j).le.2)) Yield(j)=Yield(j)+1
				MAF(j)=MAF(j)+dble(TrueSnp(k,j))
				if (ImpSnp(k,j)==TrueSnp(k,j)) Correct(j)=Correct(j)+1
			enddo

			MAF(j)=MAF(j)/dble(nInd*2)

			

			if ((Correct(j).lt.Yield(j)).and.(Yield(j).gt.0)) then
				call CalculateCorrelation(Yield(j),nInd,TrueSnp(:,j),ImpSnp(:,j),CorTrueImp)
				FinalCor(j)=CorTrueImp%Cor
			endif

			write(101,'(1i0,1x,5f10.4)') StartSnp+j-1,MAF(j),100*(dble(Yield(j))/dble(nInd)),100*(dble(Correct(j))/dble(Yield(j))),100*(dble(Yield(j)-Correct(j))/dble(Yield(j))),FinalCor(j)
		enddo
		
		deallocate(ImpSnp)
		deallocate(TrueSnp)
		deallocate(tmpSnp)
		deallocate(Yield)
		deallocate(MAF)
		deallocate(Correct)
		deallocate(FinalCor)
		
	enddo
	close(101)

	!deallocate(Id)
	
	! Calculate Stat By individual - Genotypes
	open(102, file="AlphaFamSeqStatGenoByIndividual.txt", status="unknown")
	write(102,'(1a42)') "Id Yield CorrectRate ErrorRate Correlation"


	!allocate(Id(nInd))
	allocate(ImpSnp(1,nTotSnp))
	allocate(TrueSnp(1,nTotSnp))
	allocate(tmpSnp(nTotSnp))
	allocate(Yield(nInd))
	allocate(MAF(nInd))
	allocate(Correct(nInd))
	allocate(FinalCor(nInd))

	if (nIndImp.le.nIndTrue) then

	
		do i=1,nInd
			Yield(i)=0
			Correct(i)=0
			FinalCor(i)=D_QNAN
			
			read(10,*) Id(i), ImpSnp(1,:)
			do k=1,nIndTrue
				Pos=0
				read(11,*) DumI,TrueSnp(1,:)
	 		    if (DumI==Id(i)) exit
	 		enddo
	 		rewind(11)

			do j=1,nTotSnp
				if ((ImpSnp(1,j).ge.0).and.(ImpSnp(1,j).le.2)) Yield(i)=Yield(i)+1
				if (ImpSnp(1,j)==TrueSnp(1,j)) Correct(i)=Correct(i)+1
			enddo


			call CalculateCorrelation(Yield(i),nTotSnp,TrueSnp(1,:),ImpSnp(1,:),CorTrueImp)
			FinalCor(i)=CorTrueImp%Cor
			write(102,'(1i0,1x,5f10.4)') Id(i),100*(dble(Yield(i))/dble(nTotSnp)),100*(dble(Correct(i))/dble(Yield(i))),100*(dble(Yield(i)-Correct(i))/dble(Yield(i))),FinalCor(i)

		enddo
		close(10)
		close(11)
		close(102)	

	endif


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


end module CalculateStatisticsForGenoAndPhase

!###########################################################################################################################################################

! program test

! 	use CalculateStatisticsForGenoAndPhase
! 	implicit none
	
! 	character(len=25) :: A
! 	character(len=8) :: B
! 	integer(int32)  :: nSnp
! 	integer(kind=1)  :: Geno1orPhase2
	
	
! 	A="AlphaFamSeqFinalGenos.txt"
! 	B="Geno.txt"
! 	nSnp=700000
! 	Geno1orPhase2=1

! 	call GetResultsImputation(nSnp,A,B,Geno1orPhase2)
	
! end program test

!###########################################################################################################################################################









