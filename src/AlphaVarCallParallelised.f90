module AlphaVarCallFuture
use globalGP
use ISO_Fortran_Env
implicit none

public :: AlphaVarCall

contains

    !######################################################################################################################################################

    subroutine AlphaVarCall(nAnis,nSnp,StartSnp,EndSnp,ErrorRate,Seq0Snp1Mode,ped,OutProb)

      use globalGP, only :pedigree
      use ISO_Fortran_Env
      use omp_lib

      implicit none

      integer, intent(in) :: nAnis,nSnp,StartSnp,EndSnp,Seq0Snp1Mode
      real(kind=8),intent(in) :: ErrorRate
      type(PedigreeHolder) ,target:: ped
      
      real(kind=real64),intent(inout),dimension(:,:,:)  :: OutProb(4,EndSnp-StartSnp+1,nAnis)

      integer(int64),dimension(:) :: SeqSire(nAnis),SeqDam(nAnis)
      integer(kind=2),dimension(:,:,:) :: ReadCountsTmp(1:nAnis,nSnp,2) 
      integer(kind=1),dimension(:,:) :: InputGenosTmp(1:nAnis,nSnp) 
      
      integer :: MaxFs,MaxMates,MaxReadCounts

      integer(kind=2),allocatable,dimension(:,:,:) :: ReadCounts                                         
      
      integer(kind=1),allocatable,dimension(:,:) :: InputGenos                                          

      !real(kind=8),dimension(pgg:) :: OutputMaf(EndSnp-StartSnp+1)

      real(kind=8),allocatable,dimension(:,:) :: GMatSnp 
      real(kind=8),allocatable,dimension(:,:,:) :: GMatRds

      integer :: mxeq
      INTEGER,allocatable,dimension(:) :: mate,prog,next,ifirst
      
      real(kind=8)::tstart,tend
      integer :: i

      do i=1,ped%pedigreeSize-ped%nDummys
        if ((ped%pedigree(i)%Founder)) then
          SeqSire(i)=0
          SeqDam(i)=0
        else if (.not. ped%pedigree(i)%Founder) then
          SeqSire(i)=ped%pedigree(i)%sirePointer%id
          if (SeqSire(i).gt.i) SeqSire(i)=0
          SeqDam(i)=ped%pedigree(i)%damPointer%id
          if (SeqDam(i).gt.i) SeqDam(i)=0
        endif
        print*,i,SeqSire(i),SeqDam(i)
      enddo
      ReadCountsTmp=ped%convertsequencedatatoarray()
      InputGenosTmp=ped%getgenotypesasarray()

      call GetMaxFamilySize(nAnis,SeqSire,SeqDam,MaxFs,MaxMates)
      call SetUpData(Seq0Snp1Mode,ReadCounts,InputGenos,nAnis,nSnp,EndSnp,StartSnp,ReadCountsTmp,InputGenosTmp,MaxReadCounts)
      call SetUpEquationsForSnp(Seq0Snp1Mode,GMatSnp,GMatRds,ErrorRate,MaxReadCounts)
      call CreateLinkListArrays(nAnis,SeqSire,SeqDam,mxeq,mate,ifirst,next,prog)
      
      tstart = omp_get_wtime()

      !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) PRIVATE(i) SHARED(ReadCounts,InputGenos,MaxReadCounts,GMatSnp,GMatRds,OutProb)
      do i=1,(EndSnp-StartSnp+1)
        call geneprob(i,nAnis,Seq0Snp1Mode,ReadCounts,InputGenos, &
                      maxfs,MaxMates,MaxReadCounts,GMatSnp,GMatRds,SeqSire,SeqDam, &
                      OutProb, &
                      mxeq,mate,ifirst,next,prog)
      enddo
      !$OMP END PARALLEL DO

      tend = omp_get_wtime()
      write(*,*) "Total wall time for GeneProbController is ", tend - tstart


      if (Seq0Snp1Mode==0) then
        deallocate(GMatRds)
      endif

      if (Seq0Snp1Mode==1) then
        deallocate(GMatSnp)
      endif

    deallocate(MATE,NEXT,IFIRST,prog)

 	IF( ALLOCATED(ReadCounts)) DEALLOCATE(ReadCounts) 
 	IF( ALLOCATED(InputGenos)) DEALLOCATE(InputGenos)
      
    end subroutine AlphaVarCall

    !######################################################################################################################################################
      ! BUILD A LINKLIST. EACH PARENT ANIMAL HAS A ROW, ALONG THE COLUMNS ARE
      ! OF THEIR MATES AND PROGENY. MATES CAN BE REPEATED AT SUCCESSIVE NODES
      ! REFLECT FULL SIB FAMILIES. NOTE, THE FIRST NODE FOR A PARTICULAR MATE
      ! CONTAIN THE POST(i,j) term, (the jth mate of the ith animal)
 
    subroutine CreateLinkListArrays(nAnis,SeqSire,SeqDam,mxeq,mate,ifirst,next,prog)
      use ISO_Fortran_Env
      implicit none

      integer,intent(in) :: nAnis
      integer(int64), intent(in), dimension (:) :: SeqSire(nAnis),SeqDam(nAnis)

      integer,intent(inout) :: mxeq

      INTEGER, intent(inout),allocatable,dimension(:) :: mate,prog,next,ifirst
   

      integer :: mm,nn,is,idd,ia

      mm=0
      nn=nAnis

      ALLOCATE(mate(0:2*nAnis), &
               next(2*nAnis), &
               ifirst(0:nAnis), &
               prog(0:2*nAnis))
    

      ! Create LinkList
      mate=0
      next=0
      ifirst=0
      prog=0 
      MXEQ=0 ! Added by MBattagin - it was not initialised!!
      
      do ia=1,nAnis
        is=seqsire(ia)
        idd=seqdam(ia)
        call LNKLST(is,idd,ia,1,MXEQ,nAnis,mm,nn,mate,next,ifirst,prog)  ! from MXEQ added by MBattagin
        call LNKLST(idd,is,ia,0,MXEQ,nAnis,mm,nn,mate,next,ifirst,prog)  ! from MXEQ added by MBattagin       
      end do
    end subroutine CreateLinkListArrays

    !######################################################################################################################################################

    subroutine SetUpEquationsForSnp(Seq0Snp1Mode,GMatSnp,GMatRds,ErrorRate,MaxReadCounts)
        use ISO_Fortran_Env
        implicit none


        integer,intent(in) :: Seq0Snp1Mode
        integer,intent(inout) :: MaxReadCounts
        real(kind=8),intent(in) :: ErrorRate


        real(kind=8),intent(inout),allocatable,dimension(:,:) :: GMatSnp 
        real(kind=8),intent(inout),allocatable,dimension(:,:,:) :: GMatRds


        integer :: i,k

        real(kind=8) :: ProdFactTmp !DFactorial,DFactorialInLog,
        

        if (Seq0Snp1Mode==0) then

            allocate(GMatRds(0:MaxReadCounts,3,MaxReadCounts))      
            do k=1,MaxReadCounts
                do i=0,k

                    ProdFactTmp=DFactorialInLog(k)-(DFactorialInLog(i)+DFactorialInLog(k-i)) 
                    GMatRds(i,1,k)=ProdFactTmp+(dfloat(i)*log(ErrorRate))+(dfloat(k-i)*log(1.0-ErrorRate))
                    GMatRds(i,2,k)=ProdFactTmp+(dfloat(k)*log(0.5))
                    GMatRds(i,3,k)=ProdFactTmp+(dfloat(i)*log(1.0-ErrorRate))+(dfloat(k-i)*log(ErrorRate))
                enddo
            enddo

        endif

        if (Seq0Snp1Mode==1) then
            print *, "allocating"

            allocate(GMatSnp(1:3,1:3)) ! MBattagin the first and second subscript were "0:2" instead of "1:3" - TODO: make sure it doesn't create problems

            GMatSnp(1,1)=1.0-(((2*ErrorRate)*(1.0-ErrorRate))-(ErrorRate*ErrorRate))
            GMatSnp(1,2)=(2*ErrorRate)*(1.0-ErrorRate)
            GMatSnp(1,3)=ErrorRate*ErrorRate

            GMatSnp(2,1)=ErrorRate*(1.0-ErrorRate)
            GMatSnp(2,2)=1.0-((2*ErrorRate)*(1.0-ErrorRate))
            GMatSnp(2,3)=ErrorRate*(1.0-ErrorRate)

            GMatSnp(3,1)=ErrorRate*ErrorRate
            GMatSnp(3,2)=(2*ErrorRate)*(1.0-ErrorRate)
            GMatSnp(3,3)=1.0-((2*ErrorRate)*(1.0-ErrorRate))

        endif
    end subroutine SetUpEquationsForSnp

    !######################################################################################################################################################

    subroutine SetUpData(Seq0Snp1Mode,ReadCounts,InputGenos,nAnis,nSnp,EndSnp,StartSnp,ReadCountsTmp,InputGenosTmp,MaxReadCounts)
        use ISO_Fortran_Env
        implicit none

        integer, intent(in) :: nAnis,nSnp,StartSnp,EndSnp,Seq0Snp1Mode
        integer, intent(inout) :: MaxReadCounts
        
        integer(kind=2),intent(in),dimension(:,:,:) :: ReadCountsTmp(:,:,:) !(nAnis,nSnp,2)
        integer(kind=1),intent(inout),dimension(:,:) :: InputGenosTmp(:,:)

        integer(kind=2),intent(inout),allocatable,dimension(:,:,:) :: ReadCounts
        integer(kind=1),intent(inout),allocatable,dimension(:,:) :: InputGenos

        integer :: i,j,k

        if (Seq0Snp1Mode==0) then
        	MaxReadCounts=0
            allocate(ReadCounts(nAnis,EndSnp-StartSnp+1,2)) ! Commented by MBattagin

            do i=1,nAnis
                k=1
                do j=StartSnp,EndSnp
                    ReadCounts(i,k,:)=ReadCountsTmp(i,j,:)
                    if (sum(ReadCounts(i,k,:)).gt.MaxReadCounts) MaxReadCounts=sum(ReadCounts(i,k,:))
                    k=k+1
                end do
            end do
        endif



        if (Seq0Snp1Mode==1) then

            allocate(InputGenos(nAnis,EndSnp-StartSnp+1))


            do i=1,nAnis
                do j=1,nSnp
                    if (InputGenosTmp(i,j)>2) InputGenosTmp(i,j)=9
                    if (InputGenosTmp(i,j)<0) InputGenosTmp(i,j)=9
                end do

                k=1
                do j=StartSnp,EndSnp
                    InputGenos(i,k)=InputGenosTmp(i,j)
                    k=k+1
                end do
            end do

        endif

        !call GetMaxFamilySize(nAnis,nAnis,SeqSire,SeqDam,MaxFs,MaxMates)
    end subroutine SetUpData

    !######################################################################################################################################################

    subroutine GetMaxFamilySize(nAnis,SeqSire,SeqDam,MaxFs,MaxMates) !!(maxfs, maxmates, nfamilies)
        use ISO_Fortran_Env
        implicit none

     
        integer,intent(in) :: nAnis
        integer,intent(inout) :: MaxFs,MaxMates
        integer(int64), intent(in), dimension (:) :: SeqSire(nAnis),SeqDam(nAnis)

        INTEGER (KIND= 8), allocatable :: family(:)       ! max value is 9,223,372,036,854,775,807 allowing for plenty of space
        INTEGER (KIND= 8)              :: holdfamily, multiplier
        !!INTEGER :: maxfs, maxmates, maxmates1, maxmates2, nfamilies
        INTEGER :: nFamilies,maxmates1, maxmates2
        INTEGER :: i, Noffset, Limit, Switch, fs, mates, parent1, parent2, oldparent1, oldparent2

        ALLOCATE (family(nAnis))

        !PRINT*,  ' Finding maximum FS family size and maximum mates ... '

        holdfamily=0
        multiplier = 100000000   ! allows for up to 99,999,999 in number system.

        do i=1,nAnis
         family(i) = multiplier * seqsire(i) + seqdam(i)
         !IF(family(i) /= 0) PRINT'(3i7,i15)', i, seqsire(i), seqdam(i), family(i)
        end do

          Noffset = INT(nAnis/2)
          DO WHILE (Noffset>0)
              Limit = nAnis - Noffset
              Switch=1
            DO WHILE (Switch.ne.0)
               Switch = 0
               do i = 1, Limit
                  IF (family(i).gt.family(i + Noffset)) THEN
                       holdfamily=family(i)
                       family(i)=family(i + Noffset)
                       family(i + Noffset)=holdfamily

                       Switch = i
                  endif
               enddo
               Limit = Switch - Noffset
            enddo
            Noffset = INT(Noffset/2)
          enddo

        nfamilies=0
        do i=2,nAnis
          IF(family(i) /= family(i-1)) then
            nfamilies=nfamilies+1
          endif
        end do


        MaxFs = 0
        fs=1

        do i=2,nAnis
          parent1=INT(family(i)/multiplier)
          parent2=family(i)-multiplier*parent1
         IF(parent1 /= 0 .and. parent2 /= 0) then  ! note how this is handled
          IF(family(i) == family(i-1)) then
            fs = fs + 1
          else
            IF(fs > maxfs) then
             maxfs = fs
             holdfamily=family(i-1)
            endif
            fs = 1
          endif
         endif
        end do

        if (fs > maxfs) then ! From GeneProb in AlphaImpute commited 50a88a3 by RAntolin
          maxfs = fs
          holdfamily = family(nAnis)
        end if
        
        parent1=INT(holdfamily/multiplier)
        parent2=holdfamily-multiplier*parent1
        !PRINT*, '  Maximum FS family size ... ',maxfs, ' for parents: ',parent1,' ',parent2
        !JH knocked out this line to avoid the need to allocate the vector "id"
        !PRINT*, '  Maximum FS family size ... ',maxfs, ' for parents: ',id(parent1),' ',id(parent2)

        Maxmates1 = 0
        mates=1
        oldparent1=0
        oldparent2=0

        do i=2,nAnis

          parent1=INT(family(i)/multiplier)
         IF(parent1 /= 0) then
          parent2=family(i)-multiplier*parent1
           IF(parent1 == oldparent1) then
             IF(parent2 /= oldparent2) mates=mates+1
           else
            IF(mates > maxmates1)then
             maxmates1 = mates
             Limit=oldparent1  ! to record the one with most mates
            endif
            mates=1
           endif
            If (i==nAnis .and. mates > maxmates1) Then  ! need to cover the last observation
                maxmates1 = mates
                Limit = oldparent1 ! to record the one with most mates
            End If
           oldparent1=parent1
           oldparent2=parent2
         endif
        enddo

        !JH knocked out this line to avoid the need to allocate the vector "id"
        !if (Limit > 0) then
        ! PRINT*, '  Max female mates of males=', maxmates1, ', for male ID, seqID: ', id(Limit), Limit
        !else
        ! PRINT*, '  Max female mates of males=', maxmates1
        !endif


        !Now max mates for females
        do i=1,nAnis
         family(i) = multiplier * seqdam(i) + seqsire(i)
        end do

          Noffset = INT(nAnis/2)
          DO WHILE (Noffset>0)
              Limit = nAnis - Noffset
              switch=1
            DO WHILE (Switch.ne.0)
               Switch = 0
               do i = 1, Limit
                  IF (family(i).gt.family(i + Noffset)) THEN
                       holdfamily=family(i)
                       family(i)=family(i + Noffset)
                       family(i + Noffset)=holdfamily

                       Switch = i
                  endif
               enddo
               Limit = Switch - Noffset
            enddo
            Noffset = INT(Noffset/2)
          enddo


        Maxmates2 = 0
        mates=1
        oldparent1=0
        oldparent2=0

        do i=2,nAnis
          parent1=INT(family(i)/multiplier)
         IF(parent1 /= 0) then
          parent2=family(i)-multiplier*parent1
           IF(parent1 == oldparent1) then
             IF(parent2 /= oldparent2) mates=mates+1
           else
            IF(mates > maxmates2)then
             maxmates2 = mates
             Limit=oldparent1  ! to record the one with most mates
            endif
            mates=1
           endif
            If (i==nAnis .and. mates > maxmates2) Then  ! need to cover the last observation
                maxmates2 = mates
                Limit = oldparent1 ! to record the one with most mates
            End If
           oldparent1=parent1
           oldparent2=parent2

         endif
        enddo


        !JH knocked out this line to avoid the need to allocate the vector "id"
        !if (Limit > 0) then
        ! PRINT*, '  Max male mates of females=', maxmates2, ', for female ID, seqID: ', id(Limit), Limit
        !else
        ! PRINT*, '  Max male mates of females=', maxmates2
        !endif

        maxmates = MAX(maxmates1,maxmates2)

        deallocate (family) ! added by MBattagin
    end subroutine GetMaxFamilySize

    !######################################################################################################################################################

    real(kind=8) function DFactorialInLog(n)

        implicit none
        integer,intent(in) :: n
        integer :: i
        real(kind=8) :: Ans

        Ans=0.0
        do i=1,n
            Ans=Ans+log(dfloat(i))
        enddo
        DFactorialInLog=Ans
    end function DFactorialInLog

    !######################################################################################################################################################

    real(kind=8) function DFactorial(n)

        implicit none
        integer,intent(in) :: n
        integer :: i
        real(kind=8) :: Ans

        Ans=1.0
        do i=1,n
            Ans=Ans*i
        enddo
        DFactorial=Ans
    end function DFactorial

    !######################################################################################################################################################

subroutine geneprob(currentSnp,nAnis,Seq0Snp1Mode,ReadCounts,InputGenos,maxfs,MaxMates,MaxReadCounts,GMatSnp,GMatRds,SeqSire,SeqDam,OutProb,mxeq,mate,ifirst,next,prog)
        use ISO_Fortran_Env
	      implicit none

        integer, intent(in) :: currentSnp,Seq0Snp1Mode,maxfs,MaxMates,MaxReadCounts
        integer, intent(in) :: nAnis
        integer(int64), intent(in), dimension (:) :: SeqSire(nAnis),SeqDam(nAnis)

        real(kind=8),intent(in),dimension(:,:), allocatable :: GMatSnp
        real(kind=8),intent(in),dimension(:,:,:),allocatable :: GMatRds

	      integer(kind=1),intent(in),dimension(:,:) :: InputGenos 
	      integer(kind=2),intent(in),dimension(:,:,:) :: ReadCounts 

	      real(kind=real64),intent(inout),dimension(:,:,:)  :: OutProb
        !real(kind=4),intent(inout),dimension(:,:) :: Pr00,Pr01,Pr10,Pr11 
        
	      integer,intent(in) :: mxeq

	      INTEGER, intent(inout),dimension(:) :: mate(0:2*nAnis),prog(0:2*nAnis)
	      INTEGER, intent(inout),dimension(:) :: next(2*nAnis),ifirst(0:nAnis)
	      
	      REAL (KIND=8) :: pprior, qprior, StopCrit                       ! Added by MBattagin
	      
	      INTEGER :: Imprinting,PauseAtEnd,nfreq_max                   ! Added by MBattagin
	      
	      INTEGER :: MM                                                           ! Added by MBattagin - used here and in LNKLST and ANotherOne
	      
	      REAL(KIND=8), ALLOCATABLE,dimension(:,:) :: POST                        ! Added by MBattagin - used here and in FLIPPT

	      
	      integer                 :: i, j, k, l, i2, i3,  iticks2,sumReads
	      integer                 :: ia, is, idd, ifreq_iterate, maxint, maxiter, itersused, kl, kc, kd, kj, nfams, last
	      integer                 :: nf, im, ns, mf, iaa, ii, ms, m, n, maxvalspost, ierrors, nwritten
	      integer                 :: f,ff
	      integer                 :: maxRegpoints,LeastPositive, LeastNegative, HoldInt, LimitAnimals, LimitNumber
	      integer                 :: nonzed(3,3)
	      integer                 :: ntype(3,3,3)

	      real (kind=8)                 :: tsum, prod, IMPratio, p12, p21, LeastPositiveValue, LeastNegativeValue ! ProbFit,SumFreq,
	      REAL (KIND=8)                 :: spost(3),dpost(3),fpost(3),tpost(3),temp(3),sum1(3),sum2(3),sum3(3)
	      REAL (KIND=8)                 :: pt(3,3,3)
	      REAL (KIND=8)                 :: phethw,phomhw  ! this is needed for info - or compile with dble.  Don't know why!
	      REAL (KIND=8)                 :: s0,s1,s2, areg, breg, meanX, meanY, sumY, sumXY, sumX, sumX2

	      INTEGER, allocatable,dimension(:)  :: phen,nmem,damGP ! note this is a different p1,p2 to sequence's
	      INTEGER, ALLOCATABLE,dimension(:,:):: isib

	      REAL (KIND=8), allocatable,dimension(:)       :: phom,phet,pnor,pHold,pResult,pDev
	      REAL (KIND=8), allocatable,dimension(:,:)     :: ant,term,freq
	      REAL (KIND=8), allocatable,dimension(:,:,:)   :: work

	     real (kind=8)   :: ErrorHomo,ProbHetero

       ! print*,currentSnp

	      pprior= 0.5 
	      qprior= 1-pprior
	      nfreq_max=50 
	      Imprinting=1 
	      PauseAtEnd=0 
	      StopCrit=0.0001

	      ! ----------------------------------------------------------------
	      !  P-MATRIX: PROB. OF OFFSPRING GENOTYPE GIVEN GENOTYPE OF PARENTS
	      ! ----------------------------------------------------------------

	      LimitAnimals = 0 ! 1 to invoke Limit
	      LimitNumber = 250

	      pt=0.  !  =log(1)!  the log(zero) elements should not be required ar nonzed and ntype below control addressing
	      pt(1,2,1)=log(.5)
	      pt(2,1,1)=log(.5)
	      pt(2,2,1)=log(.25)
	      pt(1,2,2)=log(.5)
	      pt(2,1,2)=log(.5)
	      pt(2,2,2)=log(.5)
	      pt(2,3,2)=log(.5)
	      pt(3,2,2)=log(.5)
	      pt(2,2,3)=log(.25)
	      pt(2,3,3)=log(.5)
	      pt(3,2,3)=log(.5)
	      nonzed(1,1)=1
	      ntype(1,1,1)=1
	      nonzed(1,2)=2
	      ntype(1,2,1)=1
	      ntype(1,2,2)=2
	      nonzed(1,3)=1
	      ntype(1,3,1)=2
	      nonzed(2,1)=2
	      ntype(2,1,1)=1
	      ntype(2,1,2)=2
	      nonzed(2,2)=3
	      ntype(2,2,1)=1
	      ntype(2,2,2)=2
	      ntype(2,2,3)=3
	      nonzed(2,3)=2
	      ntype(2,3,1)=2
	      ntype(2,3,2)=3
	      nonzed(3,1)=1
	      ntype(3,1,1)=2
	      nonzed(3,2)=2
	      ntype(3,2,1)=2
	      ntype(3,2,2)=3
	      nonzed(3,3)=1
	      ntype(3,3,1)=3

	      ALLOCATE(phen(nAnis),freq(3,0:nAnis))

	      if (Seq0Snp1Mode==0) phen=0
	      if (Seq0Snp1Mode==1) phen=9 ! covers unlisted parents

        ! THIS LOOP CONVERT THE INPUT DATA IN LOG-LIKELIHOOD
        ! Freq is the LOG-LIKELIHOOD from own information

        freq=log(1.)
        
        ErrorHomo=0.05
        ProbHetero=0.5


       ! call GetVariantErrorRate(nAnis,ReadCounts,ErrorHomo,ProbHetero,currentSnp)
       ! if ((currentSnp.ge.35001) .and. (currentSnp.le.36000)) write(*,'(1i0,1x,1f7.4)'),currentSnp,ErrorHomo

        do ia=1,nAnis 

          if (Seq0Snp1Mode==0) phen(ia) = ReadCounts(ia,currentSnp,2) !!!! was ReadCounts(i,1,2) - MBattagin
          if (Seq0Snp1Mode==1) phen(ia) = InputGenos(ia,currentSnp)

          if (Seq0Snp1Mode==0) then
              sumReads=0
              sumReads=sum(ReadCounts(ia,currentSnp,:))
              if (sumReads.eq.0) THEN
                  freq(:,ia) =log(1.)
              else
                  do i = 0, sumReads
                      IF (phen(ia).eq.i) THEN
                          IF (GMatRds(i, 1,sumReads).lt.log(.000000001))then 
                              freq(1,ia) =-9999
                          else
                              freq(1,ia) =GMatRds(i, 1,sumReads)
                          endif
                          IF(GMatRds(i, 2,sumReads).lt.log(.000000001))then
                              freq(2,ia) =-9999
                          else
                              freq(2,ia)=GMatRds(i, 2,sumReads)
                          endif
                          IF(GMatRds(i, 3,sumReads).lt.log(.000000001))then
                              freq(3,ia) =-9999
                          else
                              freq(3,ia) =GMatRds(i, 3,sumReads)
                          endif
                      endif
                  enddo
              endif
          endif

          if (Seq0Snp1Mode==1) then
              if (phen(ia).eq.9) then
                  freq(1,ia) =log(1.)
                  freq(2,ia) =log(1.)
                  freq(3,ia) =log(1.)
              else
                  do i = 1,3
                      
                      if (phen(ia).eq.i) then
                          if (GMatSnp(i,1).lt.(.000000001)) then 
                              freq(1,ia) =-9999
                          else
                              freq(1,ia) =log(GMatSnp(i, 1)) 
                          endif
                          if (GMatSnp(i,2).lt.(.000000001)) then 
                              freq(2,ia) =-9999
                          else
                              freq(2,ia) =log(GMatSnp(i, 2)) 
                          endif
                          if (GMatSnp(i,3).lt.(.000000001)) then 
                              freq(3,ia) =-9999
                          else
                              freq(3,ia) =log(GMatSnp(i, 3)) 
                          endif
                      endif
                  enddo
              endif
          endif
        end do


        ALLOCATE(post(3,0:2*nAnis), &
                ant(3,0:nAnis),phom(0:nAnis),phet(0:nAnis), &
                pnor(0:nAnis))


	      post=0. 
	      phom=0.
	      phet=0.
	      pnor=0.
	      
	      HoldInt = MAX(2*maxfs,maxmates)
	      
	      ALLOCATE(nmem(HoldInt),isib(maxmates,2*maxfs),damGP(maxmates),work(3,2*maxfs,2*maxfs),term(3,HoldInt))

	      HoldInt=0
	      isib=0


	      maxRegpoints=5

	      ifreq_iterate = -1 
	      if (nfreq_max==1) nfreq_max=2
	      ALLOCATE (pHold(0:nfreq_max), pResult(0:nfreq_max), pDev(0:nfreq_max)) 
	      
	      if (nfreq_max==0)then
	        pHold(0)=pprior  ! one hit only
	      else
	        pHold(0)=0.001 
	        pHold(1)=0.999 
	      endif

	      pDev(0)=0.
	      LeastPositiveValue =  999.
	      LeastNegativeValue = -999.
	      LeastPositive = 0
	      LeastNegative = 0


	      do WHILE (ifreq_iterate < nfreq_max -1)

	        ifreq_iterate = ifreq_iterate + 1 

	        pprior = pHold(ifreq_iterate)
	        qprior = 1-pprior

	        post=0.
	        phet=0.

	        do i=1,nAnis ! M Battagin - removed use of the prior (i.e., AlleleFrequencies)
	          ant(1,i)=log(0.5*0.5)   !log(qprior*qprior)     ! log(.000000001)! 
	          ant(2,i)=log(2*0.5*0.5) !log(2.0*pprior*qprior) ! log(.000000001)! 
	          ant(3,i)=log(0.5*0.5)   !log(pprior*pprior)     ! log(.000000001)! 
	        enddo

	        ! ----------------------------------------------
	        ! Prob(Gi) = ( Ai f(Yi|Gi) PROD - mates Pi ) / L
	        !     where L = SUM-Gi Ai f(Yi|Gi) PROD - mates Pi
	        !    Ai is the joint probability of phenotypes of members anterior to i
	        !       genotype Gi for i
	        !    PROD over mates Pi is the conditional probability of phenotypes of
	        !     posterior to i, given i has genotype Gi
	        ! ----------------------------------------------
	        ! BUILD A LINKLIST. EACH PARENT ANIMAL HAS A ROW, ALONG THE COLUMNS ARE
	        ! OF THEIR MATES AND PROGENY. MATES CAN BE REPEATED AT SUCCESSIVE NODES
	        ! REFLECT FULL SIB FAMILIES. NOTE, THE FIRST NODE FOR A PARTICULAR MATE
	        ! CONTAIN THE POST(i,j) term, (the jth mate of the ith animal)
	        ! To initialise the iterative peeling up and peeling down cycles, the an
	        ! term for founder animals is set equal to the HW probs. Post. terms for
	        ! animals and ant. terms for non founders are set to 1 (reflecting no
	        ! information)

	        maxint=9999999
	        maxiter=7
	        itersUsed=maxiter

	        do kl=1,maxiter
	          !!      IF(nfreq_max==0) print*,'  Iteration number ',kl
	          !print*,'  Iteration number ',kl
	          ! ----------------------------------------------------------------
	          ! PEEL DOWN, IE CONDENSE INFO ON MUM AND DAD ONTO PROGENY
	          ! START WITH OLDEST ANIMAL IN THE LIST. WHEN DESCENDING WE CALCULATE
	          ! ANT() TERMS FOR ALL PROGENY. PICK OUT PROGENY FROM THE ROWS OF
	          ! THE PARENT WITH ON AVERAGE THE MOST MATES OR PROGENY, USUALLY THE SIRE
	          ! IGNORE THE ROWS OF THE OTHER PARENT.
	          ! ----------------------------------------------------------------
	          

            do is=1,nAnis
              kc=ifirst(is)
	            if(kc.eq.0.or.kc.gt.nAnis) then
	              ! do nothing
	            else
	              ! load up the mates of this sire and collect posterior terms
	              idd=mate(kc)
	              spost(:)=0.0
	              nfams=0
	              last=maxint
	              do while (kc.ne.0)
	                ia=prog(kc)
	                if(idd.ne.last) then
	                  nfams=nfams+1
	                  nmem(nfams)=1
	                  damGP(nfams)=idd
	                  isib(nfams,1)=ia
	                  term(1,nfams)=post(1,kc)
	                  term(2,nfams)=post(2,kc)
	                  term(3,nfams)=post(3,kc)
	                  spost(1)=spost(1)+term(1,nfams)
	                  spost(2)=spost(2)+term(2,nfams)
	                  spost(3)=spost(3)+term(3,nfams)
	                else
	                  nmem(nfams)=nmem(nfams)+1
	                  isib(nfams,nmem(nfams))=ia
	                endif
	                last=idd
	                kc=next(kc)
	                idd=mate(kc)
	              enddo

	              do nf=1,nfams     ! GO THROUGH THROUGH THE MATES "idd" OF "is
	                idd=damGP(nf)
	                ! collect posterior term for damGP "idd" through all its mates
	                dpost(:)=0.
	                kc=ifirst(idd)
	                im=mate(kc)
	                last=maxint
	                do while (kc.ne.0)
	                  if(im.ne.is.and.im.ne.last) then
	                    dpost(1)=dpost(1)+post(1,kc)
	                    dpost(2)=dpost(2)+post(2,kc)
	                    dpost(3)=dpost(3)+post(3,kc)
	                  endif
	                  last=im
	                  kc=next(kc)
	                  im=mate(kc)
	                enddo
	                ! correct posterior prob of "is" for "idd"
	                tpost(1)=spost(1)-term(1,nf)
	                tpost(2)=spost(2)-term(2,nf)
	                tpost(3)=spost(3)-term(3,nf)
	                ! for this damGP "idd", and for each of her progeny "ia" to "is"
	                ! mark out the full sibs "iaa" to "ia" all the k mates of "iaa"
	                ! and store the term prod-k post(ia,iaa) in a work vector
	                do ns=1,nmem(nf)
	                  ia=isib(nf,ns)
	                  do mf=1,nmem(nf)
	                    iaa=isib(nf,mf)
	                    if(iaa.ne.ia) then
	                      fpost(:)=0.
	                      kc=ifirst(iaa) ! collect mates of "iaa"
	                      ms=mate(kc)
	                      last=maxint
	                      do while (kc.ne.0)
	                        if(ms.ne.last) then
	                          fpost(1)=fpost(1)+post(1,kc)
	                          fpost(2)=fpost(2)+post(2,kc)
	                          fpost(3)=fpost(3)+post(3,kc)
	                        endif
	                        last=ms
	                        kc=next(kc)
	                        ms=mate(kc)
	                      enddo
	                      do i=1,3
	                        work(i,ns,mf)=fpost(i)
	                      enddo
	                    endif
	                  enddo
	                enddo

	                !  NOW WE ARE READY To CALCULATE THE BLOODY THING
	                do ns=1,nmem(nf)
	                  ia=isib(nf,ns)
                    !if (currentSnp==1) write(*,'(2i10,3f10.5)') currentSnp,ia,ant(:,ia)
              
	                  ! Here we will get the anterior probability for the progeny in question.  For appendix equations' m, f, s and i:
	                  ! m is here is   - the sire of i as in the outer loop
	                  ! f is here imum - the dam of i
	                  ! s is here iaa  - the sibs of i
	                  ! i is here ia   - the progeny animal
	                  do i=1,3   ! MBattagin - here works only if there aren't Mendelian inconsistencies
	                    mm=0
	                    do m=1,3
	                      if (abs(i-m)/=2) then ! MBattagin - exclude opposing Homozygotes
	                        mm=mm+1
	                        ff=0
	                        do f=1,3
	                          if(i.eq.1.and.f.eq.3) then 
	                            !do nothing
	                          else if(i.eq.2.and.(m+f.lt.3.or.m+f.gt.5)) then
	                            !do nothing
	                          else if(i.eq.3.and.f.eq.1)then
	                            !do nothing
	                          else
	                            ff=ff+1
	                            prod=0.
	                            do mf=1,nmem(nf) ! go through this animal' si
	                              iaa=isib(nf,mf)
	                              if(iaa.ne.ia) then
	                                do j=1,nonzed(m,f)
	                                  k=ntype(m,f,j)
	                                  sum3(j)=pt(m,f,k)+freq(k,iaa)+work(k,ns,mf)
	                                enddo
	                                call LOGADD(sum3,nonzed(m,f))
	                                prod=prod+sum3(1)
	                              endif
	                            enddo ! line 4 calculated
	                            sum2(ff)=ant(f,idd)+freq(f,idd)+dpost(f)+pt(m,f,i)+prod
	                          endif
	                          
	                        enddo
	                        call LOGADD(sum2,ff)
	                        sum1(mm)=ant(m,is)+freq(m,is)+tpost(m)+sum2(1)
	                      endif
	                    enddo
	                    call LOGADD(sum1,mm)
	                    ant(i,ia)=sum1(1)
	                    temp(i)=sum1(1)
	                  enddo
	                  call LOGADD(temp,3)
                    do i=1,3
	                    ant(i,ia)=ant(i,ia)-temp(1)
	                  enddo
                    !write(*,'(2i10,3f10.5)') currentSnp,ia,ant(:,ia)
                    
	                enddo
	              enddo
	            endif
	            
	          enddo

	        ! ----------------------------------------------------------------
	        ! PEEL UP, IE CONDENSE INFO ON MATE AND PROGENY ONTO INDIVIDUAL
	        ! START WITH YOUNGEST IN THE LIST
	        ! ----------------------------------------------------------------
	        call flippt(mxeq,nAnis,MATE,IFIRST,NEXT,POST)
	        
	        do i=nAnis,1,-1
	          kc=ifirst(i)
	          do while(kc/=0)
	            im=mate(kc)
	            !     collect posterior probability for "im" through all mates "is" of "im"
	            kd=ifirst(im)
	            is=mate(kd)
	            spost(:)=0.
	            last=maxint
	            do while (kd.ne.0)
	              if(is.ne.last.and.is.ne.i) then
	                do ii=1,3
	                  spost(ii)=spost(ii)+post(ii,kd)
	                enddo
	              endif
	              last=is
	              kd=next(kd)
	              is=mate(kd)
	            enddo
	            !     pick up all the offspring of "i" and "im", say "idd" and
	            !     go through through the l mates "ms" of "idd", and calculate
	            !     prod-l post(idd,ms), store in work vector
	            is=im
	            kd=kc
	            k=0
	            do while (is.eq.im.and.kd.ne.0)
	              k=k+1
	              idd=prog(kd)
	              nmem(k)=idd
	              kj=ifirst(idd)
	              ms=mate(kj)
	              dpost(:)=0.
	              last=maxint
	              do while (kj.ne.0)        !Go through all mates of "idd"
	                if(ms.ne.last) then
	                  do ii=1,3
	                    dpost(ii)=dpost(ii)+post(ii,kj)
	                  enddo
	                endif
	                last=ms
	                kj=next(kj)
	                ms=mate(kj)
	              enddo
	              do ii=1,3
	                term(ii,k)=dpost(ii)
	              enddo
	              kd=next(kd)
	              is=mate(kd)
	            enddo
	          
	            !     NOW we are ready to calculate the posterior prob for match "i" and "im"
	            do ii=1,3                   !for 3 g'type of "i"
	              do j=1,3                 !for 3 g'types of "im"
	                prod=0.
	                do ns=1,k             ! through k offspring
	                  idd=nmem(ns)
	                  do m=1,nonzed(ii,j)
	                    l=ntype(ii,j,m)
	                    sum2(m)=pt(ii,j,l)+freq(l,idd)+term(l,ns)
	                  enddo
	                  call LOGADD(sum2,nonzed(ii,j))
	                  prod=prod+sum2(1)
	                enddo
	                sum1(j)=ant(j,im)+freq(j,im)+spost(j)+prod
	              enddo
	              call LOGADD(sum1,3)
	              temp(ii)=sum1(1)
	              post(ii,kc)=sum1(1)
	            enddo
	            
	            call LOGADD(temp,3)
	            do ii=1,3
	              post(ii,kc)=post(ii,kc)-temp(1)
	              ! print*,i,im,ii,post(ii,kc)
	            enddo
	            kc=kd
	          enddo
	        enddo

	        ! ------------------------------------------------------------------
	        !     NOW CALCULATE GENOTYPE PROBABILITIES
	        sum1(1)=0.
	        do i=1,nAnis
	          spost(:)=0.
	          kc=ifirst(i)
	          if (kc/=0) then
	            im=mate(kc)
	            last=maxint
	            do while (kc.ne.0)
	              if(im.ne.last) then
	                do ii=1,3
	                  spost(ii)=spost(ii)+post(ii,kc)
	                enddo
	              endif
	              last=im
	              kc=next(kc)
	              im=mate(kc)
	            enddo
	          endif
	         tsum=0.
	          do ii=1,3
	            spost(ii)=ant(ii,i)+freq(ii,i)+spost(ii)  ! just store it all in spost for convenience
	          enddo
	          maxvalspost=MAXVAL(spost)
	          do ii=1,3
	            spost(ii)=spost(ii)-maxvalspost + 20   ! bring them all up towards zero, then some (e^20 = 485 million = OK)
	          enddo
	          do ii=1,3
	            if (spost(ii).LT.-100.) then  !e^-100  =~ 10^-44
	              temp(ii)=0.
	            else
	              temp(ii)=exp(spost(ii))
	            end if
	            tsum=tsum+temp(ii)
	          enddo
	          pnor(i)=temp(1)/tsum

	          prod=temp(2)/tsum
	          sum1(1)=sum1(1)+abs(prod-phet(i))
	          phet(i)=prod
	          phom(i)=temp(3)/tsum
	        enddo


	        sum1(1)=sum1(1)/nAnis
	        !         write(*,'(f13.7)') sum1(1)
	        
	        call flippt(mxeq,nAnis,MATE,IFIRST,NEXT,POST)
	        
	        if(sum1(1).le.StopCrit) then
	          itersUsed=kl
	          exit
	        endif
	      enddo

	      s0=0.
	      s1=0.
	      s2=0.
	      n=0
	      do i=1,nAnis

	        if (ABS(1.-pnor(i)-phet(i)-phom(i)) > .0000001) then
	          PRINT*, 'Error: ',ABS(1.-pnor(i)-phet(i)-phom(i)),pnor(i),phet(i),phom(i)
	        end if

	        if(SeqSire(i).eq.0 .and. SeqDam(i).eq.0)then
	          s0=s0+pnor(i)
	          s1=s1+phet(i)
	          s2=s2+phom(i)
	          n=n+1
	        endif
	      enddo


	      phethw=2*pprior*qprior ! do these here before reset pprior so GPI comes from freq used to get probabilities
	      phomhw=pprior*pprior

	      pResult(ifreq_iterate) = (s1+2*s2)/(2*(s0+s1+s2))
	      pDev(ifreq_iterate) = pResult(ifreq_iterate) - pHold(ifreq_iterate)
	      !print*,currentSnp,ifreq_iterate,pResult(ifreq_iterate),pHold(ifreq_iterate)
	      if(pDev(ifreq_iterate) > 0.0 .and. pDev(ifreq_iterate) < LeastPositiveValue ) LeastPositive = ifreq_iterate
	      if(pDev(ifreq_iterate) < 0.0 .and. pDev(ifreq_iterate) > LeastNegativeValue ) LeastNegative = ifreq_iterate

	      !print'(2i6,3f16.8)', currentSnp,ifreq_iterate, pHold(ifreq_iterate), pResult(ifreq_iterate), pDev(ifreq_iterate)

	      if (ifreq_iterate==0) then
	        ! do nothing here
	      
	      elseif (ifreq_iterate==1) then
	        ! 2-point regression
	        breg = (pDev(1) - pDev(0)) / (pHold(1) - pHold(0))
	        areg = pDev(1) - breg*pHold(1)
	        pHold(ifreq_iterate+1)= -1*areg/breg
	        if (pDev(0)<0.  .AND. pDev(1)>0. ) then
	          nfreq_max=ifreq_iterate+1 ! Blowing out already - get sensible midpoint
	        endif
	      elseif (ifreq_iterate>1) then
	        ! j-point regression
	        j=MIN(maxRegpoints,ifreq_iterate+1)  ! 3 point seems fastest.  More points could be more robust.
	        meanX=0
	        meanY=0
	        sumY =0
	        sumXY=0
	        sumX =0
	        sumX2=0
	        do i = ifreq_iterate-(j-1), ifreq_iterate
	          meanX = meanX + pHold(i)
	          meanY = meanY + pDev(i)
	          sumY  = sumY  + pDev(i)
	          sumX  = sumX  + pHold(i)
	          sumX2 = sumX2 + pHold(i)*pHold(i)
	          sumXY = sumXY + pHold(i)*pDev(i)
	        end do
	        breg  = (sumXY - sumX*sumY/j) / (sumX2 - sumX*sumX/j)
	        meanX = meanX/j
	        meanY = meanY/j
	        areg  = meanY - breg*meanX
	        pHold(ifreq_iterate+1)= -1*areg/breg
	        !print'(2i4,3f12.8)', j, ifreq_iterate+1, areg, breg, pHold(ifreq_iterate+1)

	      end if

	      if (ifreq_iterate>=1) then

	        if(pHold(ifreq_iterate+1) > 0.9999 .OR. pHold(ifreq_iterate+1) < 0.0001) then
	          if (pHold(ifreq_iterate+1) > 0.999) pHold(ifreq_iterate+1) = 0.999 
	          if (pHold(ifreq_iterate+1) < 0.001) pHold(ifreq_iterate+1) = 0.001 
	          !nfreq_max=ifreq_iterate+1 ! and stop after getting probs from that (redundant - previous line makes a StopCrit stop)
	        endif

	      endif


	      if(nfreq_max>0) then
	        IF ( ABS(pHold(ifreq_iterate)-pHold(ifreq_iterate+1)) < StopCrit ) ifreq_iterate=nfreq_max
	      endif

	      ierrors=1
	      IF ( ifreq_iterate==nfreq_max) then
	        pnor(0)=1.- phethw - phomhw
	        phet(0)=phethw
	        phom(0)=phomhw

	        if (LimitAnimals==1 .and. nAnis>LimitNumber) then
	          PRINT*, 'Results for only', LimitNumber, ' animals will be written.'
	          nwritten=LimitNumber
	        else
	          nwritten=nAnis
	        endif

	        do i=1,nwritten
	          !call info(phet(i), phom(i), phethw, phomhw, probindex)

	         
	          if(Imprinting>0) then
	            if(phet(i)<0.0000001) then
                OutProb(1,currentSnp,i)=pnor(i)
                OutProb(2,currentSnp,i)=phet(i)
                OutProb(3,currentSnp,i)=phet(i)
                OutProb(4,currentSnp,i)=phom(i)
	            else
	              p12= (pnor(seqsire(i))+0.5*phet(seqsire(i))) * (phom( seqdam(i))+0.5*phet( seqdam(i)))  ! extra safe due to the above
	              p21= (pnor( seqdam(i))+0.5*phet( seqdam(i))) * (phom(seqsire(i))+0.5*phet(seqsire(i)))
	                if(p12+p21>0.00000001) then
	                  IMPratio= p12 / (p12+p21)
	                else
	                  IMPratio= 0.5 ! neutral but should not be invked anyway
	                endif

	                if(Imprinting==2)then
                    OutProb(1,currentSnp,i)=pnor(i)
                    OutProb(2,currentSnp,i)=phet(i)
                    OutProb(3,currentSnp,i)=(1.-2.*IMPratio)*phet(i)
                    OutProb(4,currentSnp,i)=phom(i)
	                else
                    OutProb(1,currentSnp,i)=pnor(i)
                    OutProb(2,currentSnp,i)=IMPratio*phet(i)
                    OutProb(3,currentSnp,i)=(1.-IMPratio)*phet(i)
                    OutProb(4,currentSnp,i)=phom(i)
	                endif
	            endif
	          else
              OutProb(1,currentSnp,i)=pnor(i)
              OutProb(2,currentSnp,i)=phet(i)
              OutProb(3,currentSnp,i)=phet(i)
              OutProb(4,currentSnp,i)=phom(i)
	          endif
	        enddo

	        !OutputMaf(currentSnp)=pprior
	        !print*,currentSnp,OutputMaf(currentSnp)
	        call system_clock(iticks2,i2,i3)
	      endif

	    ENDDO ! ifreq_iterate

	    !print *,'done with geneprob', currentSnp  !!*******************************************************!!

	    deALLOCATE(POST,phen,ant,phom,phet,freq,pnor)
	    deallocate(nmem)
	    deallocate(isib)
	    deallocate(damGP)
	    deallocate(work)
	    deallocate(term)
	    
	    deALLOCATE (pHold, pResult, pDev)
  end subroutine geneprob

    !######################################################################################################################################################

    SUBROUTINE LNKLST(I,J,NP,IFLAG,MXEQ,nAnis,mm,nn,mate,next,ifirst,prog)
      use ISO_Fortran_Env
      implicit none 

      integer,intent(in) :: i,j,np,iFlag,nAnis
      integer,intent(inout) :: MXEQ,MM,NN                        ! Added by MBattagin !! are MXEQ,MM,NN (inout)?
      
      INTEGER, intent(inout),dimension(:) :: MATE(0:2*nAnis),prog(0:2*nAnis)
      INTEGER, intent(inout),dimension(:) :: NEXT(2*nAnis),IFIRST(0:nAnis)  ! Added by MBattagin 
      
      integer :: IR,IP,K,NUSED
    

      IF (I.le.0.or.J.le.0) RETURN

      !     KEEP MAXIMUM EQUATION
      IF (I.GT.MXEQ)THEN
         MXEQ=I
         IF (I.GT.nAnis)STOP ' TOO MANY EQUATIONS - RECOMPILE'
      ENDIF

      IR=0
      IP=IFIRST(I)
      K=MATE(IP)
      do while(K.LE.J.AND.IP.NE.0) 
        ! IF(np.le.30) print '(6i9)',np,k,j,ip,next(ip),mate(next(ip))
         IR=IP
         IP=NEXT(IP)
         K=MATE(IP)
      enddo

      IF(IFLAG.EQ.1) THEN
         MM=MM+1
         NUSED=MM
      ELSE
         NN=NN+1
         NUSED=NN
      ENDIF

      MATE(NUSED)=J
      PROG(NUSED)=NP
      IF (IP.EQ.0) THEN
         IF (IFIRST(I).EQ.0) THEN
            IFIRST(I)=NUSED
         ELSE
            NEXT(IR)=NUSED
         ENDIF
      ELSEIF (IR.EQ.0) THEN
         NEXT(NUSED)=IFIRST(I)
         IFIRST(I)=NUSED
      ELSE
         NEXT(NUSED)=NEXT(IR)
         NEXT(IR)=NUSED
      ENDIF
    end subroutine LNKLST

    !********************************************************************

    SUBROUTINE FLIPPT(mxeq,nAnis,MATE,IFIRST,NEXT,POST)
        use ISO_Fortran_Env
        implicit none

        integer,intent(in) :: mxeq,nAnis
        integer,intent(in),dimension(:) :: MATE(0:2*nAnis) 
        
        integer,intent(inout),dimension(:) :: IFIRST(0:nAnis),NEXT(2*nAnis)
        REAL(KIND=8),intent(inout),dimension(:,:) :: POST(3,0:2*nAnis) 

        integer :: k,kc,NEWPT,IOLDPT,i

        DO K=1,MXEQ
          KC=IFIRST(K)
          NEWPT=0
          do while (kc/=0)
            IOLDPT=NEXT(KC)
            NEXT(KC)=NEWPT
            NEWPT=KC
            KC=IOLDPT
            IF(KC.NE.0.AND.MATE(KC).EQ.MATE(NEWPT)) THEN
              DO I=1,3
                POST(I,IOLDPT)=POST(I,NEWPT)
              ENDDO
            ENDIF
          enddo
          if (kc==0) then
            IFIRST(K)=NEWPT
          endif
        enddo
    end SUBROUTINE FLIPPT

    !######################################################################################################################################################

    SUBROUTINE LOGADD(summ,n)
      implicit none
      integer,intent(in) :: n
      real(kind=8),intent(inout),dimension(:) :: summ(3)
      
      integer :: i,mm,j
      real(kind=8) :: t

      if(n.eq.1) return
      if(n.eq.2) then
         summ(1)=add(summ(1),summ(2))
      else
         do i=1,2
            mm=i
            do j=i+1,3
               if(summ(j).lt.summ(mm)) then
                  t=summ(mm)
                  summ(mm)=summ(j)
                  summ(j)=t
               endif
            enddo
         enddo
         summ(1)=add(summ(1),summ(2))
         summ(1)=add(summ(1),summ(3))
      endif
      return
    end SUBROUTINE LOGADD

    !######################################################################################################################################################

    function add(x1,x2)
        implicit none
        real(kind=8) :: add,diff,expmax
        real(kind=8),intent(in) ::  x1,x2
          expmax=300.d0
          diff=x1-x2
          if(diff.gt.expmax) then
             add=x1
          elseif(-diff.gt.expmax) then
             add=x2
          else
             add=x2 + dlog(dexp(diff) + 1.d0)
          endif
          return
    end function add

    !######################################################################################################################################################

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

    subroutine ReadsLikelihood(nRef,nAlt,ErrorHomo,ProbHetero,lPr0,lPr1,lPr2)
      implicit none
      
      integer(kind=2),intent(in)  :: nRef,nAlt
      REAL (KIND=8),intent(in)    :: ErrorHomo,ProbHetero
      REAL (KIND=8),intent(inout) :: lPr0,lPr1,lPr2

          if ((nRef+nAlt)==0) then
            lPr0=log(1.)
            lPr1=log(1.)
            lPr2=log(1.)
          else
            lPr0=(dble(nRef)*log(1-ErrorHomo))+(dble(nAlt)*log(ErrorHomo))
            if (lPr0.lt.log(.000000001)) lPr0=-9999
            
            lPr1=(dble(nRef)*log(ProbHetero))+(dble(nAlt)*log(1-ProbHetero))
            if (lPr1.lt.log(.000000001)) lPr1=-9999
            
            lPr2=(dble(nAlt)*log(1-ErrorHomo))+(dble(nRef)*log(ErrorHomo))
            if (lPr2.lt.log(.000000001)) lPr2=-9999
          endif
          write(*,'(2i3,3f12.4)'),nRef,nAlt,lPr0,lPr1,lPr2
    end subroutine ReadsLikelihood

    !###########################################################################################################################################################

    subroutine GetVariantErrorRate(nInd,ReadCounts,ErrorHomo,ProbHetero,currentSnp)
      implicit none
      
      integer,intent(in)                          :: nInd,currentSnp
      integer(kind=2),intent(in),dimension(:,:,:) :: ReadCounts 
      real (kind=8), intent(inout)                :: ErrorHomo,ProbHetero

      integer :: i,nHomo,nHetero
      real (kind=8) :: lPr0,lPr1,lPr2,prob0,prob1,prob2,oldProbHetero,oldErrorHomo,cHetero,cHomo

      ErrorHomo=0.001
      ProbHetero=0.5
      oldErrorHomo=0.001
      oldProbHetero=0.5
      cHetero=1.
      cHomo=1.
      
      do while ((cHomo.gt.0.00001).or.(cHetero.gt.0.00001))

        nHomo=0
        nHetero=0
        
        oldErrorHomo=ErrorHomo
        oldProbHetero=ProbHetero

        ErrorHomo=0.001
        ProbHetero=0.5

        if (currentSnp==1) write(*,'(1i10,2f8.4)') currentSnp,oldErrorHomo, oldProbHetero

        do i=1,nInd
          call ReadsLikelihood(ReadCounts(i,currentSnp,1),ReadCounts(i,currentSnp,2),oldErrorHomo,oldProbHetero,lPr0,lPr1,lPr2) 
          
          prob0=exp(lPr0)/(exp(lPr0)+exp(lPr1)+exp(lPr2))
          prob1=exp(lPr1)/(exp(lPr0)+exp(lPr1)+exp(lPr2))
          prob2=exp(lPr2)/(exp(lPr0)+exp(lPr1)+exp(lPr2))

          if ((prob0.gt.prob1).and.(prob0.gt.prob2)) then
            nHomo=nHomo+1
            ErrorHomo=ErrorHomo+(prob0*dble(ReadCounts(i,currentSnp,2))/dble(sum(ReadCounts(i,currentSnp,:))))
!            if (currentSnp==1) write(*,'(1i10,4f10.4)') nHomo,oldErrorHomo,ErrorHomo,prob0,dble(ReadCounts(i,currentSnp,2))/dble(sum(ReadCounts(i,currentSnp,:)))
          endif

          if ((prob2.gt.prob1).and.(prob2.gt.prob0)) then 
            nHomo=nHomo+1
            ErrorHomo=ErrorHomo+(prob2*dble(ReadCounts(i,currentSnp,1))/dble(sum(ReadCounts(i,currentSnp,:))))
!            if (currentSnp==1) write(*,'(1i10,4f10.4)') nHomo,oldErrorHomo,ErrorHomo,prob0,dble(ReadCounts(i,currentSnp,1))/dble(sum(ReadCounts(i,currentSnp,:)))
          endif

          ! if ((prob1.gt.prob0).and.(prob1.gt.prob2)) then
          !   nHetero=nHetero+1
          !   ProbHetero=ProbHetero+(prob1*ReadCounts(i,currentSnp,2)/sum(ReadCounts(i,currentSnp,:)))
          ! endif
        enddo

        if (nHomo.gt.0) ErrorHomo=ErrorHomo/dble(nHomo)
        !if (nHetero.gt.0) ProbHetero=ProbHetero/dble(nHetero)

        if (nHomo.eq.0) then
          ErrorHomo=0.001
          exit
        endif
        if (nHetero.eq.0) ProbHetero=0.5

        cHomo=abs(oldErrorHomo-ErrorHomo)
        cHetero=abs(oldProbHetero-ProbHetero)

      enddo

      if (ErrorHomo.eq.0) ErrorHomo=0.0001

    end subroutine GetVariantErrorRate


end module AlphaVarCallFuture


