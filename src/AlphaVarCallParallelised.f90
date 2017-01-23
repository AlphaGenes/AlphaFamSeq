module AlphaVarCallFuture

implicit none

public :: AlphaVarCall

contains

    !######################################################################################################################################################

    subroutine AlphaVarCall(nAnis,nSnp,StartSnp,EndSnp,ErrorRate,Seq0Snp1Mode,SeqId,SeqSire,SeqDam,ReadCountsTmp,InputGenosTmp,Pr00,Pr01,Pr10,Pr11)

      use omp_lib

      implicit none

      integer, intent(in) :: nAnis,nSnp,StartSnp,EndSnp,Seq0Snp1Mode
      real(kind=8),intent(in) :: ErrorRate
      
      integer, intent(in), dimension (:) :: SeqId(nAnis),SeqSire(nAnis),SeqDam(nAnis)
      
      real(kind=8),intent(in),dimension(:,:,:) :: ReadCountsTmp 
      integer,intent(inout),dimension(:,:) :: InputGenosTmp 
      

      real(kind=4),intent(inout),dimension(:,:) :: Pr00(nAnis,EndSnp-StartSnp+1)
      real(kind=4),intent(inout),dimension(:,:) :: Pr01(nAnis,EndSnp-StartSnp+1)
      real(kind=4),intent(inout),dimension(:,:) :: Pr10(nAnis,EndSnp-StartSnp+1)
      real(kind=4),intent(inout),dimension(:,:) :: Pr11(nAnis,EndSnp-StartSnp+1)
      
      integer :: MaxFs,MaxMates,MaxReadCounts

      real(kind=8),allocatable,dimension(:) :: nReadCounts
      real(kind=8),allocatable,dimension(:,:,:) :: ReadCounts                                         
      
      integer,allocatable,dimension(:,:) :: InputGenos                                          

      real(kind=8),dimension(:) :: OutputMaf(EndSnp-StartSnp+1)

      real(kind=8),allocatable,dimension(:,:) :: GMatSnp 
      real(kind=8),allocatable,dimension(:,:,:) :: GMatRds

      integer :: mxeq
      INTEGER,allocatable,dimension(:) :: mate,prog,next,ifirst,p1,p2
      
      real(kind=8)::tstart,tend
      integer :: i ,j
      
      call GetMaxFamilySize(nAnis,SeqSire,SeqDam,MaxFs,MaxMates)
      call SetUpData(Seq0Snp1Mode,ReadCounts,InputGenos,nAnis,nSnp,EndSnp,StartSnp,nReadCounts,ReadCountsTmp,InputGenosTmp)
      call SetUpEquationsForSnp(nAnis,Seq0Snp1Mode,nReadCounts,GMatSnp,GMatRds,ErrorRate,MaxReadCounts)
      call CreateLinkListArrays(nAnis,SeqSire,SeqDam,mxeq,mate,ifirst,next,prog,p1,p2)
      tstart = omp_get_wtime()
      !!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) PRIVATE(i) SHARED(ReadCounts,nReadCounts,InputGenos,MaxReadCounts,GMatSnp,GMatRds,Pr00, Pr01, Pr10, Pr11)
      do i=1,(EndSnp-StartSnp+1)
        call geneprob(i,ErrorRate,nAnis,Seq0Snp1Mode,ReadCounts,InputGenos,nReadCounts,EndSnp,StartSnp, &
                      maxfs,MaxMates,MaxReadCounts,GMatSnp,GMatRds,OutputMaf,SeqId,SeqSire,SeqDam, &
                      Pr00,Pr01,Pr10,Pr11, &
                      mxeq,mate,ifirst,next,prog,p1,p2)
      enddo
      !!$OMP END PARALLEL DO
      
      tend = omp_get_wtime()
      write(*,*) "Total wall time for GeneProbController is ", tend - tstart


      if (Seq0Snp1Mode==0) then
        deallocate(GMatRds)
      endif

      if (Seq0Snp1Mode==1) then
        deallocate(GMatSnp)
      endif

    deallocate(MATE,NEXT,IFIRST,prog,p1,p2)
 
      
    end subroutine AlphaVarCall

    !######################################################################################################################################################
      ! BUILD A LINKLIST. EACH PARENT ANIMAL HAS A ROW, ALONG THE COLUMNS ARE
      ! OF THEIR MATES AND PROGENY. MATES CAN BE REPEATED AT SUCCESSIVE NODES
      ! REFLECT FULL SIB FAMILIES. NOTE, THE FIRST NODE FOR A PARTICULAR MATE
      ! CONTAIN THE POST(i,j) term, (the jth mate of the ith animal)
      ! To initialise the iterative peeling up and peeling down cycles, the an
      ! term for founder animals is set equal to the HW probs. Post. terms for
      ! animals and ant. terms for non founders are set to 1 (reflecting no
      ! information)

    subroutine CreateLinkListArrays(nAnis,SeqSire,SeqDam,mxeq,mate,ifirst,next,prog,p1,p2)
 
      implicit none

      integer,intent(in) :: nAnis
      integer, intent(in), dimension (:) :: SeqSire(nAnis),SeqDam(nAnis)

      integer,intent(inout) :: mxeq

      INTEGER, intent(inout),allocatable,dimension(:) :: mate,prog,next,ifirst,p1,p2 
   

      integer :: mm,nn,is,idd,ia

      mm=0
      nn=nAnis

      ALLOCATE(mate(0:2*nAnis), &
               next(2*nAnis), &
               ifirst(0:nAnis), &
               prog(0:2*nAnis),&
               p1(nAnis), &
               p2(nAnis))


      ! Create LinkList
      mate=0
      next=0
      ifirst=0
      prog=0. 
      MXEQ=0 ! Added by MBattagin - it was not initialised!!
      
      do ia=1,nAnis
        is=seqsire(ia)
        idd=seqdam(ia)
        p1(ia)=is
        p2(ia)=idd
        call LNKLST(is,idd,ia,1,MXEQ,nAnis,mm,nn,mate,next,ifirst,prog)  ! from MXEQ added by MBattagin
        call LNKLST(idd,is,ia,0,MXEQ,nAnis,mm,nn,mate,next,ifirst,prog)  ! from MXEQ added by MBattagin       
      end do
    end subroutine CreateLinkListArrays

    !######################################################################################################################################################

    subroutine SetUpEquationsForSnp(nAnis,Seq0Snp1Mode,nReadCounts,GMatSnp,GMatRds,ErrorRate,MaxReadCounts)

        implicit none


        integer,intent(in) :: nAnis,Seq0Snp1Mode
        integer,intent(inout) :: MaxReadCounts
        real(kind=8),intent(in) :: ErrorRate


        real(kind=8),intent(in),dimension(:) :: nReadCounts(nAnis)
        real(kind=8),intent(inout),allocatable,dimension(:,:) :: GMatSnp 
        real(kind=8),intent(inout),allocatable,dimension(:,:,:) :: GMatRds


        integer :: i,k

        real(kind=8) :: ProdFactTmp !DFactorial,DFactorialInLog,
        

        if (Seq0Snp1Mode==0) then

            MaxReadCounts=maxval(nReadCounts)                       

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

    subroutine SetUpData(Seq0Snp1Mode,ReadCounts,InputGenos,nAnis,nSnp,EndSnp,StartSnp,nReadCounts,ReadCountsTmp,InputGenosTmp)

        implicit none

        integer, intent(in) :: nAnis,nSnp,StartSnp,EndSnp,Seq0Snp1Mode
        
        real(kind=8),intent(in),dimension(:,:,:) :: ReadCountsTmp(:,:,:) !(nAnis,nSnp,2)
        integer,intent(inout),dimension(:,:) :: InputGenosTmp(:,:)

        real(kind=8),intent(inout),allocatable,dimension(:,:,:) :: ReadCounts
        real(kind=8),intent(inout),allocatable,dimension(:) :: nReadCounts
        integer,intent(inout),allocatable,dimension(:,:) :: InputGenos

        real(kind=8),allocatable,dimension(:) :: nReadCountsTmp

        integer :: i,j,k

        if (Seq0Snp1Mode==0) then

            allocate(ReadCounts(nAnis,EndSnp-StartSnp+1,2)) ! Commented by MBattagin
            allocate(nReadCounts(nAnis))                    ! Commented by MBattagin
            allocate(nReadCountsTmp(EndSnp-StartSnp+1))    ! For each Individual save the reads for all the SNP

            do i=1,nAnis
                nReadCountsTmp=0
                k=1
                do j=StartSnp,EndSnp
                    ReadCounts(i,k,:)=ReadCountsTmp(i,j,:)
                    !nReadCounts(i)=sum(ReadCounts(i,k,:))
                    nReadCountsTmp(k)=sum(ReadCounts(i,k,:))
                    k=k+1
                end do
                nReadCounts(i)=maxval(nReadCountsTmp(:))
            end do

            deallocate(nReadCountsTmp)

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

        implicit none

     
        integer,intent(in) :: nAnis
        integer,intent(inout) :: MaxFs,MaxMates
        integer, intent(in), dimension (:) :: SeqSire(nAnis),SeqDam(nAnis)

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
        ! IF(family(i) /= 0) PRINT'(3i7,i15)', i, seqsire(i), seqdam(i), family(i)
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

    subroutine geneprob(currentSnp,ErrorRate,nAnis,Seq0Snp1Mode,ReadCounts,InputGenos,nReadCounts,EndSnp,StartSnp, &
                        maxfs,MaxMates,MaxReadCounts,GMatSnp,GMatRds,OutputMaf,SeqId,SeqSire,SeqDam, &
                        Pr00,Pr01,Pr10,Pr11, &
                        mxeq,mate,ifirst,next,prog,p1,p2)

      implicit none

      integer, intent(in)   :: currentSnp,Seq0Snp1Mode,EndSnp,StartSnp,maxfs,MaxMates,MaxReadCounts
      integer, intent(in) :: nAnis
      integer, intent(in), dimension (:) :: SeqId(nAnis),SeqSire(nAnis),SeqDam(nAnis)

      real(kind=8),intent(in) :: ErrorRate
      real(kind=8),intent(in),dimension(:) :: nReadCounts                                          

      real(kind=8),intent(in),dimension(:,:) :: GMatSnp(1:3,1:3)
      real(kind=8),intent(in),dimension(:,:,:) :: GMatRds(0:MaxReadCounts,3,MaxReadCounts)

      integer,intent(in),dimension(:,:) :: InputGenos
      real(kind=8),intent(in),dimension(:,:,:) :: ReadCounts(nAnis,EndSnp-StartSnp+1,2)                                       

      real(kind=4),intent(inout),dimension(:,:) :: Pr00(nAnis,EndSnp-StartSnp+1),Pr01(nAnis,EndSnp-StartSnp+1),Pr10(nAnis,EndSnp-StartSnp+1),Pr11(nAnis,EndSnp-StartSnp+1)
      real(kind=8),intent(inout),dimension(:) :: OutputMaf(EndSnp-StartSnp+1)

      integer,intent(in) :: mxeq

      INTEGER, intent(inout),dimension(:) :: mate(0:2*nAnis),prog(0:2*nAnis)
      INTEGER, intent(inout),dimension(:) :: next(2*nAnis),ifirst(0:nAnis)
      INTEGER, intent(inout),dimension(:) :: p1,p2 


      REAL (KIND=8)         :: pprior_hold,qprior_hold,StopCrit_hold         ! Added by MBattagin

      INTEGER :: Imprinting_hold,PauseAtEnd_hold,nobs_hold,nfreq_max_hold     ! Added by MBattagin

      REAL (KIND=8) :: pprior, qprior, StopCrit !SeqError                      ! Added by MBattagin
      
      INTEGER :: Imprinting,PauseAtEnd,nfreq_max !phenotypes                   ! Added by MBattagin
      INTEGER, allocatable:: phenhold(:)                                      ! Added by MBattagin

      INTEGER :: MM                                                           ! Added by MBattagin - used here and in LNKLST and ANotherOne
      INTEGER :: NN                                                           ! Added by MBattagin - used here and in LNKLST   
      
      REAL(KIND=8), ALLOCATABLE :: POST(:,:)                                  ! Added by MBattagin - used here and in FLIPPT

      
      integer                 :: i, j, k, l, i2, i3,  iticks2,ifix !iticks1,
      integer                 :: ia, is, idd, ifreq_iterate, maxint, maxiter, itersused, kl, kc, kd, kj, nfams, last
      integer                 :: nf, im, ns, mf, iaa, ii, ms, m, n, maxvalspost, ierrors, iflag, nwritten
      integer                 :: f,ff
      integer                 :: maxRegpoints,LeastPositive, LeastNegative, HoldInt, LimitAnimals, LimitNumber
      integer                 :: nonzed(3,3)
      integer                 :: ntype(3,3,3)

      real (kind=8)                 :: tsum, prod, IMPratio, p12, p21, LeastPositiveValue, LeastNegativeValue ! ProbFit,SumFreq,
      REAL (KIND=8)                 :: spost(3),dpost(3),fpost(3),tpost(3),temp(3),sum1(3),sum2(3),sum3(3)
      REAL (KIND=8)                 :: pt(3,3,3)
      REAL (KIND=8)                 :: phethw,phomhw  ! this is needed for info - or compile with dble.  Don't know why!
      REAL (KIND=8)                 :: s0,s1,s2, areg, breg, meanX, meanY, sumY, sumXY, sumX, sumX2
      real                          :: LnkStart,LnkEnd

      INTEGER, allocatable,dimension(:)  :: phen,nmem,damGP ! note this is a different p1,p2 to sequence's
      INTEGER, ALLOCATABLE,dimension(:,:):: isib

      REAL (KIND=8), allocatable,dimension(:)       :: phom,phet,pnor,pHold,pResult,pDev,sumReadsCurrentSnp
      REAL (KIND=8), allocatable,dimension(:,:)     :: ant,term,freq
      REAL (KIND=8), allocatable,dimension(:,:,:)   :: work

      !JH TO FIX THESE UP LATER
      pprior_hold = 0.5
      qprior_hold = 1-pprior_hold ! MBattagin "pprior_hold" was "pprior"

      !phenotypes_hold = 3
      Imprinting_hold = 1
      PauseAtEnd_hold = 0
      nfreq_max_hold = 50
      StopCrit_hold = 0.0001
      !JH TO FIX THESE UP LATER

      pprior=pprior_hold
      qprior=qprior_hold
      nfreq_max=nfreq_max_hold
      Imprinting=Imprinting_hold
      PauseAtEnd=PauseAtEnd_hold
      StopCrit=StopCrit_hold

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

      allocate(phenhold(0:nAnis))
      ALLOCATE(post(3,0:2*nAnis),phen(nAnis), &
              ant(3,0:nAnis),phom(0:nAnis),phet(0:nAnis), &
              freq(3,0:nAnis),pnor(0:nAnis))

      if (Seq0Snp1Mode==0) then
          phenhold=0
          phen=0
      endif
      if (Seq0Snp1Mode==1) then
          phenhold=9
          phen=9 ! covers unlisted parents
      endif

      do i=1,nAnis
        if (Seq0Snp1Mode==0) then
            phenhold(i) = ReadCounts(i,currentSnp,2) !!!! was ReadCounts(i,1,2) - MBattagin
        endif
        if (Seq0Snp1Mode==1) then
            phenhold(i) = InputGenos(i,currentSnp)
        endif
      end do

      post=0. 
      phom=0.
      phet=0.
      pnor=0.
      freq=log(1.)
      HoldInt = MAX(2*maxfs,maxmates)
      
      ALLOCATE(nmem(HoldInt),isib(maxmates,2*maxfs),damGP(maxmates),work(3,2*maxfs,2*maxfs),term(3,HoldInt))

      HoldInt=0
      isib=0

      if (Seq0Snp1Mode==0) allocate(sumReadsCurrentSnp(nAnis))

      do ia=1,nAnis ! THIS LOOP CONVERT THE INPUT DATA IN LOG-LIKELIHOOD
        
        phen(ia)=phenhold(ia)   ! was phen(ia)=phenhold(passedorder(ia)) -- MBattagin
        iflag=0
        
        if (Seq0Snp1Mode==0) then
            sumReadsCurrentSnp(ia)=sum(ReadCounts(ia,currentSnp,:))
            if (sumReadsCurrentSnp(ia).eq.0) THEN
                iflag=1
                freq(:,ia) =log(1.)
            else
                do i = 0, sumReadsCurrentSnp(ia)
                    IF (phen(ia).eq.i) THEN
                        iflag=1
                        IF (GMatRds(i, 1,sumReadsCurrentSnp(ia)).lt.log(.000000001))then 
                            freq(1,ia) =-9999
                        else
                            freq(1,ia) =GMatRds(i, 1,sumReadsCurrentSnp(ia))
                        endif
                        IF(GMatRds(i, 2,sumReadsCurrentSnp(ia)).lt.log(.000000001))then
                            freq(2,ia) =-9999
                        else
                            freq(2,ia)=GMatRds(i, 2,sumReadsCurrentSnp(ia))
                        endif
                        IF(GMatRds(i, 3,sumReadsCurrentSnp(ia)).lt.log(.000000001))then
                            freq(3,ia) =-9999
                        else
                            freq(3,ia) =GMatRds(i, 3,sumReadsCurrentSnp(ia))
                        endif
                    endif
                enddo
            endif
        endif

        if (Seq0Snp1Mode==1) then
            if (phen(ia).eq.9) then
                iflag=1
                freq(1,ia) =log(1.)
                freq(2,ia) =log(1.)
                freq(3,ia) =log(1.)
            else
                do i = 0,2
                    ifix=i+1
                    if (phen(ia).eq.i) then
                        iflag=1
                        
                        if (GMatSnp(ifix,1).lt.(.000000001)) then ! MBattagin the "1" was "0", make sure it doesn't create problems
                            freq(1,ia) =-9999
                        else
                            freq(1,ia) =log(GMatSnp(ifix, 1)) ! MBattagin the "1" was "0", make sure it doesn't create problems
                        endif
                        if (GMatSnp(ifix,2).lt.(.000000001)) then ! MBattagin the "2" was "1", make sure it doesn't create problems
                            freq(2,ia) =-9999
                        else
                            freq(2,ia) =log(GMatSnp(ifix, 2)) ! MBattagin the "2" was "1", make sure it doesn't create problems
                        endif
                        if (GMatSnp(ifix,3).lt.(.000000001)) then ! MBattagin the "3" was "2", make sure it doesn't create problems
                            freq(3,ia) =-9999
                        else
                            freq(3,ia) =log(GMatSnp(ifix, 3)) ! MBattagin the "3" was "2", make sure it doesn't create problems
                        endif
                    endif
                enddo
            endif
        endif
      end do

      deallocate (phenhold)
      if (Seq0Snp1Mode==0) deallocate(sumReadsCurrentSnp)

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
        !print*,currentSnp,ifreq_iterate,pHold(ifreq_iterate)
        !   initialise
        post=0.
        phet=0.
        do i=1,nAnis
          ant(1,i)=log(.000000001)!log(qprior*qprior) 
          ant(2,i)=log(.000000001)!log(2.0*pprior*qprior) 
          ant(3,i)=log(.000000001)!log(pprior*pprior) 
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
                    ! print*, currentSnp,ia,i,ant(i,ia)
                  enddo
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

        if(p1(i).eq.0 .and. p2(i).eq.0)then
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

          iflag=0 

          !!!MBattagin - the following 8 lines in GeneProb4AlphaImpute are not commented!
          !           if (ABS(1. - pnor(i)-phet(i)-phom(i)) > .0001 ) iflag=iflag+10

          !           do j=1, phenotypes
          !               if (phen(i)==phenotype(j)) then
          !                   if (g(j,0)<0.000001 .and. pnor(i)>0.01 ) iflag=iflag+1
          !                   if (g(j,1)<0.000001 .and. phet(i)>0.01 ) iflag=iflag+1
          !                   if (g(j,2)<0.000001 .and. phom(i)>0.01 ) iflag=iflag+1
          !               endif
          !           enddo
         
          if(Imprinting>0) then
            if(phet(i)<0.0000001) then
              Pr00(i,currentSnp) = pnor(i)
              Pr01(i,currentSnp) = phet(i)
              Pr10(i,currentSnp) = phet(i)
              Pr11(i,currentSnp) = phom(i)
            else
              p12= (pnor(seqsire(i))+0.5*phet(seqsire(i))) * (phom( seqdam(i))+0.5*phet( seqdam(i)))  ! extra safe due to the above
              p21= (pnor( seqdam(i))+0.5*phet( seqdam(i))) * (phom(seqsire(i))+0.5*phet(seqsire(i)))
                if(p12+p21>0.00000001) then
                  IMPratio= p12 / (p12+p21)
                else
                  IMPratio= 0.5 ! neutral but should not be invked anyway
                endif

                if(Imprinting==2)then
                  Pr00(i,currentSnp) = pnor(i)
                  Pr01(i,currentSnp) = phet(i)
                  Pr10(i,currentSnp) = (1.-2.*IMPratio)*phet(i)
                  Pr11(i,currentSnp) = phom(i)
                else
                  Pr00(i,currentSnp) = pnor(i)
                  Pr01(i,currentSnp) = IMPratio*phet(i)
                  Pr10(i,currentSnp) = (1.-IMPratio)*phet(i)
                  Pr11(i,currentSnp) = phom(i)
                endif
            endif
          else
            Pr00(i,currentSnp) = pnor(i)
            Pr01(i,currentSnp) = phet(i)
            Pr10(i,currentSnp) = phet(i)
            Pr11(i,currentSnp) = phom(i)
          endif
        enddo

        OutputMaf(currentSnp)=pprior
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

end module AlphaVarCallFuture


