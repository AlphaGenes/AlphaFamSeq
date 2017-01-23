!###########################################################################################################################################################

module Ferdosi

	! AUTHOR : SGONEN 18.02.2016

	implicit none

	public :: FerdosiRunner

	contains

	!###########################################################################################################################################################
		subroutine FerdosiRunner(nHalfSib, nSnp, ParentGender, ParentPhase, HalfSibGeno, HalfSibPhase)

			implicit none

			integer, intent(in) :: nHalfSib, nSnp, ParentGender
			integer, intent(inout), dimension (:,:) :: HalfSibGeno(nHalfSib, nSnp), ParentPhase(nSnp, 2)
			integer, intent(inout), dimension(:,:,:) :: HalfSibPhase(nHalfSib,nSnp,2) !Phase arrays are dimension(nIndividuals, nSnp, 2) for phase 1 (paternal) and phase 2 (maternal)

			integer, allocatable, dimension(:) :: OpposingHomoPositions
			integer, allocatable, dimension(:,:) :: HalfSibRecombPos, HalfSibParentPhaseCode

			allocate(HalfSibRecombPos(nHalfsib, nSnp))
			allocate(OpposingHomoPositions(nSnp))
			allocate(HalfSibParentPhaseCode(nHalfSib, nSnp)) !1 is Pat Phase for common parent, 0 is Mat Phase for common parent

			OpposingHomoPositions = 0
			HalfSibParentPhaseCode = -999 
			HalfSibRecombPos = 0

			!Part 1 of Ferdosi: Paternal strand detection based on opposing homozygous markers in half-sib families
			
			!Part 2 of Ferdosi: Phasing and imputation of sire genotypes

			!Part 3 of Ferdosi: Phasing of genotypes of half-sib families

			call GetOpposingHomozygoteSnp(nHalfSib,nSnp,HalfSibGeno,OpposingHomoPositions)

			call DefineCommonParentPhase(nHalfSib,nSnp,ParentPhase,HalfSibGeno,OpposingHomoPositions,HalfSibRecombPos,HalfSibParentPhaseCode)

			call InferCommonParentHaplotypes(nHalfSib,nSnp,ParentPhase,HalfSibGeno,HalfSibParentPhaseCode, OpposingHomoPositions)
			
			call DefineRecombLocations(nHalfSib, nSnp, HalfSibRecombPos, ParentPhase, HalfSibGeno, HalfSibParentPhaseCode)

			call PhaseHalfSibs(nHalfSib,nSnp,ParentGender,HalfSibPhase,ParentPhase,HalfSibParentPhaseCode,HalfSibGeno)


			deallocate(HalfSibRecombPos)
			deallocate(OpposingHomoPositions)
			deallocate(HalfSibParentPhaseCode)

		end subroutine FerdosiRunner

	!###########################################################################################################################################################

		subroutine GetOpposingHomozygoteSnp(nHalfSib,nSnp,HalfSibGeno,OpposingHomoPositions)
			
			implicit none

			integer, intent(in) :: nHalfSib, nSnp
			integer, intent(inout), dimension(:) :: OpposingHomoPositions(nSnp)
			integer, intent(in), dimension(:,:,:) :: HalfSibGeno(nHalfSib,nSnp) 

			integer :: i
			integer, allocatable, dimension(:) :: TempGenos

			allocate(TempGenos(nHalfSib))

			TempGenos = 9
			!print*,"OpposingHomoPositions"
			do i=1, nSnp
				TempGenos = 9
				TempGenos = HalfSibGeno(:,i)
				if((ANY(TempGenos==0)) .and. (ANY(TempGenos==2))) then
					OpposingHomoPositions(i) = 1
					!print*,i
				endif
			enddo

			deallocate(TempGenos)

		end subroutine GetOpposingHomozygoteSnp

	!###########################################################################################################################################################

		subroutine DefineCommonParentPhase(nHalfSib,nSnp,ParentPhase,HalfSibGeno,OpposingHomoPositions,HalfSibRecombPos,HalfSibParentPhaseCode)

			implicit none

			integer, intent(in) :: nHalfSib, nSnp
			integer, intent(in), dimension(:) :: OpposingHomoPositions(nSnp)
			integer, intent(inout), dimension(:,:) :: ParentPhase(nSnp, 2), HalfSibRecombPos(nHalfsib, nSnp), HalfSibParentPhaseCode(nHalfSib,nSnp)
			integer, intent(in), dimension(:,:) :: HalfSibGeno(nHalfSib,nSnp)


			integer :: i, FirstSnp, PrevSnp
			integer, allocatable, dimension(:) :: HalfSibGenoAtSnp
			integer, allocatable, dimension(:,:) :: ParentPhaseFmv, HalfSibFmvAtSnp !1 is Pat Phase for common parent, 0 is Mat Phase for common parent

			FirstSnp = 1
			PrevSnp = 0

			allocate(ParentPhaseFmv(2,2)) !first is snp, second is phase
			allocate(HalfSibGenoAtSnp(nHalfSib))
			allocate(HalfSibFmvAtSnp(nHalfSib, 2))

			HalfSibFmvAtSnp = -999
			HalfSibParentPhaseCode = -999 
			HalfSibRecombPos = 0
			HalfSibGenoAtSnp = -999
			HalfSibFmvAtSnp = -999
			ParentPhaseFmv = -999

			do i=1, nSnp
				HalfSibGenoAtSnp = -999
				if (OpposingHomoPositions(i) == 1) then
					HalfSibGenoAtSnp(:) =  HalfSibGeno(:,i)
					call DetermineFmvAndSirePhase(nHalfSib, FirstSnp, HalfSibGenoAtSnp, HalfSibFmvAtSnp, ParentPhaseFmv)
					if ((PrevSnp/=0) .and. (ANY(ParentPhaseFmv(:,:)/=9))) then
						call DetermineHalfSibParentPhaseInhert(nHalfSib, nSnp, i, PrevSnp, HalfSibFmvAtSnp, ParentPhaseFmv, HalfSibParentPhaseCode, HalfSibRecombPos)
					endif
					ParentPhase(i,:) = ParentPhaseFmv(2,:)
					if (ANY(ParentPhaseFmv(2,:)==9)) ParentPhaseFmv(2,:) = ParentPhaseFmv(1,:)
					FirstSnp = 0
					PrevSnp = i
				endif
			enddo 

			! print*, HalfSibParentPhaseCode(18,1180:1700)
			! print*, ""

			call FillAll(nHalfSib, nSnp, HalfSibParentPhaseCode)

			! print*, HalfSibParentPhaseCode(18,1180:1700)


			! stop

			deallocate(ParentPhaseFmv)
			deallocate(HalfSibGenoAtSnp)
			deallocate(HalfSibFmvAtSnp)

		end subroutine DefineCommonParentPhase

	!###########################################################################################################################################################

		subroutine DetermineFmvAndSirePhase(nHalfSib, FirstSnp, HalfSibGenoAtSnp, HalfSibFmvAtSnp, ParentPhaseAtSnp)

			implicit none

			integer, intent(in) :: nHalfSib, FirstSnp
			integer, intent(inout), dimension(:) :: HalfSibGenoAtSnp(nHalfSib)
			integer, intent(inout), dimension(:,:) ::  ParentPhaseAtSnp(2,2), HalfSibFmvAtSnp(nHalfSib, 2) !1 is Pat Phase for common parent, 0 is Mat Phase for common parent

			integer :: i, PhaseFillZero, PhaseFillOne, FillHere, ParentAlleleFilled
			integer, allocatable, dimension(:) :: FmvCounts

			allocate(FmvCounts(4)) ! positions are counts for 0,0 ; 1,0 ; 0,1 ; 1,1 in FMV 
			FmvCounts = 0
			FillHere = 0
			PhaseFillZero = 0
			PhaseFillOne = 0
			ParentAlleleFilled = -9

			ParentPhaseAtSnp(1,:) = ParentPhaseAtSnp(2,:)
			ParentPhaseAtSnp(2,:) = 9

			if (FirstSnp == 1) then
				PhaseFillZero = 1
				PhaseFillOne = 2
				HalfSibFmvAtSnp(:,2) = HalfSibGenoAtSnp(:)
				do i=1, nHalfSib
					if (HalfSibGenoAtSnp(i) == 0) then
						FmvCounts(1) = FmvCounts(1) + 1 
					elseif (HalfSibGenoAtSnp(i) == 2) then
						FmvCounts(4) = FmvCounts(4) + 1 
					endif

				enddo
			else
				HalfSibFmvAtSnp(:,1) = HalfSibFmvAtSnp(:,2)
				HalfSibFmvAtSnp(:,2) = HalfSibGenoAtSnp(:)
				do i=1, nHalfSib
					if ((HalfSibFmvAtSnp(i,1) == 0) .and. (HalfSibFmvAtSnp(i,2) == 0)) then
						FmvCounts(1) = FmvCounts(1) + 1
					elseif ((HalfSibFmvAtSnp(i,1) == 2) .and. (HalfSibFmvAtSnp(i,2) == 0)) then
						FmvCounts(2) = FmvCounts(2) + 1
					elseif ((HalfSibFmvAtSnp(i,1) == 0) .and. (HalfSibFmvAtSnp(i,2) == 2)) then
						FmvCounts(3) = FmvCounts(3) + 1
					elseif ((HalfSibFmvAtSnp(i,1) == 2) .and. (HalfSibFmvAtSnp(i,2) == 2)) then
						FmvCounts(4) = FmvCounts(4) + 1
					endif
				enddo

				if (ParentPhaseAtSnp(1,1) ==0) then
					PhaseFillZero = 1
					PhaseFillOne = 2
				elseif (ParentPhaseAtSnp(1,2) ==0) then
					PhaseFillZero = 2
					PhaseFillOne = 1
				endif
			endif

			FillHere = PhaseFillZero
			do i=1,2
				if (FmvCounts(i) > FmvCounts(i+2)) then
					ParentPhaseAtSnp(2,FillHere) = 0
				elseif (FmvCounts(i) < FmvCounts(i+2)) then
					ParentPhaseAtSnp(2,FillHere) = 1
				endif
				FillHere = PhaseFillOne
			enddo
			!check TODO - better way?
			if ((ParentPhaseAtSnp(2,1) ==9) .and. (ParentPhaseAtSnp(2,2) /=9 )) then
				if (ParentPhaseAtSnp(2,2) == 0) then
					ParentPhaseAtSnp(2,1) = 1
				else
					ParentPhaseAtSnp(2,1) = 0
				endif
			elseif ((ParentPhaseAtSnp(2,2) ==9) .and. (ParentPhaseAtSnp(2,1) /=9 )) then
				if (ParentPhaseAtSnp(2,1) == 0) then
					ParentPhaseAtSnp(2,2) = 1
				else
					ParentPhaseAtSnp(2,2) = 0
				endif
			elseif (ParentPhaseAtSnp(2,1) == ParentPhaseAtSnp(2,2)) then
				ParentAlleleFilled = ParentPhaseAtSnp(2,1)
				if (ParentAlleleFilled == 1) then
					if (FmvCounts(3) > FmvCounts(4)) then
						ParentPhaseAtSnp(2,PhaseFillOne) = 0 
					else
						ParentPhaseAtSnp(2,PhaseFillZero) = 0
					endif
				else
					if (FmvCounts(1) > FmvCounts(2)) then
						ParentPhaseAtSnp(2,PhaseFillOne) = 1 
					else
						ParentPhaseAtSnp(2,PhaseFillZero) = 1
					endif
				endif
			endif

			deallocate(FmvCounts)

		end subroutine DetermineFmvAndSirePhase

	!###########################################################################################################################################################

		subroutine DetermineHalfSibParentPhaseInhert(nHalfSib, nSnp, CurrSnp, PrevSnp, HalfSibFmvAtSnp, ParentPhaseFmv, HalfSibParentPhaseCode, HalfSibRecombPos)

			implicit none

			integer, intent(in) :: nHalfSib, nSnp, CurrSnp, PrevSnp
			integer, intent(in), dimension(:,:) :: HalfSibFmvAtSnp(nHalfSib, 2), ParentPhaseFmv(2,2) !ParentPhaseFmv(snp,phase)
			integer, intent(inout), dimension(:,:) :: HalfSibRecombPos(nHalfSib, nSnp)
			integer, intent(inout), dimension(:,:) :: HalfSibParentPhaseCode(nHalfSib, nSnp)

			integer :: i

			do i=1, nHalfSib

				if (ALL(HalfSibFmvAtSnp(i,:) - 2*(ParentPhaseFmv(:,1)) == 0)) then ! no rec, carries paternal gamete of parent
					HalfSibParentPhaseCode(i, PrevSnp:CurrSnp) = 1 ! Set everything between two opposing homozygous snps to be Paternal gamete
				elseif (ALL(HalfSibFmvAtSnp(i,:) - 2*(ParentPhaseFmv(:,2)) == 0)) then ! no rec, carries paternal gamete of parent
					HalfSibParentPhaseCode(i, PrevSnp:CurrSnp) = 2 ! Set everything between two opposing homozygous snps to be Maternal gamete
				elseif ((HalfSibFmvAtSnp(i,1) == 2*ParentPhaseFmv(1,1)) .and. (HalfSibFmvAtSnp(i,2) == 2*ParentPhaseFmv(2,2))) then ! recombination somewhere in between - switch from parent pat to mat gamete
					HalfSibParentPhaseCode(i, PrevSnp) = 1
					HalfSibParentPhaseCode(i, CurrSnp) = 2
					HalfSibRecombPos(i, PrevSnp:CurrSnp) = 1
				elseif ((HalfSibFmvAtSnp(i,1) == 2*ParentPhaseFmv(1,2)) .and. (HalfSibFmvAtSnp(i,2) == 2*ParentPhaseFmv(2,1))) then ! recombination somewhere in between - switch from parent mat to pat gamete
					HalfSibParentPhaseCode(i, PrevSnp) = 2
					HalfSibParentPhaseCode(i, CurrSnp) = 1
					HalfSibRecombPos(i, PrevSnp:CurrSnp) = 1		
	
				endif

				
			enddo

		end subroutine DetermineHalfSibParentPhaseInhert

	!###########################################################################################################################################################

		subroutine FillAll(nHalfSib, nSnp, HalfSibParentPhaseCode)

			implicit none

			integer, intent(in) :: nHalfSib, nSnp
			integer, intent(inout), dimension(:,:) :: HalfSibParentPhaseCode(nHalfSib, nSnp)


			integer :: i, j, CurrFilled, PrevFilled

			CurrFilled = 0
			PrevFilled = 0

			do i=1, nHalfSib

				do j=1, nSnp
					if (HalfSibParentPhaseCode(i,j) /= -999) then
						PrevFilled = CurrFilled
						CurrFilled = j
						if (PrevFilled /=0) then
							if (HalfSibParentPhaseCode(i,PrevFilled) == HalfSibParentPhaseCode(i,CurrFilled)) then
								if ((CurrFilled - PrevFilled) < 100) then ! 400 is this the best way to go? 
									HalfSibParentPhaseCode(i,PrevFilled:CurrFilled) = HalfSibParentPhaseCode(i,CurrFilled)
								endif
							endif
						endif
					endif
				enddo
			enddo

		end subroutine FillAll


	!###########################################################################################################################################################

		subroutine InferCommonParentHaplotypes(nHalfsib,nSnp,ParentPhase,HalfSibGeno,HalfSibParentPhaseCode,OpposingHomoPositions)

			implicit none

			integer, intent(in) :: nHalfsib,nSnp
			integer, intent(in), dimension(:) :: OpposingHomoPositions(nSnp)
			integer, intent(in), dimension(:,:) :: HalfSibGeno(nHalfsib,nSnp)
			integer, intent(inout), dimension(:,:) :: ParentPhase(nSnp, 2), HalfSibParentPhaseCode(nHalfsib,nSnp)

			integer :: i, j
			integer :: ParentP1Count, ParentP1Indivcount, ParentP2Count, ParentP2Indivcount
			integer :: P1Average, P2Average


			do i=1, nSnp
				if (OpposingHomoPositions(i) == 0) then
					ParentP1Count = 0
					ParentP1Indivcount = 0
					ParentP2Count = 0
					ParentP2Indivcount = 0
					do j=1, nHalfsib
						if (HalfSibGeno(j,i) /= 9) then
							if (HalfSibParentPhaseCode(j,i) == 1) then
								ParentP1Count = ParentP1Count + HalfSibGeno(j,i)
								ParentP1Indivcount = ParentP1Indivcount + 1
							elseif (HalfSibParentPhaseCode(j,i) == 2) then
								ParentP2Count = ParentP2Count + HalfSibGeno(j,i)
								ParentP2Indivcount = ParentP2Indivcount + 1
							endif
						endif
					enddo

					P1Average = nint(float(ParentP1Count) / float(2*ParentP1Indivcount))
					P2Average = nint(float(ParentP2Count) / float(2*ParentP2Indivcount))


					if (P1Average < 0) then
						P1Average = 0
					endif
					if (P2Average < 0) then
						P2Average = 0
					endif
					if (ParentP1Indivcount <2) then
						P1Average = 9
					endif
					if (ParentP2Indivcount <2) then
						P2Average = 9
					endif
					ParentPhase(i,1) = P1Average !TODO checker if previously assigned genotype and this phase match (if geno avail)
					ParentPhase(i,2) = P2Average
				endif
			enddo
		end subroutine InferCommonParentHaplotypes

	!###########################################################################################################################################################
		
		subroutine DefineRecombLocations(nHalfSib, nSnp, HalfSibRecombPos, ParentPhase, HalfSibGeno, HalfSibParentPhaseCode)

			implicit none

			integer, intent(in) :: nHalfSib, nSnp
			integer, intent(in), dimension (:,:) :: HalfSibGeno(nHalfSib, nSnp), ParentPhase(nSnp, 2)
			integer, intent(inout), dimension(:,:) :: HalfSibRecombPos(nHalfSib, nSnp), HalfSibParentPhaseCode(nHalfSib, nSnp)

			integer :: i
			integer, allocatable, dimension(:) :: CurrHalfSibGeno, CurrHalfSibPhaseCode, CurrHalfSibRecombPos

			allocate(CurrHalfSibGeno(nSnp))
			allocate(CurrHalfSibPhaseCode(nSnp))
			allocate(CurrHalfSibRecombPos(nSnp))

			CurrHalfSibGeno = -9
			CurrHalfSibPhaseCode = -9
			CurrHalfSibRecombPos = -9

			do i=1, nHalfSib
				if (ANY(HalfSibRecombPos(i,:)==1)) then
					CurrHalfSibGeno = HalfSibGeno(i,:)
					CurrHalfSibPhaseCode = HalfSibParentPhaseCode(i,:)
					CurrHalfSibRecombPos = HalfSibRecombPos(i,:)
					call CheckRecombAgree(nSnp, ParentPhase, CurrHalfSibGeno, CurrHalfSibPhaseCode, CurrHalfSibRecombPos)
					HalfSibParentPhaseCode(i,:) = CurrHalfSibPhaseCode 
					HalfSibRecombPos(i,:) = CurrHalfSibRecombPos 
				endif

			enddo

		end subroutine DefineRecombLocations

	!###########################################################################################################################################################


		subroutine CheckRecombAgree(nSnp, ParentPhase, CurrHalfSibGeno, CurrHalfSibPhaseCode, CurrHalfSibRecombPos)


			implicit none

			integer, intent(in) :: nSnp
			integer, intent(in), dimension(:,:) :: ParentPhase(nSnp,2)
			integer, intent(inout), dimension(:) :: CurrHalfSibGeno(nSnp), CurrHalfSibPhaseCode(nSnp), CurrHalfSibRecombPos(nSnp)

			integer :: i, next, ParentPhaseAtPos, OtherPhase, switch, final
			logical :: skip, Compatible

			skip = .False.
			next = 0
			switch = 0
			OtherPhase = 0

			do i=1, nSnp
				if (CurrHalfSibRecombPos(i) == 1) then
					switch = 0
					if (skip) then
						cycle
					endif
					skip = .True.
					ParentPhaseAtPos = ParentPhase(i,CurrHalfSibPhaseCode(i))
					if (CurrHalfSibPhaseCode(i) == 1) then
						OtherPhase = 2
					else
						OtherPhase = 1
					endif
					next = i
					CurrHalfSibRecombPos(i) = 0
					final = i
					do 
						next = next+1

						if (next > nSnp) then
							exit
						endif
						if ((CurrHalfSibRecombPos(next) == 0)) then
							exit
						endif
						final = final+1 
					enddo
					if (final -i <=3 ) then
						switch = 3
					endif
					next = i 
					do 
						next = next+1
						if (next > nSnp) then
							exit
						endif
						if ((CurrHalfSibRecombPos(next) == 0) ) then
							exit
						endif
						Compatible = CheckGenoCompatible(ParentPhaseAtPos, CurrHalfSibGeno(next))
						if (Compatible) then
							CurrHalfSibPhaseCode(next) = CurrHalfSibPhaseCode(i)
							CurrHalfSibRecombPos(next) = 0
						else
							switch = switch + 1
							if (switch >= 3) then !not a valid recomb, set everything to missing
								CurrHalfSibRecombPos(i:final) = 9
								CurrHalfSibPhaseCode(i:final) = -999
								exit
							endif
							if (OtherPhase == 1) then
								OtherPhase = 2
							else
								OtherPhase = 1
							endif
							CurrHalfSibPhaseCode(next) = OtherPhase
						endif
					enddo
				else
					skip = .False.
				endif
			enddo

		end subroutine CheckRecombAgree

	!###########################################################################################################################################################

		subroutine PhaseHalfSibs(nHalfSib,nSnp,ParentGender,HalfSibPhase,ParentPhase,HalfSibParentPhaseCode,HalfSibGeno)

			implicit none

			integer, intent(in) :: nHalfSib,nSnp,ParentGender
			integer, intent(in), dimension(:,:) :: ParentPhase(nSnp,2), HalfSibParentPhaseCode(nHalfSib, nSnp)
			integer, intent(inout), dimension(:,:) :: HalfSibGeno(nHalfSib, nSnp)
			integer, intent(inout), dimension(:,:,:) :: HalfSibPhase(nHalfSib, nSnp, 2)

			integer :: i, j, k
			integer :: OtherParent, HSParentCode, HSGenotype, HSParentPhase, FillPhase
			logical :: Compatible

			OtherParent = 0
			HSParentCode = -999
			HSGenotype = 9
			HSParentPhase = -999
			Compatible = .False.

			if (ParentGender==1) then
				OtherParent=2
			else
				OtherParent=1
			endif


			do i=1, nHalfSib
				do j=1, nSnp
					HSParentCode = HalfSibParentPhaseCode(i, j)
					HSGenotype = HalfSibGeno(i,j)
					if (HSParentCode /= -999) then
						HSParentPhase = ParentPhase(j,HSParentCode)
						Compatible = CheckGenoCompatible(HSParentPhase, HSGenotype)
						if (Compatible) then
							HalfSibPhase(i,j,ParentGender) = HSParentPhase
							HalfSibPhase(i,j,OtherParent) = HSGenotype - HSParentPhase
						else !Genotyping error, set to missing
							HalfSibGeno(i,j) = 9
							! Genotype = CorrectedGenotype !TODO if we have both parents info...? 
						endif
					endif
					if (ANY(HalfSibPhase(i,j,:) == 9)) then
						if (HSGenotype == 2) then
							FillPhase = 1
						elseif (HSGenotype == 0) then
							FillPhase = 0
						else
							FillPhase = 9
						endif
						do k=1, 2
							HalfSibPhase(i,j,k) = FillPhase
						enddo
					endif

				enddo
			enddo

		end subroutine PhaseHalfSibs

	!###########################################################################################################################################################

		function CheckGenoCompatible(HSParentPhase, HSGenotype)

			implicit none

			integer, intent(in) :: HSParentPhase, HSGenotype
			logical :: CheckGenoCompatible

			CheckGenoCompatible = .False.

			if (HSParentPhase==9) then
				return
			endif

			if (HSGenotype == 1) then
				CheckGenoCompatible = .True.
			else
				if (HSGenotype == 2*HSParentPhase) then
					CheckGenoCompatible = .True.  
				endif
			endif

		end function CheckGenoCompatible

	!###########################################################################################################################################################

end module Ferdosi

