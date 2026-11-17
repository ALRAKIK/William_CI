!=====================================================================!
!                                                                     !
!   module eri_direct  [Nov 2026]  keyword direct                     !
!                                                                     !
!   Integral-direct HF and MP2 for the toroidal SAO path: the         !
!   two-electron integrals are never stored.                          !
!                                                                     !
!   direct_init  : once, from Tori1D/2D/3D instead of the ERI driver: !
!                  t-tables, shells, Schwarz bounds (and the          !
!                  prim_screen table), a copy of the normalized AOs.  !
!                                                                     !
!   direct_JK    : every SCF iteration. Same loop as the symmetry     !
!                  driver of ERI_toroidal.f90 (symmetry-unique shell  !
!                  quartets, same Schwarz screens), but every value   !
!                  that the driver would store in ERI(p1,p2,p3,p4)    !
!                  (p1 in the unit cell) is contracted at once:       !
!                     J(p1,p2) +=       E P(p3,p4)                    !
!                     K(p1,p4) += -1/2  E P(p2,p3)                    !
!                  (= hartree_exchange_potential on the stored array, !
!                  each stored element counted once).                 !
!                                                                     !
!   direct_slab  : MP2: all (a0 j2 | j3 , d R4) for the functions d of !
!                  one unit-cell shell and every cell R4 (no 8-fold   !
!                  symmetry), used by transform_ERI_direct_kMO.       !
!=====================================================================!

module eri_direct

      use classification_ERI, only: ERI_function
      implicit none

      private
      public :: direct_init, direct_JK, direct_slab, direct_free
      public :: d_nsh, d_nshuc, d_sh_first, d_sh_n, d_ready
      public :: d_Q, d_Qsh                    ! [Nov 2026] used by df_toroidal (density_fitting)

      logical, save :: d_ready = .false.
      integer, save :: d_ndp = 0, d_nBas = 0, d_nsh = 0, d_nshuc = 0, d_ncl(3) = 1
      type(ERI_function), allocatable, save :: d_AO(:)
      integer, allocatable, save :: d_sh_first(:), d_sh_n(:)
      double precision, allocatable, save :: d_Q(:,:), d_Qsh(:,:)
      double precision, save :: d_Qmax = 0.d0
      double precision, parameter :: schwarz_thresh = 1.0d-13

contains

      !-----------------------------------------------------------------!
      subroutine direct_init(ndp, nBas, AO)

        use keywords
        use unitcell_module
        use torus_init
        use clifford_tools, only : init_clifford, orbital_to_lxyz, setup_t_tables, bessel_nmax_estimate, &
                                   ps_noscreen, ps_qp

        integer           , intent(in) :: ndp, nBas
        type(ERI_function), intent(in) :: AO(nBas)

        integer :: i, j, a, b, nsh, lxyz(3), lt, ncomp, lmax, ncap, SI, SJ, mp, p1, p2
        double precision :: emax, amin, a_dir(3)
        logical :: same
        double precision, allocatable :: blk(:,:,:,:)
        type(ERI_function), allocatable :: A1(:), B1(:)

        d_ndp = ndp ; d_nBas = nBas
        if (allocated(d_AO)) deallocate(d_AO)
        allocate(d_AO(nBas))
        d_AO = AO

        call init_clifford()
        lmax = 0 ; emax = 0.d0
        do i = 1 , nBas
          call orbital_to_lxyz(AO(i)%orbital, lxyz)
          lmax = max(lmax, maxval(lxyz(1:ndp)))
          emax = max(emax, maxval(AO(i)%exponent))
        end do
        a_dir = (/ ax, ay, az /)
        amin  = minval(a_dir(1:ndp))
        ncap = bessel_nmax_estimate(4.d0*emax/(amin*amin)) + 2*lmax + 1
        call setup_t_tables(ncap, ndp)

        ! shells (same rule as the driver)
        if (allocated(d_sh_first)) deallocate(d_sh_first, d_sh_n)
        allocate(d_sh_first(nBas), d_sh_n(nBas))
        nsh = 0 ; i = 1
        do while (i <= nBas)
          call orbital_to_lxyz(AO(i)%orbital, lxyz)
          lt = sum(lxyz) ; ncomp = (lt+1)*(lt+2)/2
          same = (i + ncomp - 1 <= nBas)
          if (same) then
            do j = i+1 , i+ncomp-1
              call orbital_to_lxyz(AO(j)%orbital, lxyz)
              same = same .and. sum(lxyz) == lt .and. AO(j)%x == AO(i)%x .and. AO(j)%y == AO(i)%y &
                     .and. AO(j)%z == AO(i)%z .and. size(AO(j)%exponent) == size(AO(i)%exponent)
              if (same) same = all(AO(j)%exponent == AO(i)%exponent)
            end do
          end if
          if (.not. same) ncomp = 1
          nsh = nsh + 1
          d_sh_first(nsh) = i ; d_sh_n(nsh) = ncomp
          i = i + ncomp
        end do
        d_nsh = nsh

        d_ncl = (/ max(nx,1), max(ny,1), max(nz,1) /)
        d_nshuc = 0
        do SI = 1 , nsh
          if (d_sh_first(SI) <= nfuc) d_nshuc = d_nshuc + 1
        end do
        if (nsh /= d_nshuc*d_ncl(1)*d_ncl(2)*d_ncl(3) .or. nBas /= nfuc*d_ncl(1)*d_ncl(2)*d_ncl(3)) &
          error stop 'keyword direct: the shells must be the same in every cell'
        do SI = 1 , nsh - d_nshuc
          if (d_sh_n(SI+d_nshuc) /= d_sh_n(SI) .or. d_sh_first(SI+d_nshuc) /= d_sh_first(SI) + nfuc) &
            error stop 'keyword direct: the shells must be the same in every cell'
        end do

        ! Schwarz bounds (diagonal shell-pair blocks)
        if (allocated(d_Q)) deallocate(d_Q, d_Qsh)
        allocate(d_Q(nBas,nBas), d_Qsh(nsh,nsh))
        d_Q = 0.d0
        !$omp parallel do default(shared) private(SI,SJ,blk,a,b,i,j) schedule(dynamic,1)
        do SI = 1 , nsh
          do SJ = SI , nsh
            allocate(blk(d_sh_n(SI),d_sh_n(SJ),d_sh_n(SI),d_sh_n(SJ)))
            call ERI_shell_block_toroidal(ndp,d_sh_n(SI),d_AO(d_sh_first(SI):),d_sh_n(SJ),d_AO(d_sh_first(SJ):), &
                                          d_sh_n(SI),d_AO(d_sh_first(SI):),d_sh_n(SJ),d_AO(d_sh_first(SJ):),blk)
            do a = 1 , d_sh_n(SI)
              do b = 1 , d_sh_n(SJ)
                i = d_sh_first(SI)+a-1 ; j = d_sh_first(SJ)+b-1
                d_Q(i,j) = dsqrt(dabs(blk(a,b,a,b)))
                d_Q(j,i) = d_Q(i,j)
              end do
            end do
            deallocate(blk)
          end do
        end do
        !$omp end parallel do
        do SI = 1 , nsh
          do SJ = 1 , nsh
            d_Qsh(SI,SJ) = maxval(d_Q(d_sh_first(SI):d_sh_first(SI)+d_sh_n(SI)-1, &
                                      d_sh_first(SJ):d_sh_first(SJ)+d_sh_n(SJ)-1))
          end do
        end do
        d_Qmax = maxval(d_Qsh)

        ! keyword prim_screen: primitive-pair table (as in the driver), kept for all passes
        if (c_prim_screen) then
          if (allocated(ps_qp)) deallocate(ps_qp)
          mp = 0
          do SI = 1 , nsh
            mp = max(mp, size(d_AO(d_sh_first(SI))%exponent))
          end do
          allocate(ps_qp(mp,mp,d_nshuc,nsh)) ; ps_qp = 0.d0
          ps_noscreen = .true.
          !$omp parallel do default(shared) private(SI,SJ,p1,p2,A1,B1,blk,a,b) schedule(dynamic,1) collapse(2)
          do SI = 1 , d_nshuc
            do SJ = 1 , nsh
              allocate(A1(d_sh_n(SI)), B1(d_sh_n(SJ)), blk(d_sh_n(SI),d_sh_n(SJ),d_sh_n(SI),d_sh_n(SJ)))
              do p1 = 1 , size(d_AO(d_sh_first(SI))%exponent)
                do p2 = 1 , size(d_AO(d_sh_first(SJ))%exponent)
                  do a = 1 , d_sh_n(SI)
                    A1(a) = d_AO(d_sh_first(SI)+a-1)
                    A1(a)%exponent    = (/ d_AO(d_sh_first(SI)+a-1)%exponent(p1) /)
                    A1(a)%coefficient = (/ d_AO(d_sh_first(SI)+a-1)%coefficient(p1) /)
                  end do
                  do b = 1 , d_sh_n(SJ)
                    B1(b) = d_AO(d_sh_first(SJ)+b-1)
                    B1(b)%exponent    = (/ d_AO(d_sh_first(SJ)+b-1)%exponent(p2) /)
                    B1(b)%coefficient = (/ d_AO(d_sh_first(SJ)+b-1)%coefficient(p2) /)
                  end do
                  call ERI_shell_block_toroidal(ndp,d_sh_n(SI),A1,d_sh_n(SJ),B1,d_sh_n(SI),A1,d_sh_n(SJ),B1,blk)
                  do a = 1 , d_sh_n(SI)
                    do b = 1 , d_sh_n(SJ)
                      ps_qp(p1,p2,SI,SJ) = max(ps_qp(p1,p2,SI,SJ), dsqrt(dabs(blk(a,b,a,b))))
                    end do
                  end do
                end do
              end do
              deallocate(A1, B1, blk)
            end do
          end do
          !$omp end parallel do
          ps_noscreen = .false.
        end if

        d_ready = .true.
      end subroutine direct_init

      subroutine direct_free()
        use clifford_tools, only : ps_qp
        if (allocated(d_AO)) deallocate(d_AO)
        if (allocated(d_sh_first)) deallocate(d_sh_first, d_sh_n)
        if (allocated(d_Q)) deallocate(d_Q, d_Qsh)
        if (allocated(ps_qp)) deallocate(ps_qp)
        d_ready = .false.
      end subroutine direct_free

      ! primitive tables of the bra and ket pairs (keyword prim_screen)
      subroutine set_ps(SI, SJ, SK, SL)
        use clifford_tools, only : ps_qp, ps_bra, ps_ket
        integer, intent(in) :: SI, SJ, SK, SL
        integer :: ab(4)
        ps_bra => ps_qp(:,:,SI,SJ)
        ab = (/ SK, SL, SK, SK /)
        call translate_to_unit_cell(ab, d_nshuc, d_ncl)
        ps_ket => ps_qp(:,:,ab(1),ab(2))
      end subroutine set_ps

      !-----------------------------------------------------------------!
      ! J and K (rows mu <= nfuc), one pass over the unique quartets     !
      !-----------------------------------------------------------------!
      subroutine direct_JK(nBas, P, J, K)

        use keywords      , only : c_prim_screen
        use unitcell_module, only : nfuc
        use clifford_tools, only : ps_on

        integer         , intent(in)  :: nBas
        double precision, intent(in)  :: P(nBas,nBas)
        double precision, intent(out) :: J(nBas,nBas), K(nBas,nBas)

        integer :: SI, SJ, SK, SL, a, b, c, d, i, jj, kk, l, iperm, pp(4), q4(4), qc(4), nt, it, ns, h, hs, mx
        integer(8) :: key
        double precision :: val
        double precision, allocatable :: blk(:,:,:,:), Jl(:,:), Kl(:,:)
        integer(8), allocatable :: hkey(:)
        integer, allocatable :: tq(:,:), used(:)
        integer :: nsh, nshuc

        if (.not. d_ready) error stop 'direct_JK: direct_init not called'
        nsh = d_nsh ; nshuc = d_nshuc
        J = 0.d0 ; K = 0.d0
        mx = maxval(d_sh_n(1:nsh))
        hs = 1
        do while (hs < 4*8*mx**4)
          hs = 2*hs
        end do

        ! task list: (SI,SJ) bra pairs, as in the driver
        nt = 0
        allocate(tq(2, nshuc*nsh))
        do SI = 1 , nshuc
          do SJ = 1 , nsh
            if (mod(SJ-1,nshuc)+1 < SI) cycle
            if (d_Qsh(SI,SJ)*d_Qmax < schwarz_thresh) cycle
            nt = nt + 1 ; tq(:,nt) = (/ SI, SJ /)
          end do
        end do

        ps_on = c_prim_screen

        !$omp parallel default(shared) &
        !$omp private(it,SI,SJ,SK,SL,blk,a,b,c,d,i,jj,kk,l,val,q4,qc,pp,iperm,Jl,Kl,hkey,ns,h,key,used)
        allocate(Jl(nfuc,nBas), Kl(nfuc,nBas), hkey(0:hs-1), used(8*mx**4))
        hkey = -1_8
        Jl = 0.d0 ; Kl = 0.d0
        !$omp do schedule(dynamic,1)
        do it = 1 , nt
          SI = tq(1,it) ; SJ = tq(2,it)
          do SK = 1 , nsh
            if (mod(SK-1,nshuc)+1 < SI) cycle
            do SL = 1 , nsh
              if (mod(SL-1,nshuc)+1 < SI) cycle
              if (d_Qsh(SI,SJ)*d_Qsh(SK,SL) < schwarz_thresh) cycle
              q4 = (/ SI, SJ, SK, SL /)
              call canonical_quartet_toroidal(q4, nsh, nshuc, d_ncl, qc)
              if (any(qc /= q4)) cycle
              allocate(blk(d_sh_n(SI),d_sh_n(SJ),d_sh_n(SK),d_sh_n(SL)))
              if (ps_on) call set_ps(SI, SJ, SK, SL)
              call ERI_shell_block_toroidal(d_ndp,d_sh_n(SI),d_AO(d_sh_first(SI):),d_sh_n(SJ),d_AO(d_sh_first(SJ):), &
                                            d_sh_n(SK),d_AO(d_sh_first(SK):),d_sh_n(SL),d_AO(d_sh_first(SL):),blk)
              ! every stored element of this block once (the driver overwrites duplicates)
              ns = 0
              do d = 1 , d_sh_n(SL)
               do c = 1 , d_sh_n(SK)
                do b = 1 , d_sh_n(SJ)
                 do a = 1 , d_sh_n(SI)
                  i  = d_sh_first(SI)+a-1 ; jj = d_sh_first(SJ)+b-1
                  kk = d_sh_first(SK)+c-1 ; l  = d_sh_first(SL)+d-1
                  if (d_Q(i,jj)*d_Q(kk,l) < schwarz_thresh) cycle
                  val = blk(a,b,c,d)
                  do iperm = 1 , 8
                    call permute_quartet(iperm, (/ i, jj, kk, l /), pp)
                    call translate_to_unit_cell(pp, nfuc, d_ncl)
                    key = int(pp(1)-1,8) + int(nfuc,8)*(int(pp(2)-1,8) + int(nBas,8)*(int(pp(3)-1,8) &
                          + int(nBas,8)*int(pp(4)-1,8)))
                    h = int(iand(key + ishft(key,-11) + ishft(key,-23) + ishft(key,-37), int(hs-1,8)))
                    do while (hkey(h) /= -1_8 .and. hkey(h) /= key)
                      h = iand(h+1, hs-1)
                    end do
                    if (hkey(h) == key) cycle
                    hkey(h) = key ; ns = ns + 1 ; used(ns) = h
                    Jl(pp(1),pp(2)) = Jl(pp(1),pp(2)) + val*P(pp(3),pp(4))
                    Kl(pp(1),pp(4)) = Kl(pp(1),pp(4)) - 0.5d0*val*P(pp(2),pp(3))
                  end do
                 end do
                end do
               end do
              end do
              hkey(used(1:ns)) = -1_8
              deallocate(blk)
            end do
          end do
        end do
        !$omp end do
        !$omp critical (direct_jk_sum)
        J(1:nfuc,:) = J(1:nfuc,:) + Jl
        K(1:nfuc,:) = K(1:nfuc,:) + Kl
        !$omp end critical (direct_jk_sum)
        deallocate(Jl, Kl, hkey, used)
        !$omp end parallel

        ps_on = .false.
        deallocate(tq)
      end subroutine direct_JK

      !-----------------------------------------------------------------!
      ! slab(a, j2, j3, dd, R4) = (a0 j2 | j3 , d R4) , d = first of the !
      ! unit-cell shell SD + dd - 1 , a in the unit cell                 !
      !-----------------------------------------------------------------!
      subroutine direct_slab(SD, nfuc, nBas, ncell, slab)

        use keywords      , only : c_prim_screen
        use clifford_tools, only : ps_on

        integer         , intent(in)  :: SD, nfuc, nBas, ncell
        double precision, intent(out) :: slab(nfuc, nBas, nBas, d_sh_n(SD), ncell)

        integer :: SA, SJ, SK, SL, R4, a, b, c, d, i, jj, kk, l, it, nt
        integer, allocatable :: tq(:,:)
        double precision, allocatable :: blk(:,:,:,:)

        slab = 0.d0
        nt = d_nshuc*d_nsh
        allocate(tq(2,nt))
        nt = 0
        do SA = 1 , d_nshuc
          do SJ = 1 , d_nsh
            nt = nt + 1 ; tq(:,nt) = (/ SA, SJ /)
          end do
        end do
        ps_on = c_prim_screen
        !$omp parallel do default(shared) private(it,SA,SJ,SK,SL,R4,blk,a,b,c,d,i,jj,kk,l) schedule(dynamic,1)
        do it = 1 , nt
          SA = tq(1,it) ; SJ = tq(2,it)
          if (d_Qsh(SA,SJ)*d_Qmax < schwarz_thresh) cycle
          do R4 = 1 , ncell
            SL = SD + (R4-1)*d_nshuc
            do SK = 1 , d_nsh
              if (d_Qsh(SA,SJ)*d_Qsh(SK,SL) < schwarz_thresh) cycle
              allocate(blk(d_sh_n(SA),d_sh_n(SJ),d_sh_n(SK),d_sh_n(SL)))
              if (ps_on) call set_ps(SA, SJ, SK, SL)
              call ERI_shell_block_toroidal(d_ndp,d_sh_n(SA),d_AO(d_sh_first(SA):),d_sh_n(SJ),d_AO(d_sh_first(SJ):), &
                                            d_sh_n(SK),d_AO(d_sh_first(SK):),d_sh_n(SL),d_AO(d_sh_first(SL):),blk)
              do d = 1 , d_sh_n(SL)
               do c = 1 , d_sh_n(SK)
                do b = 1 , d_sh_n(SJ)
                 do a = 1 , d_sh_n(SA)
                  i  = d_sh_first(SA)+a-1 ; jj = d_sh_first(SJ)+b-1
                  kk = d_sh_first(SK)+c-1 ; l  = d_sh_first(SL)+d-1
                  if (d_Q(i,jj)*d_Q(kk,l) < schwarz_thresh) cycle
                  slab(i, jj, kk, d, R4) = blk(a,b,c,d)
                 end do
                end do
               end do
              end do
              deallocate(blk)
            end do
          end do
        end do
        !$omp end parallel do
        ps_on = .false.
        deallocate(tq)
      end subroutine direct_slab

end module eri_direct
