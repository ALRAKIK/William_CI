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
!   [Nov 2026] keyword symmetry (+ direct): skeleton J/K with the     !
!   point operations (Dupuis-King). Only the shell quartets unique    !
!   under 8 permutations x translations x point operations G are      !
!   computed, each with weight |G|/|Stab(Q)|; then                    !
!      J = 1/|G| sum_g  g.J_skel ,  (g.X)(g m, g n) = s_m s_n X(m,n)  !
!   (same for K). Exact when the density P is G-invariant: this is    !
!   checked at every call (|s P(gm,gn) - P(m,n)| <= 1e-6, else the    !
!   normal direct_JK is used). Unique quartets are found from the     !
!   orbit representatives of the unit-cell shells and the bra-pair    !
!   stabilizers, without testing the whole orbit (see direct_JK_sym). !
!   A self-check (random quartets, direct vs by symmetry, 1e-10)      !
!   switches the symmetry off if it does not hold.                    !
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
      public :: set_ps                        ! [Nov 2026] used by df_toroidal (Cholesky fitting)
      public :: direct_sym_P, d_sym           ! [Nov 2026] keyword symmetry: symmetric guess density

      logical, save :: d_ready = .false.
      integer, save :: d_ndp = 0, d_nBas = 0, d_nsh = 0, d_nshuc = 0, d_ncl(3) = 1
      type(ERI_function), allocatable, save :: d_AO(:)
      integer, allocatable, save :: d_sh_first(:), d_sh_n(:)
      double precision, allocatable, save :: d_Q(:,:), d_Qsh(:,:)
      double precision, save :: d_Qmax = 0.d0
      double precision, parameter :: schwarz_thresh = 1.0d-13

      ! [Nov 2026] keyword symmetry in the direct path (tables built in dsym_init)
      logical, save :: d_sym = .false., d_sym_warned = .false.
      integer, save :: d_ncall = 0
      integer, parameter :: d_nsym_early = 8        ! SCF iterations with averaged density if needed
      integer, save :: d_nG = 1, d_ncell = 1, d_nbra = 0
      integer(8), save :: d_nkp = 0
      integer, allocatable, save :: ds_scell(:), ds_sloc(:), ds_fcell(:), ds_floc(:), ds_relc(:,:)
      integer, allocatable, save :: ds_shmap(:,:)        ! (nsh, nG) image shell
      integer, allocatable, save :: ds_rep(:)            ! (nshuc) smallest local shell of the orbit
      integer, allocatable, save :: ds_nto(:,:)          ! (a, b) number of g with loc(g a) = b
      integer, allocatable, save :: ds_to(:,:,:)         ! (k, a, b) those g
      integer, allocatable, save :: ds_kp(:,:)           ! (2, nkp) ket shell pairs SK <= SL
      integer, allocatable, save :: ds_bra(:,:)          ! (2, nbra) bra pairs (SI rep, SJ)
      integer, allocatable, save :: ds_bgo(:), ds_bg(:)  ! stabilizer of each bra pair: ds_bg(ds_bgo(b):ds_bgo(b+1)-1)

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
        integer :: ncell_t, ct, ii, jj2
        double precision :: emax, amin, a_dir(3), qv
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
        call rys_self_check(ndp, nBas, AO)        ! [Nov 2026] keyword rys_quad (no-op otherwise)

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
        ! [Nov 2026] only SI in the unit cell is computed (nshuc x nsh pairs instead of nsh^2/2);
        ! the other pairs are the same blocks shifted by a lattice vector (as in ERI_toroidal.f90)
        if (allocated(d_Q)) deallocate(d_Q, d_Qsh)
        allocate(d_Q(nBas,nBas), d_Qsh(nsh,nsh))
        d_Q = 0.d0
        ncell_t = d_ncl(1)*d_ncl(2)*d_ncl(3)
        !$omp parallel do default(shared) private(SI,SJ,blk,a,b,i,j,ct,ii,jj2,qv) schedule(dynamic,1)
        do SI = 1 , d_nshuc
          do SJ = 1 , nsh
            allocate(blk(d_sh_n(SI),d_sh_n(SJ),d_sh_n(SI),d_sh_n(SJ)))
            call ERI_shell_block_toroidal(ndp,d_sh_n(SI),d_AO(d_sh_first(SI):),d_sh_n(SJ),d_AO(d_sh_first(SJ):), &
                                          d_sh_n(SI),d_AO(d_sh_first(SI):),d_sh_n(SJ),d_AO(d_sh_first(SJ):),blk)
            do a = 1 , d_sh_n(SI)
              do b = 1 , d_sh_n(SJ)
                i = d_sh_first(SI)+a-1 ; j = d_sh_first(SJ)+b-1
                qv = dsqrt(dabs(blk(a,b,a,b)))
                do ct = 0 , ncell_t-1                   ! shift (i,j) by every lattice vector
                  ii  = shift_cell(i, ct, nfuc)
                  jj2 = shift_cell(j, ct, nfuc)
                  d_Q(ii,jj2) = qv
                  d_Q(jj2,ii) = qv
                end do
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
        if (c_symmetry) call dsym_init()          ! [Nov 2026] keyword symmetry
      end subroutine direct_init

      ! [Nov 2026] function index m (cells x fastest, nfuc per cell) shifted by the lattice vector of cell ct
      integer function shift_cell(m, ct, nfuc)
        integer, intent(in) :: m, ct, nfuc
        integer :: c, loc, cx, cy, cz
        c   = (m-1)/nfuc ; loc = mod(m-1, nfuc)
        cx  = modulo(mod(c, d_ncl(1))          + mod(ct, d_ncl(1)),          d_ncl(1))
        cy  = modulo(mod(c/d_ncl(1), d_ncl(2)) + mod(ct/d_ncl(1), d_ncl(2)), d_ncl(2))
        cz  = modulo(c/(d_ncl(1)*d_ncl(2))     + ct/(d_ncl(1)*d_ncl(2)),     d_ncl(3))
        shift_cell = (cx + d_ncl(1)*(cy + d_ncl(2)*cz))*nfuc + loc + 1
      end function shift_cell

      subroutine direct_free()
        use clifford_tools, only : ps_qp
        if (allocated(d_AO)) deallocate(d_AO)
        if (allocated(d_sh_first)) deallocate(d_sh_first, d_sh_n)
        if (allocated(d_Q)) deallocate(d_Q, d_Qsh)
        if (allocated(ps_qp)) deallocate(ps_qp)
        call dsym_free()
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
        double precision, allocatable :: Pw(:,:)

        if (.not. d_ready) error stop 'direct_JK: direct_init not called'
        ! [Nov 2026] keyword symmetry: skeleton J/K when the density is symmetric.
        ! During the first d_nsym_early SCF iterations a non-symmetric density (partly
        ! filled degenerate level of the guess / first Fock matrices, seen for LiH 3x3x3)
        ! is replaced by its symmetry average for J/K; later the exact non-symmetric
        ! J/K is used (fallback, slow) as long as the density stays non-symmetric.
        if (d_sym) then
          d_ncall = d_ncall + 1
          if (p_invariant(nBas, P)) then
            call direct_JK_sym(nBas, P, J, K)
            return
          else if (d_ncall <= d_nsym_early) then
            allocate(Pw(nBas,nBas))
            Pw = P
            call direct_sym_P(nBas, Pw, -1)
            call direct_JK_sym(nBas, Pw, J, K)
            deallocate(Pw)
            return
          end if
        end if
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


      !=================================================================!
      ! [Nov 2026] keyword symmetry: tables for direct_JK_sym           !
      !=================================================================!
      subroutine dsym_init()
        use symmetry_clifford, only : sym_ready, sym_nop, sym_perm, sym_sign, sym_eri_note
        use unitcell_module  , only : nfuc
        use files            , only : outfile
        use omp_lib          , only : omp_get_wtime
        integer :: nsh, nshuc, SI, SJ, SK, SL, ig, i, j, a, b, c, d, c1, c2, ncl(3), k, is, ntry, nsamp, u1, uj, rj
        integer, allocatable :: ao2sh(:), gl(:)
        integer(8) :: lcg, n8
        integer :: g4(4), nb, ngl
        double precision :: dev, val2, w0, w1
        double precision, allocatable :: blk(:,:,:,:), blk2(:,:,:,:)
        logical :: ok

        call dsym_free()
        d_sym = .false. ; d_ncall = 0 ; d_sym_warned = .false.
        if (.not. sym_ready) return
        if (sym_nop <= 1) return
        if (size(sym_perm,1) /= d_nBas) return
        w0 = omp_get_wtime()
        nsh = d_nsh ; nshuc = d_nshuc ; d_nG = sym_nop
        ncl = d_ncl ; d_ncell = ncl(1)*ncl(2)*ncl(3)

        ! cell tables (cells x fastest, as super_molecule)
        allocate(ds_scell(nsh), ds_sloc(nsh), ds_fcell(d_nBas), ds_floc(d_nBas), ds_relc(0:d_ncell-1,0:d_ncell-1))
        do SI = 1 , nsh
          ds_scell(SI) = (SI-1)/nshuc ; ds_sloc(SI) = mod(SI-1,nshuc)
        end do
        do i = 1 , d_nBas
          ds_fcell(i) = (i-1)/nfuc ; ds_floc(i) = mod(i-1,nfuc)
        end do
        do c1 = 0 , d_ncell-1
          do c2 = 0 , d_ncell-1
            ds_relc(c1,c2) = modulo(mod(c2,ncl(1)) - mod(c1,ncl(1)), ncl(1)) &
                           + ncl(1)*( modulo(mod(c2/ncl(1),ncl(2)) - mod(c1/ncl(1),ncl(2)), ncl(2)) &
                           + ncl(2)*  modulo(c2/(ncl(1)*ncl(2)) - c1/(ncl(1)*ncl(2)), ncl(3)) )
          end do
        end do

        ! shell maps
        allocate(ao2sh(d_nBas), ds_shmap(nsh,d_nG))
        do SI = 1 , nsh
          ao2sh(d_sh_first(SI):d_sh_first(SI)+d_sh_n(SI)-1) = SI
        end do
        ok = .true.
        do ig = 1 , d_nG
          do SI = 1 , nsh
            i = sym_perm(d_sh_first(SI), ig)
            if (i <= 0) then
              ok = .false. ; ds_shmap(SI,ig) = SI ; cycle
            end if
            SJ = ao2sh(i) ; ds_shmap(SI,ig) = SJ
            if (d_sh_n(SJ) /= d_sh_n(SI)) ok = .false.
            do a = 1 , d_sh_n(SI)
              j = sym_perm(d_sh_first(SI)+a-1, ig)
              if (j <= 0) then
                ok = .false.
              else if (ao2sh(j) /= SJ) then
                ok = .false.
              end if
            end do
          end do
        end do
        if (.not. ok) then
          write(outfile,'(a)') ' symmetry (direct): point operations NOT used (they do not map shells on shells)'
          call dsym_free() ; deallocate(ao2sh) ; return
        end if

        ! orbit representatives of the unit-cell shells and the lists g: loc(g a) = b
        allocate(ds_rep(nshuc), ds_nto(nshuc,nshuc))
        ds_nto = 0
        do a = 1 , nshuc
          ds_rep(a) = a
          do ig = 1 , d_nG
            b = ds_sloc(ds_shmap(a,ig)) + 1
            ds_rep(a) = min(ds_rep(a), b)
            ds_nto(a,b) = ds_nto(a,b) + 1
          end do
        end do
        allocate(ds_to(max(1,maxval(ds_nto)), nshuc, nshuc))
        ds_nto = 0
        do a = 1 , nshuc
          do ig = 1 , d_nG
            b = ds_sloc(ds_shmap(a,ig)) + 1
            ds_nto(a,b) = ds_nto(a,b) + 1
            ds_to(ds_nto(a,b), a, b) = ig
          end do
        end do

        ! bra pairs: SI orbit representative, SJ not reducible by the stabilizer of SI
        allocate(gl(d_nG))
        nb = 0 ; n8 = 0
        do k = 1 , 2                                    ! k = 1 count, k = 2 fill
          if (k == 2) then
            allocate(ds_bra(2,nb), ds_bgo(nb+1), ds_bg(max(1_8,n8)))
            nb = 0 ; n8 = 0 ; ds_bgo(1) = 1
          end if
          do SI = 1 , nshuc
            if (ds_rep(SI) /= SI) cycle
            do SJ = 1 , nsh
              if (ds_rep(ds_sloc(SJ)+1) < SI) cycle
              if (d_Qsh(SI,SJ)*d_Qmax < schwarz_thresh) cycle
              ngl = 0 ; ok = .true.
              do is = 1 , ds_nto(SI,SI)
                ig = ds_to(is,SI,SI)
                u1 = ds_shmap(SI,ig) ; uj = ds_shmap(SJ,ig)
                rj = ds_relc(ds_scell(u1), ds_scell(uj))*nshuc + ds_sloc(uj) + 1
                if (rj < SJ) then
                  ok = .false. ; exit
                end if
                if (rj == SJ) then
                  ngl = ngl + 1 ; gl(ngl) = ig
                end if
              end do
              if (.not. ok) cycle
              nb = nb + 1 ; n8 = n8 + ngl
              if (k == 2) then
                ds_bra(:,nb) = (/ SI, SJ /)
                ds_bg(ds_bgo(nb):ds_bgo(nb)+ngl-1) = gl(1:ngl)
                ds_bgo(nb+1) = ds_bgo(nb) + ngl
              end if
            end do
          end do
        end do
        d_nbra = nb
        deallocate(gl)

        ! ket pairs SK <= SL above the Schwarz threshold
        n8 = 0
        do SK = 1 , nsh
          do SL = SK , nsh
            if (d_Qsh(SK,SL)*d_Qmax >= schwarz_thresh) n8 = n8 + 1
          end do
        end do
        allocate(ds_kp(2,max(1_8,n8)))
        n8 = 0
        do SK = 1 , nsh
          do SL = SK , nsh
            if (d_Qsh(SK,SL)*d_Qmax >= schwarz_thresh) then
              n8 = n8 + 1 ; ds_kp(:,n8) = (/ SK, SL /)
            end if
          end do
        end do
        d_nkp = n8

        ! self-check: random quartets computed directly and as the image of an operation
        dev = 0.d0 ; nsamp = 0
        !$omp parallel do default(shared) private(ig,is,ntry,lcg,SI,SJ,SK,SL,g4,blk,blk2,a,b,c,d,val2) &
        !$omp reduction(max:dev) reduction(+:nsamp) schedule(dynamic,1)
        do ig = 2 , d_nG
          is = 0 ; ntry = 0 ; lcg = 12345_8 + 7919_8*ig
          do while (is < merge(4, merge(2, 1, d_nG <= 48), d_nG <= 16) .and. ntry < 400)
            ntry = ntry + 1
            lcg = mod(lcg*1103515245_8 + 12345_8, 2147483648_8) ; SI = 1 + int(mod(lcg, int(nshuc,8)))
            lcg = mod(lcg*1103515245_8 + 12345_8, 2147483648_8) ; SJ = 1 + int(mod(lcg, int(nsh,8)))
            lcg = mod(lcg*1103515245_8 + 12345_8, 2147483648_8) ; SK = 1 + int(mod(lcg, int(nsh,8)))
            lcg = mod(lcg*1103515245_8 + 12345_8, 2147483648_8) ; SL = 1 + int(mod(lcg, int(nsh,8)))
            if (d_Qsh(SI,SJ)*d_Qsh(SK,SL) < 1.d-8) cycle
            is = is + 1 ; nsamp = nsamp + 1
            g4 = (/ ds_shmap(SI,ig), ds_shmap(SJ,ig), ds_shmap(SK,ig), ds_shmap(SL,ig) /)
            allocate(blk(d_sh_n(SI),d_sh_n(SJ),d_sh_n(SK),d_sh_n(SL)), blk2(d_sh_n(SI),d_sh_n(SJ),d_sh_n(SK),d_sh_n(SL)))
            call ERI_shell_block_toroidal(d_ndp,d_sh_n(SI),d_AO(d_sh_first(SI):),d_sh_n(SJ),d_AO(d_sh_first(SJ):), &
                                          d_sh_n(SK),d_AO(d_sh_first(SK):),d_sh_n(SL),d_AO(d_sh_first(SL):),blk)
            call ERI_shell_block_toroidal(d_ndp,d_sh_n(g4(1)),d_AO(d_sh_first(g4(1)):),d_sh_n(g4(2)),d_AO(d_sh_first(g4(2)):), &
                                          d_sh_n(g4(3)),d_AO(d_sh_first(g4(3)):),d_sh_n(g4(4)),d_AO(d_sh_first(g4(4)):),blk2)
            do d = 1 , d_sh_n(SL)
             do c = 1 , d_sh_n(SK)
              do b = 1 , d_sh_n(SJ)
               do a = 1 , d_sh_n(SI)
                val2 = dble(sym_sign(d_sh_first(SI)+a-1,ig)*sym_sign(d_sh_first(SJ)+b-1,ig) &
                          * sym_sign(d_sh_first(SK)+c-1,ig)*sym_sign(d_sh_first(SL)+d-1,ig))*blk(a,b,c,d)
                dev = max(dev, abs(val2 - blk2(sym_perm(d_sh_first(SI)+a-1,ig) - d_sh_first(g4(1)) + 1, &
                                               sym_perm(d_sh_first(SJ)+b-1,ig) - d_sh_first(g4(2)) + 1, &
                                               sym_perm(d_sh_first(SK)+c-1,ig) - d_sh_first(g4(3)) + 1, &
                                               sym_perm(d_sh_first(SL)+d-1,ig) - d_sh_first(g4(4)) + 1)))
               end do
              end do
             end do
            end do
            deallocate(blk, blk2)
          end do
        end do
        !$omp end parallel do
        deallocate(ao2sh)
        w1 = omp_get_wtime()
        write(outfile,'(a,i0,a,i0,a,es9.2,a,f8.2,a)') ' symmetry (direct): self-check, ', nsamp, &
             ' quartets over ', d_nG-1, ' operations: max dev ', dev, '   (', w1-w0, ' s)'
        if (dev > 1.d-10) then
          write(sym_eri_note,'(a,es9.2,a)') 'direct: point operations NOT used: self-check failed (max dev ', dev, ')'
          write(outfile,'(a)') ' symmetry (direct): '//trim(sym_eri_note(9:))
          call dsym_free() ; return
        end if
        write(outfile,'(a,i0,a,i0,a,i0,a,i0,a,i0,a)') ' symmetry (direct): ', d_nG, ' operations, ', &
             count(ds_rep == (/ (a, a = 1, nshuc) /)), ' of ', nshuc, ' unit-cell shells are orbit representatives, ', &
             d_nbra, ' bra pairs, ', d_nkp, ' ket pairs'
        write(sym_eri_note,'(a,i0,a)') 'direct J/K with point symmetry (', d_nG, ' operations, skeleton + symmetrization)'
        d_sym = .true.
        flush(outfile)
      end subroutine dsym_init

      subroutine dsym_free()
        if (allocated(ds_scell)) deallocate(ds_scell, ds_sloc, ds_fcell, ds_floc, ds_relc)
        if (allocated(ds_shmap)) deallocate(ds_shmap)
        if (allocated(ds_rep))   deallocate(ds_rep, ds_nto)
        if (allocated(ds_to))    deallocate(ds_to)
        if (allocated(ds_bra))   deallocate(ds_bra, ds_bgo, ds_bg)
        if (allocated(ds_kp))    deallocate(ds_kp)
        d_sym = .false.
      end subroutine dsym_free

      !-----------------------------------------------------------------!
      ! [Nov 2026] P <- 1/|G| sum_g g.P   (guess density only)           !
      ! A guess that fills only part of a degenerate set (e.g. a star of !
      ! k points at the Fermi level of the Huckel guess, seen for LiH    !
      ! 3x3x3) breaks the symmetry, and direct_JK_sym could not be used. !
      ! The averaged guess is symmetric; the SCF then stays symmetric    !
      ! as long as the occupied levels are not partly filled degenerate  !
      ! sets (else direct_JK falls back, with a WARNING).                !
      !-----------------------------------------------------------------!
      subroutine direct_sym_P(nBas, P, iout)
        use symmetry_clifford, only : sym_perm, sym_sign
        integer         , intent(in)    :: nBas, iout
        double precision, intent(inout) :: P(nBas,nBas)
        double precision, allocatable   :: Ps(:,:)
        integer :: ig, m, n
        double precision :: dev
        if (.not. d_sym) return
        allocate(Ps(nBas,nBas))
        Ps = 0.d0
        !$omp parallel do default(shared) private(ig,m,n) schedule(static)
        do n = 1 , nBas
          do ig = 1 , d_nG
            do m = 1 , nBas
              Ps(m,n) = Ps(m,n) + dble(sym_sign(m,ig)*sym_sign(n,ig))*P(sym_perm(m,ig),sym_perm(n,ig))
            end do
          end do
        end do
        !$omp end parallel do
        Ps = Ps/dble(d_nG)
        dev = maxval(abs(Ps - P))
        P = Ps
        deallocate(Ps)
        if (iout > 0) then
          write(iout,'(a,es9.2,a)') ' symmetry (direct): guess density averaged over the operations (max change ', dev, ')'
          flush(iout)
        end if
      end subroutine direct_sym_P

      !-----------------------------------------------------------------!
      ! is P invariant under every operation? |s_m s_n P(gm,gn)-P(m,n)| !
      ! (rows m <= nfuc). If not, the normal direct_JK is used.          !
      !-----------------------------------------------------------------!
      logical function p_invariant(nBas, P)
        use symmetry_clifford, only : sym_perm, sym_sign
        use unitcell_module  , only : nfuc
        use files            , only : outfile
        integer         , intent(in) :: nBas
        double precision, intent(in) :: P(nBas,nBas)
        integer :: ig, m, n
        double precision :: dev
        dev = 0.d0
        !$omp parallel do default(shared) private(ig,m,n) reduction(max:dev) schedule(static)
        do ig = 2 , d_nG
          do n = 1 , nBas
            do m = 1 , nfuc
              dev = max(dev, abs(dble(sym_sign(m,ig)*sym_sign(n,ig))*P(sym_perm(m,ig),sym_perm(n,ig)) - P(m,n)))
            end do
          end do
        end do
        !$omp end parallel do
        p_invariant = (dev <= 1.d-6)
        if (.not. p_invariant .and. d_ncall <= d_nsym_early) then
          write(outfile,'(a,i0,a,es9.2,a)') ' symmetry (direct): SCF call ', d_ncall, ', density not symmetric (max dev ', &
                                            dev, '): J/K from its symmetry average'
          flush(outfile)
          return
        end if
        if (.not. p_invariant .and. .not. d_sym_warned) then
          write(outfile,'(a,es9.2,a)') ' WARNING symmetry (direct): the density is not symmetric (max dev ', dev, &
                                       '): J/K without point symmetry (slow) as long as it stays so'
          flush(outfile)
          d_sym_warned = .true.
        end if
      end function p_invariant

      !-----------------------------------------------------------------!
      ! skeleton J/K over the symmetry-unique shell quartets, then      !
      ! J = 1/|G| sum_g g.J_skel (same for K)                            !
      !                                                                 !
      ! A quartet Q = (SI,SJ,SK,SL), SI in the unit cell, is unique if   !
      ! it is the smallest (lexicographic, after the translation that  !
      ! puts the first shell in the unit cell) of its orbit under       !
      ! 8 permutations x translations x G. Then:                        !
      !  - SI is the representative of its orbit, every shell of Q has  !
      !    a representative >= SI (else a smaller image exists),        !
      !  - an image can only tie or beat Q if its first shell is SI, so !
      !    only the operations g that send a shell X of Q (in the SI    !
      !    orbit) on SI are tried (lists ds_to), X = 1 only with the    !
      !    stabilizer of the bra pair (ds_bg).                          !
      ! |Stab(Q)| = number of g giving an image equal to Q.             !
      !-----------------------------------------------------------------!
      subroutine direct_JK_sym(nBas, P, J, K)

        use keywords         , only : c_prim_screen
        use unitcell_module  , only : nfuc
        use clifford_tools   , only : ps_on
        use symmetry_clifford, only : sym_perm, sym_sign

        integer         , intent(in)  :: nBas
        double precision, intent(in)  :: P(nBas,nBas)
        double precision, intent(out) :: J(nBas,nBas), K(nBas,nBas)

        integer :: ib, ik, SI, SJ, SK, SL, a, b, c, d, i, jj, kk, l, iperm, pp(4), t4(4), q(4), u(4), r(4)
        integer :: ig, is, X, ip, nst, stamp, ii, h, hs, mx, ns, nsh, nshuc, loc
        integer(8) :: key
        logical :: smaller
        double precision :: val, wgt
        double precision, allocatable :: blk(:,:,:,:), Jl(:,:), Kl(:,:), Js(:,:), Ks(:,:)
        integer(8), allocatable :: hkey(:)
        integer, allocatable :: used(:), mark(:)
        integer :: iperm8(4,8), perX(4,2,4)
        data iperm8 / 1,2,3,4, 2,1,3,4, 1,2,4,3, 2,1,4,3, 3,4,1,2, 4,3,1,2, 3,4,2,1, 4,3,2,1 /
        ! the two chemist permutations that put position X first
        data perX / 1,2,3,4, 1,2,4,3,   2,1,3,4, 2,1,4,3,   3,4,1,2, 3,4,2,1,   4,3,1,2, 4,3,2,1 /

        nsh = d_nsh ; nshuc = d_nshuc
        J = 0.d0 ; K = 0.d0
        mx = maxval(d_sh_n(1:nsh))
        hs = 1
        do while (hs < 4*8*mx**4)
          hs = 2*hs
        end do
        allocate(Js(nfuc,nBas), Ks(nfuc,nBas))
        Js = 0.d0 ; Ks = 0.d0
        ps_on = c_prim_screen

        !$omp parallel default(shared) &
        !$omp private(ib,ik,SI,SJ,SK,SL,a,b,c,d,i,jj,kk,l,iperm,pp,t4,q,u,r,ig,is,X,ip,nst,stamp,ii,h,ns,loc) &
        !$omp private(key,smaller,val,wgt,blk,Jl,Kl,hkey,used,mark)
        allocate(Jl(nfuc,nBas), Kl(nfuc,nBas), hkey(0:hs-1), used(8*mx**4), mark(d_nG))
        hkey = -1_8 ; mark = 0 ; stamp = 0
        Jl = 0.d0 ; Kl = 0.d0
        !$omp do schedule(dynamic,1)
        do ib = 1 , d_nbra
          SI = ds_bra(1,ib) ; SJ = ds_bra(2,ib)
          do ik = 1 , int(d_nkp)
            SK = ds_kp(1,ik) ; SL = ds_kp(2,ik)
            if (ds_rep(ds_sloc(SK)+1) < SI .or. ds_rep(ds_sloc(SL)+1) < SI) cycle
            if (d_Qsh(SI,SJ)*d_Qsh(SK,SL) < schwarz_thresh) cycle
            q = (/ SI, SJ, SK, SL /)
            stamp = stamp + 1 ; nst = 0 ; smaller = .false.
            ! X = 1: stabilizer of the bra pair, images (SI, SJ, g SK, g SL) and (SI, SJ, g SL, g SK)
            do is = ds_bgo(ib) , ds_bgo(ib+1) - 1
              ig = ds_bg(is)
              u(1) = ds_shmap(SI,ig) ; u(3) = ds_shmap(SK,ig) ; u(4) = ds_shmap(SL,ig)
              r(3) = ds_relc(ds_scell(u(1)), ds_scell(u(3)))*nshuc + ds_sloc(u(3)) + 1
              r(4) = ds_relc(ds_scell(u(1)), ds_scell(u(4)))*nshuc + ds_sloc(u(4)) + 1
              if (r(3) < SK .or. (r(3) == SK .and. r(4) < SL)) then
                smaller = .true. ; exit
              end if
              if (r(4) < SK .or. (r(4) == SK .and. r(3) < SL)) then
                smaller = .true. ; exit
              end if
              if ((r(3) == SK .and. r(4) == SL) .or. (r(4) == SK .and. r(3) == SL)) then
                if (mark(ig) /= stamp) then
                  mark(ig) = stamp ; nst = nst + 1
                end if
              end if
            end do
            if (smaller) cycle
            ! X = 2, 3, 4: shells of the SI orbit brought on SI
            do X = 2 , 4
              loc = ds_sloc(q(X)) + 1
              if (ds_rep(loc) /= SI) cycle
              do is = 1 , ds_nto(loc,SI)
                ig = ds_to(is,loc,SI)
                do ii = 1 , 4
                  u(ii) = ds_shmap(q(ii),ig)
                end do
                do ip = 1 , 2
                  t4 = u(perX(:,ip,X))
                  do ii = 1 , 4
                    r(ii) = ds_relc(ds_scell(t4(1)), ds_scell(t4(ii)))*nshuc + ds_sloc(t4(ii)) + 1
                  end do
                  do ii = 1 , 4
                    if (r(ii) < q(ii)) then
                      smaller = .true. ; exit
                    end if
                    if (r(ii) > q(ii)) exit
                  end do
                  if (smaller) exit
                  if (all(r == q)) then
                    if (mark(ig) /= stamp) then
                      mark(ig) = stamp ; nst = nst + 1
                    end if
                  end if
                end do
                if (smaller) exit
              end do
              if (smaller) exit
            end do
            if (smaller) cycle
            wgt = dble(d_nG)/dble(max(nst,1))

            allocate(blk(d_sh_n(SI),d_sh_n(SJ),d_sh_n(SK),d_sh_n(SL)))
            if (ps_on) call set_ps(SI, SJ, SK, SL)
            call ERI_shell_block_toroidal(d_ndp,d_sh_n(SI),d_AO(d_sh_first(SI):),d_sh_n(SJ),d_AO(d_sh_first(SJ):), &
                                          d_sh_n(SK),d_AO(d_sh_first(SK):),d_sh_n(SL),d_AO(d_sh_first(SL):),blk)
            ns = 0
            do d = 1 , d_sh_n(SL)
             do c = 1 , d_sh_n(SK)
              do b = 1 , d_sh_n(SJ)
               do a = 1 , d_sh_n(SI)
                i  = d_sh_first(SI)+a-1 ; jj = d_sh_first(SJ)+b-1
                kk = d_sh_first(SK)+c-1 ; l  = d_sh_first(SL)+d-1
                if (d_Q(i,jj)*d_Q(kk,l) < schwarz_thresh) cycle
                val = wgt*blk(a,b,c,d)
                do iperm = 1 , 8
                  t4 = (/ i, jj, kk, l /)
                  t4 = t4(iperm8(:,iperm))
                  do ii = 1 , 4
                    pp(ii) = ds_relc(ds_fcell(t4(1)), ds_fcell(t4(ii)))*nfuc + ds_floc(t4(ii)) + 1
                  end do
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
        !$omp end do
        !$omp critical (direct_jks_sum)
        Js = Js + Jl
        Ks = Ks + Kl
        !$omp end critical (direct_jks_sum)
        deallocate(hkey, used, mark)
        !$omp barrier

        ! symmetrization: J(g m, g n) += s_m s_n J_skel(m,n) / |G|  (rows brought to the unit cell)
        Jl = 0.d0 ; Kl = 0.d0
        !$omp do schedule(static)
        do ig = 1 , d_nG
          do b = 1 , nBas
            do a = 1 , nfuc
              i  = sym_perm(a,ig) ; jj = sym_perm(b,ig)
              pp(1) = ds_floc(i) + 1
              pp(2) = ds_relc(ds_fcell(i), ds_fcell(jj))*nfuc + ds_floc(jj) + 1
              val = dble(sym_sign(a,ig)*sym_sign(b,ig))
              Jl(pp(1),pp(2)) = Jl(pp(1),pp(2)) + val*Js(a,b)
              Kl(pp(1),pp(2)) = Kl(pp(1),pp(2)) + val*Ks(a,b)
            end do
          end do
        end do
        !$omp end do
        !$omp critical (direct_jks_sym)
        J(1:nfuc,:) = J(1:nfuc,:) + Jl/dble(d_nG)
        K(1:nfuc,:) = K(1:nfuc,:) + Kl/dble(d_nG)
        !$omp end critical (direct_jks_sym)
        deallocate(Jl, Kl)
        !$omp end parallel

        ps_on = .false.
        deallocate(Js, Ks)
      end subroutine direct_JK_sym

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
