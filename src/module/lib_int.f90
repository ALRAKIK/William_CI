!=======================================================================
! libcint_integrals.f90
!   All libcint wrappers, wrapped in a module so allocatable dummy
!   arguments work correctly (explicit interface is automatic).
!=======================================================================
module libcint_integrals

      use atom_basis
      use iso_c_binding, only: c_null_ptr
      implicit none

      
contains

      !---------------------------------------------------------------------
      ! build_libcint_tables : call ONCE, returns atm/bas/env/offset
      !---------------------------------------------------------------------

      subroutine build_libcint_tables(n_atoms, geometry, atoms, &
                                  natm, nbas, dim, atm, bas, env, offset)

      !=====================================================================
      ! libcint table layout
      !---------------------------------------------------------------------
      ! libcint stores all basis-set / geometry data in three integer or
      ! real arrays:
      !
      !     atm(ATM_SLOTS, natm)   -- one column per ATOM
      !     bas(BAS_SLOTS, nbas)   -- one column per SHELL
      !     env(:)                 -- flat real array holding every
      !                               number (coords, exponents,
      !                               coefficients) referenced by the
      !                               PTR_* slot indices below.
      !
      ! The number of slots per atom (ATM_SLOTS = 6) and per shell
      ! (BAS_SLOTS = 8) is FIXED by libcint's C API.  Each row has a
      ! specific, hard-coded meaning; do not reorder them.
      !
      !---------------------------------------------------------------------
      ! atm(:, iatm)  -- 6 slots per atom
      !---------------------------------------------------------------------
      !   row 1  CHARGE_OF         nuclear charge Z of this atom
      !   row 2  PTR_COORD         index into env(:) where
      !                            (x, y, z) of this atom begin;
      !                            consumes 3 real entries
      !   row 3  NUC_MOD_OF        nuclear model:
      !                              0 = point charge   (default)
      !                              1 = Gaussian charge
      !                              2 = point charge + MM point
      !   row 4  PTR_ZETA          if NUC_MOD_OF = 1, index into env
      !                            for the Gaussian exponent; else 0
      !   row 5  PTR_FRAC_CHARGE   fractional-charge pointer (QM/MM);
      !                            unused here -> 0
      !   row 6  RESERVE_ATMSLOT   reserved by libcint -- always 0
      !
      !---------------------------------------------------------------------
      ! bas(:, ish)  -- 8 slots per shell
      !---------------------------------------------------------------------
      !   row 1  ATOM_OF          which atom this shell sits on
      !                           (0-based index, so atom iatm -> iatm-1)
      !   row 2  ANG_OF           angular momentum l of the shell
      !                              0 = s  (1 cart. function)
      !                              1 = p  (3 cart. functions: x,y,z)
      !                              2 = d  (6 cart. functions: xx,xy,
      !                                     xz,yy,yz,zz)
      !                              3 = f  (10 cart. functions)
      !                              ...
      !                           NOTE:  libcint expands a single shell
      !                           entry of angular momentum l into
      !                           (l+1)(l+2)/2 Cartesian components
      !                           automatically -- do NOT create
      !                           separate bas entries for px,py,pz.
      !   row 3  NPRIM_OF         number of primitive Gaussians in this
      !                           shell (should equal the size of the
      !                           exponent array we copy into env)
      !   row 4  NCTR_OF          number of contracted functions in this
      !                           shell; 1 for STO-3G-style bases
      !   row 5  KDIM_OF          reserved for general contractions;
      !                           unused here -> 1
      !   row 6  PTR_EXP          index into env(:) where the NPRIM_OF
      !                           primitive exponents begin; consumes
      !                           NPRIM_OF real entries
      !   row 7  PTR_COEFF        index into env(:) where the NPRIM_OF
      !                           contraction coefficients begin;
      !                           consumes NPRIM_OF real entries
      !   row 8  RESERVE_BASLOT   reserved by libcint -- always 0
      !
      !---------------------------------------------------------------------
      ! Important conventions to keep in mind
      !---------------------------------------------------------------------
      ! (1) All atom indices inside bas are 0-BASED.  The shell table
      !     in this module builds them as "iatm - 1" because Fortran
      !     loop variables are 1-based.
      !
      ! (2) All shell indices passed to the libcint kernels
      !     (cint1e_ovlp_cart, cint1e_kin_cart, cint1e_nuc_cart,
      !      cint2e_cart, CINTcgto_cart, ...) are 0-BASED.  Inside the
      !     shell-pair / shell-quartet loops this is "ish-1", "jsh-1".
      !
      ! (3) PTR_* values are 0-BASED offsets into env(:).  The Fortran
      !     code stores them in an integer variable `off` and uses
      !     env(off+1 : off+n) to slice the 1-based Fortran array.
      !
      ! (4) libcint does NOT normalise the contraction coefficients
      !     you supply.  We pre-multiply each coefficient by
      !
      !           sqrt( 4*pi / (2*l + 1) )
      !
      !     so that the contracted function has unit self-overlap:
      !
      !           OV(i,i) = 1.0   for every basis function i.
      !
      !     This factor is l-dependent (fac_s, fac_p, fac_d below).
      !
      ! (5) The env(:) array is one large flat buffer holding every
      !     real value (atom coordinates, exponents, coefficients).
      !     Its size is chosen generously (20000 here).  For larger
      !     molecules grow it or compute the exact required size as
      !
      !           PTR_ENV_START
      !         + 3*natm
      !         + sum_over_shells( 2 * NPRIM_OF(ish) )
      !
      ! (6) PTR_ENV_START = 20.  Slots env(1:20) are reserved for
      !     global bookkeeping; libcint never reads them.  Any value
      !     from 0 up would work; 20 is the usual convention (adopted
      !     by PySCF and libfint).
      !
      !=====================================================================

      implicit none

      integer,          intent(in)  :: n_atoms
      double precision, intent(in)  :: geometry(n_atoms, 3)
      type(atom),       intent(in)  :: atoms(n_atoms)
      integer,          intent(out) :: natm, nbas, dim
      integer, allocatable,          intent(out) :: atm(:,:), bas(:,:), offset(:)
      double precision, allocatable, intent(out) :: env(:)

      integer, parameter :: CHARGE_OF=1, PTR_COORD=2, NUC_MOD_OF=3, PTR_ZETA=4
      integer, parameter :: ATOM_OF=1, ANG_OF=2, NPRIM_OF=3, NCTR_OF=4
      integer, parameter :: KDIM_OF=5, PTR_EXP=6, PTR_COEFF=7
      integer, parameter :: PTR_ENV_START=20

      integer, external :: CINTcgto_cart
      integer :: off, ish, iatm, ish_s, ish_p, ish_d, ish_f , ish_g, k, nprim

      double precision, parameter :: pi = 3.14159265358979323846d0
      double precision :: fac_s, fac_p, fac_d , fac_f , fac_g 
      integer          :: dim_expected


      fac_s = sqrt(4.0d0 * pi / 1.0d0)
      fac_p = sqrt(4.0d0 * pi / 3.0d0)
      fac_d = sqrt(4.0d0 * pi / 5.0d0)
      fac_f = sqrt(4.0d0 * pi / 7.0d0)
      fac_g = sqrt(4.0d0 * pi / 9.0d0)

      natm = n_atoms
      nbas = 0
      do iatm = 1, natm
         nbas = nbas + atoms(iatm)%num_s_function &
                     + atoms(iatm)%num_p_function &
                     + atoms(iatm)%num_d_function &
                     + atoms(iatm)%num_f_function &
                     + atoms(iatm)%num_g_function
      end do

      allocate(atm(6, natm), bas(8, nbas), env(20000))

      ! ---- atoms ----
      off = PTR_ENV_START
      do iatm = 1, natm
         atm(CHARGE_OF,  iatm) = atoms(iatm)%charge
         atm(PTR_COORD,  iatm) = off
         atm(NUC_MOD_OF, iatm) = 0
         atm(PTR_ZETA,   iatm) = 0
         env(off+1:off+3) = geometry(iatm, 1:3)
         off = off + 3
      end do

      ! ! ---- shells ----
      ! ish = 0
      ! do iatm = 1, natm

      !   ! s
      !   nprim = size(atoms(iatm)%exponent_s)
      !   if (nprim > 0) then
      !      do ish_s = 1, atoms(iatm)%num_s_function
      !         ish = ish + 1
      !         bas(ATOM_OF,   ish) = iatm - 1
      !         bas(ANG_OF,    ish) = 0
      !         bas(NPRIM_OF,  ish) = nprim
      !         bas(NCTR_OF,   ish) = 1
      !         bas(KDIM_OF,   ish) = 1
      !         bas(PTR_EXP, ish) = off
      !         env(off+1:off+nprim) = atoms(iatm)%exponent_s(1:nprim)
      !         off = off + nprim
      !         bas(PTR_COEFF, ish) = off
      !         do k = 1, nprim
      !            env(off+k) = atoms(iatm)%coefficient_s(k, ish_s) * fac_s
      !         end do
      !         off = off + nprim
      !      end do
      !   end if
 
      !   ! p
      !   nprim = size(atoms(iatm)%exponent_p)
      !   if (nprim > 0) then
      !      do ish_p = 1, atoms(iatm)%num_p_function
      !         ish = ish + 1
      !         bas(ATOM_OF,   ish) = iatm - 1
      !         bas(ANG_OF,    ish) = 1
      !         bas(NPRIM_OF,  ish) = nprim
      !         bas(NCTR_OF,   ish) = 1
      !         bas(KDIM_OF,   ish) = 1
      !         bas(PTR_EXP, ish) = off
      !         env(off+1:off+nprim) = atoms(iatm)%exponent_p(1:nprim)
      !         off = off + nprim
      !         bas(PTR_COEFF, ish) = off
      !         do k = 1, nprim
      !            env(off+k) = atoms(iatm)%coefficient_p(k, ish_p) * fac_p
      !         end do
      !         off = off + nprim
      !      end do
      !   end if
 
      !   ! d
      !   nprim = size(atoms(iatm)%exponent_d)
      !   if (nprim > 0) then
      !      do ish_d = 1, atoms(iatm)%num_d_function
      !         ish = ish + 1
      !         bas(ATOM_OF,   ish) = iatm - 1
      !         bas(ANG_OF,    ish) = 2
      !         bas(NPRIM_OF,  ish) = nprim
      !         bas(NCTR_OF,   ish) = 1
      !         bas(KDIM_OF,   ish) = 1
      !         bas(PTR_EXP, ish) = off
      !         env(off+1:off+nprim) = atoms(iatm)%exponent_d(1:nprim)
      !         off = off + nprim
      !         bas(PTR_COEFF, ish) = off
      !         do k = 1, nprim
      !            env(off+k) = atoms(iatm)%coefficient_d(k, ish_d) * fac_d
      !         end do
      !         off = off + nprim
      !      end do
      !   end if

      !   ! f
      !   nprim = size(atoms(iatm)%exponent_f)
      !   if (nprim > 0) then
      !      do ish_f = 1, atoms(iatm)%num_f_function
      !         ish = ish + 1
      !         bas(ATOM_OF,   ish) = iatm - 1
      !         bas(ANG_OF,    ish) = 3
      !         bas(NPRIM_OF,  ish) = nprim
      !         bas(NCTR_OF,   ish) = 1
      !         bas(KDIM_OF,   ish) = 1
      !         bas(PTR_EXP, ish) = off
      !         env(off+1:off+nprim) = atoms(iatm)%exponent_f(1:nprim)
      !         off = off + nprim
      !         bas(PTR_COEFF, ish) = off
      !         do k = 1, nprim
      !            env(off+k) = atoms(iatm)%coefficient_f(k, ish_f) * fac_f
      !         end do
      !         off = off + nprim
      !      end do
      !   end if

      !   ! g
      !   nprim = size(atoms(iatm)%exponent_g)
      !   if (nprim > 0) then
      !      do ish_g = 1, atoms(iatm)%num_g_function
      !         ish = ish + 1
      !         bas(ATOM_OF,   ish) = iatm - 1
      !         bas(ANG_OF,    ish) = 4
      !         bas(NPRIM_OF,  ish) = nprim
      !         bas(NCTR_OF,   ish) = 1
      !         bas(KDIM_OF,   ish) = 1
      !         bas(PTR_EXP, ish) = off
      !         env(off+1:off+nprim) = atoms(iatm)%exponent_g(1:nprim)
      !         off = off + nprim
      !         bas(PTR_COEFF, ish) = off
      !         do k = 1, nprim
      !            env(off+k) = atoms(iatm)%coefficient_g(k, ish_g) * fac_g
      !         end do
      !         off = off + nprim
      !      end do
      !   end if

      ! end do

      ! ! ---- shell offset table ----
      ! allocate(offset(nbas+1))
      ! offset(1) = 0
      ! do ish = 1, nbas
      !    offset(ish+1) = offset(ish) + CINTcgto_cart(ish-1, bas)
      ! end do
      ! dim = offset(nbas+1)

      ! ---- shells ----
      ish = 0
      do iatm = 1, natm

         ! --- consistency check: counters vs arrays ---
         if (atoms(iatm)%num_s_function > 0 .and. size(atoms(iatm)%exponent_s) == 0) then
            write(*,*) "ERROR: atom ", iatm, " has num_s_function > 0 but exponent_s is empty"; stop
         end if
         if (atoms(iatm)%num_p_function > 0 .and. size(atoms(iatm)%exponent_p) == 0) then
            write(*,*) "ERROR: atom ", iatm, " has num_p_function > 0 but exponent_p is empty"; stop
         end if
         if (atoms(iatm)%num_d_function > 0 .and. size(atoms(iatm)%exponent_d) == 0) then
            write(*,*) "ERROR: atom ", iatm, " has num_d_function > 0 but exponent_d is empty"; stop
         end if
         if (atoms(iatm)%num_f_function > 0 .and. size(atoms(iatm)%exponent_f) == 0) then
            write(*,*) "ERROR: atom ", iatm, " has num_f_function > 0 but exponent_f is empty"; stop
         end if
         if (atoms(iatm)%num_g_function > 0 .and. size(atoms(iatm)%exponent_g) == 0) then
            write(*,*) "ERROR: atom ", iatm, " has num_g_function > 0 but exponent_g is empty"; stop
         end if

         ! s
         nprim = size(atoms(iatm)%exponent_s)
         if (nprim > 0) then
            do ish_s = 1, atoms(iatm)%num_s_function
               ish = ish + 1
               bas(ATOM_OF,  ish) = iatm - 1
               bas(ANG_OF,   ish) = 0
               bas(NPRIM_OF, ish) = nprim
               bas(NCTR_OF,  ish) = 1
               bas(KDIM_OF,  ish) = 1
               bas(PTR_EXP,  ish) = off
               env(off+1:off+nprim) = atoms(iatm)%exponent_s(1:nprim)
               off = off + nprim
               bas(PTR_COEFF, ish) = off
               do k = 1, nprim
                  env(off+k) = atoms(iatm)%coefficient_s(k, ish_s) * fac_s
               end do
               off = off + nprim
            end do
         end if

         ! p
         nprim = size(atoms(iatm)%exponent_p)
         if (nprim > 0) then
            do ish_p = 1, atoms(iatm)%num_p_function
               ish = ish + 1
               bas(ATOM_OF,  ish) = iatm - 1
               bas(ANG_OF,   ish) = 1
               bas(NPRIM_OF, ish) = nprim
               bas(NCTR_OF,  ish) = 1
               bas(KDIM_OF,  ish) = 1
               bas(PTR_EXP,  ish) = off
               env(off+1:off+nprim) = atoms(iatm)%exponent_p(1:nprim)
               off = off + nprim
               bas(PTR_COEFF, ish) = off
               do k = 1, nprim
                  env(off+k) = atoms(iatm)%coefficient_p(k, ish_p) * fac_p
               end do
               off = off + nprim
            end do
         end if

         ! d
         nprim = size(atoms(iatm)%exponent_d)
         if (nprim > 0) then
            do ish_d = 1, atoms(iatm)%num_d_function
               ish = ish + 1
               bas(ATOM_OF,  ish) = iatm - 1
               bas(ANG_OF,   ish) = 2
               bas(NPRIM_OF, ish) = nprim
               bas(NCTR_OF,  ish) = 1
               bas(KDIM_OF,  ish) = 1
               bas(PTR_EXP,  ish) = off
               env(off+1:off+nprim) = atoms(iatm)%exponent_d(1:nprim)
               off = off + nprim
               bas(PTR_COEFF, ish) = off
               do k = 1, nprim
                  env(off+k) = atoms(iatm)%coefficient_d(k, ish_d) * fac_d
               end do
               off = off + nprim
            end do
         end if

         ! f
         nprim = size(atoms(iatm)%exponent_f)
         if (nprim > 0) then
            do ish_f = 1, atoms(iatm)%num_f_function
               ish = ish + 1
               bas(ATOM_OF,  ish) = iatm - 1
               bas(ANG_OF,   ish) = 3
               bas(NPRIM_OF, ish) = nprim
               bas(NCTR_OF,  ish) = 1
               bas(KDIM_OF,  ish) = 1
               bas(PTR_EXP,  ish) = off
               env(off+1:off+nprim) = atoms(iatm)%exponent_f(1:nprim)
               off = off + nprim
               bas(PTR_COEFF, ish) = off
               do k = 1, nprim
                  env(off+k) = atoms(iatm)%coefficient_f(k, ish_f) * fac_f
               end do
               off = off + nprim
            end do
         end if

         ! g
         nprim = size(atoms(iatm)%exponent_g)
         if (nprim > 0) then
            do ish_g = 1, atoms(iatm)%num_g_function
               ish = ish + 1
               bas(ATOM_OF,  ish) = iatm - 1
               bas(ANG_OF,   ish) = 4
               bas(NPRIM_OF, ish) = nprim
               bas(NCTR_OF,  ish) = 1
               bas(KDIM_OF,  ish) = 1
               bas(PTR_EXP,  ish) = off
               env(off+1:off+nprim) = atoms(iatm)%exponent_g(1:nprim)
               off = off + nprim
               bas(PTR_COEFF, ish) = off
               do k = 1, nprim
                  env(off+k) = atoms(iatm)%coefficient_g(k, ish_g) * fac_g
               end do
               off = off + nprim
            end do
         end if

      end do

      ! ---- post-fill sanity checks ----
      if (ish /= nbas) then
         write(*,*) "FATAL: fill loop produced ", ish, " shells but nbas = ", nbas
         write(*,*) "Some shells counted in nbas were not filled. Inspect atoms above."
         stop
      end if

      if (off > size(env)) then
         write(*,*) "FATAL: env overflow — required ", off, " but allocated ", size(env)
         stop
      end if
      write(*,'(A,I10,A,I10)') " env used = ", off, "  of ", size(env)

      ! ---- cross-check dim against direct count ----
      dim_expected = 0
      do iatm = 1, natm
         dim_expected = dim_expected &
              + 1  * atoms(iatm)%num_s_function &
              + 3  * atoms(iatm)%num_p_function &
              + 6  * atoms(iatm)%num_d_function &
              + 10 * atoms(iatm)%num_f_function &
              + 15 * atoms(iatm)%num_g_function
      end do

      ! ---- shell offset table ----
      allocate(offset(nbas+1))
      offset(1) = 0
      do ish = 1, nbas
         offset(ish+1) = offset(ish) + CINTcgto_cart(ish-1, bas)
      end do
      dim = offset(nbas+1)

      if (dim /= dim_expected) then
         write(*,*) "MISMATCH: CINTcgto_cart sum = ", dim, &
                    " vs direct count = ", dim_expected
         write(*,*) "Printing shell table to find the offending shell:"
         do ish = 1, nbas
            write(*,'(A,I4,A,I3,A,I3,A,I3)') &
                 " shell ", ish, &
                 "  ang=", bas(ANG_OF,   ish), &
                 "  nprim=", bas(NPRIM_OF, ish), &
                 "  atom=",  bas(ATOM_OF,  ish)
         end do
         stop
      end if

      write(*,*) "the table finished"

      end subroutine build_libcint_tables


      !---------------------------------------------------------------------
      ! one_electron_libcint : OV, KIN, NA — one shell-pair loop
      !---------------------------------------------------------------------

      subroutine one_electron_libcint(atm, bas, env, offset, natm, nbas, dim, &
                                  OV, KIN, NA)
      implicit none

      integer,          intent(in)  :: natm, nbas, dim
      integer,          intent(in)  :: atm(6,natm), bas(8,nbas), offset(nbas+1)
      double precision, intent(in)  :: env(:)
      double precision, intent(out) :: OV(dim,dim), KIN(dim,dim), NA(dim,dim)

      integer, external :: CINTcgto_cart
      integer, external :: cint1e_ovlp_cart, cint1e_kin_cart, cint1e_nuc_cart

      integer :: ish, jsh, di, dj, ioff, joff, not0, shls(2)
      double precision, allocatable :: buf(:,:)

      OV = 0.0d0; KIN = 0.0d0; NA = 0.0d0

      do ish = 1, nbas
         do jsh = 1, nbas
            di   = CINTcgto_cart(ish-1, bas)
            dj   = CINTcgto_cart(jsh-1, bas)
            ioff = offset(ish)
            joff = offset(jsh)
            shls = [ish-1, jsh-1]

            allocate(buf(di, dj)); buf = 0.0d0

            not0 = cint1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env)
            OV(ioff+1:ioff+di, joff+1:joff+dj) = buf

            not0 = cint1e_kin_cart(buf, shls, atm, natm, bas, nbas, env)
            KIN(ioff+1:ioff+di, joff+1:joff+dj) = buf

            not0 = cint1e_nuc_cart(buf, shls, atm, natm, bas, nbas, env)
            NA(ioff+1:ioff+di, joff+1:joff+dj) = buf

            deallocate(buf)
         end do
      end do

      end subroutine one_electron_libcint


      !-----------------------------------------------------------------!
      ! eri_libcint : full (ij|kl) tensor — shell-quartet loop
      !-----------------------------------------------------------------!

      subroutine eri_libcint_chemist(atm, bas, env, offset, natm, nbas, dim, ERI)
        use iso_c_binding, only: c_null_ptr
        implicit none
        integer,          intent(in)  :: natm, nbas, dim
        integer,          intent(in)  :: atm(6,natm), bas(8,nbas), offset(nbas+1)
        double precision, intent(in)  :: env(:)
        double precision, intent(out) :: ERI(dim,dim,dim,dim)
        integer, external             :: CINTcgto_cart, cint2e_cart

        integer :: ish, jsh, ksh, lsh, lsh_max, ipair, npairs
        integer :: di, dj, dk, dl, ioff, joff, koff, loff, not0, shls(4)
        integer :: a, b, c, d, ia, jb, kc, ld
        integer, allocatable :: pair_i(:), pair_j(:)
        double precision :: val, schwarz_thresh
        double precision, allocatable :: buf(:,:,:,:), schwarz(:,:)

        schwarz_thresh = 1.0d-12   ! tune: 1e-10 is looser but faster
        ERI = 0.0d0

        !--------------------------------------------------------------------
        ! Step 1: Schwarz bounds  Q(ish,jsh) = sqrt( max|(ij|ij)| )
        !         Only N*(N+1)/2 libcint calls — cheap
        !--------------------------------------------------------------------
        allocate(schwarz(nbas, nbas))
        schwarz = 0.0d0

        do ish = 1, nbas
          di   = CINTcgto_cart(ish-1, bas)
          do jsh = 1, ish
            dj   = CINTcgto_cart(jsh-1, bas)
            shls = [ish-1, jsh-1, ish-1, jsh-1]   ! (ij|ij)
            allocate(buf(di, dj, di, dj))
            buf  = 0.0d0
            not0 = cint2e_cart(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
            schwarz(ish,jsh) = sqrt(maxval(abs(buf)))
            schwarz(jsh,ish) = schwarz(ish,jsh)
            deallocate(buf)
          end do
        end do
     
        !--------------------------------------------------------------------
        ! Step 2: Flat pair list for clean OMP over triangular loop
        !--------------------------------------------------------------------
        npairs = nbas*(nbas+1)/2
        allocate(pair_i(npairs), pair_j(npairs))
        ipair = 0
        do ish = 1, nbas
          do jsh = 1, ish
            ipair = ipair + 1
            pair_i(ipair) = ish
            pair_j(ipair) = jsh
          end do
        end do
     
        !--------------------------------------------------------------------
        ! Step 3: Main loop — 8-fold symmetry + Schwarz screening
        !--------------------------------------------------------------------
        !$omp parallel do schedule(dynamic) &
        !$omp   private(ipair,ish,jsh,ksh,lsh,lsh_max, &
        !$omp           di,dj,dk,dl,ioff,joff,koff,loff, &
        !$omp           shls,not0,buf,a,b,c,d,ia,jb,kc,ld,val)
        do ipair = 1, npairs
          ish  = pair_i(ipair);   jsh  = pair_j(ipair)
          ioff = offset(ish);     joff = offset(jsh)
          di   = CINTcgto_cart(ish-1, bas)
          dj   = CINTcgto_cart(jsh-1, bas)
      
          do ksh = 1, ish
            lsh_max = ksh
            if (ksh == ish) lsh_max = jsh
         
            koff = offset(ksh)
            dk   = CINTcgto_cart(ksh-1, bas)
         
            do lsh = 1, lsh_max
              !--------------------------------------------------------------
              ! Schwarz screen: skip if bound < threshold
              !--------------------------------------------------------------
              if (schwarz(ish,jsh) * schwarz(ksh,lsh) < schwarz_thresh) cycle
            
              loff = offset(lsh)
              dl   = CINTcgto_cart(lsh-1, bas)
            
              shls = [ish-1, jsh-1, ksh-1, lsh-1]
              allocate(buf(di, dj, dk, dl))
              buf  = 0.0d0
              not0 = cint2e_cart(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
            
              do d = 1, dl;  ld = loff + d
                do c = 1, dk;  kc = koff + c
                  do b = 1, dj;  jb = joff + b
                    do a = 1, di;  ia = ioff + a
                      val = buf(a, b, c, d)
                      ERI(ia,jb,kc,ld) = val
                      ERI(jb,ia,kc,ld) = val
                      ERI(ia,jb,ld,kc) = val
                      ERI(jb,ia,ld,kc) = val
                      ERI(kc,ld,ia,jb) = val
                      ERI(ld,kc,ia,jb) = val
                      ERI(kc,ld,jb,ia) = val
                      ERI(ld,kc,jb,ia) = val
                    end do
                  end do
                end do
              end do
              deallocate(buf)
           
            end do   ! lsh
          end do   ! ksh
        end do   ! ipair
        !$omp end parallel do
     
        deallocate(schwarz, pair_i, pair_j)
      end subroutine eri_libcint_chemist

      !------------------------------------------------------------------
      ! Physicist notation: ERI(i,j,k,l) = <ij|kl>
      !------------------------------------------------------------------

      subroutine eri_libcint_physicist(atm, bas, env, offset, natm, nbas, dim, ERI)
        use iso_c_binding, only: c_null_ptr
        implicit none
      
        integer,          intent(in)  :: natm, nbas, dim
        integer,          intent(in)  :: atm(6,natm), bas(8,nbas), offset(nbas+1)
        double precision, intent(in)  :: env(:)
        double precision, intent(out) :: ERI(dim,dim,dim,dim)
      
        integer, external :: CINTcgto_cart, cint2e_cart
      
        integer :: ish, jsh, ksh, lsh, lsh_max, ipair, npairs
        integer :: di, dj, dk, dl, ioff, joff, koff, loff, not0, shls(4)
        integer :: a, b, c, d, ia, jb, kc, ld
        integer, allocatable :: pair_i(:), pair_j(:)
        double precision :: val, schwarz_thresh
        double precision, allocatable :: buf(:,:,:,:), schwarz(:,:)
      
        schwarz_thresh = 1.0d-12
        ERI = 0.0d0
      
        !--------------------------------------------------------------------
        ! Step 1: Schwarz bounds — identical to chemist version
        !--------------------------------------------------------------------
        allocate(schwarz(nbas, nbas))
        schwarz = 0.0d0
        do ish = 1, nbas
          di = CINTcgto_cart(ish-1, bas)
          do jsh = 1, ish
            dj   = CINTcgto_cart(jsh-1, bas)
            shls = [ish-1, jsh-1, ish-1, jsh-1]
            allocate(buf(di, dj, di, dj))
            buf  = 0.0d0
            not0 = cint2e_cart(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
            schwarz(ish,jsh) = sqrt(maxval(abs(buf)))
            schwarz(jsh,ish) = schwarz(ish,jsh)
            deallocate(buf)
          end do
        end do
     
        !--------------------------------------------------------------------
        ! Step 2: Flat pair list — identical to chemist version
        !--------------------------------------------------------------------
        npairs = nbas*(nbas+1)/2
        allocate(pair_i(npairs), pair_j(npairs))
        ipair = 0
        do ish = 1, nbas
          do jsh = 1, ish
            ipair = ipair + 1
            pair_i(ipair) = ish
            pair_j(ipair) = jsh
          end do
        end do
     
        !--------------------------------------------------------------------
        ! Step 3: Main loop — physicist version
        !
        !  3 differences from chemist:
        !  (a) shls = [ish-1, ksh-1, jsh-1, lsh-1]  ← permuted
        !  (b) allocate(buf(di, dk, dj, dl))          ← dk and dj swapped
        !  (c) val = buf(a, c, b, d)                  ← c and b swapped
        !
        !  Why: <ij|kl>_phys = (ik|jl)_chem
        !  Passing shells [ish,ksh,jsh,lsh] to libcint gives
        !  buf(a,b,c,d) = (ish_a, ksh_b | jsh_c, lsh_d)_chem
        !  so buf(a,c,b,d) = (ish_a, ksh_c | jsh_b, lsh_d)_chem
        !                  = <ish_a, jsh_b | ksh_c, lsh_d>_phys  ✓
        !
        !  8-fold physicist symmetry (real orbitals):
        !  <ij|kl> = <kj|il> = <il|kj> = <kl|ij>
        !          = <ji|lk> = <li|jk> = <jk|li> = <lk|ji>
        !--------------------------------------------------------------------
        !$omp parallel do schedule(dynamic) &
        !$omp   private(ipair,ish,jsh,ksh,lsh,lsh_max, &
        !$omp           di,dj,dk,dl,ioff,joff,koff,loff, &
        !$omp           shls,not0,buf,a,b,c,d,ia,jb,kc,ld,val)
        do ipair = 1, npairs
          ish  = pair_i(ipair);  jsh = pair_j(ipair)
          ioff = offset(ish);    joff = offset(jsh)
          di   = CINTcgto_cart(ish-1, bas)
          dj   = CINTcgto_cart(jsh-1, bas)
      
          do ksh = 1, ish
            lsh_max = ksh
            if (ksh == ish) lsh_max = jsh
            koff = offset(ksh)
            dk   = CINTcgto_cart(ksh-1, bas)
         
            do lsh = 1, lsh_max
              if (schwarz(ish,jsh) * schwarz(ksh,lsh) < schwarz_thresh) cycle
            
              loff = offset(lsh)
              dl   = CINTcgto_cart(lsh-1, bas)
            
              ! (a) permuted shell order
              shls = [ish-1, ksh-1, jsh-1, lsh-1]
              ! (b) buffer shape matches permuted shells: (di, dk, dj, dl)
              allocate(buf(di, dk, dj, dl))
              buf  = 0.0d0
              not0 = cint2e_cart(buf, shls, atm, natm, bas, nbas, env, c_null_ptr)
            
              do d = 1, dl;  ld = loff + d
                do b = 1, dj;  jb = joff + b        ! note: b loops over jsh
                  do c = 1, dk;  kc = koff + c      ! note: c loops over ksh
                    do a = 1, di;  ia = ioff + a
                      ! (c) read buf(a,c,b,d) to get <ia,jb|kc,ld>_phys
                      val = buf(a, c, b, d)
                      ! 8-fold physicist symmetry
                      ERI(ia,jb,kc,ld) = val    ! <ij|kl>
                      ERI(kc,jb,ia,ld) = val    ! <kj|il>
                      ERI(ia,ld,kc,jb) = val    ! <il|kj>
                      ERI(kc,ld,ia,jb) = val    ! <kl|ij>
                      ERI(jb,ia,ld,kc) = val    ! <ji|lk>
                      ERI(ld,ia,jb,kc) = val    ! <li|jk>
                      ERI(jb,kc,ld,ia) = val    ! <jk|li>
                      ERI(ld,kc,jb,ia) = val    ! <lk|ji>
                    end do
                  end do
                end do
              end do
              deallocate(buf)
           
            end do ! lsh
          end do ! ksh
        end do ! ipair
        !$omp end parallel do
     
        deallocate(schwarz, pair_i, pair_j)
     
      end subroutine eri_libcint_physicist

      !------------------------------------------------------------------
      ! destroy_libcint_tables : free memory
      !------------------------------------------------------------------

      subroutine destroy_libcint_tables(atm, bas, env, offset)

        implicit none
        integer, allocatable,          intent(inout) :: atm(:,:), bas(:,:), offset(:)
        double precision, allocatable, intent(inout) :: env(:)

        if (allocated(atm))    deallocate(atm)
        if (allocated(bas))    deallocate(bas)
        if (allocated(env))    deallocate(env)
        if (allocated(offset)) deallocate(offset)

      end subroutine destroy_libcint_tables

end module libcint_integrals