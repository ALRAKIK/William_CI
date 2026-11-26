! =====================================================================================================
! [Nov 2026] keyword twist_average (twist_average_M): grand-canonical twist-averaged RHF on the Clifford
! torus (1D, 2D, 3D, any number of functions per cell), run after the ordinary SAO RHF.
!
! Theory and tests: doc/metal_chain_lnN_theory.pdf, project notes metal_3D_real_test.md.
!
!   twist theta_s (periodic directions with n_s > 1): twisted boundary condition psi(x + L_s) = e^{i theta_s} psi(x)
!   translation-invariant twist gauge: every matrix element between AO mu (cell c1) and AO nu (cell c2) gets
!       PH(c1,c2) = prod_s exp(i theta_s D_s / n_s),  D = minimum-image cell displacement c2 - c1,
!       |D_s| = n_s/2 exactly: cos(theta_s/2)  (average of the two images)
!     (one-electron integrals: exact by symmetry; two-electron: exact when no pair density reaches L/2,
!      i.e. the box must be large compared with the basis functions)
!   Bloch functions: k = (2 pi m + theta)/L, block A_k(a,b) = sum_c M(a, b@c) PH(1,c) e^{i kappa.R_c}, kappa = 2 pi m/L
!   grand canonical: one common Fermi level mu for all twists, <N_el> fixed (the electron number of a twist varies)
!   Hartree: twist-averaged (neutral, real) density  Q(l,s) = <PH(l,s) D(s,l)>
!   exchange: per twist
!   twists: M midpoints per direction on (0, 2 pi), theta and -theta paired (time reversal: same energy)
!   energy per cell: E = < sum_a sum_nu D(nu,a) [2 h PH + 2 J - K](a,nu) > + E_nn/N_cell
! =====================================================================================================
subroutine twist_RHF_SAO(nBas, nO, S, T, V, ERI, ENuc, P0)

      use files
      use keywords, only : twist_M
      use constants_module
      use unitcell_module
      use omp_lib
      use torus_init, only : Lx, Ly, Lz

      implicit none

      integer         , intent(in) :: nBas, nO
      double precision, intent(in) :: S(nBas,nBas), T(nBas,nBas), V(nBas,nBas)
      double precision, intent(in) :: ERI(nfuc,nBas,nBas,nBas)
      double precision, intent(in) :: ENuc, P0(nBas,nBas)

      integer :: nn(3), ndim, dirs(3), M, ntw_all, ntw, it, i, j, k, a, b, c, c1, c2, id, mu, nu, la, sg, kb
      integer :: idx(3), pidx(3), nmo, info, lwork, iocc
      integer, allocatable :: coord(:,:), cellof(:), fofb(:), tw_idx(:,:), dcell(:,:), nocc_tw(:)
      double precision, allocatable :: tw_th(:,:), tw_w(:), eps(:,:,:), occf(:,:,:), Qbar(:,:), jrow(:,:), allE(:), allW(:), ord(:)
      double precision, allocatable :: rwork(:), sev(:), evv(:)
      complex(dpc), allocatable :: Pc(:,:), PHk(:,:), Frow(:,:,:), Kst(:,:,:), Krow(:,:), Ck(:,:,:,:), Dm(:,:)
      complex(dpc), allocatable :: Ak(:,:), Sk(:,:), Xk(:,:), Fk(:,:), work(:), Gk(:,:,:), dC(:,:,:), A2(:,:), B2(:)
      double precision :: w0, mu_F, E, Eold, Ekin, Ene, Eh, Ex, Enn, nel, dmax, beta, ph, x, sre, sim, Eseg(3), kT, ktf
      double precision, parameter :: s_min = 1.d-9
      complex(dpc) :: z, zk
      logical :: done
      character(len=16) :: tst, tst2
      integer :: stage, Rf, dxyz(3), q, nent, ios, ufile
      double precision :: Estage1, wsum, fval, Wb, Wat, Estar
      double precision, allocatable :: Fref(:,:), eref(:,:,:), occref(:,:,:), F3(:,:,:,:,:), S3(:,:,:,:,:), C3(:,:,:)
      complex(dpc), allocatable :: Frr(:,:)
      logical :: have_file
      integer :: ltst, ist, ist2

      w0 = omp_get_wtime()
      nn = (/ nx, ny, nz /)
      ndim = 0
      do id = 1, 3
        if (nn(id) > 1) then
          ndim = ndim + 1 ; dirs(ndim) = id
        end if
      end do
      M = max(1, twist_M)

      write(outfile,'(a)') ''
      write(outfile,'(a)') ' ======================================================================================'
      write(outfile,'(a)') '  Twist-averaged RHF (keyword twist_average), grand canonical'
      write(outfile,'(a)') ' ======================================================================================'
      if (ndim == 0) then
        write(outfile,'(a)') '  no periodic direction with more than one cell: nothing to twist, skipped'
        return
      end if

      ! ---- cells (x fastest, as in transform_AO_to_SAO) and AO -> cell, function ----
      allocate(coord(3,N_cell), cellof(nBas), fofb(nBas))
      do k = 0, nz-1
        do j = 0, ny-1
          do i = 0, nx-1
            c = k*nx*ny + j*nx + i + 1
            coord(:,c) = (/ i, j, k /)
          end do
        end do
      end do
      do mu = 1, nBas
        cellof(mu) = (mu-1)/nfuc + 1
        fofb(mu)   = mod(mu-1, nfuc) + 1
      end do
      ! cell index of the difference c2 - c1 (mod n)
      allocate(dcell(N_cell,N_cell))
      do c1 = 1, N_cell
        do c2 = 1, N_cell
          idx = modulo(coord(:,c2) - coord(:,c1), nn)
          dcell(c1,c2) = idx(3)*nx*ny + idx(2)*nx + idx(1) + 1
        end do
      end do

      ! ---- twists: M midpoints per twisted direction on (0, 2 pi), pairs (theta, -theta) ----
      ntw_all = M**ndim
      allocate(tw_idx(3,ntw_all), tw_th(3,ntw_all), tw_w(ntw_all))
      ntw = 0
      do i = 0, ntw_all-1
        idx = 0 ; k = i
        do id = 1, ndim
          idx(id) = mod(k, M) ; k = k/M
        end do
        pidx = 0
        pidx(1:ndim) = M - 1 - idx(1:ndim)
        ! keep the lexicographically smaller of the pair
        done = .false.
        do id = ndim, 1, -1
          if (idx(id) /= pidx(id)) then
            done = (idx(id) > pidx(id)) ; exit
          end if
        end do
        if (done) cycle
        ntw = ntw + 1
        tw_idx(:,ntw) = idx
        tw_th(:,ntw) = 0.d0
        do id = 1, ndim
          tw_th(dirs(id),ntw) = 2.d0*pi*(dble(idx(id)) + 0.5d0)/dble(M)
        end do
        tw_w(ntw) = 2.d0/dble(ntw_all)
        if (all(idx(1:ndim) == pidx(1:ndim))) tw_w(ntw) = 1.d0/dble(ntw_all)
      end do

      write(outfile,'(a,i2,a,3i4)')    '  twisted directions             :', ndim, '   cells per direction:', nn
      write(outfile,'(a,i4,a,i6,a,i6)') '  twist points per direction  M  :', M, '   total', ntw_all, '   computed (theta,-theta paired)', ntw
      write(outfile,'(a,i6,a,i4,a,i6)') '  functions                      :', nBas, '   per cell', nfuc, '   occupied (supercell)', nO

      ! ---- storage ----
      allocate(Pc(N_cell,N_cell), PHk(nBas,nBas))
      allocate(Frow(nfuc,nBas,ntw), Kst(nfuc,nBas,ntw), Krow(nfuc,nBas))
      allocate(Ck(nfuc,nfuc,N_cell,ntw), eps(nfuc,N_cell,ntw), occf(nfuc,N_cell,ntw), nocc_tw(ntw))
      allocate(Qbar(nBas,nBas), jrow(nfuc,nBas), Dm(nBas,nBas))
      allocate(Ak(nfuc,nfuc), Sk(nfuc,nfuc), Xk(nfuc,nfuc), Fk(nfuc,nfuc), sev(nfuc), evv(nfuc))
      allocate(Gk(nfuc,nfuc,N_cell), dC(nfuc,nfuc,N_cell))
      allocate(allE(nfuc*N_cell*ntw), allW(nfuc*N_cell*ntw))
      lwork = max(1, 4*nfuc)
      allocate(work(lwork), rwork(max(1,3*nfuc)))

      ! ---- initial Fock rows: core + J(P0) + K(P0) with the twist phases ----
      ! (P0 = converged RHF density of the untwisted torus; the SCF below starts with a Fermi smearing
      !  kT = 0.02 Eh that is reduced by 0.6 per iteration: avoids sticking in a higher aufbau solution)
      Qbar = 0.5d0*transpose(P0)                                       ! density coefficients D0(s,l) -> Q(l,s)
      call j_rows(Qbar, jrow)
      do i = 1, ntw
        call cell_phases(tw_th(:,i), Pc)
        call ao_phases(Pc, PHk)
        Dm = 0.5d0*P0*PHk
        call k_rows(Dm, PHk, Krow)
        do nu = 1, nBas
          do a = 1, nfuc
            Frow(a,nu,i) = (T(a,nu) + V(a,nu))*PHk(a,nu) + 2.d0*PHk(a,nu)*jrow(a,nu) - Krow(a,nu)
          end do
        end do
      end do
      kT = 0.02d0 ; ktf = 0.6d0
      call get_environment_variable('TWIST_KT0', tst2, ltst, ist2) ; if (ist2 == 0 .and. ltst > 0) read(tst2,*) kT
      call get_environment_variable('TWIST_KTF', tst2, ltst, ist2) ; if (ist2 == 0 .and. ltst > 0) read(tst2,*) ktf
      call get_environment_variable('TWIST_START', tst2, ltst, ist2)
      if (.not. (ist2 == 0 .and. trim(tst2) == 'p0')) then          ! default: core-Hamiltonian start (lowest solutions in tests)
        do i = 1, ntw
          call cell_phases(tw_th(:,i), Pc)
          call ao_phases(Pc, PHk)
          Frow(:,:,i) = (T(1:nfuc,:) + V(1:nfuc,:))*PHk(1:nfuc,:)
        end do
      end if
      ! test hook (env TWIST_TEST=core): start from the core Hamiltonian, no smearing, one iteration only
      call get_environment_variable('TWIST_TEST', tst, ltst, ist)
      if (ist == 0 .and. trim(tst) == 'core') then
        do i = 1, ntw
          call cell_phases(tw_th(:,i), Pc)
          call ao_phases(Pc, PHk)
          Frow(:,:,i) = (T(1:nfuc,:) + V(1:nfuc,:))*PHk(1:nfuc,:)
        end do
        kT = 0.d0
      end if

      allocate(Fref(nfuc,nBas), eref(nfuc,N_cell,ntw), occref(nfuc,N_cell,ntw), Frr(nfuc,nBas))
      do stage = 1, 2
      write(outfile,'(a)') ''
      if (stage == 1) then
        write(outfile,'(a)') '  Stage 1: self-consistent grand-canonical occupations (Fermi smearing annealed to 0)'
        write(outfile,'(a)') '           -> reference band and Fermi surface'
      else
        write(outfile,'(a)') '  Stage 2: occupations fixed by the reference Fermi surface (states with e_ref(k) < mu_ref)'
      end if
      write(outfile,'(a)') '   it        E/cell            mu        <N_el>/cell     max|dF|      kT (Eh)    wall (s)'
      Eold = 0.d0 ; beta = 1.d0
      if (stage == 2) then
        kT = 0.d0 ; beta = 0.6d0
      end if
      do it = 1, 200
        ! -- diagonalize every twist, every k --
        do i = 1, ntw
          call cell_phases(tw_th(:,i), Pc)
          !$omp parallel do default(shared) private(kb,Ak,Sk)
          do kb = 1, N_cell
            call k_block(Frow(:,:,i), Pc, kb, .true., Ak)
            call k_block_S(Pc, kb, Sk)
            call solve_block(Ak, Sk, Ck(:,:,kb,i), eps(:,kb,i))
          end do
          !$omp end parallel do
        end do
        ! -- common Fermi level: sum_twists w * #occupied = nO --
        k = 0
        do i = 1, ntw
          do kb = 1, N_cell
            do a = 1, nfuc
              k = k + 1 ; allE(k) = eps(a,kb,i) ; allW(k) = tw_w(i)
            end do
          end do
        end do
        if (stage == 2) then
          mu_F = Estar
        else if (ist == 0 .and. trim(tst) == 'core') then
          call get_environment_variable('TWIST_MU', tst2, ltst, ist2)
          if (ist2 == 0 .and. ltst > 0) then
            read(tst2,*) mu_F
          else
            call fermi_level(k, allE, allW, dble(nO), mu_F)
          end if
        else if (kT > 1.d-7) then
          call fermi_level_smeared(k, allE, allW, dble(nO), kT, mu_F)
        else
          kT = 0.d0
          call fermi_level(k, allE, allW, dble(nO), mu_F)
        end if
        do i = 1, ntw
          do kb = 1, N_cell
            do a = 1, nfuc
              if (stage == 2) then
                occf(a,kb,i) = occref(a,kb,i)                ! eps sorted ascending in each block
              else if (kT > 0.d0) then
                x = (eps(a,kb,i) - mu_F)/kT
                occf(a,kb,i) = 0.5d0*(1.d0 - tanh(0.5d0*x))
              else
                occf(a,kb,i) = merge(1.d0, 0.d0, eps(a,kb,i) < mu_F)
              end if
            end do
          end do
        end do
        ! -- densities, exchange, energies --
        Qbar = 0.d0 ; Ekin = 0.d0 ; Ene = 0.d0 ; Ex = 0.d0 ; nel = 0.d0
        do i = 1, ntw
          call cell_phases(tw_th(:,i), Pc)
          call ao_phases(Pc, PHk)
          call density(i, Dm, iocc)
          nocc_tw(i) = iocc
          nel = nel + tw_w(i)*2.d0*sum(occf(:,:,i))/dble(N_cell)
          Qbar = Qbar + tw_w(i)*dble(PHk*transpose(Dm))
          call k_rows(Dm, PHk, Krow)
          do nu = 1, nBas
            do a = 1, nfuc
              Ekin = Ekin + tw_w(i)*2.d0*dble(Dm(nu,a)*T(a,nu)*PHk(a,nu))
              Ene   = Ene   + tw_w(i)*2.d0*dble(Dm(nu,a)*V(a,nu)*PHk(a,nu))
              Ex   = Ex   - tw_w(i)*dble(Dm(nu,a)*Krow(a,nu))
            end do
          end do
          Kst(:,:,i) = Krow(:,:)                                        ! keep K of this twist
        end do
        call j_rows(Qbar, jrow)
        Eh = 0.d0
        do nu = 1, nBas
          do a = 1, nfuc
            Eh = Eh + 2.d0*Qbar(a,nu)*jrow(a,nu)
          end do
        end do
        Enn = ENuc/dble(N_cell)
        E = Ekin + Ene + Eh + Ex + Enn
        ! -- new Fock rows, mixed --
        dmax = 0.d0
        do i = 1, ntw
          call cell_phases(tw_th(:,i), Pc)
          call ao_phases(Pc, PHk)
          Krow(:,:) = Kst(:,:,i)
          do nu = 1, nBas
            do a = 1, nfuc
              z = (T(a,nu) + V(a,nu))*PHk(a,nu) + 2.d0*PHk(a,nu)*jrow(a,nu) - Krow(a,nu)
              dmax = max(dmax, abs(z - Frow(a,nu,i)))
              Frow(a,nu,i) = beta*z + (1.d0 - beta)*Frow(a,nu,i)
            end do
          end do
        end do
        write(outfile,'(i5,f20.12,f14.8,f14.8,es13.3,es11.2,f12.2)') it, E, mu_F, nel, dmax, kT, omp_get_wtime()-w0
        flush(outfile)
        if (kT == 0.d0 .and. abs(E - Eold) < 1.d-9 .and. dmax < 1.d-6) exit
        if (ist == 0 .and. trim(tst) == 'core') exit
        if (kT == 0.d0 .and. E > Eold + 1.d-6) beta = max(0.2d0, 0.7d0*beta)   ! oscillation: damp more
        if (it == 1) beta = 0.6d0
        kT = ktf*kT
        if (kT > 0.d0 .and. kT < 1.d-3 .and. abs(E - Eold) < 1.d-10 .and. dmax < 1.d-6) kT = 0.d0   ! annealing done
        Eold = E
      end do
      if (ist == 0 .and. trim(tst) == 'core') exit
      if (stage == 2) exit
      Estage1 = E
      write(outfile,'(a,f20.12,a)') '  Stage 1 E per cell (self-consistent occupations): ', E, ' Eh'

      ! ---------------- reference band: real-space Fock (twist phases removed, twist average) ----------------
      Fref = 0.d0
      do nu = 1, nBas
        do a = 1, nfuc
          wsum = 0.d0 ; fval = 0.d0
          do i = 1, ntw
            call cell_phases(tw_th(:,i), Pc)
            z = Pc(cellof(a), cellof(nu))
            fval = fval + tw_w(i)*dble(Frow(a,nu,i)*conjg(z))
            wsum = wsum + tw_w(i)*abs(z)**2
          end do
          Fref(a,nu) = fval/wsum
        end do
      end do
      ! file twist_reference.dat: if present, its Fock defines the Fermi surface (same surface for every n);
      ! otherwise this run writes it (copy the one of the largest cell next to the other runs)
      inquire(file='twist_reference.dat', exist=have_file)
      if (.not. have_file) then
        open(newunit=ufile, file='twist_reference.dat', status='replace', action='write')
        Rf = maxval(nn)/2
        write(ufile,'(2i6,3f16.10,a)') nfuc, Rf, Lx/dble(nx), Ly/dble(ny), Lz/dble(nz), &
             '   # nfuc, max |cell displacement|, cell lengths (bohr); then dx dy dz a b F(a,b@d) (Eh) S(a,b@d)'
        do nu = 1, nBas
          dxyz = signed_disp(cellof(nu))
          if (any(2*abs(dxyz) == nn .and. nn > 1)) cycle               ! exactly half a box: not written
          do a = 1, nfuc
            write(ufile,'(3i5,2i5,2es24.15)') dxyz, a, fofb(nu), Fref(a,nu), S(a,nu)
          end do
        end do
        close(ufile)
        write(outfile,'(a)') '  reference Fermi surface: from this run (written to twist_reference.dat)'
      else
        write(outfile,'(a)') '  reference Fermi surface: read from twist_reference.dat'
      end if
      ! [Nov 2026] the reference band is ALWAYS evaluated from the file table, with its full range of cell
      ! displacements (not truncated to this box): e_ref(k) is then the same function of k for every size n
      open(newunit=ufile, file='twist_reference.dat', status='old', action='read')
      read(ufile,*) q, Rf
      if (q /= nfuc) error stop 'twist_reference.dat: different number of functions per cell'
      allocate(F3(-Rf:Rf,-Rf:Rf,-Rf:Rf,nfuc,nfuc), S3(-Rf:Rf,-Rf:Rf,-Rf:Rf,nfuc,nfuc), C3(-Rf:Rf,-Rf:Rf,-Rf:Rf))
      F3 = 0.d0 ; S3 = 0.d0 ; C3 = 0.d0 ; nent = 0
      do
        read(ufile,*,iostat=ios) dxyz, a, b, fval, wsum
        if (ios /= 0) exit
        F3(dxyz(1),dxyz(2),dxyz(3),a,b) = fval ; S3(dxyz(1),dxyz(2),dxyz(3),a,b) = wsum
        C3(dxyz(1),dxyz(2),dxyz(3)) = 1.d0 ; nent = nent + 1
      end do
      close(ufile)
      if (nent == 0) error stop 'twist_reference.dat: no entries (old format without the S column?)'
      write(outfile,'(a,i4,a,i8)') '  reference table: cell range +-', Rf, ' , entries ', nent
      ! reference eigenvalues on every twisted grid and occupations with the exact electron number
      k = 0
      do i = 1, ntw
        !$omp parallel do default(shared) private(kb,Ak,Sk,Xk)
        do kb = 1, N_cell
          call ref_block(tw_th(:,i), kb, Ak, Sk)
          call solve_block(Ak, Sk, Xk, eref(:,kb,i))
        end do
        !$omp end parallel do
        do kb = 1, N_cell
          do a = 1, nfuc
            k = k + 1 ; allE(k) = eref(a,kb,i) ; allW(k) = tw_w(i)
          end do
        end do
      end do
      deallocate(F3, S3, C3)
      call fermi_level(k, allE, allW, dble(nO), mu_F)
      ! boundary level E* = highest level needed; states within 1e-9 of it share the remaining electrons
      Estar = -huge(1.d0)
      do q = 1, k
        if (allE(q) < mu_F) Estar = max(Estar, allE(q))
      end do
      Wb = 0.d0 ; Wat = 0.d0
      do q = 1, k
        if (allE(q) < Estar - 1.d-9) Wb = Wb + allW(q)
        if (abs(allE(q) - Estar) <= 1.d-9) Wat = Wat + allW(q)
      end do
      do i = 1, ntw
        do kb = 1, N_cell
          do a = 1, nfuc
            if (eref(a,kb,i) < Estar - 1.d-9) then
              occref(a,kb,i) = 1.d0
            else if (abs(eref(a,kb,i) - Estar) <= 1.d-9) then
              occref(a,kb,i) = (dble(nO) - Wb)/Wat
            else
              occref(a,kb,i) = 0.d0
            end if
          end do
        end do
      end do
      write(outfile,'(a,f14.8,a,f10.6)') '  reference Fermi level mu_ref : ', Estar, '   fraction on the boundary level: ', (dble(nO) - Wb)/Wat
      Estar = mu_F
      end do     ! stage

      write(outfile,'(a)') ''
      write(outfile,'(a,f20.12,a)') '  Stage 1 (self-consistent occupations) per cell : ', Estage1, ' Eh'
      write(outfile,'(a,f20.12,a)') '  Twist-averaged E_HF per cell        : ', E, ' Eh   (stage 2, reference Fermi surface)'
      write(outfile,'(a,f20.12,a)') '  Twist-averaged E_HF per supercell   : ', E*dble(N_cell), ' Eh'
      write(outfile,'(a,f20.12)')   '     kinetic                per cell : ', Ekin
      write(outfile,'(a,f20.12)')   '     nuclear attraction     per cell : ', Ene
      write(outfile,'(a,f20.12)')   '     Hartree (averaged rho) per cell : ', Eh
      write(outfile,'(a,f20.12)')   '     exchange               per cell : ', Ex
      write(outfile,'(a,f20.12)')   '     nuclear repulsion      per cell : ', Enn
      write(outfile,'(a,f20.12)')   '     electrostatic (V+J+Enn) per cell: ', Ene + Eh + Enn
      write(outfile,'(a,f14.8,a,f14.8)') '  Fermi level mu                      : ', mu_F, '   <N_el>/cell ', nel
      write(outfile,'(a,i8,a,i8)')  '  occupied states per twist  min/max  : ', minval(nocc_tw), ' / ', maxval(nocc_tw)
      write(outfile,'(a,f12.2,a)')  '  wall time                           : ', omp_get_wtime()-w0, ' s'
      write(outfile,'(a)') ' ======================================================================================'
      flush(outfile)

contains

      ! signed minimum-image cell displacement of cell c from cell 1 (exactly half a box: +n/2)
      function signed_disp(c) result(d)
        integer, intent(in) :: c
        integer :: d(3), q
        do q = 1, 3
          d(q) = modulo(coord(q,c) - coord(q,1), nn(q))
          if (2*d(q) > nn(q)) d(q) = d(q) - nn(q)
        end do
      end function signed_disp

      ! [Nov 2026] reference Bloch blocks from the file table:  A(a,b) = sum_d F(a,b@d) e^{i k.d a},  k = (2 pi m + theta)/(n a)
      subroutine ref_block(th, kb, A, B)
        double precision, intent(in) :: th(3)
        integer, intent(in) :: kb
        complex(dpc), intent(out) :: A(nfuc,nfuc), B(nfuc,nfuc)
        integer :: dx, dy, dz
        double precision :: kv(3), arg
        complex(dpc) :: z
        kv = 0.d0
        do dx = 1, 3
          if (nn(dx) > 1) kv(dx) = (2.d0*pi*dble(coord(dx,kb)) + th(dx))/dble(nn(dx))
        end do
        A = (0.d0, 0.d0) ; B = (0.d0, 0.d0)
        do dz = -Rf, Rf
          do dy = -Rf, Rf
            do dx = -Rf, Rf
              if (C3(dx,dy,dz) == 0.d0) cycle
              arg = kv(1)*dble(dx) + kv(2)*dble(dy) + kv(3)*dble(dz)
              z = exp(cmplx(0.d0, arg, kind=dpc))
              A = A + z*F3(dx,dy,dz,:,:) ; B = B + z*S3(dx,dy,dz,:,:)
            end do
          end do
        end do
        A = 0.5d0*(A + conjg(transpose(A))) ; B = 0.5d0*(B + conjg(transpose(B)))
      end subroutine ref_block

      ! reference Fock F(a, b@c) from the file table F3; half-box displacements averaged over both images
      double precision function ref_value(a, b, c)
        integer, intent(in) :: a, b, c
        integer :: d(3), sx, sy, sz, cnt, e(3)
        double precision :: acc
        d = signed_disp(c)
        acc = 0.d0 ; cnt = 0
        do sx = -1, 1, 2
          if (sx == 1 .and. .not. (nn(1) > 1 .and. 2*abs(d(1)) == nn(1))) cycle
          do sy = -1, 1, 2
            if (sy == 1 .and. .not. (nn(2) > 1 .and. 2*abs(d(2)) == nn(2))) cycle
            do sz = -1, 1, 2
              if (sz == 1 .and. .not. (nn(3) > 1 .and. 2*abs(d(3)) == nn(3))) cycle
              e = d
              if (sx == 1) e(1) = -d(1)
              if (sy == 1) e(2) = -d(2)
              if (sz == 1) e(3) = -d(3)
              cnt = cnt + 1
              if (all(abs(e) <= Rf)) then
                if (C3(e(1),e(2),e(3)) > 0.d0) acc = acc + F3(e(1),e(2),e(3),a,b)
              end if
            end do
          end do
        end do
        ref_value = acc/dble(cnt)
      end function ref_value

      ! cell-pair phases for twist th: Pc(c1,c2)
      subroutine cell_phases(th, Pc)
        double precision, intent(in) :: th(3)
        complex(dpc), intent(out) :: Pc(N_cell,N_cell)
        integer :: c1, c2, id, d
        complex(dpc) :: p
        do c2 = 1, N_cell
          do c1 = 1, N_cell
            p = (1.d0, 0.d0)
            do id = 1, 3
              if (nn(id) < 2) cycle
              d = modulo(coord(id,c2) - coord(id,c1), nn(id))
              if (2*d > nn(id)) d = d - nn(id)
              if (2*d == nn(id)) then
                p = p*cos(0.5d0*th(id))
              else
                p = p*exp(cmplx(0.d0, th(id)*dble(d)/dble(nn(id)), kind=dpc))
              end if
            end do
            Pc(c1,c2) = p
          end do
        end do
      end subroutine cell_phases

      subroutine ao_phases(Pc, PHk)
        complex(dpc), intent(in) :: Pc(N_cell,N_cell)
        complex(dpc), intent(out) :: PHk(nBas,nBas)
        integer :: mu, nu
        do nu = 1, nBas
          do mu = 1, nBas
            PHk(mu,nu) = Pc(cellof(mu), cellof(nu))
          end do
        end do
      end subroutine ao_phases

      ! e^{i kappa.R_c} for k-point kb and cell c (kappa = 2 pi m / n)
      complex(dpc) function kphase(kb, c)
        integer, intent(in) :: kb, c
        double precision :: arg
        arg = 2.d0*pi*( dble(coord(1,kb)*coord(1,c))/dble(nx) + dble(coord(2,kb)*coord(2,c))/dble(ny) &
                      + dble(coord(3,kb)*coord(3,c))/dble(nz) )
        kphase = exp(cmplx(0.d0, arg, kind=dpc))
      end function kphase

      ! Bloch block of a row matrix R(a, nu) (rows in cell 1), phases already included if withph
      subroutine k_block(R, Pc, kb, withph, A)
        complex(dpc), intent(in) :: R(nfuc,nBas), Pc(N_cell,N_cell)
        integer, intent(in) :: kb
        logical, intent(in) :: withph
        complex(dpc), intent(out) :: A(nfuc,nfuc)
        integer :: c, b
        complex(dpc) :: f
        A = (0.d0, 0.d0)
        do c = 1, N_cell
          f = kphase(kb, c)
          if (.not. withph) f = f*Pc(1,c)
          do b = 1, nfuc
            A(:,b) = A(:,b) + f*R(:, b + (c-1)*nfuc)
          end do
        end do
        A = 0.5d0*(A + conjg(transpose(A)))
      end subroutine k_block

      subroutine k_block_S(Pc, kb, A)
        complex(dpc), intent(in) :: Pc(N_cell,N_cell)
        integer, intent(in) :: kb
        complex(dpc), intent(out) :: A(nfuc,nfuc)
        integer :: c, b
        complex(dpc) :: f
        A = (0.d0, 0.d0)
        do c = 1, N_cell
          f = kphase(kb, c)*Pc(1,c)
          do b = 1, nfuc
            A(:,b) = A(:,b) + f*S(1:nfuc, b + (c-1)*nfuc)
          end do
        end do
        A = 0.5d0*(A + conjg(transpose(A)))
      end subroutine k_block_S

      ! generalized eigenproblem by canonical orthogonalization (dropped vectors -> huge energy)
      subroutine solve_block(A, Sm, C, e)
        complex(dpc), intent(in) :: A(nfuc,nfuc), Sm(nfuc,nfuc)
        complex(dpc), intent(out) :: C(nfuc,nfuc)
        double precision, intent(out) :: e(nfuc)
        complex(dpc) :: U(nfuc,nfuc), X(nfuc,nfuc), Fp(nfuc,nfuc), wk(4*nfuc)
        double precision :: sv(nfuc), rw(3*nfuc), ee(nfuc)
        integer :: q, nk, info
        U = Sm
        call zheev('V','U', nfuc, U, nfuc, sv, wk, 4*nfuc, rw, info)
        if (info /= 0) error stop 'twist_RHF_SAO: zheev failed (S_k)'
        nk = 0
        X = (0.d0, 0.d0)
        do q = 1, nfuc
          if (sv(q) > s_min*maxval(sv)) then
            nk = nk + 1 ; X(:,nk) = U(:,q)/sqrt(sv(q))
          end if
        end do
        Fp = (0.d0, 0.d0)
        Fp(1:nk,1:nk) = matmul(conjg(transpose(X(:,1:nk))), matmul(A, X(:,1:nk)))
        call zheev('V','U', nk, Fp(1:nk,1:nk), nk, ee(1:nk), wk, 4*nfuc, rw, info)
        if (info /= 0) error stop 'twist_RHF_SAO: zheev failed (F_k)'
        C = (0.d0, 0.d0) ; e = huge(1.d0)
        C(:,1:nk) = matmul(X(:,1:nk), Fp(1:nk,1:nk))
        e(1:nk) = ee(1:nk)
      end subroutine solve_block

      ! mu such that sum_states w * [e < mu] = target (midpoint between the last occupied and the next)
      subroutine fermi_level(nst, el, wl, target, mu_out)
        integer, intent(in) :: nst
        double precision, intent(in) :: el(nst), wl(nst), target
        double precision, intent(out) :: mu_out
        integer, allocatable :: perm(:)
        double precision :: acc
        integer :: q
        allocate(perm(nst))
        call sort_index(nst, el, perm)
        acc = 0.d0
        do q = 1, nst
          acc = acc + wl(perm(q))
          if (acc >= target - 1.d-9) exit
        end do
        if (q >= nst) then
          mu_out = el(perm(nst)) + 1.d0
        else
          mu_out = 0.5d0*(el(perm(q)) + el(perm(q+1)))
        end if
        deallocate(perm)
      end subroutine fermi_level

      subroutine fermi_level_smeared(nst, el, wl, target, kT, mu_out)
        integer, intent(in) :: nst
        double precision, intent(in) :: el(nst), wl(nst), target, kT
        double precision, intent(out) :: mu_out
        double precision :: lo, hi, mid, cnt
        integer :: q, iter
        lo = minval(el) - 50.d0*kT ; hi = maxval(el, mask=el < 1.d300) + 50.d0*kT
        do iter = 1, 200
          mid = 0.5d0*(lo + hi) ; cnt = 0.d0
          do q = 1, nst
            if (el(q) < 1.d300) cnt = cnt + wl(q)*0.5d0*(1.d0 - tanh(0.5d0*(el(q) - mid)/kT))
          end do
          if (cnt > target) then
            hi = mid
          else
            lo = mid
          end if
          if (hi - lo < 1.d-14) exit
        end do
        mu_out = 0.5d0*(lo + hi)
      end subroutine fermi_level_smeared

      subroutine sort_index(nst, el, perm)        ! heap sort of indices by el
        integer, intent(in) :: nst
        double precision, intent(in) :: el(nst)
        integer, intent(out) :: perm(nst)
        integer :: q, r, l, ir, ip
        double precision :: v
        do q = 1, nst
          perm(q) = q
        end do
        if (nst < 2) return
        l = nst/2 + 1 ; ir = nst
        do
          if (l > 1) then
            l = l - 1 ; ip = perm(l) ; v = el(ip)
          else
            ip = perm(ir) ; v = el(ip) ; perm(ir) = perm(1) ; ir = ir - 1
            if (ir == 1) then
              perm(1) = ip ; exit
            end if
          end if
          q = l ; r = l + l
          do while (r <= ir)
            if (r < ir) then
              if (el(perm(r)) < el(perm(r+1))) r = r + 1
            end if
            if (v < el(perm(r))) then
              perm(q) = perm(r) ; q = r ; r = r + r
            else
              r = ir + 1
            end if
          end do
          perm(q) = ip
        end do
      end subroutine sort_index

      ! TI-gauge density D'(l,s) = (1/N_cell) sum_k e^{i kappa (R_l - R_s)} G_k(a_l, a_s),  G_k = sum_{occ} C C^*
      subroutine density(i, D, nocc)
        integer, intent(in) :: i
        complex(dpc), intent(out) :: D(nBas,nBas)
        integer, intent(out) :: nocc
        integer :: kb, b, c, l, sg
        nocc = 0
        Gk = (0.d0, 0.d0)
        do kb = 1, N_cell
          do b = 1, nfuc
            if (occf(b,kb,i) > 1.d-14) then
              if (occf(b,kb,i) > 0.5d0) nocc = nocc + 1
              do sg = 1, nfuc
                Gk(:,sg,kb) = Gk(:,sg,kb) + occf(b,kb,i)*Ck(:,b,kb,i)*conjg(Ck(sg,b,kb,i))
              end do
            end if
          end do
        end do
        !$omp parallel do default(shared) private(c,kb)
        do c = 1, N_cell
          dC(:,:,c) = (0.d0, 0.d0)
          do kb = 1, N_cell
            dC(:,:,c) = dC(:,:,c) + kphase(kb, c)*Gk(:,:,kb)
          end do
          dC(:,:,c) = dC(:,:,c)/dble(N_cell)
        end do
        !$omp end parallel do
        do sg = 1, nBas
          do l = 1, nBas
            D(l,sg) = dC(fofb(l), fofb(sg), dcell(cellof(sg), cellof(l)))
          end do
        end do
      end subroutine density

      ! j(a,nu) = sum_{l,s} (a nu|l s) Q(l,s)   (Q real)
      subroutine j_rows(Q, jr)
        double precision, intent(in) :: Q(nBas,nBas)
        double precision, intent(out) :: jr(nfuc,nBas)
        integer :: a, nu, l, sg
        double precision :: acc
        !$omp parallel do default(shared) private(a,nu,l,sg,acc) collapse(2)
        do nu = 1, nBas
          do a = 1, nfuc
            acc = 0.d0
            do sg = 1, nBas
              do l = 1, nBas
                acc = acc + ERI(a,nu,l,sg)*Q(l,sg)
              end do
            end do
            jr(a,nu) = acc
          end do
        end do
        !$omp end parallel do
      end subroutine j_rows

      ! K(a,nu) = sum_{l,s} (a l|s nu) PH(a,l) PH(s,nu) D(l,s)
      subroutine k_rows(D, PH, Kr)
        complex(dpc), intent(in) :: D(nBas,nBas), PH(nBas,nBas)
        complex(dpc), intent(out) :: Kr(nfuc,nBas)
        integer :: a, nu, l, sg
        double precision :: sre, sim, tre, tim
        double precision, allocatable :: Are(:,:), Aim(:,:)
        allocate(Are(nBas,nBas), Aim(nBas,nBas))
        do a = 1, nfuc
          do sg = 1, nBas
            do l = 1, nBas
              Are(l,sg) = dble(PH(a,l)*D(l,sg)) ; Aim(l,sg) = aimag(PH(a,l)*D(l,sg))
            end do
          end do
          !$omp parallel do default(shared) private(nu,l,sg,sre,sim,tre,tim)
          do nu = 1, nBas
            tre = 0.d0 ; tim = 0.d0
            do sg = 1, nBas
              sre = 0.d0 ; sim = 0.d0
              do l = 1, nBas
                sre = sre + Are(l,sg)*ERI(a,l,sg,nu)
                sim = sim + Aim(l,sg)*ERI(a,l,sg,nu)
              end do
              tre = tre + dble(PH(sg,nu))*sre - aimag(PH(sg,nu))*sim
              tim = tim + dble(PH(sg,nu))*sim + aimag(PH(sg,nu))*sre
            end do
            Kr(a,nu) = cmplx(tre, tim, kind=dpc)
          end do
          !$omp end parallel do
        end do
        deallocate(Are, Aim)
      end subroutine k_rows

end subroutine twist_RHF_SAO
