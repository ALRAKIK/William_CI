!=====================================================================!
!                                                                     !
!          1D toroidal : NUCLEAR ATTRACTION INTEGRALS                !
!                                                                     !
!  1/|r-C|_E = 2/sqrt(pi) Int_0^inf exp(-t^2 |r-C|_E^2) dt            !
!                                                                     !
!  x (Clifford): exp(-t^2 |x-xC|_E^2) = sum_n Ĩ_n(C_t) e^{i n a(x-xC)}!
!     C_t = 2t^2/a^2 ,  Ĩ_n = e^{-C} I_n  -> precomputed table t_bessel   !
!     The pair  g_A g_B = kx a^-(i+l) sum_k d_k e^{iku} e^{A cos u},   !
!     u = a(x-xP), gives, ONCE per primitive pair (t-independent):    !
!        G_n = sum_k d_k Ĩ_{|n+k|}(A)                                 !
!     and per nucleus  F_n = Re[ e^{i n a(xP-xC)} G_n ], so that      !
!        Ix(t) = kx a^-(i+l) L [ F_0 Ĩ_0(C_t) + 2 sum_n F_n Ĩ_n(C_t) ]!
!     i.e. a dot product with the precomputed column of t.            !
!                                                                     !
!  y,z (real): Gaussian product with exp(-t^2 (y-yC)^2):              !
!     exp(-p eta yPC^2) * S_real(j,m; p+t^2, yPA-eta yPC, yPB-eta yPC) !
!     eta = t^2/(p+t^2)                                               !
!                                                                     !
!  t-integral: same nodes/weights as the ERIs (64 or 128 points).     !
!                                                                     !
!  nuclear_attraction_integral_toroidal_1D  one pair, all nuclei      !
!  nuclear_attraction_matrix_toroidal_1D    full matrix (OpenMP)      !
!=====================================================================!

subroutine nuclear_attraction_integral_toroidal_1D(l1,l2,n_atoms,geometry,atoms,r1,r2,AO1,AO2,V)

      use torus_init
      use classification_ERI
      use atom_basis
      use clifford_tools

      implicit none

      integer           ,intent(in)    :: l1(3) , l2(3)
      integer           ,intent(in)    :: n_atoms
      double precision  ,intent(in)    :: geometry(n_atoms,3)
      type(atom)        ,intent(in)    :: atoms(n_atoms)
      double precision  ,intent(in)    :: r1(3) , r2(3)
      type(ERI_function),intent(in)    :: AO1 , AO2
      double precision  ,intent(out)   :: V

      double precision, parameter      :: pi = 3.14159265358979323846d0
      integer                          :: i , j , k , n , c , M , Np , nm , kk
      double precision                 :: alpha , beta , const , pref
      double precision                 :: x1 , x2 , y1 , y2 , z1 , z2 , X , Y , Z
      double precision                 :: ax2 , inv_ax , px , xp , A , kx , ky , kz
      double precision                 :: albe , inv_albe , mu , yp , zp
      double precision                 :: ypa , ypb , zpa , zpb , ypc , zpc
      double precision                 :: phi , cn , sn , c1 , s1 , tmp
      double precision                 :: t2 , ppt , eta , sx , sy , sz , ey , ez , res
      complex(8)                       :: acc
      complex(8)      , allocatable    :: d(:) , G(:)
      double precision, allocatable    :: Iv(:) , F(:)

      !-----------------------------------------------------------------!

      x1 = r1(1) ; x2 = r2(1)
      y1 = r1(2) ; y2 = r2(2)
      z1 = r1(3) ; z2 = r2(3)
      X  = x1 - x2 ; Y = y1 - y2 ; Z = z1 - z2

      ax2    = ax * ax
      inv_ax = 1.d0 / ax
      M      = l1(1) + l2(1)

      allocate(d(-M:M))

      V = 0.d0

      do i = 1 , size(AO1%exponent)
        alpha = AO1%exponent(i)
        do j = 1 , size(AO2%exponent)
          beta  = AO2%exponent(j)
          const = AO1%coefficient(i) * AO2%coefficient(j)

          ! ---------- Clifford pair (x), t-independent ---------- !
          call bary_exponent_x        (alpha,beta,X,px)
          call bary_center_toroidal_x (alpha,beta,x1,x2,xp)

          A  = 2.d0*px/ax2
          kx = dexp(-2.d0*(alpha + beta - px)/ax2)

          call fourier_sin_product(l1(1), l2(1), ax*(xp-x1), ax*(xp-x2), d)

          Np = bessel_nmax_estimate(A)
          allocate(Iv(0:Np+M), G(0:Np), F(0:Np))
          call bessel_scaled_vector(Np+M, A, Iv)

          do n = 0 , Np
            acc = (0.d0,0.d0)
            do k = -M , M
              kk = abs(n+k)
              acc = acc + d(k) * Iv(kk)
            end do
            G(n) = acc
          end do

          ! ---------- real pair (y,z) ---------- !
          albe     = alpha + beta
          inv_albe = 1.d0 / albe
          mu       = alpha * beta * inv_albe
          yp       = (alpha*y1 + beta*y2) * inv_albe
          zp       = (alpha*z1 + beta*z2) * inv_albe
          ky       = dexp(-mu*Y*Y)
          kz       = dexp(-mu*Z*Z)
          ypa = yp - y1 ; ypb = yp - y2
          zpa = zp - z1 ; zpb = zp - z2

          pref = const * kx * ky * kz * inv_ax**M * Lx * 2.d0/dsqrt(pi)

          do c = 1 , n_atoms

            ! F_n = Re[ e^{i n phi} G_n ],  phi = a (xP - xC)
            phi = ax * (xp - geometry(c,1))
            c1 = dcos(phi) ; s1 = dsin(phi)
            cn = 1.d0      ; sn = 0.d0
            do n = 0 , Np
              if (n > 0) then
                if (mod(n,32) == 0) then
                  cn = dcos(n*phi) ; sn = dsin(n*phi)
                else
                  tmp = cn*c1 - sn*s1
                  sn  = sn*c1 + cn*s1
                  cn  = tmp
                end if
              end if
              F(n) = dble(G(n))*cn - aimag(G(n))*sn
            end do

            ypc = yp - geometry(c,2)
            zpc = zp - geometry(c,3)

            res = 0.d0
            do k = 1 , nq_t
              t2  = t_node(k) * t_node(k)
              ppt = albe + t2
              eta = t2 / ppt

              nm = min(Np, t_nmax(k,1))
              sx = F(0) * t_bessel(0,k,1)
              if (nm >= 1) sx = sx + 2.d0 * dot_product(F(1:nm), t_bessel(1:nm,k,1))

              ey = dexp(-albe*eta*ypc*ypc)
              ez = dexp(-albe*eta*zpc*zpc)
              sy = overlap_1D_real(l1(2), l2(2), ppt, ypa - eta*ypc, ypb - eta*ypc)
              sz = overlap_1D_real(l1(3), l2(3), ppt, zpa - eta*zpc, zpb - eta*zpc)

              res = res + t_weight(k) * sx * ey * sy * ez * sz
            end do

            V = V - dble(atoms(c)%charge) * pref * res

          end do

          deallocate(Iv, G, F)

        end do
      end do

      deallocate(d)

end subroutine nuclear_attraction_integral_toroidal_1D

!=====================================================================!
!  Nuclear attraction matrix: rows i = 1..nfuc (OpenMP over pairs),   !
!  the rest by fill_with_symmetry                                     !
!=====================================================================!
subroutine nuclear_attraction_matrix_toroidal_1D(number_of_atoms,nBas,geometry,atoms,AO,NA)

      use files
      use atom_basis
      use torus_init
      use classification_ERI
      use omp_lib
      use keywords
      use unitcell_module
      use clifford_tools, only : init_clifford, orbital_to_lxyz, setup_t_tables, bessel_nmax_estimate

      implicit none

      integer           ,intent(in)  :: number_of_atoms , nBas
      double precision  ,intent(in)  :: geometry(number_of_atoms,3)
      type(atom)        ,intent(in)  :: atoms(number_of_atoms)
      type(ERI_function),intent(in)  :: AO(nBas)
      double precision  ,intent(out) :: NA(nBas,nBas)

      integer                        :: i , j , ij , total_ij_pairs , lmax , ncap
      integer                        :: lxyz(3,nBas)
      integer, allocatable           :: i_index(:) , j_index(:)
      double precision               :: r1(3) , r2(3) , val , emax

      !-----------------------------------------------------------------!

      call init_clifford()

      lmax = 0 ; emax = 0.d0
      do i = 1 , nBas
        call orbital_to_lxyz(AO(i)%orbital, lxyz(:,i))
        lmax = max(lmax, lxyz(1,i))
        emax = max(emax, maxval(AO(i)%exponent))
      end do

      ! largest pair Bessel order needed by any AO pair (A <= 4 emax/ax^2)
      ncap = bessel_nmax_estimate(4.d0*emax/(ax*ax)) + 2*lmax + 1
      call setup_t_tables(ncap, 1)

      NA(:,:) = 0.d0

      total_ij_pairs = nfuc * nBas
      allocate(i_index(total_ij_pairs), j_index(total_ij_pairs))
      ij = 0
      do i = 1, nfuc
        do j = 1, nBas
          ij = ij + 1
          i_index(ij) = i
          j_index(ij) = j
        end do
      end do

      !$omp parallel do default(shared) private(ij,i,j,r1,r2,val) schedule(dynamic,4)
      do ij = 1, total_ij_pairs
        i  = i_index(ij)
        j  = j_index(ij)
        r1 = (/ AO(i)%x , AO(i)%y , AO(i)%z /)
        r2 = (/ AO(j)%x , AO(j)%y , AO(j)%z /)
        call nuclear_attraction_integral_toroidal_1D(lxyz(:,i),lxyz(:,j),number_of_atoms,geometry,atoms, &
                                                         r1,r2,AO(i),AO(j),val)
        NA(i,j) = val
      end do
      !$omp end parallel do

      deallocate(i_index, j_index)

      do i = 1 , nfuc
        do j = 1 , nBas
          if (abs(NA(i,j)) < 1e-14) NA(i,j) = 0.d0
        end do
      end do

      call fill_with_symmetry(nBas, NA)

      if (c_NA) then
        write(outfile,'(a)') ""
        write(outfile,'(a)') "-----------------------------------------------"
        write(outfile,'(a)') "The nuclear attraction matrix                 "
        write(outfile,'(a)') "-----------------------------------------------"
        write(outfile,'(a)') ""
        do i = 1 , nBas - 1
          do j = i , nBas
            write(outfile,'(i3,i3, 6x,1000(f16.12,2x))') i , j, NA(i,j)
          end do
        end do
      end if

end subroutine nuclear_attraction_matrix_toroidal_1D
