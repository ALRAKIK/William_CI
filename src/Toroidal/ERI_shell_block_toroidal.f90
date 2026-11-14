!=====================================================================!
!                                                                     !
!     Toroidal 1D/2D/3D : TWO-ELECTRON INTEGRALS, SHELL QUARTET KERNEL!
!                                                                     !
!   (ab|cd) = 2/sqrt(pi) Int_0^inf dt  Ix(t) Iy(t) Iz(t)              !
!                                                                     !
!   periodic d   : exp(-t^2|x1-x2|_E^2) = sum_n Ĩ_n(C_t) e^{in a(x1-x2)}
!       bra pair -> G_n = sum_k d_k  Ĩ_{|n+k|}(A_P)                    !
!       ket pair -> H_n = sum_k d'_k Ĩ_{|n-k|}(A_Q)                    !
!       F_n = Re[ G_n H_n e^{i n a (P-Q)} ]                           !
!       Ix(t) = [F_0 Ĩ_0 + 2 sum_n F_n Ĩ_n(C_t)]  (one dot product)    !
!   real d       : McMurchie-Davidson                                 !
!       S_N = sum_{tau+nu=N} (-1)^tau E_tau E'_nu                     !
!       Iy(t) = sum_N S_N g_N(t),  g_N = pi D e^{-wR^2} w^{N/2} H_N   !
!                                                                     !
!   All Cartesian components of a shell quartet at once: the 1D       !
!   factors depend only on the powers along each direction, so they   !
!   are computed once per power combination (ia,ib,ic,id) and then    !
!   multiplied for every component quartet.                           !
!                                                                     !
!   ERI_shell_block_toroidal (ndp,nA,SA,nB,SB,nC,SC,nD,SD,blk)       !
!       ndp          : periodic directions 1..ndp (1, 2 or 3)         !
!       SA(1:nA) ... : the Cartesian components of each shell         !
!       blk(a,b,c,d) : (ab|cd) for all components                     !
!   Needs setup_t_tables(ncap,ndp) (called by the driver).            !
!=====================================================================!

subroutine ERI_shell_block_toroidal(ndp,nA,SA,nB,SB,nC,SC,nD,SD,blk)

      use torus_init
      use classification_ERI
      use clifford_tools

      implicit none

      integer           ,intent(in)  :: ndp                 ! periodic directions: 1..ndp
      integer           ,intent(in)  :: nA , nB , nC , nD
      type(ERI_function),intent(in)  :: SA(nA) , SB(nB) , SC(nC) , SD(nD)
      double precision  ,intent(out) :: blk(nA,nB,nC,nD)

      double precision, parameter :: pi = 3.14159265358979323846d0
      double precision, parameter :: screen_tol = 2.22d-16

      integer          :: lsA, lsB, lsC, lsD, ME, NN
      integer          :: pwA(3,nA), pwB(3,nB), pwC(3,nC), pwD(3,nD)
      integer          :: n1, n2, n3, n4, nket, i, j, k, l, kl, d, n, kk, kq
      integer          :: ia, ib, ic, id, ca, cb, cc, cd, cmb, ncomb, nw, off, nm, nt, nq
      integer          :: cx, cy, cz
      double precision :: av(3), Lc(3), rA(3), rB(3), rC(3), rD(3)
      double precision :: e1, e2, esum, mu, pe, cen, Aarg, fac, fac_const, ket_max, bra_fac, cf
      double precision :: t2, D2, w, Rr, phi, c1, s1, cn, sn, tmp, gm1, g0, gp1, val, p, q, sgn
      double precision :: cAmax, cBmax, cCmax, cDmax
      complex(8)       :: acc

      ! ket storage
      double precision, allocatable :: k_fac(:), k_cen(:,:), k_q(:), k_E(:,:,:,:,:)
      integer,          allocatable :: k_N(:,:), k_off(:,:)
      complex(8),       allocatable :: Hall(:)
      ! bra / work
      double precision, allocatable :: b_cen(:), b_E(:,:,:,:), Iv(:), F(:)
      integer,          allocatable :: b_N(:)
      complex(8),       allocatable :: G(:,:,:,:), dco(:)
      double precision, allocatable :: Idir(:,:,:), SNv(:), gt(:,:), cph(:), sph(:), wq(:)
      double precision, allocatable :: coA(:,:), coB(:,:), coC(:,:), coD(:,:)

      blk = 0.d0
      if (.not. allocated(t_bessel)) error stop "ERI_shell_block_toroidal: t-table not initialized"

      do ca = 1 , nA ; call orbital_to_lxyz(SA(ca)%orbital, pwA(:,ca)) ; end do
      do cb = 1 , nB ; call orbital_to_lxyz(SB(cb)%orbital, pwB(:,cb)) ; end do
      do cc = 1 , nC ; call orbital_to_lxyz(SC(cc)%orbital, pwC(:,cc)) ; end do
      do cd = 1 , nD ; call orbital_to_lxyz(SD(cd)%orbital, pwD(:,cd)) ; end do
      lsA = sum(pwA(:,1)) ; lsB = sum(pwB(:,1)) ; lsC = sum(pwC(:,1)) ; lsD = sum(pwD(:,1))
      ME  = max(lsA+lsB, lsC+lsD)
      NN  = lsA + lsB + lsC + lsD
      ncomb = (lsA+1)*(lsB+1)*(lsC+1)*(lsD+1)

      av = (/ ax , ay , az /)
      Lc = (/ Lx , Ly , Lz /)
      rA = (/ SA(1)%x , SA(1)%y , SA(1)%z /)
      rB = (/ SB(1)%x , SB(1)%y , SB(1)%z /)
      rC = (/ SC(1)%x , SC(1)%y , SC(1)%z /)
      rD = (/ SD(1)%x , SD(1)%y , SD(1)%z /)

      fac_const = 2.d0/dsqrt(pi)
      do d = 1 , ndp
        fac_const = fac_const * Lc(d)**2
      end do

      n1 = size(SA(1)%exponent) ; n2 = size(SB(1)%exponent)
      n3 = size(SC(1)%exponent) ; n4 = size(SD(1)%exponent)
      nket = n3*n4
      nq   = nq_t
      nw   = ubound(t_bessel,1) + ME + 2

      allocate(k_fac(nket), k_cen(3,nket), k_q(nket), k_N(3,nket), k_off(3,nket))
      allocate(k_E(0:NMAX_CLIFF,0:lsC,0:lsD,3,nket))   ! hermite_E_1D fills 0:NMAX_CLIFF
      allocate(b_cen(3), b_N(3), b_E(0:NMAX_CLIFF,0:lsA,0:lsB,3))
      allocate(Iv(0:nw), F(0:nw), dco(-ME:ME), G(0:nw,0:lsA,0:lsB,ndp))
      allocate(Idir(nq,0:ncomb-1,3), SNv(0:NN), gt(0:NN,nq), cph(0:nw), sph(0:nw), wq(nq))
      wq = t_weight(1:nq)

      allocate(coA(n1,nA), coB(n2,nB), coC(n3,nC), coD(n4,nD))
      do ca = 1 , nA ; coA(:,ca) = SA(ca)%coefficient(1:n1) ; end do
      do cb = 1 , nB ; coB(:,cb) = SB(cb)%coefficient(1:n2) ; end do
      do cc = 1 , nC ; coC(:,cc) = SC(cc)%coefficient(1:n3) ; end do
      do cd = 1 , nD ; coD(:,cd) = SD(cd)%coefficient(1:n4) ; end do

      cCmax = 0.d0 ; cDmax = 0.d0
      !-------------------------- ket pairs --------------------------!
      kl = 0 ; off = 0 ; ket_max = 0.d0
      do k = 1 , n3
        do l = 1 , n4
          kl = kl + 1
          e1 = SC(1)%exponent(k) ; e2 = SD(1)%exponent(l)
          esum = e1 + e2 ; mu = e1*e2/esum
          k_q(kl) = esum
          cf = 1.d0
          do d = 1 , 3
            if (d <= ndp) then
              call bary_dir(d, e1, e2, rC(d), rD(d), pe, cen)
              k_cen(d,kl) = cen
              cf = cf * dexp(-2.d0*(esum - pe)/(av(d)*av(d)))
              k_N(d,kl)   = min(bessel_nmax_estimate(2.d0*pe/(av(d)*av(d))), nw - ME - 2)
              k_off(d,kl) = off
              off = off + (k_N(d,kl)+1)*(lsC+1)*(lsD+1)
            else
              cen = (e1*rC(d) + e2*rD(d))/esum
              k_cen(d,kl) = cen
              cf = cf * dexp(-mu*(rC(d)-rD(d))**2)
              do ic = 0 , lsC
                do id = 0 , lsD
                  call hermite_E_1D(ic, id, esum, cen - rC(d), cen - rD(d), k_E(:,ic,id,d,kl))
                end do
              end do
            end if
          end do
          k_fac(kl) = cf                     ! coefficient-free
          ket_max = max(ket_max, dabs(cf) * maxval(abs(coC(k,:))) * maxval(abs(coD(l,:))))
        end do
      end do

      allocate(Hall(0:max(off,1)-1))
      kl = 0
      do k = 1 , n3
        do l = 1 , n4
          kl = kl + 1
          e1 = SC(1)%exponent(k) ; e2 = SD(1)%exponent(l)
          do d = 1 , ndp
            call bary_dir(d, e1, e2, rC(d), rD(d), pe, cen)
            Aarg = 2.d0*pe/(av(d)*av(d))
            call bessel_scaled_vector_fast(k_N(d,kl)+lsC+lsD, Aarg, Iv(0:k_N(d,kl)+lsC+lsD))
            do ic = 0 , lsC
              do id = 0 , lsD
                call fourier_sin_product(ic, id, av(d)*(cen-rC(d)), av(d)*(cen-rD(d)), dco(-(ic+id):(ic+id)))
                do n = 0 , k_N(d,kl)
                  acc = (0.d0,0.d0)
                  do kk = -(ic+id) , ic+id
                    acc = acc + dco(kk) * Iv(abs(n-kk))
                  end do
                  Hall(k_off(d,kl) + n + (k_N(d,kl)+1)*(ic + (lsC+1)*id)) = acc
                end do
              end do
            end do
          end do
        end do
      end do

      !------------------------ bra pairs / quartets -----------------!
      do i = 1 , n1
        do j = 1 , n2
          e1 = SA(1)%exponent(i) ; e2 = SB(1)%exponent(j)
          esum = e1 + e2 ; mu = e1*e2/esum
          p  = esum
          cAmax = maxval(abs(coA(i,:)))
          cBmax = maxval(abs(coB(j,:)))
          bra_fac = fac_const
          do d = 1 , 3
            if (d <= ndp) then
              call bary_dir(d, e1, e2, rA(d), rB(d), pe, cen)
              b_cen(d) = cen
              bra_fac  = bra_fac * dexp(-2.d0*(esum - pe)/(av(d)*av(d)))
              Aarg     = 2.d0*pe/(av(d)*av(d))
              b_N(d)   = min(bessel_nmax_estimate(Aarg), nw - ME - 2)
            else
              cen = (e1*rA(d) + e2*rB(d))/esum
              b_cen(d) = cen
              bra_fac  = bra_fac * dexp(-mu*(rA(d)-rB(d))**2)
            end if
          end do
          if (dabs(bra_fac)*cAmax*cBmax*ket_max <= screen_tol) cycle

          do d = 1 , 3
            if (d <= ndp) then
              call bary_dir(d, e1, e2, rA(d), rB(d), pe, cen)
              Aarg = 2.d0*pe/(av(d)*av(d))
              call bessel_scaled_vector_fast(b_N(d)+lsA+lsB, Aarg, Iv(0:b_N(d)+lsA+lsB))
              do ia = 0 , lsA
                do ib = 0 , lsB
                  call fourier_sin_product(ia, ib, av(d)*(cen-rA(d)), av(d)*(cen-rB(d)), dco(-(ia+ib):(ia+ib)))
                  do n = 0 , b_N(d)
                    acc = (0.d0,0.d0)
                    do kk = -(ia+ib) , ia+ib
                      acc = acc + dco(kk) * Iv(abs(n+kk))
                    end do
                    G(n,ia,ib,d) = acc
                  end do
                end do
              end do
            else
              do ia = 0 , lsA
                do ib = 0 , lsB
                  call hermite_E_1D(ia, ib, p, b_cen(d) - rA(d), b_cen(d) - rB(d), b_E(:,ia,ib,d))
                end do
              end do
            end if
          end do

          kl = 0
          do k = 1 , n3
            do l = 1 , n4
              kl = kl + 1
              fac = bra_fac * k_fac(kl)
              if (dabs(fac)*cAmax*cBmax*maxval(abs(coC(k,:)))*maxval(abs(coD(l,:))) <= screen_tol) cycle
              q = k_q(kl)

              !---- 1D factors for every power combination ----!
              do d = 1 , 3
                if (d <= ndp) then
                  nm  = min(b_N(d), k_N(d,kl))
                  phi = av(d)*(b_cen(d) - k_cen(d,kl))
                  c1 = dcos(phi) ; s1 = dsin(phi)
                  cn = 1.d0      ; sn = 0.d0
                  do n = 0 , nm
                    if (n > 0) then
                      if (mod(n,32) == 0) then
                        cn = dcos(n*phi) ; sn = dsin(n*phi)
                      else
                        tmp = cn*c1 - sn*s1
                        sn  = sn*c1 + cn*s1
                        cn  = tmp
                      end if
                    end if
                    cph(n) = cn ; sph(n) = sn
                  end do
                  do id = 0 , lsD
                   do ic = 0 , lsC
                    do ib = 0 , lsB
                     do ia = 0 , lsA
                      cmb = ia + (lsA+1)*(ib + (lsB+1)*(ic + (lsC+1)*id))
                      do n = 0 , nm
                        acc  = G(n,ia,ib,d) * Hall(k_off(d,kl) + n + (k_N(d,kl)+1)*(ic + (lsC+1)*id))
                        F(n) = dble(acc)*cph(n) - aimag(acc)*sph(n)
                      end do
                      sgn = av(d)**(-(ia+ib+ic+id))
                      do kq = 1 , nq
                        nt = min(nm, t_nmax(kq,d))
                        val = F(0)*t_bessel(0,kq,d)
                        if (nt >= 1) val = val + 2.d0*dot_product(F(1:nt), t_bessel(1:nt,kq,d))
                        Idir(kq,cmb,d) = val * sgn
                      end do
                     end do
                    end do
                   end do
                  end do
                else
                  Rr = b_cen(d) - k_cen(d,kl)
                  do kq = 1 , nq
                    t2 = t_node(kq)*t_node(kq)
                    D2 = 1.d0/(p*q + (p+q)*t2)
                    w  = p*q*t2*D2
                    gt(0,kq) = pi*dsqrt(D2)*dexp(-w*Rr*Rr)
                    gm1 = 0.d0 ; g0 = gt(0,kq)
                    do n = 0 , NN - 1
                      gp1 = 2.d0*w*Rr*g0 - 2.d0*dble(n)*w*gm1
                      gt(n+1,kq) = gp1
                      gm1 = g0 ; g0 = gp1
                    end do
                  end do
                  do id = 0 , lsD
                   do ic = 0 , lsC
                    do ib = 0 , lsB
                     do ia = 0 , lsA
                      cmb = ia + (lsA+1)*(ib + (lsB+1)*(ic + (lsC+1)*id))
                      SNv = 0.d0
                      do n = 0 , ia+ib
                        do kk = 0 , ic+id
                          SNv(n+kk) = SNv(n+kk) + merge(-1.d0,1.d0,mod(n,2)==1) * b_E(n,ia,ib,d) * k_E(kk,ic,id,d,kl)
                        end do
                      end do
                      do kq = 1 , nq
                        Idir(kq,cmb,d) = dot_product(SNv(0:ia+ib+ic+id), gt(0:ia+ib+ic+id,kq))
                      end do
                     end do
                    end do
                   end do
                  end do
                end if
              end do
              ! fold the weights into x
              do cmb = 0 , ncomb-1
                Idir(:,cmb,1) = Idir(:,cmb,1) * wq
              end do

              !---- all component quartets ----!
              do cd = 1 , nD
               do cc = 1 , nC
                do cb = 1 , nB
                 do ca = 1 , nA
                  cx = pwA(1,ca) + (lsA+1)*(pwB(1,cb) + (lsB+1)*(pwC(1,cc) + (lsC+1)*pwD(1,cd)))
                  cy = pwA(2,ca) + (lsA+1)*(pwB(2,cb) + (lsB+1)*(pwC(2,cc) + (lsC+1)*pwD(2,cd)))
                  cz = pwA(3,ca) + (lsA+1)*(pwB(3,cb) + (lsB+1)*(pwC(3,cc) + (lsC+1)*pwD(3,cd)))
                  val = 0.d0
                  do kq = 1 , nq
                    val = val + Idir(kq,cx,1)*Idir(kq,cy,2)*Idir(kq,cz,3)
                  end do
                  blk(ca,cb,cc,cd) = blk(ca,cb,cc,cd) + fac * coA(i,ca)*coB(j,cb)*coC(k,cc)*coD(l,cd) * val
                 end do
                end do
               end do
              end do

            end do
          end do
        end do
      end do

      deallocate(k_fac, k_cen, k_q, k_N, k_off, k_E, b_cen, b_N, b_E, Iv, F, dco, G, Idir, SNv, gt, cph, sph, wq, Hall, coA, coB, coC, coD)

end subroutine ERI_shell_block_toroidal
