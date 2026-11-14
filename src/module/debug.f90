module debug_tools
      implicit none
      private
      public :: check_symmetric, check_scf_state
      public :: check_eri_symmetry 

contains

      subroutine check_symmetric(label, A, unit, tol)
          character(*),     intent(in) :: label
          double precision, intent(in) :: A(:,:)
          integer,          intent(in) :: unit
          double precision, intent(in), optional :: tol
          double precision :: d, t
          integer :: loc(2)
          t = 1d-10 ; if (present(tol)) t = tol
          d   = maxval(abs(A - transpose(A)))
          loc = maxloc(abs(A - transpose(A)))
          write(unit,'(6x,a,es12.4,a,2i5)') '|'//label//' - '//label//'^T| = ', d, '   at', loc
          if (d > t) write(unit,'(1x,a)') '  WARNING: '//label//' is not symmetric'
          write(unit,'(a)') ""
      end subroutine

      subroutine check_scf_state(n, F, P, S, X, unit)
          integer,          intent(in) :: n, unit
          double precision, intent(in) :: F(n,n), P(n,n), S(n,n), X(n,n)
          double precision :: M(n,n)
          integer :: i
          call check_symmetric('F', F, unit)
          write(unit,'(6x,a,es12.4)') '|PSP - 2P|    = ', &
                maxval(abs(matmul(P,matmul(S,P)) - 2d0*P))
          M = matmul(transpose(X), matmul(S,X))
          do i = 1, n
            M(i,i) = M(i,i) - 1d0
          end do
          write(unit,'(6x,a,es12.4)') '|X^T S X - I| = ', maxval(abs(M))
      end subroutine

      subroutine check_eri_symmetry(n, ERI, unit)
          use unitcell_module
          integer,          intent(in) :: n, unit
          double precision, intent(in) :: ERI(n,n,n,n)
          double precision :: dmax
          integer :: p, q, r, s
          dmax = 0d0
          do s=1,nfuc; do r=1,n; do q=1,n; do p=1,n     ! column-major order
            dmax = max(dmax, abs(ERI(p,q,r,s)-ERI(q,p,r,s)), &
                             abs(ERI(p,q,r,s)-ERI(p,q,s,r)), &
                             abs(ERI(p,q,r,s)-ERI(r,s,p,q)))
          end do; end do; end do; end do
          write(unit,'(1x,a,es12.4)') 'max ERI symmetry violation = ', dmax
      end subroutine

end module debug_tools
