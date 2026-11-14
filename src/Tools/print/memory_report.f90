!=====================================================================!
!                                                                     !
!              MEMORY REPORT (printed before the calculation)          !
!                                                                     !
!   Estimates the memory of the large arrays of every step that the  !
!   keywords switch on, the peak, and compares it with the physical  !
!   memory of the machine.                                            !
!                                                                     !
!   N = nBas , f = nfuc (functions per unit cell) , C = N_cell        !
!   d = 8 bytes (real) , z = 16 bytes (complex)                       !
!                                                                     !
!   normal path                                                       !
!     integrals   ERI(N,N,N,N) (Tori: rows > nfuc by translation,     !
!                 in place)                          d N^4            !
!     RHF         ERI + P,J,K,F,Fp,error + DIIS      d N^4 + ...      !
!     UHF         ERI + ~14 matrices                 d N^4 + 14 d N^2 !
!     AO -> MO    IN PLACE (MP2, MCDE)               d N^4            !
!                 + per thread 2 pb N^2 (pb = 32)                     !
!     MP2/MCDE    ERI_MO (= the same array)          d N^4            !
!                                                                     !
!   SAO path (Tori1D/2D/3D)                                           !
!     integrals   ERI(f,N,N,N)                       d f N^3          !
!     RHF_SAO     ERI + 6 real + 5 complex N x N     d f N^3 + ...    !
!     AO -> k-MO  ERI + ERI_kMO + X1 + nb X2         d f N^3 + out    !
!                 (one pass, k4 by k4, no ERI_SAO)   + 16(1+nb)fN^3/C !
!                 + per thread f^4 slabs                              !
!                 out: CCSD  all (pq|rs)             z f^4 C^3        !
!                      MP2   only (ov|ov)            z o^2 v^2 C^3    !
!     k-MP2       out                                                 !
!     k-CCSD      g + t2,r2,Wvoov,Wvovo (o o v v)    z f^4 C^3        !
!                 + Woooo (o^4); Wvvvv kept in g     + 4 z o^2v^2 C^3 !
!                 + DIIS vectors (n = n1 + n2)       ~ 248 n + 28 n2  !
!                   n1 = C o v , n2 = C^3 o^2 v^2 , o = nO/C, v = f-o  !
!     (T)         g + t2 + triples tables                             !
!                                                                     !
!   (f N^3 = f^4 C^3 = number of entries of g = ERI_kMO)              !
!=====================================================================!

subroutine memory_report(calculation_type, nBas, nfuc_in, N_cell_in, nO)

      use files
      use keywords
      use HF_parameter , only : max_diis
      use omp_lib      , only : omp_get_max_threads

      implicit none

      ! - input - !

      character(len=*) , intent(in) :: calculation_type
      integer          , intent(in) :: nBas , nfuc_in , N_cell_in , nO

      ! - local - !

      integer         , parameter   :: max_step = 16
      integer         , parameter   :: ndiis_cc = 6           ! ndiis in ccsd_call.f90
      double precision, parameter   :: d = 8.d0 , z = 16.d0

      character(len=40)             :: step_name(max_step)
      character(len=60)             :: step_info(max_step)
      double precision              :: step_mem (max_step)
      integer                       :: ns , i , ipeak
      logical                       :: tori
      double precision              :: N , f , C , N2 , N4 , B3
      double precision              :: o , v , n1 , n2cc , ncc , thr , nb , out , slabs
      double precision              :: mat , peak , phys
      character(len=40)             :: lab

      !-----------------------------------------------------------------!

      tori = any(calculation_type == ["Tori1D","Tori2D","Tori3D"])

      N  = dble(nBas)
      f  = dble(nfuc_in)
      C  = dble(max(N_cell_in,1))
      N2 = N*N
      N4 = N2*N2
      B3 = f*N*N*N                                       ! = f^4 C^3

      ! matrices always allocated in main: S,T,V,Hc,X,coeff (real), coeff_SAO (complex)
      mat = 6.d0*d*N2 + z*N2

      ns = 0

      if (c_SAO) then

        ! ------------------------- SAO path ------------------------- !

        call add("Integrals (SAO)",  "ERI(nfuc,nBas,nBas,nBas)", d*B3 + d*N2)

        if (c_UHF) then
          call add("UHF (SAO path)", "WARNING: UHF expects ERI(nBas^4)", d*B3 + 14.d0*d*N2)
        else
          call add("RHF (SAO)", "ERI + 6 real + 5 complex nBas^2 + DIIS", &
                   d*B3 + 6.d0*d*N2 + 5.d0*z*N2 + 2.d0*d*N2*max_diis)
        end if

        if (c_MP2) then
          o    = dble(ceiling(dble(nO)/C))               ! occupied bands per k (estimate)
          v    = max(f - o, 0.d0)
          thr  = dble(omp_get_max_threads())
          nb   = min(C, max(1.d0, dble(ceiling(4.d0*thr/(C*C)))))   ! k4 batch (AO_to_MO_SAO.f90)
          if (c_ccsdt) then
            out   = z*B3                                                  ! all (pq|rs)
            slabs = thr*z*(4.d0*f**4)
          else
            out   = z*o*o*v*v*C**3                                        ! (ov|ov) only
            slabs = thr*z*(f**4 + f**3*v + f*f*o*v + f*o*v*v)
          end if
          call add("AO -> k-MO transformation", "ERI + ERI_kMO + X1 + X2 (one pass)", &
                   d*B3 + out + z*(1.d0+nb)*B3/C + slabs)
          call add("k-MP2", "ERI_kMO", out)

          if (c_ccsdt) then
            n1   = C*o*v
            n2cc = C**3 * o*o * v*v
            ncc  = n1 + n2cc
            call add("k-CCSD", "g + t2,r2,W (occ/vir blocks) + DIIS", &
                     z*B3 + 4.d0*z*n2cc + z*o**4*C**3 + (3.d0*z + d + 2.d0*z*ndiis_cc)*ncc &
                     + 12.d0*n1 + 28.d0*n2cc)
            call add("(T)", "g + t2 + triples tables", &
                     z*B3 + z*n2cc + 20.d0*dble(nO)**3*C + 18.d0*z*dble(nO)**3)
          end if
        end if

      else

        ! ------------------------ normal path ----------------------- !

        if (tori) then
          call add("Integrals", "ERI(nBas^4)", d*N4 + d*N2)
        else
          call add("Integrals (libcint)", "ERI(nBas^4)", d*N4)
        end if

        if (c_UHF) then
          call add("UHF", "ERI + ~14 nBas^2 matrices", d*N4 + 14.d0*d*N2)
        else
          call add("RHF", "ERI + 6 nBas^2 + DIIS", d*N4 + 6.d0*d*N2 + 2.d0*d*N2*max_diis)
        end if

        if (c_MP2 .or. c_MCDE) then
          call add("AO -> MO transformation", "in place: ERI(nBas^4) + thread scratch", &
                   d*N4 + dble(omp_get_max_threads())*d*(2.d0*32.d0+3.d0)*N2)
        end if

        if (c_MP2)  call add("MP2",  "ERI_MO(nBas^4)", d*N4)
        if (c_MCDE) call add("MCDE", "ERI_MO(nBas^4) (+ MCDE work arrays)", d*N4)

      end if

      ! ------------------------- peak / machine ----------------------- !

      peak = 0.d0 ; ipeak = 1
      do i = 1 , ns
        step_mem(i) = step_mem(i) + mat
        if (step_mem(i) > peak) then
          peak = step_mem(i) ; ipeak = i
        end if
      end do

      call physical_memory(phys)

      ! ------------------------------ print --------------------------- !

      call header("Memory estimate",-1)

      write(outfile,'(6x,a,i10,a,i10,a,i10,a,i8)') "nBas = ", nBas, "    nfuc = ", nfuc_in, &
                                                  "    N_cell = ", max(N_cell_in,1), "    nO = ", nO
      if (c_SAO) then
        write(outfile,'(6x,a)') "path : SAO  (ERI stored as ERI(nfuc,nBas,nBas,nBas))"
      else
        write(outfile,'(6x,a)') "path : normal  (ERI stored as ERI(nBas,nBas,nBas,nBas))"
      end if
      write(outfile,*)
      write(outfile,'(6x,a)') "-----------------------------------------------------------------------------------------"
      write(outfile,'(6x,a40,a40,a12)') "step                                    ", &
                                        "largest arrays                          ", "   memory"
      write(outfile,'(6x,a)') "-----------------------------------------------------------------------------------------"
      do i = 1 , ns
        write(outfile,'(6x,a40,a40,a12)') step_name(i), step_info(i)(1:40), human(step_mem(i))
      end do
      write(outfile,'(6x,a)') "-----------------------------------------------------------------------------------------"
      lab = "PEAK  ("//trim(step_name(ipeak))//")"
      write(outfile,'(6x,a40,a40,a12)') lab, " ", human(peak)
      if (phys > 0.d0) then
        lab = "physical memory of this machine"
        write(outfile,'(6x,a40,a40,a12)') lab, " ", human(phys)
      end if
      write(outfile,'(6x,a)') "-----------------------------------------------------------------------------------------"
      write(outfile,'(6x,a)') "(memory allocated at each step, incl. the nBas x nBas matrices of main; estimate)"

      if (phys > 0.d0 .and. peak > phys) then
        write(outfile,*)
        write(outfile,'(6x,a)') "WARNING: the estimated peak memory is LARGER than the physical memory."
        write(outfile,'(6x,a)') "         The calculation will swap heavily or be killed."
        write(*,'(a)')          " WARNING: estimated peak memory ("//trim(adjustl(human(peak)))// &
                                ") > physical memory ("//trim(adjustl(human(phys)))//")"
      else if (phys > 0.d0 .and. peak > 0.8d0*phys) then
        write(outfile,*)
        write(outfile,'(6x,a)') "WARNING: the estimated peak memory is above 80% of the physical memory."
        write(*,'(a)')          " WARNING: estimated peak memory ("//trim(adjustl(human(peak)))// &
                                ") above 80% of physical memory ("//trim(adjustl(human(phys)))//")"
      end if
      write(outfile,*)
      flush(outfile)

      write(*,'(a,a,a,a,a)') " Memory estimate: peak ", trim(adjustl(human(peak))), &
                             "  (", trim(step_name(ipeak)), ")"

contains

      subroutine add(name, info, mem)
        character(len=*), intent(in) :: name, info
        double precision, intent(in) :: mem
        ns = ns + 1
        step_name(ns) = name
        step_info(ns) = info
        step_mem (ns) = mem
      end subroutine add

      ! bytes -> "12.34 GB"
      character(len=12) function human(bytes)
        double precision, intent(in) :: bytes
        character(len=2), parameter  :: unit(5) = ["B ","KB","MB","GB","TB"]
        double precision             :: x
        integer                      :: k
        x = bytes ; k = 1
        do while (x >= 1024.d0 .and. k < 5)
          x = x / 1024.d0 ; k = k + 1
        end do
        write(human,'(f9.2,1x,a2)') x, unit(k)
      end function human

      ! physical memory in bytes (macOS: sysctl, Linux: /proc/meminfo), 0 if unknown
      subroutine physical_memory(bytes)
        double precision, intent(out) :: bytes
        character(len=*), parameter   :: fname = ".memory_report_tmp"
        integer                       :: u, ios
        integer(8)                    :: b
        bytes = 0.d0
        call execute_command_line("(sysctl -n hw.memsize 2>/dev/null || " // &
             "awk '/MemTotal/{printf ""%.0f\n"", $2*1024}' /proc/meminfo) > " // fname // " 2>/dev/null", &
             exitstat=ios)
        open(newunit=u, file=fname, status='old', action='read', iostat=ios)
        if (ios /= 0) return
        read(u,*,iostat=ios) b
        if (ios == 0 .and. b > 0) bytes = dble(b)
        close(u, status='delete')
      end subroutine physical_memory

end subroutine memory_report
