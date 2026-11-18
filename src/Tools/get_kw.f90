subroutine  get_keywords(keyword)

      use keywords

      implicit none 

      character*(10) ,intent(in)  ::   keyword(100)

      c_integral    = any(keyword == 'Integrals' )
      c_read        = any(keyword == 'Read'      )
      c_trexio      = any(keyword == 'Trexio'    )
      c_Angstrom    = any(keyword == 'Angstrom'  )
      c_plot        = any(keyword == 'Plot'      )     ! 5
      c_details     = any(keyword == 'Details'   )
      c_MO          = any(keyword == 'MO'        )
      c_UHF         = any(keyword == 'UHF'       )
      c_Huckel      = any(keyword == 'Huckel'    )
      c_MP2         = any(keyword == 'MP2'       )     ! 10
      c_OV          = any(keyword == 'OV'        )
      c_K           = any(keyword == 'K'         )
      c_NA          = any(keyword == 'NA'        )
      c_ERI         = any(keyword == 'ERI'       )
      c_One         = any(keyword == 'One'       )     ! 15
      c_Orbitals    = any(keyword == 'Orbitals'  )
      c_SAO         = any(keyword == 'SAO'       )
      c_save        = any(keyword == 'Save'      )
      c_DIIS        = any(keyword == 'DIIS'      )     
      c_MCDE        = any(keyword == 'MCDE'      )     ! 20
      c_nquad_128   = any(keyword == 'nquad_128' )
      c_fast_quad   = any(keyword == 'fast_quad' )
      c_prim_screen = any(keyword == 'primscreen')     ! input word: prim_screen
      c_physicist   = any(keyword == 'physicist' )
      c_ccsdt       = any(keyword == 'ccsdt'     )
      c_debug       = any(keyword == 'debug'     )     ! 25 
      c_fast_ccsdt  = any(keyword == 'fast_ccsdt')    
      c_fast_real   = any(keyword == 'fcc_real'  )     
      c_frozen_core = any(keyword == 'frozcore'  )     
      c_checkpoint  = any(keyword == 'checkpoint')     
      c_restart     = any(keyword == 'restart'   )     ! 30
      if (c_restart) c_checkpoint = .true.
      if (c_checkpoint) call set_ckpt_dir()

      ! [Nov 2026] keyword cholesky: Cholesky-decomposed integrals in fast k-CCSD(T)
      c_cholesky    = any(keyword == 'cholesky' ) .or. any(keyword == 'cholesky_4') .or. any(keyword == 'cholesky_8')
      cholesky_tol  = 1.0d-6
      if (any(keyword == 'cholesky_4')) cholesky_tol = 1.0d-4
      if (any(keyword == 'cholesky_8')) cholesky_tol = 1.0d-8
      if (c_cholesky .and. .not. c_fast_real) c_fast_ccsdt = .true.
      ! [Nov 2026] keyword CC_to_disk: vvvv k-blocks on disk, read by block in the ladder
      c_cc_disk     = any(keyword == 'cc_disk'  )
      if (c_cc_disk .and. c_cholesky) then
        write(*,'(a)') ' WARNING: CC_to_disk and cholesky together: cholesky is used (no vvvv integrals to store)'
        c_cc_disk = .false.
      end if
      if (c_cc_disk .and. .not. c_fast_real) c_fast_ccsdt = .true.
      ! [Nov 2026] keyword direct: integral-direct RHF (SAO) and k-MP2, no stored ERIs
      c_direct      = any(keyword == 'direct'   )
      if (c_direct) then
        if (.not. c_SAO) error stop 'keyword direct needs SAO (toroidal RHF)'
        if (c_ccsdt .or. c_cholesky .or. c_cc_disk .or. c_MCDE .or. c_UHF) &
          error stop 'keyword direct: only RHF (SAO) and MP2 (no ccsdt, fast_ccsdt, cholesky, CC_to_disk, MCDE, UHF)'
      end if
      ! [Nov 2026] keyword density_fitting: fitted J/K and (ov|ov), no 4-index ERIs stored
      !            (uses the direct-mode plumbing: no ERI array, AO shells and Schwarz bounds)
      c_df          = any(keyword == 'df'       )
      c_df_etb      = any(keyword == 'df_etb'   )      ! [Nov 2026]
      df_cd_tol     = 1.0d-8
      if (any(keyword == 'df_cd_6' )) df_cd_tol = 1.0d-6
      if (any(keyword == 'df_cd_10')) df_cd_tol = 1.0d-10
      if (c_df) then
        if (.not. c_SAO) error stop 'keyword density_fitting needs SAO (toroidal RHF)'
        if (c_cholesky .or. c_cc_disk .or. c_MCDE .or. c_UHF) &
          error stop 'keyword density_fitting: not with cholesky, CC_to_disk, MCDE or UHF (DF replaces cholesky in the CC)'
        ! [Nov 2026] with ccsdt / fast_ccsdt / fast_ccsdt_real: DF-CCSD(T), the DF vectors
        ! take the place of the Cholesky vectors (same CC path as keyword cholesky)
        if (c_ccsdt .or. c_fast_ccsdt .or. c_fast_real) then
          c_cholesky = .true.
          if (.not. c_fast_real) c_fast_ccsdt = .true.
        end if
        c_direct = .true.
      end if

      if (c_fast_ccsdt) c_ccsdt = .true.
      if (c_fast_real) then
        c_fast_ccsdt = .true. ; c_ccsdt = .true.
      end if

      if (c_checkpoint .and. .not. c_fast_ccsdt) then
        write(*,'(a)') ' WARNING: keywords checkpoint/restart need fast_ccsdt or fast_ccsdt_real: ignored'
        c_checkpoint = .false. ; c_restart = .false.
      end if

      if (c_frozen_core .and. .not. c_fast_ccsdt) then
        write(*,'(a)') ' WARNING: keyword frozen_core needs fast_ccsdt or fast_ccsdt_real: ignored'
        c_frozen_core = .false.
      end if

      
      if (c_MCDE) c_physicist = .true.

      if (c_SAO .and. c_MCDE) then
        write(*,'(a)') ' ERROR: MCDE is not available with SAO (remove one of the two keywords)'
        stop
      end if
      if (c_SAO .and. c_UHF) then
        write(*,'(a)') ' ERROR: UHF is not available with SAO (remove one of the two keywords)'
        stop
      end if
      if (c_SAO .and. c_physicist) then
        write(*,'(a)') ' WARNING: keyword physicist ignored with SAO (the SAO path uses chemist order)'
        c_physicist = .false.
      end if

end subroutine get_keywords

! [Nov 2026] directory of the restart files:
!   $CIPOP_CHECKPOINT_DIR  >  $SLURM_SUBMIT_DIR/checkpoint  >  ./checkpoint
subroutine set_ckpt_dir()
      use keywords, only: ckpt_dir
      implicit none
      character(len=512) :: v
      integer :: l, st
      call get_environment_variable('CIPOP_CHECKPOINT_DIR', v, l, st)
      if (st == 0 .and. l > 0) then
        ckpt_dir = v
      else
        call get_environment_variable('SLURM_SUBMIT_DIR', v, l, st)
        if (st == 0 .and. l > 0) then
          ckpt_dir = trim(v)//'/checkpoint'
        else
          ckpt_dir = './checkpoint'
        end if
      end if
      call execute_command_line('mkdir -p "'//trim(ckpt_dir)//'"', wait=.true.)
end subroutine set_ckpt_dir
