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
