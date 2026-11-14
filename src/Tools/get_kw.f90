subroutine  get_keywords(keyword)

      use keywords

      implicit none 

      character*(10) ,intent(in)  ::   keyword(100)

      c_integral    = any(keyword == 'Integrals')
      c_read        = any(keyword == 'Read'     )
      c_trexio      = any(keyword == 'Trexio'   )
      c_Angstrom    = any(keyword == 'Angstrom' )
      c_plot        = any(keyword == 'Plot'     )     ! 5
      c_details     = any(keyword == 'Details'  )
      c_MO          = any(keyword == 'MO'       )
      c_UHF         = any(keyword == 'UHF'      )
      c_Huckel      = any(keyword == 'Huckel'   )
      c_MP2         = any(keyword == 'MP2'      )     ! 10
      c_OV          = any(keyword == 'OV'       )
      c_K           = any(keyword == 'K'        )
      c_NA          = any(keyword == 'NA'       )
      c_ERI         = any(keyword == 'ERI'      )
      c_One         = any(keyword == 'One'      )     ! 15
      c_Orbitals    = any(keyword == 'Orbitals' )
      c_SAO         = any(keyword == 'SAO'      )
      c_save        = any(keyword == 'Save'     )
      c_DIIS        = any(keyword == 'DIIS'     )     
      c_MCDE        = any(keyword == 'MCDE'     )     ! 20
      c_nquad_128   = any(keyword == 'nquad_128')
      c_physicist   = any(keyword == 'physicist')
      c_ccsdt       = any(keyword == 'ccsdt'    )
      c_debug       = any(keyword == 'debug'    )     ! 25 

      ! [Nov 2026] ONE rule for the order in which the AO two-electron
      ! integrals are stored (used by libcint, Tori integrals, J/K, AO->MO):
      !   c_physicist = .true.  ->  ERI(i,k,j,l) = (ij|kl)  (physicist)
      ! MCDE reads physicist-ordered integrals, so MCDE switches it on.
      ! (before: libcint used MCDE, the Tori code used physicist and J/K
      !  used "physicist or MCDE" -> wrong Fock matrix for some combinations)
      if (c_MCDE) c_physicist = .true.

      ! [Nov 2026] the SAO path (k-space transformation, J/K) works with
      ! chemist-ordered integrals only; physicist order is only needed by
      ! MCDE, which does not run in the SAO path
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