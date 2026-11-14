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

end subroutine get_keywords