module keywords
  
      implicit none

      public  :: c_One , c_OV , c_K , c_NA , c_ERI
      public  :: c_read, c_Integral, c_trexio
      public  :: c_Angstrom, c_plot, c_details
      public  :: c_MO      , c_UHF , c_Huckel
      public  :: c_MP2
      public  :: c_SAO
      public  :: c_save
      public  :: c_DIIS
      public  :: c_MCDE
      public  :: c_nquad_128
      public  :: c_fast_quad
      public  :: c_prim_screen
      public  :: c_cholesky, cholesky_tol
      public  :: c_cc_disk
      public  :: c_direct, c_df
      public  :: c_symmetry                         ! [Nov 2026] keyword symmetry
      public  :: c_df_etb, df_cd_tol, s_eig_min   ! [Nov 2026] density fitting: aux set, Cholesky tolerance
      public  :: c_df_auxfile, c_df_chol          ! [Nov 2026] density fitting: aux set from file / Cholesky set
      public  :: c_fs_corr                        ! [Nov 2026] finite-size correction (electrostatics)
      public  :: c_physicist
      public  :: c_ccsdt
      public  :: c_fast_ccsdt
      public  :: c_fast_real      
      public  :: c_frozen_core    
      public  :: c_checkpoint, c_restart, ckpt_dir, ckpt_orb_ok
      public  :: c_debug


      public  :: c_Spin_adapt
      public  :: c_Lanczos
      public  :: c_screen
      public  :: c_block_lanczos              ! [Nov 2026] MCDE: block Lanczos spectral function


      logical :: c_One , c_OV , c_K , c_NA , c_ERI
      logical :: c_read, c_Integral, c_trexio
      logical :: c_Angstrom, c_plot, c_details
      logical :: c_MO      , c_UHF , c_Huckel
      logical :: c_MP2     , c_Orbitals
      logical :: c_SAO
      logical :: c_save
      logical :: c_DIIS
      logical :: c_MCDE
      logical :: c_nquad_128
      logical :: c_fast_quad   = .true.
      logical :: c_prim_screen = .true.
      logical :: c_cholesky    = .false.
      
      logical :: c_cc_disk  = .false.                                    ! [Nov 2026] keyword CC_to_disk
      logical :: c_direct   = .false.                                    ! [Nov 2026] keyword direct
      logical :: c_df       = .false.                                    ! [Nov 2026] keyword density_fitting
      logical :: c_symmetry = .false.                                    ! [Nov 2026] keyword symmetry: point symmetry used in V, ERIs, direct J/K, DF
      logical :: c_df_etb   = .false.                                    ! [Nov 2026] keyword df_etb: even-tempered aux basis (PySCF rule)
      logical :: c_df_auxfile = .false.                                  ! [Nov 2026] keyword df_auxfile: aux basis from the file AuxBasis
      logical :: c_fs_corr  = .false.                                    ! [Nov 2026] keyword fs_correction
      logical :: c_df_chol  = .false.                                    ! [Nov 2026] keywords df_threshold_6/8/10: Cholesky fitting set
      logical :: c_physicist
      logical :: c_ccsdt
      logical :: c_fast_ccsdt 
      logical :: c_fast_real  
      logical :: c_frozen_core
      logical :: c_checkpoint = .false., c_restart = .false.
      logical :: ckpt_orb_ok  = .false.
      
      logical :: c_debug



      logical :: c_Spin_adapt
      logical :: c_Lanczos
      logical :: c_screen
      logical :: c_block_lanczos = .false.                               ! [Nov 2026] keyword block_lanczos

      double precision   :: cholesky_tol = 1.0d-6
      double precision   :: df_cd_tol = 1.0d-8                           ! [Nov 2026] Cholesky fitting tolerance (df_threshold_6/8/10)
      double precision   :: s_eig_min = 1.0d0                            ! [Nov 2026] smallest overlap eigenvalue (set in check_the_overlap)
      character(len=512) :: ckpt_dir = './checkpoint'


end module keywords