!---------------------------------------------------------------------!
! [Nov 2026] keywords checkpoint / restart (fast_ccsdt path)           !
! HF orbitals: written once after HF (binary, stream access).          !
! With restart: if the file matches (sizes, mesh) and the orbital      !
! energies agree to 1e-6, the saved orbitals replace the new ones, so  !
! that the saved CC amplitudes keep their meaning (same phases).       !
!---------------------------------------------------------------------!
subroutine ckpt_orbitals(nBas, nO, coeff_SAO, e, blk)

      use files
      use constants_module, only: dp, dpc
      use unitcell_module , only: nfuc, N_cell
      use k_blocks        , only: K_block, build_k_blocks
      use keywords        , only: c_restart, ckpt_dir, ckpt_orb_ok

      implicit none

      integer      , intent(in)    :: nBas, nO
      complex(dpc) , intent(inout) :: coeff_SAO(nBas,nBas)
      real(dp)     , intent(inout) :: e(nBas)
      type(K_block), intent(inout) :: blk(N_cell)

      integer, parameter :: magic = 20261115
      integer      :: u, ios, h(5), st
      logical      :: ex
      real(dp)     :: de
      complex(dpc), allocatable :: cs(:,:)
      real(dp)    , allocatable :: es(:)
      character(len=600) :: f

      f = trim(ckpt_dir)//'/hf_orbitals.bin'
      ckpt_orb_ok = .false.

      if (c_restart) then
        inquire(file=trim(f), exist=ex)
        if (ex) then
          open(newunit=u, file=trim(f), form='unformatted', access='stream', status='old', iostat=ios)
          read(u, iostat=ios) h
          if (ios == 0 .and. all(h == [magic, nBas, nO, N_cell, nfuc])) then
            allocate(cs(nBas,nBas), es(nBas))
            read(u, iostat=ios) es, cs
            if (ios == 0) then
              de = maxval(abs(es - e))
              if (de < 1.0d-6) then
                coeff_SAO = cs ; e = es
                call build_k_blocks(nBas, nO, N_cell, coeff_SAO, e, blk)
                ckpt_orb_ok = .true.
                write(outfile,'(a,es9.2,a)') '  [restart] HF orbitals read from '//trim(f)//' (max |de| =', de, ')'
              else
                write(outfile,'(a,es9.2,a)') '  [restart] WARNING: orbital energies differ (', de, &
                                             '): saved orbitals and amplitudes NOT used'
              end if
            end if
            deallocate(cs, es)
          else
            write(outfile,'(a)') '  [restart] WARNING: '//trim(f)//' does not match this system: not used'
          end if
          close(u)
        else
          write(outfile,'(a)') '  [restart] no file '//trim(f)//': normal start'
        end if
      end if

      if (.not. ckpt_orb_ok) then
        open(newunit=u, file=trim(f)//'.tmp', form='unformatted', access='stream', status='replace')
        write(u) [magic, nBas, nO, N_cell, nfuc]
        write(u) e, coeff_SAO
        close(u)
        call rename(trim(f)//'.tmp', trim(f), st)
        write(outfile,'(a)') '  [checkpoint] HF orbitals written to '//trim(f)
      end if
      flush(outfile)

end subroutine ckpt_orbitals
