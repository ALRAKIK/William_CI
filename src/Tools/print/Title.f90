subroutine Title()

      use files 

      implicit none 

      character(len=8)  :: date
      character(len=10) :: time

      ! --------------------------------------------------------------- !

      write(outfile,*) 
      write(outfile,*)                                        
                                   
      write(outfile,"(6x,a)") "__| |________________________________________________________________________________________| |__"
      write(outfile,"(6x,a)") "__   ________________________________________________________________________________________   __"                                                    
      write(outfile,"(6x,a)") "  | |      ___                       ___       ___                   ___            ___      | |  "
      write(outfile,"(6x,a)") "  | |     /\__\          ___        /\__\     /\__\      ___        /\  \          /\__\     | |  "
      write(outfile,"(6x,a)") "  | |    /:/ _/_        /\  \      /:/  /    /:/  /     /\  \      /::\  \        /::|  |    | |  "
      write(outfile,"(6x,a)") "  | |   /:/ /\__\       \:\  \    /:/  /    /:/  /      \:\  \    /:/\:\  \      /:|:|  |    | |  "
      write(outfile,"(6x,a)") "  | |  /:/ /:/ _/_      /::\__\  /:/  /    /:/  /       /::\__\  /::\~\:\  \    /:/|:|__|__  | |  "
      write(outfile,"(6x,a)") "  | | /:/_/:/ /\__\  __/:/\/__/ /:/__/    /:/__/     __/:/\/__/ /:/\:\ \:\__\  /:/ |::::\__\ | |  "
      write(outfile,"(6x,a)") "  | | \:\/:/ /:/  / /\/:/  /    \:\  \    \:\  \    /\/:/  /    \/__\:\/:/  /  \/__/~~/:/  / | |  "
      write(outfile,"(6x,a)") "  | |  \::/_/:/  /  \::/__/      \:\  \    \:\  \   \::/__/          \::/  /         /:/  /  | |  "
      write(outfile,"(6x,a)") "  | |   \:\/:/  /    \:\__\       \:\  \    \:\  \   \:\__\          /:/  /         /:/  /   | |  "
      write(outfile,"(6x,a)") "  | |    \::/  /      \/__/        \:\__\    \:\__\   \/__/         /:/  /         /:/  /    | |  "
      write(outfile,"(6x,a)") "  | |     \/__/                     \/__/     \/__/                 \/__/          \/__/     | |  "
      write(outfile,"(6x,a)") "__| |_______________________________________________________________________________________ | |__"
      write(outfile,"(6x,a)") "__   _______________________________________________________________________________________    __"
      write(outfile,"(6x,a)") "  | |                                                                                        | |  "                                                  

      write(outfile,'(a)')
      write(outfile,'(a)')
      write(outfile,'(a)')  "       William_CI: Calculate Integrals, is a program to calculate   the" 
      write(outfile,'(a)')  "                   molecular integrals using Clifford Toroidal Gaussian"
      write(outfile,'(a)')  "                   functions for advanced quantum chemistry computations"
      write(outfile,'(a)')  ""
      write(outfile,'(2a)') "      ",repeat("=",67)
      write(outfile,'(a)')  ""
      write(outfile,'(2a)') "      ", "Amer Alrakik        Email  : alrakikamer@gmail.com"
      write(outfile,'(2a)') "      ", "                    Github : https://github.com/ALRAKIK/William_CI"
      write(outfile,'(a)')  ""
      write(outfile,'(2a)') "      ", "Contributors        Arjan Berger"
      write(outfile,'(2a)') "      ", "                    Stefano Evangelisti"
      write(outfile,'(a)')  ""
      write(outfile,'(2a)') "      ",repeat("=",67)
      write(outfile,'(a)')  "                          Version: 0.4.1 | Quantum Chemistry | 2026"
      write(outfile,'(2a)') "      ",repeat("=",67)
      write(outfile,'(a)')  ""


      call date_and_time(date, time)

      write(outfile,'(6x,10a)') "Executing Date : ", date(7:8), "/", date(5:6), "/", date(1:4)
      write(outfile,'(6x,10a)') "Executing time : ", time(1:2), ":", time(3:4), ":", time(5:6)
      write(outfile,'(a)')  ""

      FLUSH(outfile)

      ! --------------------------------------------------------------- !

end subroutine Title