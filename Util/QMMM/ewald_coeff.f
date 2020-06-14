
      subroutine ewald_coeff(rcorte,ewald_alpha,kewald_cutoff)

      use precision, only: dp

      implicit none

      real(dp) rcorte
      real(dp) ewald_alpha, kewald_cutoff

      real(dp) pi

      real(dp) sfactor
  
      sfactor=2.6d0

      ewald_alpha=(sfactor/rcorte)**2

      pi=ACOS(-1.0d0)
      kewald_cutoff=2.0d0*sfactor*sqrt(ewald_alpha)

      end
