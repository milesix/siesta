       MODULE linres_matrices
C      --------------------------------------------------------------
C      Linres module to store the matrices which are non self-consistent
C      in the Linres part. As these matrices can be computed
C      only once but are used in every scf call they are declared in
C      this module and resetFirstmatrices() is called once the perturbed
C      density matrix is converged. L.Riches and S.Illera (2016).
C      --------------------------------------------------------------
C real(dp) dRhoatm(3,ntpl)        : Derivative of the atomic density in the mesh points.
C                                   rho^atm* der(phi_mu)*phi_nu 
C                                      + rho^atm * phi_mu * der(phi_nu)
C real(dp) dRhoscf0(3,ntpl,nspin) : Non-scf part of the perturbed density corresponding
C                                   to:  rho^GD_mu,nu* der(phi_mu)*phi_nu 
C                                      + rho^GD_mu,nu * phi_mu * der(phi_nu)
C                                   where rho^GD corresponds to the
C				    ground state density matrix
C real (dp) drhopcc(1:ntlp,1:3)   : gradient/derivative of the core correction density

C NOTE: the first time that dRhoatm and dRhoscf0 and drhopcc are used, they store the gradient 
C instead of the radial derivative.

       use precision, only: grid_p
       real(grid_p), pointer :: dRhoatm(:,:) => null()
       real(grid_p), pointer :: dRhoscf0(:,:,:) => null()
       real(grid_p), pointer :: dRhopcc(:,:) => null()

       CONTAINS

       subroutine resetFirstmatrices()
       use alloc, only: de_alloc
       implicit none
       call de_alloc(dRhoatm, 'dRhoatm', 'linres_matrices')
       call de_alloc(dRhoscf0, 'dRhoscf0', 'linres_matrices')
       call de_alloc(drhopcc, 'drhopcc', 'linres_matrices')
       end subroutine resetFirstmatrices

       END MODULE linres_matrices

