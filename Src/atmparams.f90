!     This file is part of the SIESTA package.
!     
!     Copyright (c) Fundacion General Universidad Autonoma de Madrid:
!     E.Artacho, J.Gale, A.Garcia, J.Junquera, P.Ordejon, D.Sanchez-Portal
!     and J.M.Soler, 1996-2006.
!     
!     Use of this software constitutes agreement with the full conditions
!     given in the SIESTA license, as signed by all legitimate users.
!     
      module atmparams

      implicit none 
!
!    Hard-wired parameters to be eliminated in the future
!

! INTEGER  NZETMX   : Maximum number of PAOs or polarization orbitals
!                     with the same angular  momentum and 
!                     for the same species.       

         integer, parameter  :: nzetmx =    3  

! INTEGER  NKBMX    : Maximum number of Kleinman-Bylander projectors
!                     for each angular momentum

         integer, parameter  :: nkbmx  =    2  

! INTEGER  NSMX    : Maximum number of semicore shells for each angular
!                    momentum present in the atom ( for normal atom nsmx=0)

         integer, parameter  :: nsmx  =    1  
         integer, parameter  :: nsemx = 1 + nsmx  



! INTEGER LMAXD     : Maximum angular momentum for both orbitals and 
!                      projectors.

         integer, parameter  :: lmaxd  =    4  
         integer, parameter  :: lmx2   = (lmaxd+1)*(lmaxd+1)  

! INTEGER  NRMAX    : Maximum number of points in the functions read
!                     from file '.vps' or '.psatom.data' (this number is
!                     determined by the parameter nrmax in the
!                     program atm, which generates the files with
!                     the pseudopotential information).

         integer, parameter  :: nrmax  = 6000  

         integer, parameter      :: maxos=2*nzetmx*lmx2*nsemx

      end module atmparams
