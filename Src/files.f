! ---
! Copyright (C) 1996-2014	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
! 
!
      module files
!
!     Contains the short system label, used to generate file names
!     slabel is currently set in reinit.
!
      integer, parameter, public                  :: label_length = 60
      character(len=label_length), save, public   :: slabel

      private

      end module files
