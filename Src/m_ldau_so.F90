
! Copyright (C) 1996-2016       The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt.
! See Docs/Contributors.txt for a list of contributors.
!

!> \brief General purpose of the m_ldau_so module:
!! Compute the LSDA+U potential in the case of non-collinear magnetism
!! or spin orbit coupling.
!!
!! Starting from Eq. (5) of the paper by Liechtenstein {\it et al.} 
!! \cite Liechtenstein-95,
!!
!! \f{eqnarray*}{
!!   V_{m m^\prime}^{\sigma} = \sum_{m^{\prime \prime} m^{\prime\prime\prime}} &
!!   \left[ \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime},m^{\prime\prime\prime} \rangle
!!          n_{m^{\prime \prime}m^{\prime\prime\prime}}^{-\sigma} \right.  +
!!   \nonumber \\
!!        & \left. \left(
!!        \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime},m^{\prime\prime\prime} \rangle -
!!        \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime\prime\prime}, m^{\prime} \rangle
!!        \right) n_{m^{\prime \prime}m^{\prime\prime\prime}}^{\sigma}
!!   \right]
!!   \nonumber \\
!!   & - U \left( n - \frac{1}{2} \right) + 
!!   J \left( n^{\sigma} - \frac{1}{2} \right),
!! \f}
!! where
!! \f{eqnarray*}{
!!   n^{\sigma} & = \mathrm{Tr} \left( n_{m m^{\prime}}^{\sigma} \right)
!!   \nonumber \\
!!   & = \sum_{m m^{\prime}} n_{m m^{\prime}}^{\sigma} \delta_{m m^{\prime}} ,
!! \f}
!! and
!! \f{eqnarray*}{
!!  n = n^{\uparrow} + n^{\downarrow}
!! \f}
!! Replacing the two last Equations in the double counting terms 
!! of the first Equation,
!! \f{eqnarray*}{
!!   V_{m m^\prime}^{\sigma} = \sum_{m^{\prime \prime} m^{\prime\prime\prime}} &
!!   \left[ \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime},m^{\prime\prime\prime} \rangle
!!          n_{m^{\prime \prime}m^{\prime\prime\prime}}^{-\sigma} \right.  +
!!   \nonumber \\
!!        & \left. \left(
!!        \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime},m^{\prime\prime\prime} \rangle -
!!        \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime\prime\prime}, m^{\prime} \rangle
!!        \right) n_{m^{\prime \prime}m^{\prime\prime\prime}}^{\sigma}
!!   \right]
!!   \nonumber \\
!!   & - U \left[ \sum_{m^{\prime\prime} m^{\prime \prime \prime}} \left(
!!       n_{m^{\prime\prime} m^{\prime \prime \prime}}^{\sigma} \delta_{m^{\prime\prime} m^{\prime \prime \prime}} +
!!       n_{m^{\prime\prime} m^{\prime \prime \prime}}^{-\sigma} \delta_{m^{\prime\prime} m^{\prime \prime \prime}}\right) \right] + U \frac{1}{2}
!!   \nonumber \\
!!   & + J \left[ \sum_{m^{\prime\prime} m^{\prime \prime \prime}} \left(  n_{m^{\prime\prime} m^{\prime \prime \prime}}^{\sigma} \delta_{m^{\prime\prime} m^{\prime \prime \prime}} \right) \right] - J \frac{1}{2}
!!   \nonumber \\
!!   = \sum_{m^{\prime \prime} m^{\prime\prime\prime}} &
!!   \left[ \left( \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime},m^{\prime\prime\prime} \rangle - U \delta_{m^{\prime\prime} m^{\prime \prime \prime}}   \right)
!!          n_{m^{\prime \prime}m^{\prime\prime\prime}}^{-\sigma} \right.  +
!!   \nonumber \\
!!        & \left. \left(
!!        \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime},m^{\prime\prime\prime} \rangle -
!!        \langle m, m^{\prime \prime} \vert
!!          V_{ee} \vert m^{\prime\prime\prime}, m^{\prime} \rangle
!!        - (U - J) \delta_{m^{\prime\prime} m^{\prime \prime \prime}}
!!        \right) n_{m^{\prime \prime}m^{\prime\prime\prime}}^{\sigma}
!!   \right]
!!   \nonumber \\
!!   & + \frac{1}{2} \left( U - J \right)
!! \f}
!! If we replace
!! \f{eqnarray*}{
!!   m & \rightarrow 1,
!!   \nonumber \\
!!   m^{\prime} & \rightarrow 2,
!!   \nonumber \\
!!   m^{\prime \prime} & \rightarrow 3,
!!   \nonumber \\
!!   m^{\prime \prime \prime} & \rightarrow 4,
!! \f}
!! in the former Equation, we arrive to Eq.(3) in the paper
!! by Bousquet and Spaldin \cite Bousquet-10,
!! \f{eqnarray*}{
!!    V_{1, 2}^{\sigma} = \sum_{3,4} &
!!   \left[ \left( \langle 1, 3 \vert
!!          V_{ee} \vert 2, 4 \rangle - U \delta_{3,4}   \right)
!!          n_{3,4}^{-\sigma} \right.  +
!!   \nonumber \\
!!        & \left. \left(
!!        \langle 1, 3 \vert
!!          V_{ee} \vert 2, 4 \rangle -
!!        \langle 1, 3 \vert
!!          V_{ee} \vert 4, 2 \rangle
!!        - (U - J) \delta_{3,4}
!!        \right) n_{3,4}^{\sigma}
!!   \right]
!!   \nonumber \\
!!   & + \frac{1}{2} \left( U - J \right)
!! \f}



module m_ldau_so

  use precision,       only: dp    ! Double precision
  use m_vee_integrals, only 

  implicit none

  private


  CONTAINS 

end module m_ldau_so
