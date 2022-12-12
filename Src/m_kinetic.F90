module m_kinetic

  use precision, only: dp
	implicit none

  public

  real(dp):: vn      = 0.0_dp   ! Velocity (time derivative) of the Nose thermostat
  real(dp):: vn_n(1,2) = 0.0_dp !Velocities of the Nose thermostats
  real(dp):: vpr     = 0.0_dp   ! Velocity (time derivative) of the PR variables

  real(dp):: tempion = 0.0_dp   ! Ionic temperature
  real(dp):: tempion_n(1,2)= 0.0_dp ! Ionic temperature Nose baths

  real(dp):: kn       = 0.0_dp  ! Kinetic energy of the Nose' thermostat
  real(dp):: kn_n(1,2)  = 0.0_dp        ! Kinetic energies of the Nose'
!thermostats
  real(dp):: kpr      = 0.0_dp  ! Kinetic energy of the Parrinello-Rahman variables

end module m_kinetic



