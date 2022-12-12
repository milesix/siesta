module m_steps
  implicit none
  public
  integer:: inicoor           ! Initial step in geommetry iteration
  integer:: fincoor           ! Final step in geommetry iteration
  integer:: inicoor_bath        !Initial and final steps in the geommetry iteration
  integer:: inicoor_Verlet      !for the Nose and Verlet Simulations (idyn=9)
  integer:: fincoor_bath
  integer:: fincoor_Verlet
  character(len=150) :: sys_prep 
  logical :: TwoBaths
  integer:: istp              ! Geommetry iteration step starting in istp=1
  logical:: final=.false.     ! Last geometry step?
  integer :: natoms_1

end module m_steps
