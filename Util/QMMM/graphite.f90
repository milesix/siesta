
subroutine find_graphite_layers(nac,attype,ng1,bondxat,graphite_layer_no)

  implicit none

  integer, intent(in) :: nac, ng1(nac,6), bondxat(nac)
  character(len=4), dimension(nac), intent(in) :: attype(nac)
  integer, dimension(nac),  intent(inout) :: graphite_layer_no

  integer :: i, num_of_graphite_layers

  num_of_graphite_layers=0
  do i=1,nac
     if (graphite_layer_no(i)/=0 .or. attype(i)/='C#') cycle
     num_of_graphite_layers=num_of_graphite_layers+1
     graphite_layer_no(i)=num_of_graphite_layers
     call check_graphite_connectivity(i,nac,attype,ng1,bondxat,graphite_layer_no)
  enddo
  
end subroutine find_graphite_layers

recursive subroutine check_graphite_connectivity(i,nac,attype,ng1,bondxat,graphite_layer_no)

  implicit none

  integer, intent(in) :: i, nac, ng1(nac,6), bondxat(nac)
  character(len=4), dimension(nac), intent(in) :: attype(nac)
  integer, dimension(nac),  intent(inout) :: graphite_layer_no

  integer :: j

  do j=1,bondxat(i)
     if (graphite_layer_no(ng1(i,j))/=0 .or. attype(i)/='C#') cycle
     graphite_layer_no(ng1(i,j))=graphite_layer_no(i)
     call check_graphite_connectivity(ng1(i,j),nac,attype,ng1,bondxat,graphite_layer_no)
  enddo

end subroutine check_graphite_connectivity
