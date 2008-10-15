module compare_m
use flib_dom
use compare_tol_m, only: findTol, sp, dp
use m_strings
use corresponding_node
use string_utilities

implicit none



public compare
private

logical, public :: compare_debug = .false.
logical :: compare_node_debug = .false.
logical :: compare_values_debug = .false.
logical :: matching_debug = .false.
logical :: compare_node_list_debug = .false.
logical, public :: STOP_ON_ERROR = .false.
integer, public :: MAX_NUMBER_OF_ERRORS = huge(1)
integer         :: n_errors = 0

!This is crude ....
integer ::max_char_len = 100
contains

!-------------------------------------------------------------------------

recursive subroutine compare(reference,modified)
use flib_dom
! Given two nodes lists this subroutine compares all the nodes, one by one.
type(fnode), pointer :: reference
type(fnode), pointer :: modified

!Internal vars.
type(fnodelist),pointer  :: sub_list_ref, sub_list_mod, scf_ref,scf_mod
integer                  :: md_steps, md_ref,md_mod,i
type(fnode),pointer      :: ref,mod !, ref_tmp, ref_mod

!Parameters List
sub_list_ref => getElementsByTagName(reference,"parameterList")
sub_list_mod => getElementsByTagName(modified,"parameterList")
call compareList(sub_list_ref,sub_list_mod)

sub_list_ref => getElementsByTagAttrName(reference,"module","title","Initial System")
sub_list_mod => getElementsByTagAttrName(modified,"module","title","Initial System")
call compareList(sub_list_ref,sub_list_mod)

!kpoints
sub_list_ref => getElementsByTagAttrName(reference,"propertyList","title","k-points")
sub_list_mod => getElementsByTagAttrName(modified,"propertyList","title","k-points")
call compareList(sub_list_ref,sub_list_mod)

!kpoints
sub_list_ref => getElementsByTagAttrName(reference,"property","dictRef","siesta:kscell")
sub_list_mod => getElementsByTagAttrName(modified,"property","dictRef","siesta:kscell")
call compareList(sub_list_ref,sub_list_mod)

!more kpoints?
sub_list_ref => getElementsByTagAttrName(reference,"property","dictRef","siesta:kdispl")
sub_list_mod => getElementsByTagAttrName(modified,"property","dictRef","siesta:kdispl")
call compareList(sub_list_ref,sub_list_mod)

!MD steps
!Get the list of MD steps
sub_list_ref => getElementsByTagAttrName(reference,"module","dictRef","MD")
sub_list_mod => getElementsByTagAttrName(modified,"module","dictRef","MD")

md_ref = getLength(sub_list_ref)
md_mod = getLength(sub_list_mod)

if (md_ref /= md_mod) print *, "Warning: the number of MD steps differs: ref, mod=",md_ref,md_mod
md_steps = min(md_ref,md_mod)

!For each MD
do i=0,md_steps-1
!Get step i
   !print *,"starting md"
   
   ref => item(sub_list_ref,i)
   mod => item(sub_list_mod,i)

   !Geometry check
   call checkgeometry(ref,mod)
   
   !Get all the scf steps of this MD step
   scf_ref => getElementsByTagAttrName(ref,"module","dictRef","SCF")
   scf_mod => getElementsByTagAttrName(mod,"module","dictRef","SCF")
   !Compare them
   call compareList(scf_ref,scf_mod)

   scf_ref => getElementsByTagAttrName(reference,"module","title","SCF Finalization")
   scf_mod => getElementsByTagAttrName(modified,"module","title","SCF Finalization")
   call compareList(scf_ref,scf_mod)
end do

!Final things missing.
!Forces, stress, pressure, finalization block.
sub_list_ref => getElementsByTagAttrName(reference,"module","title","Finalization")
sub_list_mod => getElementsByTagAttrName(modified,"module","title","Finalization")

call compareList(sub_list_ref,sub_list_mod)

end subroutine compare


!--------------------------------------------------------------------------

recursive subroutine compareList(referenceList,modifiedList)
use flib_dom
! Given two nodes lists this subroutine compares all the nodes, one by one.
type(fnodeList), pointer :: referenceList
type(fnodeList), pointer :: modifiedList

!Internal vars.
integer                  :: i,length_ref,last_found
type(fnode),pointer      :: reference,modified
type(string)             :: sn_ref,sn_mod,sv_ref, sv_mod


length_ref = getLength(referenceList)

if(compare_node_list_debug) print *, "List begin"
if(compare_node_list_debug) print *, "---List: number of elements in reference list=",length_ref
if(compare_node_list_debug) print *, "---List: number of elements in modified  list=",getLength(modifiedList)


last_found = 0
do i=0,length_ref-1
   if(compare_node_list_debug) print *, "------List element:",i
   reference => item(referenceList,i)

   modified => Null()
   modified => findMatchingNodeinList(reference,modifiedList,last_found)

   if(.not. associated(modified))then
      if(compare_node_list_debug) then
         sn_ref = getNodeName(reference)  
         print *, "--------List:compare_node_list: node: ",char(sn_ref), &
           "--------has no corresponding node in the modified list!"
      endif
      cycle
   endif

   if (compare_node_list_debug) then
      !print *, "    compare_node_list: begin element=",i,"^^^^^^^^^^^^^^^^^^^"
      sn_ref = getNodeName(reference)      
      sn_mod = getNodeName(modified)
      sv_ref = getNodevalue(reference)
      sv_mod = getNodevalue(modified)
      !print *, "--------List: names=",char(sn_ref),"|",char(sn_mod)

      if (len(sv_ref) > 0 .and. len(sv_mod) > 0) then
      !   print *, "--------list: values=",char(sv_ref),"|",char(sv_mod)
      endif
   endif

   call compareNode(reference,modified) 

   if (compare_node_list_debug) then
      print *, "---List element:",i," end"
   endif

enddo

if (compare_node_list_debug) then
   print *, "---list: end"
   print *, "            "
   print *, "            "
endif

end subroutine compareList

!--------------------------------------------------------------------------

function findMatchingNodeinList(ref,modList,last) result(mod)
!Given a reference node this sub finds the corresponding node in modList.
  use flib_dom

  type(fnode), pointer        :: ref
  type(fnodeList), pointer    :: modList
  type(fnode), pointer        :: mod
  integer, intent(inout)      :: last

  integer               :: i,j,k,length_ref,length_mod
  type(string)          :: sn_ref, sn_mod,name_attr_ref,name_attr_mod,val_attr_ref,val_attr_mod
  type(fnode), pointer  :: modified 
  type(fnamedNodeMap),pointer :: attr_ref,attr_mod
  type(fnode), pointer :: attr_ref_j, attr_mod_j
  logical              :: error = .false.


  mod => Null()
  sn_ref = getNodeName(ref)  

  if ( matching_debug) print *,"********Match of:",char(sn_ref)
    
  do i=last,getLength(modList)-1
     if ( matching_debug)  print *, "*********Match: mod i=",i
     modified => Item(modList,i)
     sn_mod = getNodeName(modified)

      if (sn_ref == sn_mod) then
         if(.not. hasAttributes(ref) .and. .not. hasAttributes(modified))then
            mod => Item(modList,i)
            if ( matching_debug)  print *, "**********Match: found (no attr)"
            last = i
            exit !We found the node
        
         else !Check Attributes
            attr_ref => getAttributes(ref)
            attr_mod => getAttributes(modified)

            length_ref = getlength(attr_ref)
            length_mod = getlength(attr_mod)
 
            if (length_ref .ne. length_mod) then
               error = .true.
               cycle !If the number of attrs is different, they are different!!!
            else !Check the attrs, one by one.
               error = .false.
               do j=0,length_ref-1
                  !print *, "     Attrs j=",j
                  attr_ref_j => Item(attr_ref,j)
                  name_attr_ref = getNodeName(attr_ref_j)
                  val_attr_ref =  getNodeValue(attr_ref_j)

                  attr_mod_j => Item(attr_mod,j)
                  name_attr_mod = getNodeName(attr_mod_j)   
                  val_attr_mod =  getNodeValue(attr_mod_j)  
   

                  if ( matching_debug) then
                     print *, "************* Match: Attrs name:",trim(char(name_attr_ref)),"|",trim(char(name_attr_mod)),"|"

                     print *, "************* Match: Attrs value:",trim(char(val_attr_ref)),"|",trim(char(val_attr_mod)),"|"
                  endif

                  if(name_attr_ref .ne. name_attr_mod .or. val_attr_ref .ne. val_attr_mod)then
                     error = .true.
                     exit
                  endif

               enddo
               if (.not. error) then
                  mod =>Item(modList,i) !If we reached this point all the attrs and name are the same  
                  if ( matching_debug)  print *, "*********Match: Found corresponding node!"
                  last = i
                  exit
               endif                             
            endif
         endif
      endif
   enddo
   if ( matching_debug .and. error)  print *, "*********Match: Cant find corresponding node!"
  
end function findMatchingNodeinList

!--------------------------------------------------------------------------

recursive subroutine compareNode(ref,mod)
  use flib_dom

  type(fnode),pointer :: ref
  type(fnode),pointer :: mod
  

  integer :: i,n_type_ref, n_type_mod !,length_ref,length_mod , length
  integer :: n_children_ref, n_children_mod, n_children
  type(string) :: sn_ref,sn_mod,sv_ref, sv_mod
  character(len=50) ::  label,cn_parent_ref
  !type(fnamedNodeMap),pointer :: attr_ref,attr_mod
  !type(fnode), pointer :: c_node_ref, c_node_mod
  type(fnode), pointer  :: sub_ref,sub_mod
  type(fnodelist),pointer  :: sub_list_ref, sub_list_mod 
  logical :: same_values = .true.
  

  !Check name
  sn_ref= getNodeName(ref)
  sn_mod= getNodeName(mod)

  if (char(sn_ref) /= "" .and. char(sn_mod) /= "") then
     sn_ref = clean_string(sn_ref)
     sn_mod = clean_string(sn_mod)
  endif
  if (compare_node_debug) then
     !print *,"          1++++++++++++++++++++++++++++++++++++++++++++++++++"
     print *,"+++++++++++++Node: Names=","|",char(sn_ref),"|", char(sn_mod),"|"
  endif

  if (char(sn_ref) /= char(sn_mod)) call dump_error(ref,mod,name=sn_ref)

  !Check node type
  n_type_ref = getNodeType(ref)
  n_type_mod = getNodeType(mod)

  !if (n_type_ref /= n_type_mod ) call dump_error(ref,mod,type=n_type_ref)

  !Check node values. This part is the most complicate.
  sv_ref = getNodevalue(ref)
  sv_mod = getNodevalue(mod)

 

  !Check values
  if (len(sv_ref) == 0) then
  
     if (compare_node_debug) print *, "+++++++++++++ Node:No values, so no problem"
  
  else 
  
     sv_ref = clean_string(sv_ref)
     sv_mod = clean_string(sv_mod)

     if (compare_node_debug) print *,"+++++++++++++ Node: Value=","|",char(sv_ref),"|", char(sv_mod),"|"
    


     !If the name is empty then look for parents name.
     if (sn_ref == "" .or. sn_ref=="#text")then
        call getParentNodeProperties(ref,cn_parent_ref)
        label = trim(cn_parent_ref)
     else
        label = trim(char(sn_ref))
     endif

     !Main comparision
     same_values = compare_values(ref,mod,label)

     if ( .not. same_values) then
        call dump_error(ref,mod,value=sv_ref)
        return
     else if(compare_node_debug)then
         print *,"++++++++++No errors, they are identical!"
     endif
  endif
  
  !Check children
  if ( hasChildNodes(ref) .and. hasChildNodes(mod)) then
   sub_list_ref   => getchildNodes(ref)
   n_children_ref =  getLength(sub_list_ref)

   sub_list_mod   => getchildNodes(mod)

   
   if(compare_node_debug) print *, "++++++++++++ Node: compare children"
   call compareList(sub_list_ref,sub_list_mod)
   if(compare_node_debug) print *, "++++++++++++ Node: end compare children"

endif

 

endsubroutine compareNode

!-----------------------------------------------------------

recursive subroutine getParentNodeproperties(node,name)
type(fnode), pointer                   :: node
character(len=50)                      :: name  !Name

!Internal vars
integer              :: i,length
type(fnode), pointer :: parent,attr
type(string)         :: n_attr
type(fnamedNodeMap),pointer :: attributes
character(len=50), save     :: oldName

parent => NULL()
parent => getParentNode(node)

if (.not.associated(parent)) then
   name=""
else
   name = trim(getNodeName(parent))
   
   !Not useful go up
   if (name == "" .or. name == "#text" .or. name == "scalar" .or. name == "array" &
        .or.name == "matrix")then ! .or. name == "property") then 
      call getParentNodeProperties(parent,name)
   elseif(name == "parameter" .or. name == "property")then

      attributes => NULL()
      if (associated(parent)) attributes => getAttributes(parent)

      if (.not. associated(attributes)) then
         call getParentNodeProperties(parent,name)
      else
         !Loop over the attributes and find the one with the name.
         length = getLength(attributes)
         do i=0,length-1
            attr => item(attributes,i)
            n_attr = getNodeName(attr)
            if ( n_attr /= "" .or. n_attr /= "#text" .or. n_attr /= "scalar" &
                 .or. n_attr == "title" .or. n_attr == "dictRef")then
               name = getNodeValue(attr)
               exit
            endif
         enddo
      endif
   endif
endif

if(name=="")then
   name=oldName
else
   oldName=name
endif

end subroutine getParentNodeproperties

!-----------------------------------------------

subroutine dump_error(ref,mod,name,value,attr)
type(fnode), pointer :: ref,mod

type(string), intent(in), optional ::  name, value
!integer, intent(in), optional :: type
type(fnamedNodeMap),intent(in),optional :: attr

!Internal vars.
!integer :: n_type_ref, n_type_mod
type(string) :: sn_ref,sn_mod,sv_ref, sv_mod

logical ::  name_e = .false., value_e = .false., attr_e= .false.
character(len=50) :: n_parent_ref,n_parent_mod 
type(fnamedNodeMap),pointer :: attr_ref,attr_mod


!Find the type of error

!Find if the error is in the name:
if (present(name)) name_e = .true.

!Find if the error is in the value:
if (present(value)) value_e = .true.

!Find if the errror is in the attr.
if (present(attr)) attr_e = .true.


!Find all the info to make the error report as
!complete as possible.

!Find parent properties.
call getParentNodeProperties(ref,n_parent_ref)
call getParentNodeProperties(mod,n_parent_mod)

!Name 
sn_ref = getNodeName(ref)
sn_mod = getNodeValue(mod)

!If the name is empty then use the parents name.
if (len(sn_ref) == 0) then
   sn_ref = n_parent_ref
else
    if(sn_ref=="#text") then
       sn_ref = n_parent_ref
    else
       sn_ref = clean_string(sn_ref)
       sn_mod = clean_string(sn_mod)
    endif
endif

!Type
!n_type_ref = getNodeType(ref)
!n_type_mod = getNodeType(mod)

!Values
sv_ref = getNodevalue(ref)
sv_mod = getNodevalue(mod)

if (len(sv_ref) > 0 .and. len(sv_mod) > 0)then
   sv_ref = clean_string(sv_ref)
   sv_mod = clean_string(sv_mod)
endif

!Attr.
attr_ref => getAttributes(ref)
attr_mod => getAttributes(mod)

if ( len(trim(sn_ref)) < 3 ) then
   print *,"parent:",n_parent_ref
   call dump_error_heading(sn_ref,parent=n_parent_ref)
else 
   call dump_error_heading(sn_ref)
endif

if (name_e) then  
   print *,"Different names:"
   if (len(sn_ref) >0) then
      print *,"  Ref: ",char(sn_ref)
   else
      print *,"  Ref: no name"
   endif
   if (len(sn_mod)>0) then
      print *,"  Mod: ",char(sn_mod)
   else
      print *,"  Mod: no name"
   endif
elseif(value_e)then
   print *, "Different values:"
   if (len(sv_ref) >0) then
      print *,"  Ref: ","|",char(sv_ref),"|"
   else
      print *,"  Ref: no value"
   endif
   if (len(sv_mod)>0) then
      print *,"  Mod: ","|",char(sv_mod),"|"
   else
      print *,"  Mod: no value"
   endif
elseif(attr_e)then
   print *, "Different attr:"
   !print *, char(sn_ref),char(sn_mod)
endif

call handle_error()
print*, "---------------------------------------------------------"
end subroutine dump_error

!---------------------------------------------------

subroutine dump_error_heading (str,parent)

type(string), intent(in) :: str
character(len=*), intent(in), optional :: parent
print*, "---------------------------------------------------------"
print *,"There is an error in node: ","|",trim(char(str)),"|"
if (present(parent)) print *,"   which is a son of node: ", parent
end subroutine dump_error_heading

!---------------------------------------------------

function compare_values(ref,mod,label)
  ! Function wich compares two scalars (introduced in string format)
  ! If they do not differ then compare_scalar = .true.
  type(fnode), pointer            :: ref
  type(fnode), pointer            :: mod
  logical                         :: compare_values
  character(len=*), intent(in)    :: label


  !Internal vars.
  type(string) ::  n_ref, v_ref, n_mod, v_mod
  real(dp) :: tol

  n_ref = getNodeName(ref)
  v_ref = getNodeValue(ref)

  n_mod = getNodeName(mod)
  v_mod = getNodeValue(mod)

  if (len(v_ref) > 0 .and. len(v_mod) > 0)then
        v_ref = clean_string(v_ref)
        v_mod = clean_string(v_mod)
  endif

  tol = findTol(trim(label))

  !check if there are numbers in the string.
  if ( only_numbers(v_ref) .and. only_numbers(v_mod) )then 
     if (compare_values_debug) print "(7a,f10.6)", "            compare_values: label,ref,mod,tol= ",&
          trim(label),", |",char(v_ref),"|",char(v_mod),"|, ",tol
     compare_values = compare_only_numbers(v_ref,v_mod,label,tol)
  else
     compare_values = compare_alpha(v_ref,v_mod)
  endif
 
end function compare_values

!---------------------------------------------------

function should_be_checked(node) result(should)
  type (fnode), pointer :: node
  logical :: should 


  character(len=400)       :: tag
  type(string)             :: name
 
  should = .true.

  !Look for the properties of this node.
  if(compare_debug) print *,"------------------------------------------"

  tag = getTagName(node)

  if(compare_debug)then
     print *,"Should be checked: Tag: ", "|",trim(tag),"|"
     name = getNodeName(node)
     print *,"Shoulb be checked: Name:",char(name)
  endif

  if (tag == "" .or. trim(tag) == "metadata") then
     should = .false.
  endif
  
  if(compare_debug)then
     print *,"Should_be_checked: should?",should
     print *,"----------------------------------"
  endif

end function should_be_checked

subroutine handle_error()

if (STOP_ON_ERROR) then
   STOP
else
   n_errors = n_errors + 1
   if (n_errors == MAX_NUMBER_OF_ERRORS) then
      STOP
   else
      continue
   endif
endif
end subroutine handle_error

!---------------------------------------------------------------------------------

subroutine checkGeometry(ref, mod)
!For a given MD step this subroutine checks that both geometries are equivalent.
type(fnode),pointer      :: ref,mod

type(fnodeList), pointer :: ref_list,mod_list

ref_list => getElementsByTagName(ref,"molecule")
mod_list => getElementsByTagName(mod,"molecule")

call compareList(ref_list,mod_list)

ref_list => getElementsByTagName(ref,"crystal")
mod_list => getElementsByTagName(mod,"crystal")

call compareList(ref_list,mod_list)

end subroutine



end module compare_m
