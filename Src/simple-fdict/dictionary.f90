! @LICENSE@, see README.md
!
! Generated from master sources (https://github.com/zerothi/fdict)
! at commit 2c4534efacf
!
! Generic purpose dictionary as in any scripting language
! It has the power to contain any data by using the variable type.
module dictionary
  !! A key-value dictionary module to contain _any_ data in fortran.
  !!
  !! This module implements a generic dictionary-type (`type(dictionary_t)`)
  !! which may contain _any_ data-type using the `variable_t` data-type defined
  !! in `variable`.
  !!
  !! Example:
  !!
  !!```fortran
  !! real :: r
  !! real :: ra(10)
  !! real, target :: rb(10)
  !! type(dictionary_t) :: dict
  !! dict = ('Value'.kv.r)
  !!```
  !!
  !! Since it is f90 standard there are some limitations to the implementation
  !! that is perhaps not standard in usage.
  !! When concatenating two dictionaries one should suspect that the
  !! conactenated dictionaries are corrupted and should be nullifed after:
  !!
  !! ```fortran
  !! type(dictionary_t) :: d1, d2, d3
  !! ...
  !! d1 = d2
  !! call nullify(d2)
  !! call nullify(d3)
  !! ```
  !! One should _not_ operate on dictionaries which has
  !! been concatenated.
  !! Also, any keys that are in both will delete the *left* one.
  !! So if you have a pointer stored and you want to retain the values
  !! you will have to store the pointer somewhere else and nullify the key
  !! first.
  use, intrinsic :: iso_c_binding
  use variable
  implicit none
  private
  integer, parameter :: ih = selected_int_kind(4)
  integer, parameter :: is = selected_int_kind(9)
  integer, parameter :: il = selected_int_kind(18)
  integer, parameter :: sp = selected_real_kind(p=6)
  integer, parameter :: dp = selected_real_kind(p=15)
  ! Internal variables for determining the maximum size of the dictionaries.
  ! We could consider changing this to a variable size string
  ! However, that will increase the dependencies and will most likely not yield
  ! a better interface.
  !> Maximum character length of the keys in the dictionary, no
  !! index/key can be longer than this.
  integer, parameter :: DICTIONARY_KEY_LENGTH = 48
  public :: DICTIONARY_KEY_LENGTH
  ! A parameter returned if not found.
  character(len=DICTIONARY_KEY_LENGTH), parameter :: DICTIONARY_NOT_FOUND = 'ERROR: key not found'
  public :: DICTIONARY_NOT_FOUND
  !> The dictionary container it-self
  !!
  !! All contained variables are private.
  type :: dictionary_t
    ! We will keep the dictionary private so that any coding
    ! has to use .KEY. and .VAL. etc.
    type(dictionary_entry_t), pointer :: first => null()
    integer :: len = 0
  end type dictionary_t
  public :: dictionary_t
  ! We need to create a linked list to create arbitrarily long dictionaries...
  ! The dictionary entry is not visible outside.
  type :: dictionary_entry_t
    character(len=DICTIONARY_KEY_LENGTH) :: key = ' '
    ! in order to extend the dictionary to contain a dictionary
    ! we simply need to add the dictionary type to the variable
    ! library.
    type(variable_t) :: value
    integer :: hash = 0
    type(dictionary_entry_t), pointer :: next => null()
  end type dictionary_entry_t
  !> Return the length of a dictionary, by internal counting algorithms
  interface len
    module procedure len_
  end interface
  public :: LEN
  !> Actually count number of elements in the dictionary by forcing traversing the linked-list
  interface llen
    module procedure llen_
  end interface
  public :: LLEN
  !> Print out all keys and which data-type it contains as well as the hash-number
  interface print
    module procedure print_
  end interface
  public :: print
  ! Concatenate dicts or list of dicts to list of dicts
  !> Concatenate, or extend, dictionaries, this can
  !! be done on it-self `dic = dic // ('key'.kv.1)
  interface operator( // )
    module procedure d_cat_d
  end interface
  public :: operator( // )
  ! Retrieve the key from a dictionary (unary)
  !> Returns the key of the current _top_ entry,
  interface operator( .KEY. )
    module procedure key
  end interface operator( .KEY. )
  public :: operator(.KEY.)
  ! check whether key exists in dictionary
  !> Returns .true. if the key exists in the dictionary, else returns false.
  interface operator( .IN. )
    module procedure in
  end interface operator( .IN. )
  public :: operator(.IN.)
  ! check whether key not exists in dictionary
  !> Returns .not. ('key' .in. dict)
  interface operator( .NIN. )
    module procedure nin
  end interface operator( .NIN. )
  public :: operator(.NIN.)
  ! Retrieve the value from a dictionary (unary)
  !> Returns the value from a dictionary by copy
  interface operator( .VAL. )
    module procedure value
  end interface operator( .VAL. )
  public :: operator(.VAL.)
  !> Returns the value from a dictionary by pointer
  interface operator( .VALP. )
    module procedure value_p
  end interface operator( .VALP. )
  public :: operator(.VALP.)
  ! Retrieve the hash value from a dictionary entry (unary)
  interface operator( .HASH. )
    module procedure hash_
  end interface operator( .HASH. )
  public :: operator(.HASH.)
  interface hash
    module procedure hash_
  end interface hash
  public :: hash
  ! Checks for two dicts have all the same keys
  !> Checks whether all keys are the same in two dictionaries.
  interface operator( .EQ. )
    module procedure d_eq_d
  end interface operator( .EQ. )
  public :: operator(.EQ.) ! Overloaded
  ! Checks for two dicts do not share any common keys
  !> Checks whether not all keys are the same in two dictionaries.
  interface operator( .NE. )
    module procedure d_ne_d
  end interface operator( .NE. )
  public :: operator(.NE.) ! Overloaded
  ! Steps one time in the dictionary (unary)
  !> Looping construct.
  interface operator( .NEXT. )
    module procedure d_next
  end interface operator( .NEXT. )
  public :: operator(.NEXT.)
  interface next
    module procedure d_next
  end interface next
  public :: next
  ! Retrieve the first of a dictionary (unary)
  !> Returns the first entry
  interface operator( .FIRST. )
    module procedure d_first
  end interface operator( .FIRST. )
  public :: operator(.FIRST.)
  interface first
    module procedure d_first
  end interface first
  public :: first
  ! Check whether the dictionary is empty (unary)
  !> Checks if it is an empty dictionary, i.e. no keys exist
  interface operator( .EMPTY. )
    module procedure d_empty
    module procedure e_empty
  end interface operator( .EMPTY. )
  public :: operator(.EMPTY.)
  interface empty
    module procedure d_empty
    module procedure e_empty
  end interface empty
  public :: empty
  interface hash_coll
    module procedure hash_coll_
  end interface hash_coll
  public :: hash_coll
  interface delete
    module procedure delete_
  end interface delete
  public :: delete
  interface pop
    module procedure pop_
  end interface pop
  public :: pop
  interface copy
    module procedure copy_
  end interface copy
  public :: copy
  interface nullify
    module procedure nullify_
    module procedure nullify_key_
  end interface nullify
  public :: nullify
  interface extend
    module procedure sub_d_cat_d
  end interface extend
  public :: extend
  interface which
    module procedure d_key_which
  end interface which
  public :: which
  public :: assign, associate
interface operator(.KV.)
module procedure d_kv_a0_0
module procedure d_kv_var
module procedure d_kv_a1
module procedure d_kv_s0
module procedure d_kv_s1
module procedure d_kv_s2
module procedure d_kv_s3
module procedure d_kv_d0
module procedure d_kv_d1
module procedure d_kv_d2
module procedure d_kv_d3
module procedure d_kv_c0
module procedure d_kv_c1
module procedure d_kv_c2
module procedure d_kv_c3
module procedure d_kv_z0
module procedure d_kv_z1
module procedure d_kv_z2
module procedure d_kv_z3
module procedure d_kv_b0
module procedure d_kv_b1
module procedure d_kv_b2
module procedure d_kv_b3
module procedure d_kv_h0
module procedure d_kv_h1
module procedure d_kv_h2
module procedure d_kv_h3
module procedure d_kv_i0
module procedure d_kv_i1
module procedure d_kv_i2
module procedure d_kv_i3
module procedure d_kv_l0
module procedure d_kv_l1
module procedure d_kv_l2
module procedure d_kv_l3
module procedure d_kv_cp0
module procedure d_kv_cp1
module procedure d_kv_fp0
module procedure d_kv_fp1
end interface
interface operator(.KVP.)
module procedure d_kvp_var
module procedure d_kvp_dict
module procedure d_kvp_a1
module procedure d_kvp_s0
module procedure d_kvp_s1
module procedure d_kvp_s2
module procedure d_kvp_s3
module procedure d_kvp_d0
module procedure d_kvp_d1
module procedure d_kvp_d2
module procedure d_kvp_d3
module procedure d_kvp_c0
module procedure d_kvp_c1
module procedure d_kvp_c2
module procedure d_kvp_c3
module procedure d_kvp_z0
module procedure d_kvp_z1
module procedure d_kvp_z2
module procedure d_kvp_z3
module procedure d_kvp_b0
module procedure d_kvp_b1
module procedure d_kvp_b2
module procedure d_kvp_b3
module procedure d_kvp_h0
module procedure d_kvp_h1
module procedure d_kvp_h2
module procedure d_kvp_h3
module procedure d_kvp_i0
module procedure d_kvp_i1
module procedure d_kvp_i2
module procedure d_kvp_i3
module procedure d_kvp_l0
module procedure d_kvp_l1
module procedure d_kvp_l2
module procedure d_kvp_l3
module procedure d_kvp_cp0
module procedure d_kvp_cp1
module procedure d_kvp_fp0
module procedure d_kvp_fp1
end interface
interface assign
module procedure d_get_val
module procedure d_get_val_a_
module procedure d_get_val_a1
module procedure d_get_val_first_a1
module procedure d_get_val_s0
module procedure d_get_val_first_s0
module procedure d_get_val_s1
module procedure d_get_val_first_s1
module procedure d_get_val_s2
module procedure d_get_val_first_s2
module procedure d_get_val_s3
module procedure d_get_val_first_s3
module procedure d_get_val_d0
module procedure d_get_val_first_d0
module procedure d_get_val_d1
module procedure d_get_val_first_d1
module procedure d_get_val_d2
module procedure d_get_val_first_d2
module procedure d_get_val_d3
module procedure d_get_val_first_d3
module procedure d_get_val_c0
module procedure d_get_val_first_c0
module procedure d_get_val_c1
module procedure d_get_val_first_c1
module procedure d_get_val_c2
module procedure d_get_val_first_c2
module procedure d_get_val_c3
module procedure d_get_val_first_c3
module procedure d_get_val_z0
module procedure d_get_val_first_z0
module procedure d_get_val_z1
module procedure d_get_val_first_z1
module procedure d_get_val_z2
module procedure d_get_val_first_z2
module procedure d_get_val_z3
module procedure d_get_val_first_z3
module procedure d_get_val_b0
module procedure d_get_val_first_b0
module procedure d_get_val_b1
module procedure d_get_val_first_b1
module procedure d_get_val_b2
module procedure d_get_val_first_b2
module procedure d_get_val_b3
module procedure d_get_val_first_b3
module procedure d_get_val_h0
module procedure d_get_val_first_h0
module procedure d_get_val_h1
module procedure d_get_val_first_h1
module procedure d_get_val_h2
module procedure d_get_val_first_h2
module procedure d_get_val_h3
module procedure d_get_val_first_h3
module procedure d_get_val_i0
module procedure d_get_val_first_i0
module procedure d_get_val_i1
module procedure d_get_val_first_i1
module procedure d_get_val_i2
module procedure d_get_val_first_i2
module procedure d_get_val_i3
module procedure d_get_val_first_i3
module procedure d_get_val_l0
module procedure d_get_val_first_l0
module procedure d_get_val_l1
module procedure d_get_val_first_l1
module procedure d_get_val_l2
module procedure d_get_val_first_l2
module procedure d_get_val_l3
module procedure d_get_val_first_l3
module procedure d_get_val_cp0
module procedure d_get_val_first_cp0
module procedure d_get_val_cp1
module procedure d_get_val_first_cp1
module procedure d_get_val_fp0
module procedure d_get_val_first_fp0
module procedure d_get_val_fp1
module procedure d_get_val_first_fp1
end interface
interface associate
module procedure d_get_p_val
module procedure d_get_p_dict
module procedure d_get_p_a1
module procedure d_get_p_first_a1
module procedure d_get_p_s0
module procedure d_get_p_first_s0
module procedure d_get_p_s1
module procedure d_get_p_first_s1
module procedure d_get_p_s2
module procedure d_get_p_first_s2
module procedure d_get_p_s3
module procedure d_get_p_first_s3
module procedure d_get_p_d0
module procedure d_get_p_first_d0
module procedure d_get_p_d1
module procedure d_get_p_first_d1
module procedure d_get_p_d2
module procedure d_get_p_first_d2
module procedure d_get_p_d3
module procedure d_get_p_first_d3
module procedure d_get_p_c0
module procedure d_get_p_first_c0
module procedure d_get_p_c1
module procedure d_get_p_first_c1
module procedure d_get_p_c2
module procedure d_get_p_first_c2
module procedure d_get_p_c3
module procedure d_get_p_first_c3
module procedure d_get_p_z0
module procedure d_get_p_first_z0
module procedure d_get_p_z1
module procedure d_get_p_first_z1
module procedure d_get_p_z2
module procedure d_get_p_first_z2
module procedure d_get_p_z3
module procedure d_get_p_first_z3
module procedure d_get_p_b0
module procedure d_get_p_first_b0
module procedure d_get_p_b1
module procedure d_get_p_first_b1
module procedure d_get_p_b2
module procedure d_get_p_first_b2
module procedure d_get_p_b3
module procedure d_get_p_first_b3
module procedure d_get_p_h0
module procedure d_get_p_first_h0
module procedure d_get_p_h1
module procedure d_get_p_first_h1
module procedure d_get_p_h2
module procedure d_get_p_first_h2
module procedure d_get_p_h3
module procedure d_get_p_first_h3
module procedure d_get_p_i0
module procedure d_get_p_first_i0
module procedure d_get_p_i1
module procedure d_get_p_first_i1
module procedure d_get_p_i2
module procedure d_get_p_first_i2
module procedure d_get_p_i3
module procedure d_get_p_first_i3
module procedure d_get_p_l0
module procedure d_get_p_first_l0
module procedure d_get_p_l1
module procedure d_get_p_first_l1
module procedure d_get_p_l2
module procedure d_get_p_first_l2
module procedure d_get_p_l3
module procedure d_get_p_first_l3
module procedure d_get_p_cp0
module procedure d_get_p_first_cp0
module procedure d_get_p_cp1
module procedure d_get_p_first_cp1
module procedure d_get_p_fp0
module procedure d_get_p_first_fp0
module procedure d_get_p_fp1
module procedure d_get_p_first_fp1
end interface
  ! Create a dict type: 'key' .KV. 'val'
  public :: operator(.KV.)
  ! Create a dict type: 'key' .KVP. 'pointer'
  public :: operator(.KVP.)
contains
  pure function hash_val(key) result(val)
    character(len=*), intent(in) :: key
    integer :: val
    integer :: i
    ! This is 32-bit integers, hence a 32-bit hash
    integer, parameter :: FNV_OFF = 28491 ! see crt_hash_basis.f90
    integer, parameter :: FNV_PRIME = 16777619
    integer, parameter :: MAX_32 = huge(1)
    ! Initialize by the FNV_OFF hash for 32 bit
    val = FNV_OFF
    do i = 1 , min(DICTIONARY_KEY_LENGTH,len_trim(key))
      val = ieor(val,iachar(key(i:i)))
      val = mod(val * FNV_PRIME, MAX_32)
    end do
  end function hash_val
  pure &
      function new_d_key(key) result(d)
    character(len=*), intent(in) :: key
    type(dictionary_t) :: d
    allocate(d%first)
    if ( len_trim(key) > DICTIONARY_KEY_LENGTH ) then
      d%first%key = key(1:DICTIONARY_KEY_LENGTH)
    else
      d%first%key = trim(key)
    end if
    d%first%hash = hash_val(key)
    d%len = 1
    nullify(d%first%next)
  end function new_d_key
  ! Retrieves the key value in a dictionary type (or a list)
  ! We expect that the key will only be called on single element dictionaries...
  pure function key(d)
    type(dictionary_t), intent(in) :: d
    character(len=DICTIONARY_KEY_LENGTH) :: key
    key = d%first%key
  end function key
  ! Retrieves the value value in a dictionary type (or a list)
  function value(d)
    type(dictionary_t), intent(in) :: d
    type(variable_t) :: value
    call assign(value,d%first%value)
  end function value
  function value_p(d)
    type(dictionary_t), intent(in) :: d
    type(variable_t) :: value_p
    call associate(value_p,d%first%value)
  end function value_p
  ! Returns the hash value of the dictionary first item...
  pure function hash_(d)
    type(dictionary_t), intent(in) :: d
    integer :: hash_
    hash_ = d%first%hash
  end function hash_
  ! Returns number of collisions in the hash-table
  ! The optional keyword 'max' can be used to
  ! extract the maximum number of collisions for
  ! one hash-value (i.e. not total collisions).
  function hash_coll_(this,max) result(col)
    type(dictionary_t), intent(inout) :: this
    logical, intent(in), optional :: max
    integer :: col
    integer :: chash, max_now, same
    type(dictionary_entry_t), pointer :: ld
    col = 0
    if ( .empty. this ) return
    ! Initialize
    same = 0
    max_now = 0
    ld => this%first
    chash = ld%hash
    do while ( associated(ld) )
      if ( chash == ld%hash ) then
        ! total collisions
        col = col + 1
        ! count total current collisions
        max_now = max_now + 1
      else
        chash = ld%hash
        if ( max_now > same ) then
          same = max_now
        end if
        max_now = 0
      end if
      ld => ld%next
    end do
    ! If the user requests maximum collisions
    ! for any given hash value
    if ( present(max) ) then
      if ( max ) col = same
    end if
    ! return col
  end function hash_coll_
  function in(key,d)
    character(len=*), intent(in) :: key
    type(dictionary_t), intent(in) :: d
    type(dictionary_t) :: ld
    integer :: hash, lhash
    logical :: in
    hash = hash_val(key)
    ld = .first. d
    search: do while ( .not. (.empty. ld) )
      lhash = .hash. ld
      if ( hash > lhash ) then
        ! skip to next search
      else if ( hash < lhash ) then
        exit search
      else if ( hash == lhash ) then
        if ( key .eq. (.KEY. ld) ) then
          in = .true.
          return
        end if
      end if
      ld = .next. ld
    end do search
    in = .false.
  end function in
  function nin(key,d)
    character(len=*), intent(in) :: key
    type(dictionary_t), intent(in) :: d
    logical :: nin
    nin = .not. in(key,d)
  end function nin
  ! Compares two dict types against each other
  ! Will do comparison by hash.
  function d_eq_d(d1,d2) result(bool)
    type(dictionary_t), intent(in) :: d1,d2
    logical :: bool
    type(dictionary_t) :: tmp1, tmp2
    bool = len(d1) == len(d2)
    if ( .not. bool ) return
    bool = (.hash. d1) == (.hash. d2)
    if ( .not. bool ) return
    ! if all the keys are going to be the same
    ! the we know that the hash-tags are going to
    ! be the same... :)
    tmp1 = .first. d1
    tmp2 = .first. d2
    do while ( .not. (.empty. tmp1) )
      bool = (.hash. tmp1) == (.hash. tmp2)
      if ( .not. bool ) return
      tmp1 = .next. tmp1
      tmp2 = .next. tmp2
    end do
  end function d_eq_d
  ! Compares two dict types against each other
  ! not necessarily the negative of .eq.
  function d_ne_d(d1,d2) result(bool)
    type(dictionary_t), intent(in) :: d1,d2
    logical :: bool
    type(dictionary_t) :: tmp1, tmp2
    tmp1 = .first. d1
    do while ( .not. (.empty. tmp1) )
      tmp2 = .first. d2
      do while ( .not. (.empty. tmp2) )
        bool = (.hash. tmp1) == (.hash. tmp2)
        if ( bool ) then
          bool = .false.
          return
        end if
        tmp2 = .next. tmp2
      end do
      tmp1 = .next. tmp1
    end do
  end function d_ne_d
  ! Concatenate two dictionaries to one dictionary...
  ! it does not work with elemental as the
  function d_cat_d(d1,d2) result(d)
    type(dictionary_t), intent(in) :: d1,d2
    type(dictionary_t) :: d
    if ( .empty. d1 ) then
      if ( .empty. d2 ) return
      call copy_assign(d2,d)
      return
    end if
    call copy_assign(d1,d)
    call sub_d_cat_d(d,d2)
  end function d_cat_d
  ! Concatenate two dictionaries to one dictionary...
  ! it does not work with elemental as the
  subroutine sub_d_cat_d(d,d2)
    type(dictionary_t), intent(inout) :: d
    type(dictionary_t), intent(in) :: d2
    type(dictionary_entry_t), pointer :: ladd, lnext
    type(dictionary_t) :: fd
    integer :: kh
    if ( .empty. d ) then
      if ( .empty. d2 ) return
      call copy_assign(d2,d)
      return
    end if
    if ( .empty. d2 ) return
    ladd => d2%first
    fd%len = 0
    fd%first => d%first
    do
      ! step ...
      lnext => ladd%next ! we need to get the next
      kh = fd%first%hash
      call d_insert(fd,ladd)
      if ( .not. associated(lnext) ) exit
      ladd => lnext
    end do
    ! Update the first entry
    d%first => fd%first
    d%len = d%len + fd%len
  end subroutine sub_d_cat_d
  subroutine d_insert(d,entry)
    type(dictionary_t), intent(inout) :: d
    type(dictionary_entry_t), intent(inout), pointer :: entry
    type(dictionary_entry_t), pointer :: search, prev
    ! if the dictionary is empty
    ! simply put it first
    if ( .not. associated(d%first) ) then
      d%first => entry
      d%len = 1
      return
    end if
    nullify(prev)
    search => d%first
    ! Matching the first entry
    if ( search%hash > entry%hash ) then
      ! The added hash is smaller than the first hash
      ! of the dictionary, so we get a new *first* entry
      entry%next => d%first
      d%first => entry
      d%len = d%len + 1
      return
    else if ( search%hash == entry%hash ) then
      ! If the key already exists we will simply overwrite
      if ( search%key == entry%key ) then
        ! deletion
        call delete(search%value)
        ! key and hash are the same, no need to transfer those
        entry%next => search%next
        d%first => entry
        deallocate(search)
        return
      end if
    end if
    search_loop: do
      ! step...
      prev => search
      ! step...
      search => prev%next
      if ( .not. associated(search) ) exit search_loop
      if ( search%hash > entry%hash ) then
        prev%next => entry
        entry%next => search
        d%len = d%len + 1
        return
      else if ( search%hash == entry%hash ) then
        ! If the key already exists we will simply overwrite
        ! If not, we will put it *after* the same hash
        if ( search%key == entry%key ) then
          call delete(search%value)
          prev%next => entry
          entry%next => search%next
          deallocate(search)
          return
        end if
      end if
    end do search_loop
    prev%next => entry
    ! Increment length of the dictionary...
    d%len = d%len + 1
    ! As we could insert from a dictionary we have to reset, to not do endless loops...
    nullify(entry%next)
  end subroutine d_insert
  !> Generate the copy routine
  subroutine copy_(from, to)
    type(dictionary_t), intent(in) :: from
    type(dictionary_t), intent(inout) :: to
    type(dictionary_entry_t), pointer :: d
    type(variable_t) :: v
    ! Delete the dictionary
    call delete(to)
    d => from%first
    do while ( associated(d) )
      ! Associate data...
      call associate(v, d%value)
      ! Copy data, hence .kv.
      to = to
      d => d%next
    end do
    ! Clean up pointers...
    call nullify(v)
    nullify(d)
  end subroutine copy_
  ! Retrieve the length of the dictionary...
  pure function len_(d)
    type(dictionary_t), intent(in) :: d
    integer :: len_
    len_ = d%len
  end function len_
  function llen_(this)
    type(dictionary_t), intent(inout) :: this
    type(dictionary_entry_t), pointer :: d
    integer :: llen_
    llen_ = 0
    d => this%first
    do while ( associated(d) )
      llen_ = llen_ + 1
      d => d%next
    end do
  end function llen_
  function d_next(d)
    type(dictionary_t), intent(in) :: d
    type(dictionary_t) :: d_next
    d_next%first => d%first%next
    d_next%len = d%len - 1
  end function d_next
  pure function d_empty(d)
    type(dictionary_t), intent(in) :: d
    logical :: d_empty
    d_empty = .not. associated(d%first)
  end function d_empty
  pure function e_empty(e)
    type(dictionary_entry_t), intent(in) :: e
    logical :: e_empty
    e_empty = empty(e%value)
  end function e_empty
  function d_first(d)
    type(dictionary_t), intent(in) :: d
    type(dictionary_t) :: d_first
    call copy_assign(d,d_first)
  end function d_first
  subroutine copy_assign(din,dcopy)
    type(dictionary_t), intent(in) :: din
    type(dictionary_t), intent(inout) :: dcopy
    dcopy%first => din%first
    dcopy%len = din%len
  end subroutine copy_assign
  subroutine print_(d)
    type(dictionary_t), intent(in) :: d
    type(dictionary_t) :: ld
    ld = .first. d
    do while ( .not. (.empty. ld) )
      write(*,'(t2,a,tr1,a,i0,a)') trim(.key. ld), &
          '['//trim(ld%first%value%t)//'] (',.hash. ld,')'
      ld = .next. ld
    end do
  end subroutine print_
  subroutine delete_(this,key,dealloc)
    type(dictionary_t), intent(inout) :: this
    character(len=*), intent(in), optional :: key
    logical, intent(in), optional :: dealloc
    type(dictionary_entry_t), pointer :: de, pr
    logical :: ldealloc
    integer :: kh
    ! We default to de-allocation of everything
    ldealloc = .true.
    if ( present(dealloc) ) ldealloc = dealloc
    ! if no keys are present, simply return
    if ( .not. associated(this%first) ) then
      this%len = 0
      return
    end if
    if ( present(key) ) then
      ! we only need to delete the one key
      kh = hash_val(key)
      pr => this%first
      if ( kh == pr%hash ) then
        if ( key == pr%key ) then
          this%first => pr%next
          this%len = this%len - 1
          call delete(pr%value,dealloc=ldealloc)
          nullify(pr%next)
          deallocate(pr)
          nullify(pr)
          return
        end if
      end if
      ! more complicated case
      de => pr%next
      do while ( associated(de) )
        ! We know it is sorted with hash-tags.
        ! So if we are beyond the hash, we just quit.
        if ( kh < de%hash ) exit ! it does not exist
        if ( de%hash == kh ) then
          if ( de%key == key ) then
            pr%next => de%next
            call delete(de%value,dealloc=ldealloc)
            nullify(de%next)
            deallocate(de)
            this%len = this%len - 1
            exit
          end if
        end if
        pr => de
        de => de%next
      end do
      return
    end if
    ! delete the entire entry-tree
    call del_dictionary_entry_t_tree(this%first,dealloc=ldealloc)
    call delete(this%first%value,dealloc=ldealloc)
    deallocate(this%first)
    nullify(this%first)
    this%len = 0
  contains
    recursive subroutine del_dictionary_entry_t_tree(d,dealloc)
      type(dictionary_entry_t), pointer :: d
      logical, intent(in) :: dealloc
      if ( associated(d) ) then
        if ( associated(d%next) ) then
          call del_dictionary_entry_t_tree(d%next,dealloc)
          call delete(d%next%value,dealloc=dealloc)
          deallocate(d%next)
          nullify(d%next)
        end if
      end if
    end subroutine del_dictionary_entry_t_tree
  end subroutine delete_
  subroutine pop_(val,this,key,dealloc)
    type(variable_t), intent(inout) :: val
    type(dictionary_t), intent(inout) :: this
    character(len=*), intent(in) :: key
    logical, intent(in), optional :: dealloc
    type(dictionary_entry_t), pointer :: de, pr
    ! Here the default is to de-allocate
    ! even though we use the association feature
    ! Hence, we need a variable here
    logical :: ldealloc
    integer :: kh
    ldealloc = .true.
    if ( present(dealloc) ) ldealloc = dealloc
    ! if no keys are present, simply return
    if ( .not. associated(this%first) ) then
      this%len = 0
      call val_delete_request(val,dealloc=ldealloc)
      return
    end if
    pr => this%first
    if ( pr%key == key ) then
      this%first => pr%next
      call associate(val,pr%value,dealloc=ldealloc)
      ! Ensures that the encoding gets removed
      call nullify(pr%value)
      deallocate(pr)
      this%len = this%len - 1
      return
    end if
    kh = hash_val(key)
    de => pr%next
    do while ( associated(de) )
      ! Check if even exists
      if ( kh < de%hash ) exit
      if ( kh == de%hash ) then
        if ( de%key == key ) then
          pr%next => de%next
          call associate(val,de%value,dealloc=ldealloc)
          ! Ensures that the encoding gets removed
          call nullify(de%value)
          deallocate(de)
          this%len = this%len - 1
          exit
        end if
      end if
      pr => de
      de => de%next
    end do
  end subroutine pop_
  elemental subroutine nullify_key_(this,key)
    type(dictionary_t), intent(inout) :: this
    character(len=*), intent(in) :: key
    type(dictionary_entry_t), pointer :: de, pr
    integer :: kh
    ! if no keys are present, simply return
    if ( .not. associated(this%first) ) then
      this%len = 0
      return
    end if
    pr => this%first
    if ( pr%key == key ) then
      this%first => pr%next
      ! Ensures that the encoding gets removed
      call nullify(pr%value)
      deallocate(pr)
      this%len = this%len - 1
      return
    end if
    kh = hash_val(key)
    de => pr%next
    do while ( associated(de) )
      ! Check if even exists
      if ( kh < de%hash ) exit
      if ( kh == de%hash ) then
        if ( de%key == key ) then
          pr%next => de%next
          ! Ensures that the encoding gets removed
          call nullify(de%value)
          deallocate(de)
          this%len = this%len - 1
          exit
        end if
      end if
      pr => de
      de => de%next
    end do
  end subroutine nullify_key_
  elemental subroutine nullify_(this)
    type(dictionary_t), intent(inout) :: this
    ! This will simply nullify the dictionary, thereby
    ! remove all references to all objects.
    nullify(this%first)
    this%len = 0
  end subroutine nullify_
  subroutine d_get_val(val,d,key,dealloc)
    type(variable_t), intent(inout) :: val
    type(dictionary_t), intent(inout) :: d
    character(len=*), intent(in), optional :: key
    logical, intent(in), optional :: dealloc
    type(dictionary_t) :: ld
    integer :: hash, lhash
    if ( .not. present(key) ) then
      if ( .not. (.empty. d) ) then
        call assign(val,d%first%value,dealloc=dealloc)
      else
        call val_delete_request(val,dealloc=dealloc)
      end if
      return
    end if
    hash = hash_val(key)
    ld = .first. d
    search: do while ( .not. (.empty. ld) )
      lhash = .hash. ld
      if ( hash > lhash ) then
        ! skip to next search
      else if ( hash < lhash ) then
        ! the key does not exist, delete if requested, else clean it
        call val_delete_request(val,dealloc=dealloc)
        exit search
      else if ( hash == lhash ) then
        if ( key .eq. (.KEY. ld) ) then
          call assign(val,ld%first%value,dealloc=dealloc)
          return
        end if
      end if
      ld = .next. ld
    end do search
  end subroutine d_get_val
  subroutine d_get_p_val(val,d,key,dealloc)
    type(variable_t), intent(inout) :: val
    type(dictionary_t), intent(inout) :: d
    character(len=*), intent(in), optional :: key
    logical, intent(in), optional :: dealloc
    type(dictionary_t) :: ld
    integer :: hash, lhash
    if ( .not. present(key) ) then
      if ( .not. (.empty. d) ) then
        call associate(val,d%first%value,dealloc=dealloc)
      else
        call val_delete_request(val,dealloc=dealloc)
      end if
      return
    end if
    hash = hash_val(key)
    ld = .first. d
    search: do while ( .not. (.empty. ld) )
      lhash = .hash. ld
      if ( hash > lhash ) then
        ! skip to next search
      else if ( hash < lhash ) then
        ! the key does not exist, delete if requested, else clean it
        call val_delete_request(val,dealloc=dealloc)
        exit search
      else if ( hash == lhash ) then
        if ( key .eq. (.KEY. ld) ) then
          call associate(val,ld%first%value,dealloc=dealloc)
          return
        end if
      end if
      ld = .next. ld
    end do search
  end subroutine d_get_p_val
  subroutine d_get_val_a_(val,d,key,dealloc)
    character(len=*), intent(out) :: val
    type(dictionary_t), intent(inout) :: d
    character(len=*), intent(in), optional :: key
    logical, intent(in), optional :: dealloc
    type(variable_t) :: v
    type(dictionary_t) :: ld
    integer :: hash, lhash
    val = ' '
    if ( .not. present(key) ) then
      if ( .not. (.empty. d) ) then
        call associate(v,d%first%value)
      end if
      return
    end if
    hash = hash_val(key)
    ld = .first. d
    search: do while ( .not. (.empty. ld) )
      lhash = .hash. ld
      if ( hash > lhash ) then
        ! skip to next search
      else if ( hash < lhash ) then
        exit search
      else if ( hash == lhash ) then
        if ( key .eq. (.KEY. ld) ) then
          call assign(val, ld%first%value)
          return
        end if
      end if
      ld = .next. ld
    end do search
  end subroutine d_get_val_a_
  function d_kv_a0_0(key,val) result(this)
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: val
    type(dictionary_t) :: this
    this = new_d_key(key)
    call assign(this%first%value,val)
  end function d_kv_a0_0
  function d_kv_var(key,val) result(this)
    character(len=*), intent(in) :: key
    type(variable_t), intent(in) :: val
    type(dictionary_t) :: this
    this = new_d_key(key)
    call assign(this%first%value,val)
  end function d_kv_var
  function d_kvp_var(key,val) result(this)
    character(len=*), intent(in) :: key
    type(variable_t), intent(in) :: val
    type(dictionary_t) :: this
    this = new_d_key(key)
    call associate(this%first%value,val)
  end function d_kvp_var
  function d_key_which(this,key) result(t)
    type(dictionary_t), intent(in) :: this
    character(len=*), optional, intent(in) :: key
    character(len=VARIABLE_TYPE_LENGTH) :: t
    type(dictionary_t) :: ld
    integer :: hash, lhash
    if ( present(key) ) then
      hash = hash_val(key)
      ld = .first. this
      search: do while ( .not. (.empty. ld) )
        lhash = .hash. ld
        if ( hash > lhash ) then
          ! skip to next search
        else if ( hash < lhash ) then
          t = '  '
          exit search
        else if ( hash == lhash ) then
          if ( key .eq. (.KEY. ld) ) then
            t = which(ld%first%value)
            return
          end if
        end if
        ld = .next. ld
      end do search
    else
      t = which(this%first%value)
    end if
  end function d_key_which
function d_kv_a1(key,val) result(this)
character(len=*), intent(in) :: key
character(len=1), intent(in), dimension(:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_a1
function d_kvp_a1(key, val) result(this)
character(len=*), intent(in) :: key
character(len=1), intent(in), dimension(:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_a1
subroutine d_get_val_a1(val,this,key,success)
character(len=1), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_a1
subroutine d_get_val_first_a1(val,this,success)
character(len=1), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_a1
subroutine d_get_p_a1(val,this,key,success)
character(len=1), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_a1
subroutine d_get_p_first_a1(val,this,success)
character(len=1), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_a1
function d_kv_s0(key,val) result(this)
character(len=*), intent(in) :: key
real(sp), intent(in) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_s0
function d_kvp_s0(key, val) result(this)
character(len=*), intent(in) :: key
real(sp), intent(in), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_s0
subroutine d_get_val_s0(val,this,key,success)
real(sp), intent(out) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_s0
subroutine d_get_val_first_s0(val,this,success)
real(sp), intent(out) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_s0
subroutine d_get_p_s0(val,this,key,success)
real(sp), pointer :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_s0
subroutine d_get_p_first_s0(val,this,success)
real(sp), pointer :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_s0
function d_kv_s1(key,val) result(this)
character(len=*), intent(in) :: key
real(sp), intent(in), dimension(:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_s1
function d_kvp_s1(key, val) result(this)
character(len=*), intent(in) :: key
real(sp), intent(in), dimension(:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_s1
subroutine d_get_val_s1(val,this,key,success)
real(sp), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_s1
subroutine d_get_val_first_s1(val,this,success)
real(sp), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_s1
subroutine d_get_p_s1(val,this,key,success)
real(sp), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_s1
subroutine d_get_p_first_s1(val,this,success)
real(sp), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_s1
function d_kv_s2(key,val) result(this)
character(len=*), intent(in) :: key
real(sp), intent(in), dimension(:,:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_s2
function d_kvp_s2(key, val) result(this)
character(len=*), intent(in) :: key
real(sp), intent(in), dimension(:,:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_s2
subroutine d_get_val_s2(val,this,key,success)
real(sp), intent(out), dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_s2
subroutine d_get_val_first_s2(val,this,success)
real(sp), intent(out), dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_s2
subroutine d_get_p_s2(val,this,key,success)
real(sp), pointer , dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_s2
subroutine d_get_p_first_s2(val,this,success)
real(sp), pointer , dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_s2
function d_kv_s3(key,val) result(this)
character(len=*), intent(in) :: key
real(sp), intent(in), dimension(:,:,:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_s3
function d_kvp_s3(key, val) result(this)
character(len=*), intent(in) :: key
real(sp), intent(in), dimension(:,:,:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_s3
subroutine d_get_val_s3(val,this,key,success)
real(sp), intent(out), dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_s3
subroutine d_get_val_first_s3(val,this,success)
real(sp), intent(out), dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_s3
subroutine d_get_p_s3(val,this,key,success)
real(sp), pointer , dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_s3
subroutine d_get_p_first_s3(val,this,success)
real(sp), pointer , dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_s3
function d_kv_d0(key,val) result(this)
character(len=*), intent(in) :: key
real(dp), intent(in) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_d0
function d_kvp_d0(key, val) result(this)
character(len=*), intent(in) :: key
real(dp), intent(in), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_d0
subroutine d_get_val_d0(val,this,key,success)
real(dp), intent(out) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_d0
subroutine d_get_val_first_d0(val,this,success)
real(dp), intent(out) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_d0
subroutine d_get_p_d0(val,this,key,success)
real(dp), pointer :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_d0
subroutine d_get_p_first_d0(val,this,success)
real(dp), pointer :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_d0
function d_kv_d1(key,val) result(this)
character(len=*), intent(in) :: key
real(dp), intent(in), dimension(:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_d1
function d_kvp_d1(key, val) result(this)
character(len=*), intent(in) :: key
real(dp), intent(in), dimension(:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_d1
subroutine d_get_val_d1(val,this,key,success)
real(dp), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_d1
subroutine d_get_val_first_d1(val,this,success)
real(dp), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_d1
subroutine d_get_p_d1(val,this,key,success)
real(dp), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_d1
subroutine d_get_p_first_d1(val,this,success)
real(dp), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_d1
function d_kv_d2(key,val) result(this)
character(len=*), intent(in) :: key
real(dp), intent(in), dimension(:,:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_d2
function d_kvp_d2(key, val) result(this)
character(len=*), intent(in) :: key
real(dp), intent(in), dimension(:,:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_d2
subroutine d_get_val_d2(val,this,key,success)
real(dp), intent(out), dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_d2
subroutine d_get_val_first_d2(val,this,success)
real(dp), intent(out), dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_d2
subroutine d_get_p_d2(val,this,key,success)
real(dp), pointer , dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_d2
subroutine d_get_p_first_d2(val,this,success)
real(dp), pointer , dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_d2
function d_kv_d3(key,val) result(this)
character(len=*), intent(in) :: key
real(dp), intent(in), dimension(:,:,:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_d3
function d_kvp_d3(key, val) result(this)
character(len=*), intent(in) :: key
real(dp), intent(in), dimension(:,:,:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_d3
subroutine d_get_val_d3(val,this,key,success)
real(dp), intent(out), dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_d3
subroutine d_get_val_first_d3(val,this,success)
real(dp), intent(out), dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_d3
subroutine d_get_p_d3(val,this,key,success)
real(dp), pointer , dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_d3
subroutine d_get_p_first_d3(val,this,success)
real(dp), pointer , dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_d3
function d_kv_c0(key,val) result(this)
character(len=*), intent(in) :: key
complex(sp), intent(in) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_c0
function d_kvp_c0(key, val) result(this)
character(len=*), intent(in) :: key
complex(sp), intent(in), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_c0
subroutine d_get_val_c0(val,this,key,success)
complex(sp), intent(out) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_c0
subroutine d_get_val_first_c0(val,this,success)
complex(sp), intent(out) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_c0
subroutine d_get_p_c0(val,this,key,success)
complex(sp), pointer :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_c0
subroutine d_get_p_first_c0(val,this,success)
complex(sp), pointer :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_c0
function d_kv_c1(key,val) result(this)
character(len=*), intent(in) :: key
complex(sp), intent(in), dimension(:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_c1
function d_kvp_c1(key, val) result(this)
character(len=*), intent(in) :: key
complex(sp), intent(in), dimension(:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_c1
subroutine d_get_val_c1(val,this,key,success)
complex(sp), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_c1
subroutine d_get_val_first_c1(val,this,success)
complex(sp), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_c1
subroutine d_get_p_c1(val,this,key,success)
complex(sp), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_c1
subroutine d_get_p_first_c1(val,this,success)
complex(sp), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_c1
function d_kv_c2(key,val) result(this)
character(len=*), intent(in) :: key
complex(sp), intent(in), dimension(:,:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_c2
function d_kvp_c2(key, val) result(this)
character(len=*), intent(in) :: key
complex(sp), intent(in), dimension(:,:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_c2
subroutine d_get_val_c2(val,this,key,success)
complex(sp), intent(out), dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_c2
subroutine d_get_val_first_c2(val,this,success)
complex(sp), intent(out), dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_c2
subroutine d_get_p_c2(val,this,key,success)
complex(sp), pointer , dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_c2
subroutine d_get_p_first_c2(val,this,success)
complex(sp), pointer , dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_c2
function d_kv_c3(key,val) result(this)
character(len=*), intent(in) :: key
complex(sp), intent(in), dimension(:,:,:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_c3
function d_kvp_c3(key, val) result(this)
character(len=*), intent(in) :: key
complex(sp), intent(in), dimension(:,:,:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_c3
subroutine d_get_val_c3(val,this,key,success)
complex(sp), intent(out), dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_c3
subroutine d_get_val_first_c3(val,this,success)
complex(sp), intent(out), dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_c3
subroutine d_get_p_c3(val,this,key,success)
complex(sp), pointer , dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_c3
subroutine d_get_p_first_c3(val,this,success)
complex(sp), pointer , dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_c3
function d_kv_z0(key,val) result(this)
character(len=*), intent(in) :: key
complex(dp), intent(in) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_z0
function d_kvp_z0(key, val) result(this)
character(len=*), intent(in) :: key
complex(dp), intent(in), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_z0
subroutine d_get_val_z0(val,this,key,success)
complex(dp), intent(out) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_z0
subroutine d_get_val_first_z0(val,this,success)
complex(dp), intent(out) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_z0
subroutine d_get_p_z0(val,this,key,success)
complex(dp), pointer :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_z0
subroutine d_get_p_first_z0(val,this,success)
complex(dp), pointer :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_z0
function d_kv_z1(key,val) result(this)
character(len=*), intent(in) :: key
complex(dp), intent(in), dimension(:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_z1
function d_kvp_z1(key, val) result(this)
character(len=*), intent(in) :: key
complex(dp), intent(in), dimension(:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_z1
subroutine d_get_val_z1(val,this,key,success)
complex(dp), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_z1
subroutine d_get_val_first_z1(val,this,success)
complex(dp), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_z1
subroutine d_get_p_z1(val,this,key,success)
complex(dp), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_z1
subroutine d_get_p_first_z1(val,this,success)
complex(dp), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_z1
function d_kv_z2(key,val) result(this)
character(len=*), intent(in) :: key
complex(dp), intent(in), dimension(:,:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_z2
function d_kvp_z2(key, val) result(this)
character(len=*), intent(in) :: key
complex(dp), intent(in), dimension(:,:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_z2
subroutine d_get_val_z2(val,this,key,success)
complex(dp), intent(out), dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_z2
subroutine d_get_val_first_z2(val,this,success)
complex(dp), intent(out), dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_z2
subroutine d_get_p_z2(val,this,key,success)
complex(dp), pointer , dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_z2
subroutine d_get_p_first_z2(val,this,success)
complex(dp), pointer , dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_z2
function d_kv_z3(key,val) result(this)
character(len=*), intent(in) :: key
complex(dp), intent(in), dimension(:,:,:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_z3
function d_kvp_z3(key, val) result(this)
character(len=*), intent(in) :: key
complex(dp), intent(in), dimension(:,:,:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_z3
subroutine d_get_val_z3(val,this,key,success)
complex(dp), intent(out), dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_z3
subroutine d_get_val_first_z3(val,this,success)
complex(dp), intent(out), dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_z3
subroutine d_get_p_z3(val,this,key,success)
complex(dp), pointer , dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_z3
subroutine d_get_p_first_z3(val,this,success)
complex(dp), pointer , dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_z3
function d_kv_b0(key,val) result(this)
character(len=*), intent(in) :: key
logical, intent(in) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_b0
function d_kvp_b0(key, val) result(this)
character(len=*), intent(in) :: key
logical, intent(in), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_b0
subroutine d_get_val_b0(val,this,key,success)
logical, intent(out) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_b0
subroutine d_get_val_first_b0(val,this,success)
logical, intent(out) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_b0
subroutine d_get_p_b0(val,this,key,success)
logical, pointer :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_b0
subroutine d_get_p_first_b0(val,this,success)
logical, pointer :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_b0
function d_kv_b1(key,val) result(this)
character(len=*), intent(in) :: key
logical, intent(in), dimension(:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_b1
function d_kvp_b1(key, val) result(this)
character(len=*), intent(in) :: key
logical, intent(in), dimension(:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_b1
subroutine d_get_val_b1(val,this,key,success)
logical, intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_b1
subroutine d_get_val_first_b1(val,this,success)
logical, intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_b1
subroutine d_get_p_b1(val,this,key,success)
logical, pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_b1
subroutine d_get_p_first_b1(val,this,success)
logical, pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_b1
function d_kv_b2(key,val) result(this)
character(len=*), intent(in) :: key
logical, intent(in), dimension(:,:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_b2
function d_kvp_b2(key, val) result(this)
character(len=*), intent(in) :: key
logical, intent(in), dimension(:,:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_b2
subroutine d_get_val_b2(val,this,key,success)
logical, intent(out), dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_b2
subroutine d_get_val_first_b2(val,this,success)
logical, intent(out), dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_b2
subroutine d_get_p_b2(val,this,key,success)
logical, pointer , dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_b2
subroutine d_get_p_first_b2(val,this,success)
logical, pointer , dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_b2
function d_kv_b3(key,val) result(this)
character(len=*), intent(in) :: key
logical, intent(in), dimension(:,:,:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_b3
function d_kvp_b3(key, val) result(this)
character(len=*), intent(in) :: key
logical, intent(in), dimension(:,:,:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_b3
subroutine d_get_val_b3(val,this,key,success)
logical, intent(out), dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_b3
subroutine d_get_val_first_b3(val,this,success)
logical, intent(out), dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_b3
subroutine d_get_p_b3(val,this,key,success)
logical, pointer , dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_b3
subroutine d_get_p_first_b3(val,this,success)
logical, pointer , dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_b3
function d_kv_h0(key,val) result(this)
character(len=*), intent(in) :: key
integer(ih), intent(in) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_h0
function d_kvp_h0(key, val) result(this)
character(len=*), intent(in) :: key
integer(ih), intent(in), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_h0
subroutine d_get_val_h0(val,this,key,success)
integer(ih), intent(out) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_h0
subroutine d_get_val_first_h0(val,this,success)
integer(ih), intent(out) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_h0
subroutine d_get_p_h0(val,this,key,success)
integer(ih), pointer :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_h0
subroutine d_get_p_first_h0(val,this,success)
integer(ih), pointer :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_h0
function d_kv_h1(key,val) result(this)
character(len=*), intent(in) :: key
integer(ih), intent(in), dimension(:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_h1
function d_kvp_h1(key, val) result(this)
character(len=*), intent(in) :: key
integer(ih), intent(in), dimension(:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_h1
subroutine d_get_val_h1(val,this,key,success)
integer(ih), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_h1
subroutine d_get_val_first_h1(val,this,success)
integer(ih), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_h1
subroutine d_get_p_h1(val,this,key,success)
integer(ih), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_h1
subroutine d_get_p_first_h1(val,this,success)
integer(ih), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_h1
function d_kv_h2(key,val) result(this)
character(len=*), intent(in) :: key
integer(ih), intent(in), dimension(:,:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_h2
function d_kvp_h2(key, val) result(this)
character(len=*), intent(in) :: key
integer(ih), intent(in), dimension(:,:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_h2
subroutine d_get_val_h2(val,this,key,success)
integer(ih), intent(out), dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_h2
subroutine d_get_val_first_h2(val,this,success)
integer(ih), intent(out), dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_h2
subroutine d_get_p_h2(val,this,key,success)
integer(ih), pointer , dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_h2
subroutine d_get_p_first_h2(val,this,success)
integer(ih), pointer , dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_h2
function d_kv_h3(key,val) result(this)
character(len=*), intent(in) :: key
integer(ih), intent(in), dimension(:,:,:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_h3
function d_kvp_h3(key, val) result(this)
character(len=*), intent(in) :: key
integer(ih), intent(in), dimension(:,:,:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_h3
subroutine d_get_val_h3(val,this,key,success)
integer(ih), intent(out), dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_h3
subroutine d_get_val_first_h3(val,this,success)
integer(ih), intent(out), dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_h3
subroutine d_get_p_h3(val,this,key,success)
integer(ih), pointer , dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_h3
subroutine d_get_p_first_h3(val,this,success)
integer(ih), pointer , dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_h3
function d_kv_i0(key,val) result(this)
character(len=*), intent(in) :: key
integer(is), intent(in) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_i0
function d_kvp_i0(key, val) result(this)
character(len=*), intent(in) :: key
integer(is), intent(in), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_i0
subroutine d_get_val_i0(val,this,key,success)
integer(is), intent(out) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_i0
subroutine d_get_val_first_i0(val,this,success)
integer(is), intent(out) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_i0
subroutine d_get_p_i0(val,this,key,success)
integer(is), pointer :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_i0
subroutine d_get_p_first_i0(val,this,success)
integer(is), pointer :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_i0
function d_kv_i1(key,val) result(this)
character(len=*), intent(in) :: key
integer(is), intent(in), dimension(:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_i1
function d_kvp_i1(key, val) result(this)
character(len=*), intent(in) :: key
integer(is), intent(in), dimension(:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_i1
subroutine d_get_val_i1(val,this,key,success)
integer(is), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_i1
subroutine d_get_val_first_i1(val,this,success)
integer(is), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_i1
subroutine d_get_p_i1(val,this,key,success)
integer(is), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_i1
subroutine d_get_p_first_i1(val,this,success)
integer(is), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_i1
function d_kv_i2(key,val) result(this)
character(len=*), intent(in) :: key
integer(is), intent(in), dimension(:,:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_i2
function d_kvp_i2(key, val) result(this)
character(len=*), intent(in) :: key
integer(is), intent(in), dimension(:,:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_i2
subroutine d_get_val_i2(val,this,key,success)
integer(is), intent(out), dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_i2
subroutine d_get_val_first_i2(val,this,success)
integer(is), intent(out), dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_i2
subroutine d_get_p_i2(val,this,key,success)
integer(is), pointer , dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_i2
subroutine d_get_p_first_i2(val,this,success)
integer(is), pointer , dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_i2
function d_kv_i3(key,val) result(this)
character(len=*), intent(in) :: key
integer(is), intent(in), dimension(:,:,:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_i3
function d_kvp_i3(key, val) result(this)
character(len=*), intent(in) :: key
integer(is), intent(in), dimension(:,:,:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_i3
subroutine d_get_val_i3(val,this,key,success)
integer(is), intent(out), dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_i3
subroutine d_get_val_first_i3(val,this,success)
integer(is), intent(out), dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_i3
subroutine d_get_p_i3(val,this,key,success)
integer(is), pointer , dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_i3
subroutine d_get_p_first_i3(val,this,success)
integer(is), pointer , dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_i3
function d_kv_l0(key,val) result(this)
character(len=*), intent(in) :: key
integer(il), intent(in) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_l0
function d_kvp_l0(key, val) result(this)
character(len=*), intent(in) :: key
integer(il), intent(in), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_l0
subroutine d_get_val_l0(val,this,key,success)
integer(il), intent(out) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_l0
subroutine d_get_val_first_l0(val,this,success)
integer(il), intent(out) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_l0
subroutine d_get_p_l0(val,this,key,success)
integer(il), pointer :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_l0
subroutine d_get_p_first_l0(val,this,success)
integer(il), pointer :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_l0
function d_kv_l1(key,val) result(this)
character(len=*), intent(in) :: key
integer(il), intent(in), dimension(:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_l1
function d_kvp_l1(key, val) result(this)
character(len=*), intent(in) :: key
integer(il), intent(in), dimension(:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_l1
subroutine d_get_val_l1(val,this,key,success)
integer(il), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_l1
subroutine d_get_val_first_l1(val,this,success)
integer(il), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_l1
subroutine d_get_p_l1(val,this,key,success)
integer(il), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_l1
subroutine d_get_p_first_l1(val,this,success)
integer(il), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_l1
function d_kv_l2(key,val) result(this)
character(len=*), intent(in) :: key
integer(il), intent(in), dimension(:,:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_l2
function d_kvp_l2(key, val) result(this)
character(len=*), intent(in) :: key
integer(il), intent(in), dimension(:,:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_l2
subroutine d_get_val_l2(val,this,key,success)
integer(il), intent(out), dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_l2
subroutine d_get_val_first_l2(val,this,success)
integer(il), intent(out), dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_l2
subroutine d_get_p_l2(val,this,key,success)
integer(il), pointer , dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_l2
subroutine d_get_p_first_l2(val,this,success)
integer(il), pointer , dimension(:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_l2
function d_kv_l3(key,val) result(this)
character(len=*), intent(in) :: key
integer(il), intent(in), dimension(:,:,:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_l3
function d_kvp_l3(key, val) result(this)
character(len=*), intent(in) :: key
integer(il), intent(in), dimension(:,:,:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_l3
subroutine d_get_val_l3(val,this,key,success)
integer(il), intent(out), dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_l3
subroutine d_get_val_first_l3(val,this,success)
integer(il), intent(out), dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_l3
subroutine d_get_p_l3(val,this,key,success)
integer(il), pointer , dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_l3
subroutine d_get_p_first_l3(val,this,success)
integer(il), pointer , dimension(:,:,:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_l3
function d_kv_cp0(key,val) result(this)
character(len=*), intent(in) :: key
type(c_ptr), intent(in) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_cp0
function d_kvp_cp0(key, val) result(this)
character(len=*), intent(in) :: key
type(c_ptr), intent(in), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_cp0
subroutine d_get_val_cp0(val,this,key,success)
type(c_ptr), intent(out) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_cp0
subroutine d_get_val_first_cp0(val,this,success)
type(c_ptr), intent(out) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_cp0
subroutine d_get_p_cp0(val,this,key,success)
type(c_ptr), pointer :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_cp0
subroutine d_get_p_first_cp0(val,this,success)
type(c_ptr), pointer :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_cp0
function d_kv_cp1(key,val) result(this)
character(len=*), intent(in) :: key
type(c_ptr), intent(in), dimension(:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_cp1
function d_kvp_cp1(key, val) result(this)
character(len=*), intent(in) :: key
type(c_ptr), intent(in), dimension(:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_cp1
subroutine d_get_val_cp1(val,this,key,success)
type(c_ptr), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_cp1
subroutine d_get_val_first_cp1(val,this,success)
type(c_ptr), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_cp1
subroutine d_get_p_cp1(val,this,key,success)
type(c_ptr), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_cp1
subroutine d_get_p_first_cp1(val,this,success)
type(c_ptr), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_cp1
function d_kv_fp0(key,val) result(this)
character(len=*), intent(in) :: key
type(c_funptr), intent(in) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_fp0
function d_kvp_fp0(key, val) result(this)
character(len=*), intent(in) :: key
type(c_funptr), intent(in), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_fp0
subroutine d_get_val_fp0(val,this,key,success)
type(c_funptr), intent(out) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_fp0
subroutine d_get_val_first_fp0(val,this,success)
type(c_funptr), intent(out) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_fp0
subroutine d_get_p_fp0(val,this,key,success)
type(c_funptr), pointer :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_fp0
subroutine d_get_p_first_fp0(val,this,success)
type(c_funptr), pointer :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_fp0
function d_kv_fp1(key,val) result(this)
character(len=*), intent(in) :: key
type(c_funptr), intent(in), dimension(:) :: val
type(dictionary_t) :: this
this = new_d_key(key)
call assign(this%first%value,val)
end function d_kv_fp1
function d_kvp_fp1(key, val) result(this)
character(len=*), intent(in) :: key
type(c_funptr), intent(in), dimension(:), target :: val
type(dictionary_t) :: this
this = new_d_key(key)
call associate(this%first%value,val)
end function d_kvp_fp1
subroutine d_get_val_fp1(val,this,key,success)
type(c_funptr), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call assign(val,v,success=success)
call nullify(v)
end subroutine d_get_val_fp1
subroutine d_get_val_first_fp1(val,this,success)
type(c_funptr), intent(out), dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call assign(val,this%first%value,success=success)
end subroutine d_get_val_first_fp1
subroutine d_get_p_fp1(val,this,key,success)
type(c_funptr), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
character(len=*), intent(in) :: key
logical, intent(out), optional :: success
type(variable_t) :: v
call associate(v,this,key=key)
call associate(val,v,success=success)
call nullify(v)
end subroutine d_get_p_fp1
subroutine d_get_p_first_fp1(val,this,success)
type(c_funptr), pointer , dimension(:) :: val
type(dictionary_t), intent(inout) :: this
logical, intent(out), optional :: success
call associate(val,this%first%value,success=success)
end subroutine d_get_p_first_fp1
  ! helper routines for often used stuff
  subroutine val_delete_request(val,dealloc)
    type(variable_t), intent(inout) :: val
    logical, intent(in), optional :: dealloc
    if ( present(dealloc) ) then
      if ( dealloc ) call delete(val)
    end if
    call nullify(val)
  end subroutine val_delete_request
  ! Create a routine for making the dictionary point to the data
  ! key.
  function d_kvp_dict(key,dic) result(this)
    character(len=*), intent(in) :: key
    type(dictionary_t), intent(in) :: dic
    type(dictionary_t) :: this
    type :: pdictionary_entry_t
      type(dictionary_entry_t), pointer :: d => null()
    end type pdictionary_entry_t
    type(pdictionary_entry_t) :: pd
    type(variable_t) :: v
    character(len=1) :: c(1)
    pd%d => dic%first
    call associate_type(v,transfer(pd,c), which="dict")
    this = (key.kvp.v)
    call nullify(v)
  end function d_kvp_dict
  ! In case the value of the dictionary is a dictionary we can request that
  ! dictionary directly
  subroutine d_key2dict(dic,d,key,dealloc)
    type(dictionary_t), intent(inout) :: dic
    type(dictionary_t), intent(inout) :: d
    character(len=*), intent(in), optional :: key
    logical, intent(in), optional :: dealloc
    ! Retrieving a dictionary will NEVER
    ! be copying the entire dictionary.
    call d_get_p_dict(dic,d,key=key,dealloc=dealloc)
  end subroutine d_key2dict
  subroutine d_get_p_dict(dic,d,key,dealloc)
    type(dictionary_t), intent(inout) :: dic
    type(dictionary_t), intent(inout) :: d
    character(len=*), intent(in), optional :: key
    logical, intent(in), optional :: dealloc
    ! Instead of saving the data-type dict
    ! we save the first pointer.
    ! This will allow greater flexibility as the
    ! parent container can then be re-used with out
    ! worries.
    ! I.e.
    ! if one uses :
    ! type :: pdict
    ! type(dictionary_t), pointer :: d
    ! end type
    ! then the address of the "parenting" dictionary is saved,
    ! And hence, doing:
    ! dic1 = ('a'.kv.1)
    ! dic2 = ('dic1'.kvp.dic1)
    ! call nullify(dic1)
    ! dic1 = ('b'.kv.1)
    ! will make dic1 in dic2 contain ('b'.kv.1)
    ! Specifically because the address of the dic1 does not change.
    ! However, the dictionary_entry_t pointer is irrespective of parent locality.
    type :: pdictionary_entry_t
      type(dictionary_entry_t), pointer :: d => null()
    end type pdictionary_entry_t
    type(pdictionary_entry_t) :: pd
    type(dictionary_t) :: ld
    type(variable_t) :: v
    character(len=1), allocatable :: c(:)
    integer :: i
    logical :: ldealloc
    ldealloc = .false.
    if ( present(dealloc) ) ldealloc = dealloc
    if ( ldealloc ) then
      call delete(dic)
    else
      call nullify(dic)
    end if
    ! Retrieve the dictionary key
    call associate(v,d,key=key)
    if ( v%t .eq. '    ' ) then
      call nullify(v)
      return
    end if
    i = size_enc(v)
    allocate(c(i))
    call enc(v,c)
    pd = transfer(c,pd)
    deallocate(c)
    dic%first => pd%d
    call nullify(v)
    ! we need to re-count the number of entries in
    ! the dictionary_entry_t tree.
    ! Sadly, this is because we contain the dictionary_entry_t
    ! type, and NOT the dict type :(
    ! However, it makes the programming style more
    ! intuitive (dependent on how you look at it)
    ld = .first. dic
    dic%len = 0
    do while ( .not. (.empty. ld) )
      dic%len = dic%len + 1
      ld = .next. ld
    end do
  end subroutine d_get_p_dict
end module dictionary
