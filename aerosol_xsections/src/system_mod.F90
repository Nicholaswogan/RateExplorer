module system_mod

implicit none
public

! set numerical precision
integer,parameter :: f = selected_real_kind(12) ! 8 byte real
!integer,parameter :: f = selected_real_kind(8) ! 8 byte real

! Return Codes
integer, public, parameter :: RC_OK             = 0   !! Success
integer, public, parameter :: RC_ERROR          = -1  !! Failure
integer, public, parameter :: RC_WARNING        = 1   !! Warning
integer, public, parameter :: RC_WARNING_RETRY  = 2   !! Warning, Retry Suggested



end module