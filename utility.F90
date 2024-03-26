subroutine locateval(xx,n,x,ii)
!
 use mod_post, only: mykind
 implicit none
!
 integer, intent(in) :: n
 integer, intent(out) :: ii
 real(mykind), dimension(1:n), intent(in) :: xx
 real(mykind), intent(in) :: x
!
 if (x>xx(n)) then
  ii = n
 else
  ii = 0
  do while (xx(ii+1)<x)
   ii = ii+1
  enddo
 endif
!
end subroutine locateval


