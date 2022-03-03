program main

    use open_linalg_m, only: inv, operator(.i.)

    real :: x(2, 2) = reshape([1, 2, 3, 4], [2, 2])

    x = .i.x

    do i = 1, 2
        print 100, x(i, :)
    end do
    
100 format(*(f6.2))

end program main

!@note: Fortran array column first
! -2.00  1.50
!  1.00 -0.50
  