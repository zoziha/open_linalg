program main

    use open_linalg_m, only: operator(.x.), matmul
    real :: x(2, 1) = reshape([1, 2], [2, 1]), &
            y(1, 2) = reshape([3, 4], [1, 2])
    real :: z(2, 2)

    z = matmul(x, y)
    do i = 1, 2
        print 100, z(i, :)
    end do

100 format(*(f6.2))

end program main

!@note: Fortran array column first
!  3.00  4.00
!  6.00  8.00
