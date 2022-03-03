program main

    use open_linalg_m, only: det
    real :: x(2, 2) = reshape([1, 2, 3, 4], [2, 2])

    print *, det(x)

end program main

!@note: Fortran array column first
! -2.00000000