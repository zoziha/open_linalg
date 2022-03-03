program main

    use open_linalg_m, only: cross => cross_product
    real :: x(3) = [1, 0, 0], y(3) = [0, 1, 0]
    print *, cross(x, y)

end program main

!    0.00000000       0.00000000       1.00000000
