program main

    use open_linalg_m, only: solve

    real :: a(2, 2) = reshape([1.0, 3.0, 2.0, 4.0], [2, 2]), &
            b(2, 1) = reshape([7.0, 15.0], [2, 1])

    b = solve(a, b)
    print 100, b

100 format(*(f6.2))

end program main

!   1.00  3.00
