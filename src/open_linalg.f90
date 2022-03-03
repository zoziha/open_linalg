module open_linalg_m

    use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64, stderr => error_unit
    implicit none
    private

    public :: cross_product, det, inv, matmul, operator(.i.), operator(.x.), solve

    interface operator(.x.)
        module procedure matmul_sp, matmul_dp
    end interface

    interface matmul
        module procedure matmul_sp, matmul_dp
    end interface matmul

    interface det
        module procedure det_sp, det_dp
    end interface det

    interface operator(.i.)
        module procedure inv_sp, inv_dp
    end interface

    interface inv
        module procedure inv_sp, inv_dp
    end interface inv

    interface cross_product
        module procedure cross_product_sp, cross_product_dp
    end interface cross_product

    interface solve
        module procedure solve_sp, solve_dp
    end interface solve

contains

    pure function cross_product_sp(a, b) result(c)
    !! calculate single precision cross product
        real(sp), intent(in) :: a(3), b(3)
        real(sp) c(3)

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)

    end function cross_product_sp

    pure function cross_product_dp(a, b) result(c)
    !! calculate double precision cross product
        real(dp), intent(in) :: a(3), b(3)
        real(dp) c(3)

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)

    end function cross_product_dp

    function inv_sp(a) result(b)
    !! calculate inverse of a single precision matrix
        real(sp), intent(in) :: a(:, :)
        real(sp) b(size(a, 2), size(a, 1))
        integer ipiv(size(a, 1)), info
        real(sp) work(size(a, 2))

        b = a
        ! http://www.netlib.org/lapack/explore-html/d8/ddc/group__real_g_ecomputational_ga8d99c11b94db3d5eac75cac46a0f2e17.html
        call sgetrf(size(a, 1), size(a, 2), b, size(a, 1), ipiv, info)
        if (info < 0) then
            write (stderr, *) 'sgetrf: illegal value in argument ', info
            error stop
        else if (info > 0) then
            write (stderr, *) 'sgetrf: singular matrix, U(i,i) is exactly zero, info = ', info
            error stop
        end if

        ! http://www.netlib.org/lapack/explore-html/d8/ddc/group__real_g_ecomputational_ga1af62182327d0be67b1717db399d7d83.html
        call sgetri(size(a, 2), b, size(a, 1), ipiv, work, size(a, 2), info)
        if (info < 0) then
            write (stderr, *) 'sgetri: illegal value in argument ', info
            error stop
        else if (info > 0) then
            write (stderr, *) 'sgetri: singular matrix, U(i,i) is exactly zero, info = ', info
            error stop
        end if

    end function inv_sp

    function inv_dp(a) result(b)
    !! calculate inverse of a double precision matrix
        real(dp), intent(in) :: a(:, :)
        real(dp) b(size(a, 2), size(a, 1))
        integer ipiv(size(a, 1)), info
        real(dp) work(size(a, 2))

        b = a
        ! http://www.netlib.org/lapack/explore-html/d8/ddc/group__real_g_ecomputational_ga8d99c11b94db3d5eac75cac46a0f2e17.html
        call dgetrf(size(a, 1), size(a, 2), b, size(a, 1), ipiv, info)
        if (info < 0) then
            write (stderr, *) 'dgetrf: illegal value in argument ', info
            error stop
        else if (info > 0) then
            write (stderr, *) 'dgetrf: singular matrix, U(i,i) is exactly zero, info = ', info
            error stop
        end if

        ! http://www.netlib.org/lapack/explore-html/d8/ddc/group__real_g_ecomputational_ga1af62182327d0be67b1717db399d7d83.html
        call dgetri(size(a, 2), b, size(a, 1), ipiv, work, size(a, 2), info)
        if (info < 0) then
            write (stderr, *) 'dgetri: illegal value in argument ', info
            error stop
        else if (info > 0) then
            write (stderr, *) 'dgetri: singular matrix, U(i,i) is exactly zero, info = ', info
            error stop
        end if

    end function inv_dp

    function solve_sp(a, b) result(x)
    !! solve linear system of single precision
        real(sp), intent(in) :: a(:, :), b(:, :)
        real(sp) x(size(b, 1), size(b, 2))

        real(sp) a_(size(a, 1), size(a, 2))
        integer ipiv(size(a, 1))
        integer info

        a_ = a; x = b
        ! http://www.netlib.org/lapack/explore-html/d0/db8/group__real_g_esolve_ga3b05fb3999b3d7351cb3101a1fd28e78.html
        call sgesv(size(a, 1), size(b, 2), a_, size(a, 1), ipiv, x, size(b, 1), info)

        if (info < 0) then
            write (stderr, *) 'sgesv: illegal value in argument ', info
            error stop
        else if (info > 0) then
            write (stderr, *) 'sgesv: singular matrix, U(i,i) is exactly zero, info = ', info
            error stop
        end if

    end function solve_sp

    function solve_dp(a, b) result(x)
    !! solve linear system of double precision
        real(dp), intent(in) :: a(:, :), b(:, :)
        real(dp) x(size(b, 1), size(b, 2))

        real(dp) a_(size(a, 1), size(a, 2))
        integer ipiv(size(a, 1))
        integer info

        a_ = a; x = b
        ! http://www.netlib.org/lapack/explore-html/d7/d3b/group__double_g_esolve_ga5ee879032a8365897c3ba91e3dc8d512.html
        call dgesv(size(a, 1), size(b, 2), a_, size(a, 1), ipiv, x, size(b, 1), info)

        if (info < 0) then
            write (stderr, *) 'dgesv: illegal value in argument ', info
            error stop
        else if (info > 0) then
            write (stderr, *) 'dgesv: singular matrix, U(i,i) is exactly zero, info = ', info
            error stop
        end if

    end function solve_dp

    function det_sp(a) result(d)
    !! calculate determinant of a single precision matrix
        real(sp), intent(in) :: a(:, :)
        real(sp) d
        real(sp) a_(size(a, 1), size(a, 2))
        integer ipiv(size(a, 1)), info, i

        a_ = a
        ! http://www.netlib.org/lapack/explore-html/d8/ddc/group__real_g_ecomputational_ga8d99c11b94db3d5eac75cac46a0f2e17.html
        call sgetrf(size(a, 1), size(a, 2), a_, size(a, 1), ipiv, info)
        if (info < 0) then
            write (stderr, *) 'sgetrf: illegal value in argument ', info
            error stop
        else if (info > 0) then
            write (stderr, *) 'sgetrf: singular matrix, U(i,i) is exactly zero, info = ', info
            error stop
        end if

        d = 1.0_sp
        do i = 1, size(a, 2)
            if (ipiv(i) /= i) then
                d = -d*a_(i, i)
            else
                d = d*a_(i, i)
            end if
        end do

    end function det_sp

    function det_dp(a) result(d)
    !! calculate determinant of a double precision matrix
        real(dp), intent(in) :: a(:, :)
        real(dp) d
        real(dp) a_(size(a, 1), size(a, 2))
        integer ipiv(size(a, 1)), info, i

        a_ = a
        ! http://www.netlib.org/lapack/explore-html/d8/ddc/group__real_g_ecomputational_ga8d99c11b94db3d5eac75cac46a0f2e17.html
        call dgetrf(size(a, 1), size(a, 2), a_, size(a, 1), ipiv, info)
        if (info < 0) then
            write (stderr, *) 'dgetrf: illegal value in argument ', info
            error stop
        else if (info > 0) then
            write (stderr, *) 'dgetrf: singular matrix, U(i,i) is exactly zero, info = ', info
            error stop
        end if

        d = 1.0_dp
        do i = 1, size(a, 2)
            if (ipiv(i) /= i) then
                d = -d*a_(i, i)
            else
                d = d*a_(i, i)
            end if
        end do

    end function det_dp

    function matmul_sp(a, b) result(c)
    !! matrix multiplication of single precision matrices
        real(sp), intent(in) :: a(:, :), b(:, :)
        real(sp) c(size(a, 1), size(b, 2))

        ! http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
        call sgemm( &
            'n', 'n', &
            size(a, 1), size(b, 2), size(a, 2), &
            1.0_sp, a, size(a, 1), b, size(b, 1), &
            0.0_sp, c, size(a, 1) &
            )

    end function matmul_sp

    function matmul_dp(a, b) result(c)
    !! matrix multiplication of double precision matrices
        real(dp), intent(in) :: a(:, :), b(:, :)
        real(dp) c(size(a, 1), size(b, 2))

        ! http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
        call dgemm( &
            'n', 'n', &
            size(a, 1), size(b, 2), size(a, 2), &
            1.0_dp, a, size(a, 1), b, size(b, 1), &
            0.0_dp, c, size(a, 1) &
            )

    end function matmul_dp

end module open_linalg_m
