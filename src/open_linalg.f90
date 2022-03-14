module open_linalg_m

    use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64, stderr => error_unit
    implicit none
    private

    public :: det, inv, matmul, operator(.i.), operator(.x.), solve

    interface operator(.x.)
        module procedure matmul_sp, matmul_dp, cmatmul_sp, cmatmul_dp, &
            rcmatmul_sp, rcmatmul_dp, crmatmul_sp, crmatmul_dp
    end interface

    interface matmul
        module procedure matmul_sp, matmul_dp, cmatmul_sp, cmatmul_dp, &
            rcmatmul_sp, rcmatmul_dp, crmatmul_sp, crmatmul_dp
    end interface matmul

    interface det
        module procedure det_sp, det_dp
    end interface det

    interface operator(.i.)
        module procedure inv_sp, inv_dp, cinv_sp, cinv_dp
    end interface

    interface inv
        module procedure inv_sp, inv_dp, cinv_sp, cinv_dp
    end interface inv

    interface solve
        module procedure solve_sp, solve_dp, csolve_sp, csolve_dp
    end interface solve

contains

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

    function cinv_sp(a) result(b)
    !! calculate inverse of a single precision matrix
        complex(sp), intent(in) :: a(:, :)
        complex(sp) b(size(a, 2), size(a, 1))
        integer ipiv(size(a, 1)), info
        complex(sp) work(size(a, 2))

        b = a
        ! http://www.netlib.org/lapack/explore-html/d8/ddc/group__real_g_ecomputational_ga8d99c11b94db3d5eac75cac46a0f2e17.html
        call cgetrf(size(a, 1), size(a, 2), b, size(a, 1), ipiv, info)
        if (info < 0) then
            write (stderr, *) 'cgetrf: illegal value in argument ', info
            error stop
        else if (info > 0) then
            write (stderr, *) 'cgetrf: singular matrix, U(i,i) is exactly zero, info = ', info
            error stop
        end if

        ! http://www.netlib.org/lapack/explore-html/d8/ddc/group__real_g_ecomputational_ga1af62182327d0be67b1717db399d7d83.html
        call cgetri(size(a, 2), b, size(a, 1), ipiv, work, size(a, 2), info)
        if (info < 0) then
            write (stderr, *) 'cgetri: illegal value in argument ', info
            error stop
        else if (info > 0) then
            write (stderr, *) 'cgetri: singular matrix, U(i,i) is exactly zero, info = ', info
            error stop
        end if

    end function cinv_sp

    function cinv_dp(a) result(b)
    !! calculate inverse of a double precision matrix
        complex(dp), intent(in) :: a(:, :)
        complex(dp) b(size(a, 2), size(a, 1))
        integer ipiv(size(a, 1)), info
        complex(dp) work(size(a, 2))

        b = a
        ! http://www.netlib.org/lapack/explore-html/d8/ddc/group__real_g_ecomputational_ga8d99c11b94db3d5eac75cac46a0f2e17.html
        call zgetrf(size(a, 1), size(a, 2), b, size(a, 1), ipiv, info)
        if (info < 0) then
            write (stderr, *) 'zgetrf: illegal value in argument ', info
            error stop
        else if (info > 0) then
            write (stderr, *) 'zgetrf: singular matrix, U(i,i) is exactly zero, info = ', info
            error stop
        end if

        ! http://www.netlib.org/lapack/explore-html/d8/ddc/group__real_g_ecomputational_ga1af62182327d0be67b1717db399d7d83.html
        call zgetri(size(a, 2), b, size(a, 1), ipiv, work, size(a, 2), info)
        if (info < 0) then
            write (stderr, *) 'zgetri: illegal value in argument ', info
            error stop
        else if (info > 0) then
            write (stderr, *) 'zgetri: singular matrix, U(i,i) is exactly zero, info = ', info
            error stop
        end if

    end function cinv_dp

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
    
    function csolve_sp(a, b) result(x)
    !! solve linear system of single precision
        complex(sp), intent(in) :: a(:, :), b(:, :)
        complex(sp) x(size(b, 1), size(b, 2))

        complex(sp) a_(size(a, 1), size(a, 2))
        integer ipiv(size(a, 1))
        integer info

        a_ = a; x = b
        ! http://www.netlib.org/lapack/explore-html/d0/db8/group__real_g_esolve_ga3b05fb3999b3d7351cb3101a1fd28e78.html
        call cgesv(size(a, 1), size(b, 2), a_, size(a, 1), ipiv, x, size(b, 1), info)

        if (info < 0) then
            write (stderr, *) 'cgesv: illegal value in argument ', info
            error stop
        else if (info > 0) then
            write (stderr, *) 'cgesv: singular matrix, U(i,i) is exactly zero, info = ', info
            error stop
        end if

    end function csolve_sp
    
    function csolve_dp(a, b) result(x)
    !! solve linear system of double precision
        complex(dp), intent(in) :: a(:, :), b(:, :)
        complex(dp) x(size(b, 1), size(b, 2))

        complex(dp) a_(size(a, 1), size(a, 2))
        integer ipiv(size(a, 1))
        integer info

        a_ = a; x = b
        ! http://www.netlib.org/lapack/explore-html/d7/d3b/group__double_g_esolve_ga5ee879032a8365897c3ba91e3dc8d512.html
        call zgesv(size(a, 1), size(b, 2), a_, size(a, 1), ipiv, x, size(b, 1), info)

        if (info < 0) then
            write (stderr, *) 'zgesv: illegal value in argument ', info
            error stop
        else if (info > 0) then
            write (stderr, *) 'zgesv: singular matrix, U(i,i) is exactly zero, info = ', info
            error stop
        end if

    end function csolve_dp

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
        integer m, n, k
        real(sp) c(size(a, 1), size(b, 2))

        m = size(a, 1)
        n = size(b, 2)
        k = size(a, 2)
        ! http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
        call sgemm( &
            'N', 'N', &
            m, n, k, &
            1.0_sp, a, m, b, k, &
            0.0_sp, c, m)

    end function matmul_sp

    function matmul_dp(a, b) result(c)
    !! matrix multiplication of double precision matrices
        real(dp), intent(in) :: a(:, :), b(:, :)
        real(dp) c(size(a, 1), size(b, 2))
        integer m, n, k

        m = size(a, 1)
        n = size(b, 2)
        k = size(a, 2)
        ! http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
        call dgemm( &
            'N', 'N', &
            m, n, k, &
            1.0_dp, a, m, b, k, &
            0.0_dp, c, m)

    end function matmul_dp

    function cmatmul_sp(a, b) result(c)
    !! matrix multiplication of single precision complex matrices
        complex(sp), intent(in) :: a(:, :), b(:, :)
        complex(sp) c(size(a, 1), size(b, 2))
        integer m, n, k

        m = size(a, 1)
        n = size(b, 2)
        k = size(a, 2)
        ! http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
        call cgemm( &
            'N', 'N', &
            m, n, k, &
            (1.0_sp, 0.0_sp), a, m, b, k, &
            (0.0_sp, 0.0_sp), c, m)

    end function cmatmul_sp

    function cmatmul_dp(a, b) result(c)
    !! matrix multiplication of double precision complex matrices
        complex(dp), intent(in) :: a(:, :), b(:, :)
        complex(dp) c(size(a, 1), size(b, 2))
        integer m, n, k

        m = size(a, 1)
        n = size(b, 2)
        k = size(a, 2)
        ! http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
        call zgemm( &
            'N', 'N', &
            m, n, k, &
            (1.0_dp, 0.0_dp), a, m, b, k, &
            (0.0_dp, 0.0_dp), c, m)

    end function cmatmul_dp

    function rcmatmul_sp(a, b) result(c)
    !! matrix multiplication of single precision complex matrices
        real(sp), intent(in) :: a(:, :)
        complex(sp), intent(in) :: b(:, :)
        complex(sp) c(size(a, 1), size(b, 2))
        integer m, n, k

        m = size(a, 1)
        n = size(b, 2)
        k = size(a, 2)
        ! http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
        call cgemm( &
            'N', 'N', &
            m, n, k, &
            (1.0_sp, 0.0_sp), cmplx(a, kind=sp), m, b, k, &
            (0.0_sp, 0.0_sp), c, m)

    end function rcmatmul_sp

    function rcmatmul_dp(a, b) result(c)
    !! matrix multiplication of double precision complex matrices
        real(dp), intent(in) :: a(:, :)
        complex(dp), intent(in) :: b(:, :)
        complex(dp) c(size(a, 1), size(b, 2))
        integer m, n, k

        m = size(a, 1)
        n = size(b, 2)
        k = size(a, 2)
        ! http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
        call zgemm( &
            'N', 'N', &
            m, n, k, &
            (1.0_dp, 0.0_dp), cmplx(a, kind=dp), m, b, k, &
            (0.0_dp, 0.0_dp), c, m)

    end function rcmatmul_dp

    function crmatmul_sp(a, b) result(c)
    !! matrix multiplication of single precision complex matrices
        complex(sp), intent(in) :: a(:, :)
        real(sp), intent(in) :: b(:, :)
        complex(sp) c(size(a, 1), size(b, 2))
        integer m, n, k

        m = size(a, 1)
        n = size(b, 2)
        k = size(a, 2)
        ! http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
        call cgemm( &
            'N', 'N', &
            m, n, k, &
            (1.0_sp, 0.0_sp), a, m, cmplx(b, kind=sp), k, &
            (0.0_sp, 0.0_sp), c, m)

    end function crmatmul_sp

    function crmatmul_dp(a, b) result(c)
    !! matrix multiplication of double precision complex matrices
        complex(dp), intent(in) :: a(:, :)
        real(dp), intent(in) :: b(:, :)
        complex(dp) c(size(a, 1), size(b, 2))
        integer m, n, k

        m = size(a, 1)
        n = size(b, 2)
        k = size(a, 2)
        ! http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_gaeda3cbd99c8fb834a60a6412878226e1.html
        call zgemm( &
            'N', 'N', &
            m, n, k, &
            (1.0_dp, 0.0_dp), a, m, cmplx(b, kind=dp), k, &
            (0.0_dp, 0.0_dp), c, m)

    end function crmatmul_dp

end module open_linalg_m
