module poly
    use mpmod
    implicit none

    contains


        subroutine left_shift(P)
            real(mp), dimension(0:) :: P
            integer :: n, i

            n = size(P) - 1

            do i = n, 1, -1
                P(i) = P(i-1)
            enddo

            P(0) = 0.0_mp
        end subroutine left_shift



        subroutine right_shift(P)
            real(mp), dimension(0:) :: P
            integer :: n, i

            n = size(P) - 1

            do i = 0, n-1
                P(i) = P(i+1)
            enddo

            P(n) = 0.0_mp
        end subroutine right_shift

        subroutine bernulli(x, Pn, n)
            real(mp), dimension (0:) :: Pn
            real(mp), dimension (0:n) :: y
            real(mp) :: x, xp, eps
            integer :: n, i, k
            eps = 1e-8_mp
            x = 0.5_mp
            xp = 0
            forall(i = 0:n)
                y(i) = i
            end forall


            do while (abs(x-xp) > eps)
                !y = 2*y/y(n)
                xp = x
                call right_shift(y)
                y(n) = -dot_product(Pn(0:n), y(0:n))/Pn(n)
                x = y(n)/y(n-1)
            enddo
        end subroutine bernulli


        subroutine gorner(P, x, n)
            real(mp), dimension (0:) :: P
            real(mp) :: x
            integer :: n, i

            do i = n-1, 0, -1
                P(i) = P(i+1)*x + P(i)
            enddo

            call right_shift(P)

        end subroutine gorner
end module poly
