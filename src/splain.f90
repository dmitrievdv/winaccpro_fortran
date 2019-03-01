module splain
    use mymod


    implicit none
    integer :: spl_n
    contains

        subroutine fivepoints(BQBt, D, S) !������������ ��������
            real(mp), dimension (5, spl_n) :: BQBt
            real(mp), dimension (spl_n) :: D, alpha, beta, a
            real(mp), dimension (spl_n) :: S
            real(mp), dimension (0:spl_n-1) :: b
            real(mp), dimension (-1:spl_n-2) :: c
            real(mp), dimension (-1:spl_n) :: p, q
            real(mp), dimension (0:spl_n) :: r
            integer :: i

            forall(i = 1:spl_n)
                c(i-2) = BQBt(1, i)
                b(i-1) = BQBt(2, i)
                a(i) = BQBt(3, i)
            end forall

            p = 0
            q = 0
            r = 0
            alpha = 0
            beta = 0


            do i = 1, spl_n, 1
                beta(i) = b(i-1) - p(i-2)*c(i-2)
                alpha(i) = a(i) - p(i-1)*beta(i) - q(i-2)*c(i-2)
                r(i) = (D(i) - r(i-1)*beta(i) - r(i-2)*c(i-2))/alpha(i)
                p(i) = (b(i) - q(i-1)*beta(i))/alpha(i)
                q(i) = c(i)/alpha(i)
            enddo

            S(spl_n) = r(spl_n)
            S(spl_n-1) = r(spl_n-1) - p(spl_n-1)*S(spl_n)

            do i = spl_n-2, 1, -1
                S(i) = r(i) - p(i)*S(i+1) - q(i)*S(i+2)
            enddo

        end subroutine fivepoints


        subroutine init_mat(dat, A, B) ! ������������� ������ A � �
            real(mp), dimension (2, 0:spl_n+1) :: dat
            real(mp), dimension (3, spl_n) :: A, B
            integer :: i

            A = 0; B = 0

            forall(i = 1:spl_n)
                A(1, i) = dat(1, i) - dat(1, i-1)
                A(3, i) = dat(1, i+1) - dat(1, i)
                A(2, i) = 2*(dat(1, i+1) - dat(1, i-1))
                B(1, i) = 1/(dat(1, i) - dat(1, i-1))
                B(3, i) = 1/(dat(1, i+1) - dat(1, i))
                B(2, i) = -(B(1, i) + B(3, i))
            end forall

            A(3, spl_n-1) = 0
            A(1, 2) = 0
            A(3, 1) = 0
            A(1, spl_n) = 0

            B(:, 1) = 0
            B(:, spl_n) = 0

        end subroutine init_mat

        real(mp) function splain_value(x, xn, spl_cf)
            real(mp), dimension (0:3, spl_n-1) :: spl_cf
            real(mp), dimension (spl_n) :: xn
            real(mp) :: x

            integer, dimension(1) :: loc
            integer :: k
            real(mp) :: t

            loc = minloc(abs(x - xn))
            k = loc(1)

            if(x < xn(k) .and. k /= 1) then; k = k-1; endif
            if(k == spl_n) then; k = k-1; endif

            t = (x - xn(k))/(xn(k+1) - xn(k))

            splain_value = spl_cf(0, k) + t*spl_cf(1, k) + t**2*spl_cf(2, k) + t**3*spl_cf(3, k)

        end function splain_value

        subroutine get_splain(dat, spl_cf)

            real(mp), dimension (2, 0:spl_n+1) :: dat
            real(mp), dimension (0:3, spl_n-1) :: spl_cf
            real(mp), dimension (3, spl_n) :: A1, B
            real(mp), dimension (5, spl_n) :: A
            real(mp), dimension (-1:spl_n+2) :: S
            real(mp), dimension (spl_n) :: R, D
            real(mp) :: h, x1, x2
            integer :: i, j, k
            S = 0

            call init_mat(dat, A1, B)

            forall(i = 1:spl_n, j = 2:4) A(j, i) = A1(j-1, i)

            forall(i = 1:spl_n) D(i) = 6*(dot_product(B(:, i), dat(2, i-1:i+1)))

            call fivepoints(A, D, S(1:spl_n))

            forall(i = 1:spl_n) R(i) = dat(2, i)

            !call make_unif_net(xmin, xmax, m, spl(1, :))


            do k = 1, spl_n-1
                x2 = dat(1, k+1)
                x1 = dat(1, k)
                h = x2-x1
                spl_cf(0, k) = R(k)
                spl_cf(1, k) = (R(k+1) - R(k) - (h**2)/6*(2*S(k) + S(k+1)))!/h
                spl_cf(2, k) = S(k)/2*h**2
                spl_cf(3, k) = (S(k+1) - S(k))/6*h**2!/h
            enddo


        end subroutine get_splain

        subroutine init_derivates(x, p, m)
            real(mp), dimension(spl_n) :: x, p, m

            integer :: k
            real(mp) :: c1, c2, c3, dp0, dp1, dx0, dx1, c

            do k=2, spl_n-1
                c = 0.9
                dp0 = p(k) - p(k-1)
                dp1 = p(k+1) - p(k)
                dx0 = x(k) - x(k-1)
                dx1 = x(k+1) - x(k)
                m(k) = tan((atan(dp0/dx0) + atan(dp1/dx1))/2)
                !m(k) = (dp0/dx0 + dp1/dx1)/2
                !m(k) = (1-c)*(dp1 + dp0)/(dx1+dx0)
            enddo

            m(1) = (p(2) - p(1))/(x(2) - x(1))
            m(spl_n) = (p(spl_n) - p(spl_n-1))/(x(spl_n) - x(spl_n-1))

            !c1 = p(1)/(x(1)-x(2))/(x(1) - x(3))
            !c2 = p(2)/(x(2)-x(1))/(x(2) - x(3))
            !c3 = p(3)/(x(3)-x(1))/(x(3) - x(2))

            !m(1) = c1*(2*x(1)-x(2)-x(3)) + c2*(x(1)-x(3)) + c3*(x(1)-x(2))

            !c1 = p(spl_n)/(x(spl_n)-x(spl_n-1))/(x(spl_n) - x(spl_n-2))
            !c2 = p(spl_n-1)/(x(spl_n-1)-x(spl_n))/(x(spl_n-1) - x(spl_n-2))
            !c3 = p(spl_n-2)/(x(spl_n-2)-x(spl_n))/(x(spl_n-2) - x(spl_n-1))

            !m(spl_n) = c1*(2*x(spl_n)-x(spl_n-1)-x(spl_n-2)) + c2*(x(spl_n)-x(spl_n-2)) + c3*(x(spl_n)-x(1))

        end subroutine init_derivates

        subroutine get_hermite(x, p, spl_cf)
            real(mp), dimension(spl_n) :: x, p
            real(mp), dimension(0:3, spl_n-1) :: spl_cf

            real(mp), dimension(spl_n) :: m
            real(mp), dimension(0:3) :: h00, h01, h10, h11
            real(mp) :: h

            integer :: k

            h00(0:) = (/1.0_mp, 0.0_mp, -3.0_mp, 2.0_mp/)
            h10(0:) = (/0.0_mp, 1.0_mp, -2.0_mp, 1.0_mp/)
            h01(0:) = (/0.0_mp, 0.0_mp, 3.0_mp, -2.0_mp/)
            h11(0:) = (/0.0_mp, 0.0_mp, -1.0_mp, 1.0_mp/)

            call init_derivates(x, p, m)

            do k=1, spl_n-1
                h = x(k+1) - x(k)
                spl_cf(:, k) = p(k)*h00 + h*m(k)*h10 + p(k+1)*h01 + h*m(k+1)*h11
            enddo

        end subroutine get_hermite
end module
