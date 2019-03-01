module modmod
    use mymod
    use mpmod
    implicit none

    contains

        subroutine lagran(xn, yn, inter, a, b)
            real(mp), dimension(:) :: xn, yn
            real(mp), dimension(:) :: inter
            real(mp) :: a, b
            integer :: n, n_int

            real(mp) :: step, x
            real(mp), dimension(size(xn)) :: l, t1, t2
            integer :: i, j

            n = size(xn)
            n_int = size(inter)

            step = (b-a)/n_int

            do i=1, n_int
                x = a + step*i - step/2.0_mp
                do j=1, n
                    t1 = x - xn
                    t2 = xn(j) - xn
                    l(j) = product(t1, MASK = (t2 /= 0.0_mp))/product(t2, MASK = (t2 /= 0.0_mp))
                enddo

                inter(i) = dot_product(l, yn)
            enddo

        end subroutine lagran

        subroutine quad(a, b, c, x_1, x_2)
            real(mp) :: a, b, c
            real(mp) :: x_1, x_2

            real(mp) :: D

            D = b**2 - 4.0_mp*a*c

            if(D >= 0.0_mp) then
                x_1 = (-b + sqrt(D))/(2.0_mp*a)
                x_2 = (-b - sqrt(D))/(2.0_mp*a)
                if(abs(D) < 1.0e-15_mp) then; x_1 = imp_r; x_2 = imp_r; endif
                if(abs(a) < 1.0e-15_mp) then
                    x_1 = -c/b
                    x_2 = imp_r
                endif
            else
                x_1 = imp_r
                x_2 = imp_r
            endif

            if(abs(a) < 1.0e-15_mp .and. abs(b) < 1.0e-15_mp) then
                x_1 = 3.0_mp
                x_2 = 3.0_mp
            endif

        end subroutine quad

        subroutine sort(s)
            real(mp), dimension(:) :: s
            real(mp) :: b
            integer :: i, j, n

            n = size(s)

            do i = 1, n
                do j = 1, n-i
                    if(s(j)>s(j+1)) then
                        b = s(j)
                        s(j) = s(j+1)
                        s(j+1) = b
                    endif
                enddo
            enddo

        end subroutine

        function sgn(a)
            real(mp) :: sgn, a

            if(a < 0) then
                sgn = -1.0_mp
            else
                if (a == 0) then
                    sgn = 0.0_mp
                else
                    sgn = 1.0_mp
                endif
            endif
        end function sgn

        function check(x, y, z, theta, phi_1, phi_2)
            real(mp) :: theta, phi_1, phi_2
            real(mp) :: x, y, z
            logical :: check

            real(mp) :: a, r, r_xy

            r = sqrt(x**2 + y**2 + z**2)
            a = abs(sin(theta)*y + cos(theta)*z)/r
            r_xy = sqrt(x**2 + y**2)

            check = a <= cos(phi_1) .and. a >= cos(phi_2)! .and. r > 1.0_mp .and. ((z < 0 .and. r_xy > 1.0_mp) .or. (z >= 0))

        end function check

        ! subroutine right_shift(P)
        !     real(mp), dimension(0:) :: P
        !     integer :: n, i

        !     n = size(P) - 1

        !     do i = 0, n-1
        !         P(i) = P(i+1)
        !     enddo

        !     P(n) = 0.0_mp
        ! end subroutine right_shift

        subroutine derivate_pol(dP, P, n)
            real(mp), dimension (0:) :: dP, P
            integer :: n

            real(mp), dimension (0:n-1) :: cf
            integer :: i
            dP = 0d0

            do i=0, n-1
                cf(i) = real(i+1)
            enddo

            dP(0:n-1) = P(1:n)
            dP(0:n-1) = cf*dP(0:n-1)
        end subroutine derivate_pol

        function calc_pol(Pn, x, n)
            real(mp), dimension (0:n) :: Pn
            real(mp) :: x, calc_pol
            integer :: n
            
            real(mp), dimension (0:n) :: xn
            integer :: i

            do i=0, n
                xn(i) = x**i
            enddo

            calc_pol = sum(Pn*xn)
        end function calc_pol

        subroutine bisection_pol(x, Pn, a, b, n)
            real(mp), dimension(0:n) :: Pn
            real(mp) :: x, a, b
            integer :: n

            real(mp) :: h, aa, bb, pa, pb, px

            aa = a; bb = b; h = bb-aa

            do while (h > 1d-10)
                x = aa + h/2
                pa = calc_pol(Pn, aa, n)
                pb = calc_pol(Pn, bb, n)
                px = calc_pol(Pn, x, n)
                !write(*, *) h
                !write(*, *) pa, px, pb
                !write(*, *) aa, x, bb
                if(pa == 0d0) then; x = aa; exit; endif
                if(pb == 0d0) then; x = bb; exit; endif
                if(pa*px < 0d0) then; bb=x; h=h/2
                else; aa=x; h=h/2
                endif
            enddo
        end subroutine bisection_pol

        subroutine bisection_pol_all(roots, Pn, m_int, n, m, nr)
            real(mp), dimension (0:n) :: Pn
            real(mp), dimension (n) :: roots
            real(mp), dimension (m) :: m_int
            integer :: n, m, nr

            integer :: i, nm
            real(mp) :: a, b, h, pa, pb
            
            nr = 0
            !write(*, *) n, m
            if (n == 1) then
                a = m_int(1); b = m_int(2); h = b-a
                if(h > 1d-11) then
                    roots(1) = abs(calc_pol(Pn,a,n))/(abs(calc_pol(Pn,a,n))&
                            + abs(calc_pol(Pn,b,n)))*h + a
                    !write(*, *) calc_pol(Pn,a,n), calc_pol(Pn,b,n), a ,b
                    nr = nr+1
                endif
            else
                do i = 1, m-1
                    a = m_int(i); b = m_int(i+1)
                    pa = calc_pol(Pn, a, n)
                    pb = calc_pol(Pn, b, n)
                    !write(*, *) pa, pb, a, b
                    if(b - a > 1d-11 .and. pa*pb < 0d0) then
                        call bisection_pol(roots(i), Pn, a, b, n)
                        nr = nr+1
                        !write(*, *) roots(i)                      
                    endif
                enddo   
            endif
            !write(*, *) nr
            m = 2+nr
        end subroutine bisection_pol_all

        subroutine pol_roots(x, Pn, n, a, b)
            real(mp), dimension (0:) :: Pn
            real(mp), dimension (:) :: x
            real(mp) :: a, b
            integer :: n

            real(mp), dimension(0:n,0:n-1) :: dPn ! производные Pn
            real(mp), dimension(n+2) :: m_int ! промежутки монотонности
            real(mp), dimension(n) :: roots
            real(mp) :: s
            integer :: i, m, mm, nr ! mm - количество промежутков монотонности + 1
            !write(*, *) 'pfpfp'
            dPn(:, 0) = Pn


            do i=1, n-1
                call derivate_pol(dPn(:, i), dPn(:, i-1), n)
            enddo

            !do i=0, n-1
            !    write(*, '(6f7.2)') dPn(:, i)
            !enddo

            m_int = 0d0
            m_int(1:2) = (/a, b/)
            mm = 2
            m = 0
            

            do i=n-1, 0, -1
                m = m+1
                roots = b
                !write(*, *) i
                call bisection_pol_all(roots(:m), dPn(0:m, i), m_int(:mm), m, mm, nr)
                !write(*, '(8f10.2)') m_int
                m_int(1:2) = (/a, b/)
                !write(*, '(6f10.2)') dPn(:, i)
                m_int(3:m+2) = roots(:m)
                call sort(m_int(:m+2)) 
                !write(*, '(8f10.2)') m_int
                !read(*, *) 
            enddo
            s = 0
            do i = 2, nr+1
                s = s + calc_pol(Pn, m_int(i), n)
            enddo
            !write(*, *) s, nr
            !write(*, *) m_int
            x = m_int(2:size(x)+2)
            !write(*, *) x
            !read(*, *) 
                       
        end subroutine pol_roots

end module modmod
