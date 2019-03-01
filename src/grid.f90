module grid
    use mymod
    implicit none

    contains

        subroutine grid_func(i, j, x, y, m_0, dS, field_parameters)
            real(mp) :: x, y, dS
            integer :: i, j, m_0

            real(mp), dimension(0:4) :: field_parameters
            real(mp) :: h_x1, h_x2, x1, x2
            real(mp) :: x1_s, x1_f, x2_s, x2_f

            x1_s = field_parameters(1)
            x1_f = field_parameters(2)

            if(field_parameters(0) < 5d-1) then
                x2_s = field_parameters(3)/1.8d2*pi
                x2_f = field_parameters(4)/1.8d2*pi
                h_x1 = (x1_f - x1_s)/m_0
                h_x2 = (x2_f - x2_s)/(m_0-1)

                x1 = x1_s + h_x1*i - 0.5_mp*h_x1
                x2 = x2_s + h_x2*(j-1)
                x = x1*cos(x2)
                y = x1*sin(x2)

                dS = x1*h_x1*h_x2
            else
                x2_s = field_parameters(3)
                x2_f = field_parameters(4)
                h_x1 = (x1_f - x1_s)/m_0
                h_x2 = (x2_f - x2_s)/m_0

                x1 = x1_s + h_x1*i - 0.5_mp*h_x1
                x2 = x2_s + h_x2*j - 0.5_mp*h_x2
                x = x1
                y = x2

                dS = h_x1*h_x2
            endif

        end subroutine grid_func


end module grid
