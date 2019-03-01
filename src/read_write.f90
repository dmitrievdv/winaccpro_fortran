module read_write
    use mymod
    use mpmod
    implicit none

    contains

       ! function calcSource(r, n1, r_0) !CHECK THE STRING WHERE TAU IS CALCULATED
       !      real(mp) :: r, n1, calcSource

       !      real(mp) :: beta, tau, psi, mu, omega, omega_s, tet, r_0
       !      real(mp) :: v_0, v_inf
       !      real(mp) :: alpha, nu, k
       !      integer :: i, n

       !      n=5000

       !      tet = sqrt(1d0-(R_s/r)**(2))


       !      v_0 = 1.0d6
       !      v_inf = 5.0d7
       !      alpha = 1.8d0
       !      k = 7.8d-13

       !      omega = 0d0
       !      omega_s = 0d0

       !      !OPEN(11, FILE = 'output/debug3', STATUS = 'REPLACE')

       !      do i=1, n
       !          !mu = i*2d0/n - 1d0/n - 1d0
       !          mu = i*1d0/n - 0.5d0/n
       !          psi = r_0*alpha*(v_inf-v_0)*(1d0-r_0/r)**(alpha-1d0)/r**2
       !          !write(*, *) (1d0-mu**2), R_s*r
       !          psi = psi*mu**2 + (v_0 + (v_inf-v_0)*(1d0 - r_0/r)**alpha)/r*(1d0-mu**2)
       !          !tau = n1*!ut/psi*k*sqrt(pi)*ut*2.8d14/v_c !!!!this must be uncommented if you need to use calcSource
       !          beta = (1d0 - exp(-tau))/tau 
       !          !write(*, *) mu, pi/n*i - pi/2d0/n
       !          !write(*, *) psi, tau
       !          !write(*, *) beta, r, 1/r
       !          !write(*, *) tet, W, n1
       !          !write(*, *)
       !          !write(*, *)
       !          !write(11, *) r, tet, mu, psi, tau, beta
       !          omega = omega + beta/n
       !          if(mu > tet) then; omega_s = omega_s + beta/n; endif
       !      enddo

       !      calcSource = omega_s/omega/2d0

       !      !CLOSE(11)

       !  end function calcSource

        subroutine read_nlnu(r, nl, nu, u, l)
            real(mp), dimension(:) :: r, nl, nu
            integer :: l, u
            integer :: n

            real(mp), dimension(8) :: file_data
            integer :: i, input_file

            n = size(r)
            input_file = 1090

            OPEN(input_file, FILE = "in_data/level_data")

            i=1

            read(input_file, *)
            read(input_file, *)

            do i=1, n
                read(input_file, *) file_data
             !write(*, *) file_data
                r(i) = file_data(1)
                nl(i) = file_data(l+2)
                nu(i) = file_data(u+2)
        !write(*, *) u, l
            enddo

            r = r/R_s


        end subroutine read_nlnu

        subroutine read_field(T_acr, i, psi, alpha, phi_1, phi_2, m_0, m_z)
            real(mp) :: T_acr, i, psi, alpha, phi_1, phi_2
            integer :: m_0, field_in, m_z

            field_in = 9876

            OPEN(field_in, FILE = 'in_data/field.dat')

            read(field_in, *) T_acr, i, psi, alpha, phi_1, phi_2, r_end, m_0, m_z, r_0
            i = i/180.0_mp*pi
            psi = psi/180.0_mp*pi
            alpha = alpha/180.0_mp*pi
            phi_1 = phi_1/180.0_mp*pi
            phi_2 = phi_2/180.0_mp*pi

            imp_r = r_end + 1.0_mp

            CLOSE(field_in)

        end subroutine read_field

        subroutine read_line(l, u, f_lu, n, nu, trig)
            real(mp) :: f_lu, nu
            integer :: n, l, u, trig

            real(mp), dimension(8, 7) :: A
            integer :: line_in

            line_in = 4567

            A(:, 1) = (/0.00d0, 4.67d8, 5.54d7, 1.27d7, 4.10d6, 1.64d6, 7.53d5, 3.85d5/)
            A(:, 2) = (/0.00d0, 0.00d0, 4.39d7, 8.37d6, 2.52d6, 9.68d5, 4.37d5, 2.20d5/)
            A(:, 3) = (/0.00d0, 0.00d0, 0.00d0, 8.94d6, 2.19d6, 7.74d5, 3.34d5, 1.64d5/)
            A(:, 4) = (/0.00d0, 0.00d0, 0.00d0, 0.00d0, 2.68d6, 7.67d5, 3.03d5, 1.42d5/)
            A(:, 5) = (/0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 1.02d6, 3.24d5, 1.38d5/)
            A(:, 6) = (/0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 4.50d5, 1.55d5/)
            A(:, 7) = (/0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0, 2.26d5/)

            OPEN(line_in, FILE = 'in_data/line.dat')
        
        write(*, *) trig
        
            if(trig == 0) then; read(line_in, *) l, u, n
        else; read(line_in, *) n, n, n;
        endif

            nu = v_c*109677.583407_mp*(1.0_mp/l**2 - 1.0_mp/u**2)
            !nu = 2.8e14

            write(*, *) m_e, v_c, q_e, pi

            f_lu = m_e*v_c**3/8d0/pi**2/q_e**2/nu**2*u**2*1d0/l**2*A(u, l)
            !f_lu =v_c**2*A(u, l)/8d0/pi/nu**2*real(u)**2/real(l)**2

            write(*, *) A(u, l), f_lu

            CLOSE(line_in)
        end subroutine read_line

        subroutine read_star(T_eff, R_s, v_esc, omega_m)
            real(mp) :: R_s, v_esc, omega_m, T_eff

            integer :: star_in

            star_in = 7876

            OPEN(star_in, FILE = 'in_data/star.dat')

            read(star_in, *) T_eff, R_s, v_esc, omega_m
            v_esc = sqrt(v_esc/R_s)*618_mp
            R_s = R_s*6.96e5_mp
            R_s = 3.5d6 
            omega_m = omega_m/R_s
            R_s = R_s*1.0e5_mp
            v_esc = v_esc*1.0e5_mp

            CLOSE(star_in)

        end subroutine read_star

        subroutine write_prof(prof, emm, absorb, nu_0)
            real(mp), allocatable, dimension(:) :: prof, emm, absorb
            integer :: n

            real(mp) :: nu_s, nu_step, nu_0
            integer :: i, absorb_out, emm_out, prof_out
            n = size(emm)


            nu_s = v_esc*sqrt(1.0_mp - 1.0_mp/r_0)/v_c*nu_0*1.1 !here nu_s is the dopler shift 
            !                                                    on the maximum of the accretion velocity
            !                                                    multiplied by 1.1
            nu_step = nu_s/n*2.0_mp
            nu_s = nu_0-nu_s !                                   actual nu_s

            !write(*, '(a6, e10.3)') 'x_0 = ', x_0*ut/1d5

            emm_out = 2
            absorb_out = 1
            prof_out = 3

            OPEN(1, FILE = 'output/absorb.dat', STATUS = 'REPLACE')

            OPEN(2, FILE = 'output/emm.dat', STATUS = 'REPLACE')

            OPEN(3, FILE = 'output/prof.dat', STATUS = 'REPLACE')

            do i=1, n
                write(absorb_out, *) -(nu_s+i*nu_step-0.5_mp*nu_step-nu_0)/nu_0*v_c/1d5, absorb(i), i
                write(emm_out, *) -(nu_s+i*nu_step-0.5_mp*nu_step-nu_0)/nu_0*v_c/1d5, emm(i), i
                write(prof_out, *) -(nu_s+i*nu_step-0.5_mp*nu_step-nu_0)/nu_0*v_c/1d5, prof(i), i
            enddo

            CLOSE(absorb_out)
            CLOSE(emm_out)

        end subroutine write_prof

        subroutine make_plotmap(S_max, n, nu_0, nu_s, nu_step)
            real(mp) :: S_max, nu_0, nu_s, nu_step
            integer :: n

            integer :: i

            OPEN(1, FILE = 'output/plotmap.gnu', STATUS = 'REPLACE')

            write(1, *) "set pm3d map"
            write(1, *) "set size ratio -1"
            !write(1, *) "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb"black""
            write(1, *) "set palette defined (1 1 1 1, 1 0 0 0)"
            write(1, *)
            write(1, *) "set term eps size 3, 3"
            write(1, *) 'set output "map.eps"'
            write(1, *) "splot 'map.dat' u 1:2:3"
            write(1, *)
            write(1, *) "set term gif animate size 500, 500"
            write(1, *) 'set output "map.gif"'
            !write(1, *) "set logscale cb"
            write(1, *) "set cbrange [0:", S_max, "]"

            do i=1, n
                write(1, *)
                write(1, '(a11, f10.3, a1)') 'set title "', -(nu_s+i*nu_step-0.5_mp*nu_step-nu_0)/nu_0*v_c/1d5, '"'
                write(1, '(a19, i3, a13)') "splot 'map_data/map", i+100, ".dat' u 1:2:3"
                write(1, *)
            enddo

            CLOSE(1)

        end subroutine make_plotmap

        subroutine read_arguments(argument_list)
            real(mp), dimension(0:, :) :: argument_list
            integer :: n

            character(10) :: argument
            integer :: i

            n = size(argument_list(:, 1))*size(argument_list(0, :))
            argument_list = 0d0

            do i=1, n
                call getarg(i, argument)

                if(len_trim(argument) == 0) then
                  exit
                endif
                if(argument_list(0, 9) < 5d-1 .and. &
                    ( argument(1:1) == '+' .or. argument(1:1) == '-' &
                    .or. argument(1:1) == ':' .or. argument(1:1) == '?')) then

                    call process_argument(trim(argument), argument_list, i)
                endif
            enddo

            if(argument_list(0, 1) == 1d0) then
                argument_list(0, 2) = 0d0
            endif

        end subroutine read_arguments

        subroutine process_argument(argument, argument_list, argument_number)
            real(mp), dimension(0:, :) :: argument_list
            integer :: n, argument_number
            character(*) :: argument

            character(20) :: arg
            integer :: i


            n = size(argument_list)

            select case(argument)
                case("-main")
                    argument_list(0, 1) = 1d0
                case("+map")
                    argument_list(0, 2) = 1d0
                case("+z_data")
                    argument_list(0, 3) = 1d0
                    do i=1, 3
                        call getarg(argument_number+i, arg)
                        read(arg, *) argument_list(i, 3)
                    enddo
                case(":cone")
                    if(argument_list(0,4)<0.5d0) then
                        stop("Err:cone: Undefined axis!")
                    endif
                    if(argument_list(0,4)>1.5d0) then
                        stop("Err:cone: Field is defined already")
                    endif
                    argument_list(0, 4) = 2d0
                    do i=4, 5
                        call getarg(argument_number+i-3, arg)
                        read(arg, *) argument_list(i, 4)
                    enddo
                case(":dipl")
                    if(argument_list(0,4)<0.5d0) then
                        stop("Err:dip: Undefined axis!")
                    endif
                    if(argument_list(0,4)>1.5d0) then
                        stop("Err:dip: Field is defined already")
                    endif
                    argument_list(0, 4) = 3d0
                    do i=4, 5
                        call getarg(argument_number+i-3, arg)
                        read(arg, *) argument_list(i, 4)
                    enddo
                case(":axis")
                    argument_list(0, 4) = 1d0
                    do i=1, 3
                        call getarg(argument_number+i, arg)
                        read(arg, *) argument_list(i, 4)
                    enddo
                case(":acc")
                    argument_list(0, 5) = 1d0
                    do i=1, 3
                        call getarg(argument_number+i, arg)
                        read(arg, *) argument_list(i, 5)
                    enddo
                    call getarg(argument_number+4, arg)
                    argument_list(4, 5) = 0d0
                    select case(arg)
                    case('aft')
                        argument_list(4, 5) = 2d0
                    case('bef')
                        argument_list(4, 5) = 1d0
                    end select
                case(":wind")
                    argument_list(0, 6) = 1d0
                    call getarg(argument_number+1, arg)
                    read(arg, *) argument_list(1, 6)
                    call getarg(argument_number+2, arg)
                    read(arg, *) argument_list(2, 6)
                    argument_list(2, 6) = 1d5*argument_list(2, 6)
                case(":accr")
                    argument_list(0, 6) = 1d0
                    argument_list(1, 6) = -1d0
                    call getarg(argument_number+1, arg)
                    read(arg, *) argument_list(2, 6)
                case(":dom")
                    argument_list(0, 7) = 1d0
                    do i=1, 6
                        call getarg(argument_number+i, arg)
                        if(arg(1:1) == 'd') then
                            argument_list(i, 7) = 0d0
                        else
                            argument_list(i, 7) = 1d0
                            read(arg, *) argument_list(i, 8)
                        endif
                    enddo
                case(":grid")
                    call getarg(argument_number+1, arg)
                    select case(arg)
                        case("dec")
                            argument_list(0, 8) = 1d0
                        case("pol")
                            argument_list(0, 8) = 0d0
                    end select
                case (":line")
                        argument_list(0, 10) = 1d0
                        i=3
                        call getarg(argument_number+1, arg)
                        if(trim(arg)=='2d') then
                          argument_list(i, 10) = 1d0
                        endif
                        do i=1, 2
                            call getarg(argument_number+i+1, arg)
                            read(arg, *) argument_list(i, 10)
                        enddo
                        i=3
                case (":star")
                    argument_list(0, 11) = 1d0
                    do i=1, 4
                        call getarg(argument_number+i, arg)
                        read(arg, *) argument_list(i, 11)
                    enddo
                case(":rad")    
                    argument_list(0, 12) = 1d0
                    call getarg(argument_number+1, arg)
                    read(arg, *) argument_list(1, 12)
                case (":temp")
                    argument_list(0, 13) = 1d0
                    call getarg(argument_number+1, arg)
                    read(arg, *) argument_list(1, 13)
                    call getarg(argument_number+2, arg)
                    select case(arg)
                    case('iso')
                        argument_list(2, 13) = 1d0
                    case('pow')
                        argument_list(2, 13) = 2d0
                    case('mag')
                        argument_list(2, 13) = 3d0
                    end select
                case (":stark")
                    argument_list(0,14) = 1d0
                case ("+hotspot")
                    argument_list(0, 15) = 1d0
                    call getarg(argument_number+1, arg)
                    read(arg, *) argument_list(1, 15)
                case("?")
                    call getarg(argument_number+1, arg)
                    select case(arg)
                        case("-main")
                            write(*, *) "Switch off profile calculation"
                        case("+map")
                            write(*, *) "Add map output and generate gnuplot file to plot maps"
                            write(*, *) "Generate two maps:"
                            write(*, *) "1. Column density"
                            write(*, *) "2. Emission map on each frequency unit"
!                             write(*, *)
!                             write(*, *) "Colunn density map output data locate in 'map.dat' file. ", &
!                                         "File 'map.dat' cotains 3 columns and m*(m+1) (grid accuracy) lines, each m+1 line is empty. ", &
!                                         "First column is x coordinate, second is y and third is column density of gas in (x, y). ", &
!                                         "Each block of m+1 lines have constant phi coordinate (polar angle)."
!                             write(*, *) "cos(phi) = x/sqrt(x**2+y**2)"
!                             write(*, *) "sin(phi) = y/sqrt(x**2+y**2)"
!                             write(*, *) 
!                             write(*, *) "Emission map output data locate in 'map_data' directory. ", &
!                                         "This directory contains n(number of frequency units) files, ", &
!                                         "each file has same format, as 'map.dat' file, but instead of column density third column contains " 
                        case("+z_data")
                            write(*, *) "Add output of dependency of emission parameters", &
                            "(on given frequency) from z coordinate on given x and y coordinates"
                            write(*, *)
                            write(*, *) "Syntax: +z_data [x] [y] [v_z]"
                            write(*, *) " [x]: x coordinate on tangent plane (in star radii)"
                            write(*, *) " [y]: y coordinate on tangent plane (in star radii)"
                            write(*, *) " [v_z_nu]: radial velocity in km/s (need to set emission frequency)"
                            write(*, *)
                            write(*, *) "Output format:"
                            write(*, *) " File: z_data"
                            write(*, *) " Contains 10 columns and m_z (number of integration steps) lines"
                            write(*, *) " Columns:"
                            write(*, *) " 1. z: z coordinate [star radii]"
                            write(*, *) " 2. dz: Integration step [SGS]"
                            write(*, *) " 3. v_z: Radial velocity [km/s]"
                            write(*, *) " 4. n_u: Upper level population [SGS]"
                            write(*, *) " 5. n_l: Lower level population [SGS]"
                            write(*, *) " 6. v_z-v_z_nu [km/s]"
                            write(*, *) " 7. alpha: extinction profile [SGS]"
                            write(*, *) " 8. k_lu: extinction coefficient [SGS]"
                            write(*, *) " 9. tau: optical depth"
                            write(*, *) " 10. S_ul: source function [SGS]"
                            write(*, *)
                            write(*, *) "Example: "
                            write(*, *) " +z_data 0 1.2 0"
                        case (":cone")
                            write(*, *) "Set emission domain parameters"
                            write(*, *)
                            write(*, *) "Syntax: :cone [i] [psi] [alpha] [phi_1] [phi_2]"
                            write(*, *) " [i]: angle between star spin axis and view axis [degrees]"
                            write(*, *) " [psi]: rotation phase (zero means that", &
                             " view axis, cone axis and rotation axis lean in one plane", &
                             " and cone axis leans between view and rotation axis) [degrees]"
                            write(*, *) " [alpha]: angle between rotation axis and cone axis [degrees]"
                            write(*, *) " [phi_1]: inner cone angle [degrees]"
                            write(*, *) " [phi_2]: outer cone angle [degrees]"
                            write(*, *)
                            write(*, *) "Example: "
                            write(*, *) " :cone 0 0 0 0 90 - sphere"
                        case (":acc")
                            write(*, *) "Set accuracy"
                            write(*, *)
                            write(*, *) "Syntax: :acc [m] [m_z] [n]"
                            write(*, *) " [m]: accuracy of grid on tangent plane"
                            write(*, *) " [m_z]: z axis domain accuracy"
                            write(*, *) " [n]: number of frequency units"
                            write(*, *) "Example: "
                            write(*, *) " :acc 100 100 100"
                        case (":wind")
                            write(*, *) "Calculate profile for wind instead of accretion"
                            write(*, *)
                            write(*, *) "Syntax: :wind [beta] [v_inf]"
                            write(*, *) " [beta]: wind parameter. Wind velocity calculation formula:"
                            write(*, *) "          thermal_velocity + v_inf*(1 - 1/r)^beta."
                            write(*, *) " [v_inf]: infinity velocity [km/s]"
                            write(*, *)
                            write(*, *) "Example: "
                            write(*, *) " :wind 0.5"
                        case (":dom")
                            write(*, *) "Set emission domain truncation"
                            write(*, *)
                            write(*, *) "Syntax: :dom [coord1_s] [coord1_f] [coord2_s] [coord2_f] [z_s] [z_f]"
                            write(*, *) "Program will calculate profile only in domain from", &
                              " coord1_s to coord1_f and coord2_s to coord2_f on tangent plane and", &
                              " from z_s to z_f on z axis. Use 'd' instead of _s and _f", &
                              " number if you don't want to truncate some coordinate"
                            write(*, *)
                            write(*, *) "Example: "
                            write(*, *) " :dom d d d d 0 2 - truncate only z axis"
                        case (":grid")
                            write(*, *) "Choose grid for tangent plane. ", &
                                        "Use :grid dec for decart or :grid pol for polar coordinats"
                        case default
                            write(*, *) "This program calculate profile of emission line in tTau star spectrum", &
                            " that formed in case of accretion or wind."
                            write(*, *)
                            write(*, *) "Input files:"
                            write(*, *)
                            write(*, *) " field.dat - domain parameters"
                            write(*, *) "  Lines:"
                            write(*, *) "  1. [T]: Gas temperature [K]"
                            write(*, *) "  2. [i]: angle between star spin axis and view axis [degrees]"
                            write(*, *) "  3. [psi]: rotation phase (zero means that", &
                             " view axis, cone axis and rotation axis lean in one plane", &
                             " and cone axis leans between view and rotation axis) [degrees]"
                            write(*, *) "  4. [alpha]: angle between rotation axis and cone axis [degrees]"
                            write(*, *) "  5. [phi_1]: inner cone angle [degrees]"
                            write(*, *) "  6. [phi_2]: outer cone angle [degrees]"
                            write(*, *) "  7. [r_end]: domain radius [star radii]"
                            write(*, *) "  8. [m]: accuracy of polar grid on tangent plane"
                            write(*, *) "  9. [m_z]: z axis domain accuracy"
                            write(*, *) "  10. [r_0]: accretion start radius [star radii]"
                            write(*, *)
                            write(*, *) " line.dat - line parameters"
                            write(*, *) "  Lines:"
                            write(*, *) "  1. [l]: lower level"
                            write(*, *) "  2. [u]: upper level"
                            write(*, *) "  3. [n]: number of frequency units"
                            write(*, *)
                            write(*, *) " star.dat - star parameters"
                            write(*, *) "  Lines: "
                            write(*, *) "  1. [T_s]: star temperature [K]"
                            write(*, *) "  2. [R_s]: star radius [Sun radii]"
                            write(*, *) "  3. [M_s]: star mass [Sun masses]"
                            write(*, *) "  4. [omega]: rotation speed on equator [km/s]"
                            write(*, *)
                            write(*, *) "Possible arguments:"
                            write(*, *) "1. -main"
                            write(*, *) "2. +map"
                            write(*, *) "3. +z_data"
                            write(*, *) "4. :cone"
                            write(*, *) "5. :acc"
                            write(*, *) "6. :wind"
                            write(*, *) "7. :dom"
                            write(*, *) "8. :grid"
                            write(*, *) 'To learn about each argument use "? argument_name"'
                        end select
                        argument_list(0, 9) = 1d0
                                        
                                                                
            end select
        end subroutine process_argument
end module read_write
