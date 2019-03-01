program kurs
    use mymod
    use mpmod
    use mainmodule
    use read_write
    use init
    use splain
    !use mainmod
    implicit none
    
    !параметры области и точность задания сетки
    !массив границ области на заданной сетке
    real(mp) :: theta, phi_1, phi_2
    real(mp) :: i, psi, alpha
    real(mp), allocatable, dimension(:, :, :) :: borders
    integer :: m, m_z
    integer :: n
    
    integer :: l, u
    real(mp) :: f_lu, nu
    real(mp), allocatable, dimension(:) :: rul
    real(mp), allocatable, dimension(:,:) :: spl_nl, spl_nu

    real(mp), dimension(3) :: omega
    real(mp) :: omega_m
    real(mp) :: T_acr, T_s, T_hot, Mdot

    real(mp), dimension(0:6, 15) :: argument_list

    real(mp), allocatable, dimension(:) :: prof, emm, absorb
    real(mp) :: h1, h2

    character(2) :: sym

    logical :: d2, nostark

    integer :: field_type

    argument_list = 0d0
 
    call read_arguments(argument_list)


    nostark = argument_list(0,14) < 5d-1

    if(argument_list(0,9) < 5d-1) then
        call read_field(T_acr, i, psi, alpha, phi_1, phi_2, m, m_z) !считывание параметров области
        field_type = 2
        if(argument_list(0, 4) > 5d-1) then
            i = argument_list(1, 4)/18d1*pi
            psi = argument_list(2, 4)/18d1*pi
            alpha = argument_list(3, 4)/18d1*pi
            phi_1 = argument_list(4, 4)/18d1*pi
            phi_2 = argument_list(5, 4)/18d1*pi
            field_type = int(argument_list(0,4))
        endif

        theta = acos(cos(i)*cos(alpha) + sin(i)*sin(alpha)*cos(psi))

        omega(1) = sin(i)*sin(alpha)*sin(psi)/sin(theta)
        if(sin(theta) == 0.0_mp) then; omega(1) = 0.0_mp; endif
        omega(2) = sin(i)*(cos(alpha)*sin(i) - cos(i)*sin(alpha)*cos(psi))/sin(theta)
        omega(3) = cos(i)
        
    
    if(argument_list(0,12)>5d-1) then; r_end = argument_list(1, 12); endif
    
    if(argument_list(0,13)>5d-1) then; T_acr = argument_list(1, 13); endif
        
    if(argument_list(0, 10) > 5d-1) then
        u = int(argument_list(1, 10))
        l = int(argument_list(2, 10))
            
    endif

    call read_line(l, u, f_lu, n, nu, int(argument_list(0, 10)))

    !nu = 2.8d14
    
    write(*, '(2i4)') u, l, int(argument_list(0, 10))

        if(argument_list(0, 5) > 5d-1) then
            m = int(argument_list(1, 5))
            m_z = int(argument_list(2, 5))
            n = int(argument_list(3, 5))
            write(*, *) m, m_z, n
        endif

                
    
        allocate(borders(-3:8, m, m))
        write(*, *) "yes"

        write(*,*) phi_1*180/pi, phi_2*180/pi
        
        call read_star(T_s, R_s, v_esc, omega_m)
        
        if(argument_list(0, 11)>5d-1) then
            T_s = argument_list(1, 11)
            R_s = argument_list(2, 11)
            v_esc =  argument_list(3 ,11)
            omega_m = argument_list(4, 11)
            v_esc = sqrt(v_esc/R_s)*618_mp
            R_s = R_s*6.96e5_mp
            !R_s = 3.5d6
            omega_m = omega_m/R_s
            R_s = R_s*1.0e5_mp
            v_esc = v_esc*1.0e5_mp
        endif
        
        select case(field_type)
        case(2)
            h1 = phi_1
            h2 = phi_2
        case(3)
            h1 = asin(sqrt(1/(phi_2/pi*180)))
            h2 = asin(sqrt(1/(phi_1/pi*180)))
        end select

        if(argument_list(0,15) > 5d-1) then
            Mdot = 1d1**argument_list(1,15)
            write(*,*) h1,h2,abs(cos(h2)-cos(h1)), (1d0 - 2d0/(1d0/sin(h1)**2+1d0/sin(h2)**2))
            T_hot = (v_esc**2/2d0*Mdot*2d33/365.25d0/86400d0*&
                (1d0 - 2d0/(1d0/sin(h1)**2+1d0/sin(h2)**2))/R_s**2/4d0/pi/abs(cos(h2)-cos(h1))/5.67d-5)**0.25d0
            if(T_hot < T_s) then; write(*,*) "Hotspot isn't hot"; endif
        else
            T_hot = T_s
        endif

        write(*,*) "T_hot = ", T_hot


        call init_field(theta, phi_1, phi_2, borders, m, argument_list(:, 7:8), field_type)
        
        if(argument_list(0, 6) > 5d-1) then
            r_0 = argument_list(2, 6)
            !if(argument_list(1, 6) > 0d0) then
            !    v_esc = argument_list(2, 6)
            !endif
        endif

        allocate(prof(n), emm(n), absorb(n))

        OPEN(10, FILE = 'in_data/level_data')

        read(10, *) sym, spl_n

        ! write(*, *) spl_n
        CLOSE(10)

        allocate(spl_nu(0:3, spl_n-1), spl_nl(0:3, spl_n-1), rul(spl_n))

        

        omega = omega_m*omega
        !ut = sqrt(2.0_mp*k_b*T_acr/4d0/aem)
        d2 = argument_list(3,10) > 5d-1
        if(.not. d2) then
            call init_spl(rul, spl_nl, spl_nu, u, l)
        else
            call init_2d_interpolation(u,l,phi_1/pi*180, phi_2/pi*180)
        endif


        call main_subroutine(prof, emm, absorb, borders, omega, T_s, T_acr, nu, f_lu, l, u, spl_nl, spl_nu, rul,&
         m_z, argument_list, theta, d2, nostark, sin(h1), sin(h2), T_hot)

        call write_prof(prof, emm, absorb, nu)

        deallocate(spl_nu, spl_nl, rul)
        deallocate(borders, emm, absorb)
    endif
end program kurs
