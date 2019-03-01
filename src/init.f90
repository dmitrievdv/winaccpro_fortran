module init
  use mymod
  use mpmod
  use modmod
  use grid
  use read_write
  use splain
  use inter2d_mod
  implicit none

  real(8), allocatable, dimension(:,:) :: popul_tab_u, popul_tab_l, stark_tab, Ttab
  real(8), allocatable, dimension(:) :: rm_grid, d_grid

  private :: popul_tab_u, popul_tab_l, Ttab

contains
  
  subroutine init_T_hot(T_hot)
    real(8) :: T_hot

    T_hot = Ttab(1,1)
  end subroutine

  function theta_to_normed_d(theta, rm) result(d)
    real(8) :: theta, rm, d

    real(8) :: d0, dk

    ! write(*,*) rm, sqrt(1d0/rm)

    d0 = d_th(asin(sqrt(1d0/rm)))
    dk = d_th(pi/2d0)

    ! write(*,*) d0,dk

    d = (d_th(theta)-d0)/(dk-d0)
  
  contains
    function d_th(th)
      real(8) :: d_th, th

      d_th = rm/2d0*(log(1d0 + sqrt(2d0)) + sqrt(2d0) - &
           log(cos(th) + sqrt(1d0 + cos(th)**2)) - cos(th)*sqrt(1d0 + cos(th)**2))

      ! write(*,*) th, log(1d0 + sqrt(2d0)) + sqrt(2d0), log(cos(th)+sqrt(1d0+cos(th)**2)) + &
      !                      cos(th)*sqrt(1d0+cos(th)**2)
    end function d_th
  end function theta_to_normed_d

  subroutine init_n(rul, n_l, n_u, stark, ut, x, y0, z0, theta, spl_nl, spl_nu, d2)
    real(mp), dimension(:, :) :: spl_nl, spl_nu
    real(mp), dimension(:) :: z0, n_l, n_u, stark, ut
    real(mp), dimension(:) :: rul
    logical :: d2
    real(8) :: x, y0, theta

    integer :: n, i, d2file=1090
    real(8) :: phi, rm, rmi, rmo, dt
    real(mp), dimension(size(z0)) :: r, z, y

    n = size(z0)
    y = y0
    z = z0
    r = sqrt(x**2 + y**2 + z**2)


    ! OPEN(d2file, FILE = 'test2dinter.dat')
    ! do i=1,size(popul_tab_l(1,:))*100
    !   rm = 2.2d0
    !   phi = (pi/2d0-asin(sqrt(1d0/rm)))/(size(popul_tab_l(1,:))*100)*i + asin(sqrt(1d0/rm)) 
    !   ! call inter2d_init_something(y_start = asin(sqrt(1d0/rm)))
    !   dt = theta_to_normed_d(phi, rm)
    !   write(d2file,*) rm*sin(phi)**2, exp(inter2d_interpolate(rm,dt,log(popul_tab_u))), &
    !                              exp(inter2d_interpolate(rm,dt,log(popul_tab_l))), &
    !                              inter2d_interpolate(rm,dt,stark_tab)
    ! enddo
    ! CLOSE(d2file)

    if(.not. d2) then
      stark = -1d0
      do i=1, n
        n_l(i) = splain_value(r(i), rul, spl_nl)
        n_u(i) = splain_value(r(i), rul, spl_nu)
      enddo
    else 
      ! write(*,*) 'int'

      y = -z0*sin(theta) + y0*cos(theta)
      z = z0*cos(theta) + y0*sin(theta)
      do i=1, n
        phi = asin(sqrt(y(i)**2+x**2)/r(i))
        rm = r(i)**3/(y(i)**2+x**2)
        ! call inter2d_init_something(y_start = asin(sqrt(1d0/rm)))
        ! call inter2d_get_something(x_start = rmi, x_end = rmo)
        ! write(*,*) r(i), z(i), phi
        ! write(*,*) rm
        ! if(rm >= rmi .and. rm <= rmo) then
          ! write(*,*) 'yes'
          ! write(*,*) rm_grid
          ! write(*,*) d_grid
          dt = theta_to_normed_d(phi, rm)
          ! write(*,*) 'yes', dt
          n_l(i) = exp(inter2d_interpolate(rm, dt, log(popul_tab_l)))
          n_u(i) = exp(inter2d_interpolate(rm, dt, log(popul_tab_u)))
          ut(i) = exp(inter2d_interpolate(rm, dt, log(Ttab)))
          stark(i) = exp(inter2d_interpolate(rm, dt, log(stark_tab)))
          ! write(*,*) 'yes', n_l(i), n_u(i), ut(i)
        ! endif
      enddo
    endif
  end subroutine init_n

  subroutine init_2d_interpolation(u,l,rmi,rmo)
    integer :: u,l
    real(8) :: rmi, rmo

    integer :: data = 9809
    integer :: nr, ntet, i, k
    character(1) :: s
    real(8), dimension(12) :: dat
    real(8) :: C_rad, C_vdW, C_stark, nh, ne, T

    C_rad = 6.5d-12; C_vdW = 4.4d-12; C_stark = 1.174d-11

    OPEN(data, FILE = 'in_data/level_data')
    read(data,*) s,nr,ntet

    allocate(popul_tab_l(nr,ntet), popul_tab_u(nr,ntet), Ttab(nr,ntet), stark_tab(nr,ntet))
    allocate(rm_grid(nr), d_grid(ntet))

    read(data,*)

    do i=1, nr 
      do k=1, ntet
        read(data, *) dat
        rm_grid(i) = dat(1)
        d_grid(k) = theta_to_normed_d(dat(2), dat(1))
        popul_tab_l(i,k) = dat(l+4)
        popul_tab_u(i,k) = dat(u+4)
        T = dat(3)
        Ttab(i,k) = T
        nh = dat(4); ne = dat(12)
        stark_tab(i,k) = C_rad + C_stark*(ne/1d12)**(2d0/3d0) + C_vdW*(nh/1d16)*(T/5d3)**3d-1
        write(*,*) stark_tab(i,k)
      enddo
    enddo



    call inter2d_init(rm_grid, d_grid)

  end subroutine init_2d_interpolation

  subroutine init_spl(r, spl_nl, spl_nu, u, l)
    real(mp), dimension(:, :) :: spl_nl, spl_nu
    real(mp), dimension(:) :: r
    integer :: i, n, u, l

    real(mp), dimension(size(r)) :: nl, nu
    real(mp) :: step, x

    n = size(r)

    call read_nlnu(r, nl, nu, u, l)
    !call read_nlnu(r, nl, nu, "brams.fs", 1, 3, 4, 4)

    call get_hermite(r, nl, spl_nl)
    call get_hermite(r, nu, spl_nu)

    OPEN(1, FILE = 'output/debug1', STATUS = 'REPLACE')
    OPEN(2, FILE = 'output/debug2', STATUS = 'REPLACE')

    step = 1d0/(n*100)

    do i = 1, n*100
      x = 1d0 + step*i - step/2
      write(2, *) x, splain_value(x, r, spl_nu), splain_value(x, r, spl_nl)
    enddo

    do i = 1, n
      write(1, *) r(i), nu(i), nl(i)

    enddo

    CLOSE(1)
    CLOSE(2)

  end subroutine init_spl

  subroutine init_ut_and_dv_d(ut, dv_d, r, T, nu_0, T_mod)
    real(mp), dimension(:) :: ut, dv_d, r
    real(mp) :: T, nu_0, T_mod

    ! write(*,*) 'kek', T_mod
    if(T_mod < 1.5d0) then; ut = sqrt(2.0_mp*k_b*T/aem)        ! T = const
    else if(T_mod < 2.5d0) then; ut = sqrt(2.0_mp*k_b*T*r**(-1d0/3d0)/aem)            ! T = T_e*(r/R_*)^(-1/3)
    else if(T_mod < 3.5d0) then
      ! write(*,*) ut
      ut = sqrt(2d0*k_b*ut/aem)
      ! write(*,*) ut 
    endif

    dv_d = nu_0*ut/v_c
  end subroutine init_ut_and_dv_d

  subroutine init_dipl_cosv_z(v_z, rs, x, y_0, z_0, theta)
    real(mp), dimension(:) :: v_z, z_0, rs
    real(mp) :: x, y_0, theta

    real(mp), dimension(size(z_0)) :: eta, z, y, r, phi, r2_0
    real(mp), dimension(size(z_0)) :: v_xy, v_x, v_y, v_y0
    integer :: i

    y = -z_0*sin(theta) + y_0*cos(theta)
    z = z_0*cos(theta) + y_0*sin(theta)

    eta = atan(y/x)
    if(x < 0d0) then; eta = eta + pi
    endif
    if(x == 0d0) then
      eta = sign(pi/2, y)
    endif

    r = sqrt(x**2 + y**2)
    phi = acos(z/sqrt(r**2+z**2))
    r2_0 = (r**2 + z**2)/sin(phi)**4

    v_xy = (6d0*z*(z**2 + r**2)**2)*sign(1d0, z)
    v_z = -(6d0*r*(z**2 + r**2)**2 - 4d0*r2_0*r**3)*sign(1d0, z)

    ! v_z = (6d0*z*(z**2 + r**2)**2)*sign(1d0, z)
    ! v_xy = (6d0*r*(z**2 + r**2)**2 - 4d0*r2_0*r**3)*sign(1d0, z)
    
    v_x = v_xy*cos(eta)/sqrt(v_z**2 + v_xy**2)
    v_y = v_xy*sin(eta)/sqrt(v_z**2 + v_xy**2)
    v_z = v_z/sqrt(v_z**2 + v_xy**2)

    v_y0 = v_z*sin(theta) + v_y*cos(theta)
    v_z = v_z*cos(theta) - v_y*sin(theta)

    y = eta*cos(theta) + phi*sin(theta)
    z = -eta*sin(theta) + phi*cos(theta)
    
    do i=1, size(z)
      if(abs(sqrt(r2_0(i))-3d0) < 1d-2) then
        if(abs(sqrt(r2_0(i-1))-3d0) >= 1d-2) then
          write(1,*) x, y_0, z_0(i), -v_x(i), -v_y0(i), -v_z(i)
        endif
      endif
      !write(1,*) r2_0(i), y(i), z(i), v_x(i), v_y0(i), v_z(i)
    enddo

    rs = sqrt(r2_0)
    ! read(*, *)
  end subroutine init_dipl_cosv_z

  subroutine init_v_z(v_z, r, z, ox, oy, omega, beta, v_inf, field_type, theta)
    real(mp), dimension(:) :: v_z, z, r
    real(mp), dimension(:) :: omega
    real(mp) :: ox, oy, theta
    real(mp) :: v_0, v_inf, beta
    integer :: field_type

    real(mp), dimension(size(z)) :: cosv_z, rs

    !write(*, *) field_type

    select case(field_type)
    case(2)
      cosv_z = z/r
      rs = r_0
    case(3)
      call init_dipl_cosv_z(cosv_z, rs, ox, oy, z, theta)

    case default
      stop('Err init_v_z: Invalid field type!')
    end select

    v_0 = 1d6
    if(v_inf < 5d-1) then
      v_inf = v_esc
    endif
    !beta = 1.0d0
    !x = 0.5_mp*v_esc*z/r/ut + R_s/ut*(omega(1)*oy + omega(2)*ox)
    !x = 0.5_mp*v_esc*sqrt(1.0_mp/(r_end - r + 1.0_mp) - 1.0_mp/r_0)*z/r/ut  + R_s/ut*(omega(1)*oy + omega(2)*ox)

    if(beta > 0) then
      v_z = (v_0 + (v_inf - v_0)*(1d0 - 3d0/r)**beta)*cosv_z + R_s*(omega(1)*oy - omega(2)*ox)
    else
      v_z = -v_esc*sqrt(1.0_mp/r - 1.0_mp/rs)*cosv_z  + R_s*(omega(1)*oy - omega(2)*ox)
    endif

    !x = 3d7*z/r/ut
      !x = 0.0_mp

  end subroutine init_v_z

  subroutine init_n_appr(n_l, n_u, mu_l, mu_u, r)
    real(mp), dimension(:) :: n_l, n_u, r
    real(mp) :: mu_l, mu_u

    n_l = 1.0e3_mp*r**(-mu_l)
    n_u = 1.0e3_mp*r**(-mu_u)

      !n_l = 5.0_mp
      !n_u = 5.0_mp
  end subroutine init_n_appr

  subroutine init_z(z_coord, ddz, border)
    real(mp), dimension(:) :: z_coord, ddz
    real(mp), dimension(:) :: border
    integer :: m_z

    real(mp), dimension(4) :: d, z0, dz
    integer :: i, n, nn

    d = 0.0_mp
    z0 = 0.0_mp
    dz = 0.0_mp
    m_z = size(z_coord)

    n = 4

    do i=1, 4
      d(i) = border(2*i) - border(2*i-1)
      if(d(i) < 1.0e-8_mp) then
        n = n-1
      else
        z0(i) = border(2*i-1)
      endif
    enddo

    if(n /= 0) then
      nn = m_z/n
      dz = d/nn

      do i=0, m_z-1

        z_coord(i+1) = z0(i/nn+1) + (mod(i, nn)+1)*dz(i/nn+1) - dz(i/nn+1)*0.5_mp
        ddz(i+1) = dz(i/nn+1)
      enddo
    else
      ddz = 0.0_mp
        !z_coord = 3.5_mp
    endif

  end subroutine init_z

  !todo - add z_max to check function

  subroutine cone_border(theta, phi_1, phi_2, border, borders, r_start, z_max, z_min, z_in, x, y, r)
    real(mp) :: theta, phi_1, phi_2
    real(mp), dimension(-3:) :: border
    real(mp), dimension(-3:) :: borders
    real(mp) :: x, y, r, z_max, z_min, z_in, r_start

    real(mp) :: a, b, c
    integer :: n, k

    n = 0

    border(1:) = imp_r

    if(r < r_end) then
    
      z_in = sqrt(r_start**2 - r**2)

      if(r <= 1.0_mp .and. z_min < sqrt(1d0 - r**2)) then
        z_min = sqrt(1.0_mp - r**2)
      endif

      if(r <= 1.0_mp .and. z_max < sqrt(1d0 - r**2)) then
        z_max = sqrt(1.0_mp - r**2)
      endif

      if(r > 1.0_mp .and. z_max > sqrt(r_end**2 - r**2)) then
        z_max = sqrt(r_end**2 - r**2)
      endif

      if(r > 1.0_mp .and. z_min < -sqrt(r_end**2 - r**2)) then
        z_min = -sqrt(r_end**2 - r**2)
      endif

      a = cos(theta)**2*tan(phi_1)**2 - sin(theta)**2
      b = 2.0_mp*y*cos(theta)*sin(theta)*(tan(phi_1)**2+1.0_mp)
      c = y**2*(sin(theta)**2*tan(phi_1)**2 - cos(theta)**2) - x**2
      call quad(a, b, c, border(1), border(2))

      a = cos(theta)**2*tan(phi_2)**2 - sin(theta)**2
      b = 2.0_mp*y*cos(theta)*sin(theta)*(tan(phi_2)**2+1.0_mp)
      c = y**2*(sin(theta)**2*tan(phi_2)**2 - cos(theta)**2) - x**2
      call quad(a, b, c, border(3), border(4))

      do k=1, 4
        if(border(k) > z_max .and. border(k) /= imp_r) then
          border(k) = imp_r
        endif
        if(border(k) < z_min) then
          border(k) = imp_r
        endif
      enddo

      border(5:8) = imp_r  
      do k=1, 6
            
        if(border(k) /= imp_r) then; n = n+1; endif
      enddo

      call sort(border(1:))

      if((z_max - z_min) > 1d-8 .and. border(1) /= z_min .and. check(x, y, z_min, theta, phi_1, phi_2)) then
        ! borders(0:) = border(0:)
        ! border(2:) = borders(1:7)
        border(2:) = border(1:7)
        border(1) = z_min
        n = n+1
      endif

      if((z_max - z_min) > 1d-8 .and. border(n) /= z_max .and. check(x, y, z_max, theta, phi_1, phi_2)) then
        border(n+1) = z_max
        n = n+1
      endif 

      do k=1, 6
        if(border(k) /= imp_r .and. abs(border(k)-border(k+1))<1d-8) then
          border(k) = imp_r
          border(k+1) = imp_r
        endif
      enddo

      call sort(border(1:))

      if(r < r_start) then
        !write(*, *) r, z_in
        !write(*, '(9f6.2)') border(0:)
        do k=1,8
          if(border(k) > -z_in .and. border(k) < z_in) then
            if(mod(k,2) == 1) then
              if(border(k+1) < z_in) then; border(k) = imp_r
              else; border(k) = z_in
              endif
            else
              if(border(k-1) > -z_in) then; border(k) = imp_r
              else; border(k) = -z_in
              endif
            endif
          else
            if(mod(k,2) == 1 .and. border(k) < -z_in .and. border(k+1) > z_in) then
              border(7) = -z_in
              border(8) = z_in
            endif
          endif
        enddo
        !write(*, *) i, j 
        !read(*, *)
      endif
        call sort(border(1:))

        do k=1, 8, 2
          border(0) = border(0) + (border(k+1) - border(k))
        enddo
            
    endif
  end subroutine cone_border

  subroutine dipl_border(theta, r_1, r_2, border, z_max, z_min, x, y, r)
    real(mp) :: theta, r_1, r_2
    real(mp), dimension(-3:8) :: border
    real(mp) :: x, y, r, z_max, z_min

    real(mp), dimension(0:6) :: Pn
    real(mp) :: c,s,a
    integer :: i, n

    !x = 0.5d0; y=0.0d0; r = sqrt(x**2+y**2)
    border(1:8) = imp_r

    if(r < r_end) then 
      
      c = cos(theta)
      s = sin(theta)
      a = x**2 + y**2*c**2
      n = 6

      !write(*, *) r_1, r_2, c, s, a

      Pn(0) = r**6 - r_1**2*a**2
      Pn(1) = 4d0*a*c*r_1**2*s*y
      Pn(2) = 3d0*r**4 - 4d0*c**2*s**2*r_1**2*y**2 - 2*a*r_1**2*s**2
      Pn(3) = 4d0*r_1**2*s**3*c*y
      Pn(4) = 3d0*r**2 - r_1**2*s**4
      Pn(5) = 0d0
      Pn(6) = 1d0
      
      call pol_roots(border(1:4), Pn, n, -imp_r, imp_r)


      Pn(0) = r**6 - r_2**2*a**2
      Pn(1) = 4d0*a*c*r_2**2*s*y
      Pn(2) = 3d0*r**4 - 4d0*c**2*s**2*r_2**2*y**2 - 2*a*r_2**2*s**2
      Pn(3) = 4d0*r_2**2*s**3*c*y
      Pn(4) = 3d0*r**2 - r_2**2*s**4
      Pn(5) = 0d0
      Pn(6) = 1d0

      call pol_roots(border(5:8), Pn, n, -imp_r, imp_r)

      call sort(border(1:8))

      !write(*, *) border(1:8)

      if(r <= 1d0) then
        do i=8, 1, -1
          if(border(i) <= sqrt(1d0 - r**2)) then
            if(mod(i, 2) == 0) then 
              border(1:i) = imp_r
              exit
            else
              border(i) = sqrt(1d0 - r**2)
              border(1:i-1) = imp_r
            endif
          endif

        enddo
      endif

      call sort(border(1:8))

      do i=1, 8, 2
        border(0) = border(0) + (border(i+1) - border(i))
      enddo
      !write(*, *) border(1:8)
      !read(*, *) 

    endif
  end subroutine dipl_border

  subroutine init_field(theta, phi_1, phi_2, borders, m_0, field_arguments, field_type)
    real(mp) :: theta, phi_1, phi_2
    real(mp), dimension(-3:, :, :) :: borders
    real(mp), dimension(-3:8) :: border
    real(mp), dimension(0:6, 0:1) :: field_arguments
    integer :: m_0, field_type

    real(mp) :: x, y, r, z_max, z_min, z_in, r_start
    integer :: j, i
    logical :: clipped_z
    real(mp), dimension(4) :: default_arguments

    imp_r = r_end+1d0

    r_start = 1.30000d12/3.5d11
    r_start = 1d0
    r_start = 3d0

    OPEN(1, FILE = 'output/map.dat', STATUS = 'REPLACE')

    clipped_z = field_arguments(5, 0) > 5d-1 .and. field_arguments(6, 0) > 5d-1

    if(field_arguments(0, 1) < 5d-1) then
      default_arguments = (/0d0, r_end, 0d0, 3.6d2/)
    else
      default_arguments = (/-r_end, r_end, -r_end, r_end/)
    endif

    do i=1, 4
      if(field_arguments(i, 0) < 5d-1) then
        field_arguments(i, 1) = default_arguments(i)
      endif
    enddo

    write(*, *) field_arguments(1:, 1)
    write(*, *) clipped_z
    write(*,*) phi_1/pi*180, phi_2/pi*180
    !read(*, *)
    !h_r = r_end/m_0
    !h_phi = 2.0_mp*pi/m_0

    OPEN(10, FILE = 'output/debug3', STATUS = 'REPLACE')

    if(clipped_z) then
      z_max = field_arguments(6, 1)
      z_min = field_arguments(5, 1)
    endif
    
    do i=1, m_0
      do j=1, m_0
        border(0) = 0.0_mp
        call grid_func(i, j, x, y, m_0, borders(-3, i, j), field_arguments(0:4, 1))
        r = sqrt(x**2 + y**2)

        if(.not. clipped_z) then
          z_max = sqrt(r_end**2 - r**2)
          z_min = -z_max
        endif

          select case(field_type)
          case(2)
            call cone_border(theta, phi_1, phi_2, border, borders(:,i,j), r_start, z_max, z_min, z_in, x, y, r)
          case(3)
            call dipl_border(theta, phi_1/pi*180d0, phi_2/pi*180d0, border, z_max, z_min, x, y, r)
          case default
            stop("Err init_field: Invalid field type!")
          end select
        
          !if(border(0) > 0.0_mp) then
          borders(0:, i, j) = border(0:)
          borders(-2:-1, i, j) = (/x, y /)


          !endif
          write(1, *) x, y, border(0), border(1:8)


          if(j == m_0/2) then; write(10, *) r, z_min, border(0); endif

        enddo
        write(1, *)
      enddo

      CLOSE(1)
      CLOSE(10)
    end subroutine init_field

    subroutine init_flux(flux, T, nu_0, nu_s, nu_step)
      real(mp) :: nu_0, nu_s, nu_step
      real(mp), dimension(:) :: flux

      integer :: n, i
      real(mp) :: nu, T

      n = size(flux)
      !T = ut**2*aem/k_b/3d0

      !T = 5d3
      write(*, *) h*nu_0/k_b
      write(*, *) T, 1d0/(exp(h*nu_0/k_b/T) - 1d0)

      do i=1, n
        nu = (nu_s + nu_step*i - 0.5_mp*nu_step)

        flux(i) = 2d0*h*nu**3/v_c**2*(1d0/(exp(h*nu/k_b/T) - 1d0)) !+ 1d-1/(exp(h*nu/k_b/8d3) - 1d0))
          !write(*, *) flux(i)
      enddo
    end subroutine init_flux
  end module init

