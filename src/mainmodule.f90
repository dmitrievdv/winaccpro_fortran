! module gauss
!   implicit none
!   private

!   real(8), dimension(200) :: A, x

!   public :: init, A, x

! contains
  
!   subroutine init
!     OPEN(1234, FILE = 'gauss.dat')
!     read(1234,*) x
!     read(1234,*) A
!     CLOSE(1234)
!   end subroutine init

! end module gauss

! module gauss_mod
!   use gauss, only: &
!     gauss_init => init, &
!     gauss_A => A, &
!     gauss_x => x
! end module gauss_mod


module mainmodule
  use mymod
  use mpmod
  use read_write
  use init
  ! use gauss_mod

  implicit none

contains

  subroutine main_subroutine&
    (prof, emm, absorb, borders, omega, T_s, T_acr, nu_0, f_lu, l, u, spl_nl, spl_nu, rul,&
     m_z, argument_list, theta, d2, nostark, sin_h1, sin_h2, T_hot)
    real(mp), dimension(3) :: omega
    real(mp), dimension(-3:, :, :) :: borders
    real(mp), dimension(:) :: prof, emm, absorb
    real(mp), dimension(:, :) :: spl_nl, spl_nu
    real(mp), dimension(:) :: rul
    real(mp) :: nu_0, f_lu, theta, sin_h2, sin_h1
    integer :: l, u
    integer :: m, n, m_z
    real(mp), dimension(0:, :) :: argument_list
    logical :: d2, nostark

    real(mp), dimension(size(emm)) :: flux, star_flux, hot_flux, cont
    real(mp), dimension(m_z) :: z_coord, dz, n_l, n_u, stark, v_z, r, dx
    !real(mp), dimension(m_z) :: alpha, k_lu, tau, S_ul
    real(mp) :: alpha, k_lu, tau, S_ul
    real(mp), dimension(size(emm)) :: emm_z, abs_z
    real(mp) :: nu_s, nu_step !start frequency, frequency step
    real(mp) :: r_xy, x_cord, y_cord
    real(mp) :: I_s
    real(mp) :: S_max, S_min
    real(mp) :: dS, S_s
    real(mp), dimension(m_z) :: dv_d, ut
    real(mp) :: T_s, T_acr, T_hot                    !temperature of the star and temperature of the gas
    real(mp) :: time1, time2
    real(mp) :: k_lu_0
    real(mp) :: g_l, g_u
    character(26) :: filename
    integer :: i, j, k

    logical :: main, gif, z_data                          !аргументы
    real(mp) :: z_data_x, z_data_y, z_data_v

    integer, dimension(2) :: loc

    call cpu_time(time1)

    ! call gauss_init

    main = (argument_list(0, 1) < 5d-1)
    gif = (argument_list(0, 2) > 5d-1)
    z_data = (argument_list(0, 3) > 5d-1)

    z_data_x = argument_list(1, 3)
    z_data_y = argument_list(2, 3)
    z_data_v = argument_list(3, 3)


    m = size(borders(1, 1, :))
    n = size(emm)

    spl_n = size(rul)

    !dv_d = ut/v_c*nu_0
    I_s = 1.0_mp

    S_s = 0d0

    g_l = 2d0*l**2
    g_u = 2d0*u**2

    k_lu_0 = f_lu*par
    !k_lu_0 = 7.8d-13*sqrt(pi)*dv_d
    write(*,*) 'theta = ', theta/pi*180d0
    write(*, '(a7, e16.10)') ' lam = ', v_c/nu_0*1.0e8
    write(*, '(a7, e16.10)') ' nu = ', nu_0
    write(*, '(a7, e16.10, f10.5)') 'k_lu = ', k_lu_0, f_lu



    !write(*, *) f_lu*par/dv_d/sqrt(pi)

    k = 1

    nu_s = v_esc*sqrt(1.0_mp - 1.0_mp/r_0)/v_c*nu_0*1.1 !here nu_s is the dopler shift 
    !                                                    on the maximum of the accretion velocity
    !                                                    multiplied by 1.1
    nu_step = nu_s/n*2.0_mp
    nu_s = nu_0-nu_s !                                   actual nu_s

    call init_flux(star_flux, T_s, nu_0, nu_s, nu_step)  !calculation of the flux of the star (Plank)
    call init_flux(hot_flux, T_hot, nu_0, nu_s, nu_step)
    write(*,*) T_hot, hot_flux(1)/star_flux(1)
    cont = 0d0
    write(*, '(2(a6, e10.3))') 'l_s = ', v_c/nu_s*1d8, 'R_s = ', R_s
    ! OPEN(1, FILE = 'diplv.dat', STATUS = 'REPLACE')
    ! do i=1, 1000
    !     y = rand()*8d0 - 4d0
    !     x = rand()*8d0 - 4d0
    !     z_coord(1:8000) = (/ (-4d0 + j*8d0/8000d0 - 8d0/16000d0, j = 1,8000) /)
    !     call init_dipl_cosv_z(v_z(1:8000), dz, x, y, z_coord(1:8000), theta)
    ! enddo
    ! CLOSE(1)
    if(main) then

      
      write(*, *) 'main'
      emm = 0.0_mp
      absorb = 0.0_mp

      !u=1; l=1
      S_max = 0.0_mp; S_min = 0.0_mp

      if(gif) then
        do i=1, n
          write(filename, '(a19, i3, a4)') "output/map_data/map", 100+i, ".dat"
          OPEN(100+i, FILE = filename, STATUS = 'REPLACE')
          write(100+i, '(a8,f10.3)') '# v_z = ', -(nu_s + nu_step*i - 0.5_mp*nu_step-nu_0)/nu_0*v_c/1d5
        enddo
      endif
      ! borders(0, :,:) = 0d0
      !n_u = 1d-2; n_l = 1d-2
      OPEN(123456, FILE = 'output/star_map.dat', STATUS = 'REPLACE')
      do i = 1, m                 ! i - phi, j - r; i - x, j - y
        write(*,*) i*1d0/m*100, "%"
        do j = 1, m
          x_cord = borders(-2, j, i); y_cord = borders(-1, j, i)
          r_xy = sqrt(x_cord**2 + y_cord**2)
          dS = borders(-3, j, i)*R_s**2
          emm_z = 0d0
          abs_z = 0d0
          flux = choose_flux(x_cord, y_cord)
          cont = cont + flux*dS
          write(123456, *) x_cord, y_cord, flux(1)/star_flux(1)
          !if(r_xy <= 1d0) then; S_min = S_min + dS; endif
          if(borders(0, j, i) > 0.0_mp) then
            !write(*, * ) borders(0, j, i)
            call init_z(z_coord, dz, borders(1:8, j, i))
            dz = dz*R_s
            r = sqrt(x_cord**2 + y_cord**2 + z_coord**2)

            call init_n(rul, n_l, n_u, stark, ut, x_cord, y_cord, z_coord, theta, spl_nl, spl_nu, d2)
            if(u /= 3 .or. l /= 2) then; stark = -1d0; endif
            call init_v_z(v_z, r, z_coord, x_cord, y_cord, omega,&
                argument_list(1, 6), argument_list(2, 6), int(argument_list(0, 4)), theta)
            call init_ut_and_dv_d(ut, dv_d, r, T_acr, nu_0, argument_list(2, 13))
            do k = 1, n
              call emission_on_ray(nu_s + nu_step*k - 0.5_mp*nu_step, emm_z(k), abs_z(k))

              if(gif) then

                if(emm_z(k) > S_max) then
                  S_max = emm_z(k)
                endif

              endif
              emm(k) = emm(k) + emm_z(k)*dS

              if(r_xy < 1d0) then
                absorb(k) = absorb(k) + flux(k)*abs_z(k)*dS
                if(k==1) then; S_s = S_s + dS; endif
                !write(*, *) abs_z(k)/flux(k)
                !if(abs_z(k)/flux(k) > 1d0) then; read(*, *); endif  
              endif
              
            enddo

            

          else if(r_xy < 1d0) then
            absorb = absorb + flux*dS
            S_s = S_s + dS
          endif
          if(gif) then
            do k=1, n
              write(k + 100, *) x_cord, y_cord, emm_z(k)
            enddo
          endif
        enddo
        write(123456, *) 
        if(gif) then
          do k=1, n
            write(k + 100, *)
          enddo
        endif
      enddo

      !do i=1, n
      !  write(*, *) i, (x_0 + x_step*i - 0.5_mp*x_step)*ut/1d5, emm(i), absorb(i), flux(i)*pi*R_s**2
      !enddo

      !write(*, *) x_0*ut/1d5, S_min
      write(*, *)  S_s/pi/R_s**2
      write(*, *) hot_flux(1)/star_flux(1), cont(1)/star_flux(1)/S_s 
      emm = emm/cont
      absorb = absorb/cont
      prof = emm + absorb
      emm = emm + 1d0
      CLOSE(123456)

      if(gif) then
        do j = 1, n
          CLOSE(100+j)
        enddo
        call make_plotmap(S_max, n, nu_0, nu_s, nu_step)
      endif



    endif

    if(z_data) then
      loc = minloc(abs(borders(-2, :, :)-z_data_x) + abs(borders(-1, :, :)-z_data_y))
      i = loc(2)
      j = loc(1)
      write(*, *)
      write(*, '(a6, f6.3, 2(a7, f6.3))') 'x_z = ', x_cord, ' y_z = ', y_cord, ' d_z = ', borders(0, j, i)

      if(borders(0, j, i) > 1d-8) then
        call init_z(z_coord, dz, borders(1:8, j, i))
        dz = dz*R_s
        r = sqrt(x_cord**2 + y_cord**2 + z_coord**2)

        r_xy = sqrt(x_cord**2 + y_cord**2)
        dS = borders(-3, j, i)*R_s**2
        call init_n(rul, n_l, n_u, stark, ut, x_cord, y_cord, z_coord, theta, spl_nl, spl_nu, d2)
        if(u /= 3 .or. l /= 2 .or. nostark) then; stark = -1d0; endif
        call init_v_z(v_z, r, z_coord, x_cord, y_cord, omega,&
         argument_list(1, 6), argument_list(2, 6),int(argument_list(0, 4)), theta)
        call init_ut_and_dv_d(ut, dv_d, r, T_acr, nu_0, argument_list(2,13))
        write(*, *) maxloc(n_u), 'asdasd'
        k = 0
        write(*, *) k, m_z
        call emission_on_ray(nu_0-z_data_v/v_c*nu_0*1d5, emm_z(1), abs_z(1))

        write(*, *) dz(1)

      endif
    endif

    call cpu_time(time2)

    write(*, '(3(a7, e10.3))') 'time = ', time2-time1, ' v_0 = ', v_esc*sqrt(1.0_mp - 1.0_mp/r_0)*1.0e-5, ' v_e = ', v_esc*1.0e-5

    write(*, *) sum(prof-1d0)



  contains

    function choose_flux(x,y) result(flux)
      real(8) :: x, y
      real(8), dimension(size(emm)) :: flux

      real(8) :: sin_h_squared
      
      flux = 0d0

      if(x**2+y**2 <= 1d0) then
        sin_h_squared = sin(theta)**2 + y**2*cos(2*theta) + x**2*cos(theta)**2 - y*sqrt(1-x**2-y**2)*sin(2*theta)
        ! write(*,*) asin(sqrt(sin_h_squared))/pi*180, asin(sin_h1)/pi*180, asin(sin_h2)/pi*180
        if(sin_h1**2 <= sin_h_squared .and. sin_h_squared <= sin_h2**2) then
          flux = hot_flux
        else
          flux = star_flux
        endif
      endif
    end function choose_flux

    function voigt(a, v) result(W_real)

      implicit none

      real    (kind = 8), intent(in)  :: a, v
      complex (kind = 8) :: W
      complex (kind = 8)              :: z, u
      real    (kind = 8)              :: s, W_real

      z = cmplx(a, -v)
      s = abs(v) + a

      if (s >= 15.0) then

      !* --- Approximation in region I --               -------------- *!

      W = (z * 0.5641896) / (0.5 + (z * z))
      else if (s >= 5.5) then

      !* --- Approximation in region II --              -------------- *!
     
      u = z * z
      W = (z * (1.410474 + u*0.5641896)) / (0.75 + (u*(3.0 + u)))
      else if (a >= 0.195*abs(v) - 0.176) then

      !* --- Approximation in region III --             -------------- *!

        W = (16.4955 + z*(20.20933 + z*(11.96482 + z*(3.778987 + &
          0.5642236*z)))) / &
          (16.4955 + z*(38.82363 + z*(39.27121 + z*(21.69274 + &
          z*(6.699398 + z)))))
      else
      !* --- Approximation in region IV --              -------------- *!

        u = z * z
        W = exp(u) - (z*(36183.31 - u*(3321.99 - u*(1540.787 - &
          u*(219.031 - u*(35.7668 - u*(1.320522 - u*0.56419)))))) / &
          (32066.6 - u*(24322.84 - u*(9022.228 - u*(2186.181 - &
          u*(364.2191 - u*(61.57037 - u*(1.841439 - u))))))))
      endif

      W_real = real(W) 

    end function voigt

    subroutine emission_on_ray(nu, emm_ray, abs_ray)
      real(mp) :: nu, emm_ray, abs_ray

      integer :: p
      real(mp) :: eps
      !real(mp) :: nu

      emm_ray = 0d0
      abs_ray = 0d0

      

      tau = 0d0
      k_lu = 0d0
      S_ul = 0d0

      !write(*, *) k, m_z
      ! OPEN(23, FILE = 'voigt_test.dat', STATUS = 'REPLACE')
      if(k==0) then
        OPEN(1, FILE = 'output/z_data', STATUS = 'REPLACE')
      endif

      do p = m_z, 1, -1
        eps = 5d-6/sqrt(pi)/dv_d(p)
        
        if(stark(p) < 0d0) then
          alpha = exp(-((nu-nu_0)/dv_d(p) - v_z(p)/ut(p))**2)/sqrt(pi)/dv_d(p)
        else
          alpha = voigt(stark(p)/(dv_d(p)*v_c/nu_0**2), (nu-nu_0)/dv_d(p) - v_z(p)/ut(p))/sqrt(pi)/dv_d(p)
        endif
        ! write(23,*) z_coord(p), alpha, exp(-((nu-nu_0)/dv_d(p) - v_z(p)/ut(p))**2)/sqrt(pi)/dv_d(p),&
        !          stark(p)/(dv_d(p)*v_c/nu_0**2)
        !S_ul = flux(k)*n_u(p)
        !if(k == 0) then; s_ul = flux(1)*n_u(p); endif
        if(alpha > eps) then
          ! write(*, *) 'yes'
          !write(*, *) nu
          !alpha = alpha
          !write(*, *) alpha
          k_lu = k_lu_0*n_l(p)*alpha*(1d0 - n_u(p)/n_l(p)*g_l/g_u)
          S_ul = 2d0*h*nu_0**3/v_c**2/(n_l(p)/n_u(p)*g_u/g_l - 1d0)
          !write(*, *) S_ul
          if(argument_list(4, 5) < 0.5d0) then
            emm_ray = emm_ray + S_ul*exp(-tau)*(1.0_mp-exp(-k_lu*dz(p)))
            tau = tau + k_lu*dz(p)
          else if(argument_list(4, 5) < 1.5d0) then 
            emm_ray = emm_ray + S_ul*exp(-tau)*k_lu*dz(p)
            tau = tau + k_lu*dz(p)
          else if(argument_list(4, 5) < 2.5d0) then
            tau = tau + k_lu*dz(p)
            emm_ray = emm_ray + S_ul*exp(-tau)*k_lu*dz(p)
          endif

        endif
       
        if(k == 0) then
          k_lu = k_lu_0*n_l(p)*alpha*(1d0 - n_u(p)/n_l(p)*g_l/g_u)
          S_ul = 2d0*h*nu_0**3/v_c**2/(n_l(p)/n_u(p)*g_u/g_l - 1d0)
          if(abs(abs(z_coord(p)-z_coord(mod(p,m_z)+1))-dz(p)/R_s) > 1d-5) then; write(1,*) 
          endif
          write(1, *) z_coord(p), r(p), -v_z(p)*1d-5, n_u(p), n_l(p), dv_d(p), ut(p)/1d5, &
                      alpha, k_lu, tau, S_ul, emm_ray/flux(1), S_ul/flux(1), &
                      dz(p), S_ul*exp(-tau+k_lu*dz(p))*(1.0_mp-exp(-k_lu*dz(p)))/flux(1), &
                      exp(-((nu-nu_0)/dv_d(p) - v_z(p)/ut(p))**2)/sqrt(pi)/dv_d(p), &
                      stark(p)/(dv_d(p)*v_c/nu_0**2), k_lu*dz(p), tau/0.7/R_s
        endif
      enddo
      abs_ray = exp(-tau)

      ! CLOSE(23)
      ! write(*,*) i,j, m_z
      ! read(*,*)

      if(k == 0) then
        CLOSE(1)
      endif

    end subroutine emission_on_ray

  end subroutine main_subroutine

end module mainmodule
