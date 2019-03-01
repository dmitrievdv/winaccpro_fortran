module inter2d
  implicit none

  private

  integer :: NX, NY
  real(8) :: XS, XF, YS, YF
  real(8), allocatable, dimension(:) :: X_GRID
  real(8), allocatable, dimension(:) :: Y_GRID

  public :: init, interpolate, deinit

contains

  subroutine init(grid_x, grid_y)
    real(8), dimension(:) :: grid_x, grid_y

    real(8), dimension(1) :: maxmin

    maxmin = maxval(grid_x); XF = maxmin(1)
    maxmin = maxval(grid_y); YF = maxmin(1)
    maxmin = minval(grid_x); XS = maxmin(1)
    maxmin = minval(grid_y); YS = maxmin(1)
    NX = size(grid_x); NY = size(grid_y)

    allocate(X_GRID(0:NX+1), Y_GRID(0:NY+1))
    
    X_GRID(1:NX) = grid_x; Y_GRID(1:NY) = grid_y
    X_GRID(0) = 2d0*X_GRID(1) - X_GRID(2)
    X_GRID(NX+1) = 2d0*X_GRID(NX) - X_GRID(NX-1)
    Y_GRID(0) = 2d0*Y_GRID(1) - X_GRID(2)
    Y_GRID(NY+1) = 2d0*Y_GRID(NY) - X_GRID(NY-1) 

  end subroutine init

  subroutine deinit
    deallocate(X_GRID, Y_GRID)
  end subroutine deinit

  function interpolate_one_spline(x,data)
    real(8), dimension(2,4) :: data
    real(8) :: x,y
    real(8) :: interpolate_one_spline

    real(8), dimension(2) :: m,p
    real(8) :: dp1,dx1,dp2,dx2,h,t
    integer :: i
    real(8), dimension(0:3) :: h00, h01, h10, h11
    real(8), dimension(0:3) :: spl


    do i=1,2
      dp1 = data(2,i+1)-data(2,i); dp2 = data(2,i+2)-data(2,i+1)
      dx1 = data(1,i+1)-data(1,i); dx2 = data(1,i+2)-data(1,i+1)
      ! m(i) = tan((atan(dp1/dx1)+atan(dp2/dx2))/2d0)
      m(i) = (dp1/dx1+dp2/dx2)/2d0
    enddo

    ! write(*,'(4f9.3)') data(1,:)
    ! write(*,'(4f9.3)') data(2,:)
    ! write(*,'(2f9.3)') m
    ! pause

    p(1) = data(2,2); p(2) = data(2,3)

    h00(0:) = (/1d0, 0d0, -3d0, 2d0/)
    h10(0:) = (/0d0, 1d0, -2d0, 1d0/)
    h01(0:) = (/0d0, 0d0, 3d0, -2d0/)
    h11(0:) = (/0d0, 0d0, -1d0, 1d0/)

    h = data(1,3)-data(1,2)
    spl = p(1)*h00 + h*m(1)*h10 + p(2)*h01 + h*m(2)*h11
    !write(*,*) spl
    t = (x-data(1,2))/h
    y = spl(0) + spl(1)*t+spl(2)*t**2 + spl(3)*t**3
    interpolate_one_spline = y

  end function interpolate_one_spline

  function interpolate_segment(segment, grid_x, grid_y, coords) result(int)
    real(8), dimension(-1:2,-1:2) :: segment
    real(8), dimension(4) :: grid_x, grid_y
    real(8), dimension(2) :: coords
    real(8) :: int

    real(8), dimension(2,4) :: data
    real(8), dimension(-1:2) :: y
    integer :: i

    data(1,:) = grid_y

    do i=-1,2
      data(2,:) = segment(i,:)
      y(i) = interpolate_one_spline(coords(2), data)
    enddo

    data(1,:) = grid_x
    data(2,:) = y
    int = interpolate_one_spline(coords(1), data)

  end function interpolate_segment

  function interpolate(x, y, func) result(int)
    real(8) :: x,y,int
    real(8), dimension(NX, NY) :: func

    real(8), dimension(0:NX+1, 0:NY+1) :: func_plus
    real(8), dimension(4,4) :: segment
    real(8), dimension(2) :: seg_coords
    integer :: i,k,j
    real(8) :: hx,hy

    func_plus(1:NX, 1:NY) = func
    func_plus(1:NX, 0) = func(:,1) + 2d0*(func(:,1) - func(:,2))
    func_plus(1:NX, NY+1) = func(:,NY) + 2d0*(func(:,NY) - func(:,NY-1))
    func_plus(0, :) = func_plus(1,:) + 2d0*(func_plus(1,:) - func_plus(2,:))
    func_plus(NX+1, :) = func_plus(NX,:) + 2d0*(func_plus(NX,:) - func_plus(NX-1,:))

    i=1
    do while(X_GRID(i+1) < x)
      if(i > NX) then; exit; endif
      i=i+1
    enddo

    k=1
    do while(Y_GRID(k+1) < y)
      if(i > NX) then; exit; endif
      k=k+1
    enddo

    ! if(i >= NX) then; i = NX-1; endif
    ! if(k >= NY) then; k = NY-1; endif
    ! if(i < 1) then; i = 1; endif
    ! if(k < 1) then; k = 1; endif

    segment(:,:) = func_plus(i-1:i+2, k-1:k+2)

    ! do j=1,4
    !   write(*,'(4f8.5)') segment(j,:)
    ! enddo
    seg_coords = [x, y]
    ! write(*,*) func(:,GRIDY)
    ! write(*,*) segment
    int = interpolate_segment(segment, X_GRID(i-1:i+2), Y_GRID(k-1:k+2), seg_coords)
    ! int = exp(int)
    ! write(*,*) int
    ! write(*,*) int

  end function interpolate  

end module inter2d

module inter2d_mod
  use inter2d, only: &
    inter2d_interpolate => interpolate, &
    inter2d_init => init, &
    inter2d_deinit => deinit

end module inter2d_mod


! program test_int
!   use inter2d_mod
!   implicit none

!   real(8), dimension(3) :: gridx, gridy
!   real(8), dimension(3,3) :: func
!   real(8) :: x,y,f
!   integer :: i,j
!   character(1) :: s

  
!   read(*,*) s,gridx
!   read(*,*) s,gridy
!   read(*,*) s,func(1,:)
!   read(*,*) s,func(2,:)
!   read(*,*) s,func(3,:)

!   ! write(*,*) gridx, gridy

!   call inter2d_init(gridx, gridy)

!   OPEN(1, FILE='tab.dat', STATUS = 'REPLACE')
!   do i=1,3
!     do j=1,3
!       write(1,*) gridx(i), gridy(j), func(i,j)
!     enddo
!   enddo
!   CLOSE(1)

!   do i=-1,11
!     do j=-1,11
!       x = gridx(1) + i*(gridx(3)-gridx(1))/10
!       y = gridy(1) + j*(gridy(3)-gridy(1))/10
!       f = inter2d_interpolate(x, y, func)
!       write(*,*) x,y,f
!       ! pause
!     enddo
!     write(*,*) 
!   enddo

!   call inter2d_deinit
! end program test_int