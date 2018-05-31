!
! read_routines.f90
! Copyright (C) 2017 dlilien <dlilien@berens>
!
! Distributed under terms of the MIT license.
!

      module read_routines
          contains
          subroutine get_dat(fn, x, y, z, dat)
              implicit None
              character (len=128):: fn
              integer i, r, np
              real, dimension(:), allocatable :: x, y, z, dat
              open(10, access='DIRECT', recl=4, file=fn)
              read(10, rec=1) np
              allocate (x(np), y(np), z(np), dat(np))
              r=2
              do i=1,np
                  read(10, rec=r) x(i)
                  r = r+1
              end do
              do i=1,np
                  read(10, rec=r) y(i)
                  r = r+1
              end do
              do i=1,np
                  read(10, rec=r) z(i)
                  r = r+1
              end do
              do i=1,np
                  read(10, rec=r) dat(i)
                  r = r+1
              end do
          end

          subroutine get_z(fn, x, y, z)
              implicit None
              character (len=128):: fn
              integer i, r, np
              real, dimension(:), allocatable :: x, y, z
              open(10, access='DIRECT', recl=4, file=fn)
              read(10, rec=1) np
              allocate (x(np), y(np), z(np))
              r=2
              do i=1,np
                  read(10, rec=r) x(i)
                  r = r+1
              end do
              do i=1,np
                  read(10, rec=r) y(i)
                  r = r+1
              end do
              do i=1,np
                  read(10, rec=r) z(i)
                  r = r+1
              end do
          end

          subroutine get_threed_grid(fn, x, y, z)
            implicit None
            character (len=128):: fn
            integer nx, ny, nz, r, i, j, k
            real, dimension(:), allocatable :: x, y
            real, dimension(:, :, :), allocatable :: z
            open(10, access='DIRECT', recl=4, file=fn)
            read(10, rec=1) ny
            read(10, rec=3) nx
            read(10, rec=5) nz
            r = 7
            allocate (z(nx, ny, nz), x(nx), y(ny))
            do i=1,nx
                read(10, rec=r) x(i)
                r = r + 1
            end do
            do i=1,ny
                read(10, rec=r) y(i)
                r = r + 1
            end do
            do k=1,nz
                do j=1,ny
                    do i=1,nx
                        read(10, rec=r) z(i, j, k)
                        r = r + 1
                    end do
                end do
            end do
         return
         end subroutine get_threed_grid

         subroutine get_twod_grid(fn, x, y, z)
             implicit None
             character (len=128), INTENT(IN):: fn
             integer r, i, j
             integer*4 nx, ny
             real, dimension(:), allocatable, INTENT(OUT) :: x, y
             real, dimension(:, :), allocatable, INTENT(OUT) :: z
             open(10, access='DIRECT', recl=4, file=fn)
             read(10, rec=1) nx
             read(10, rec=2) ny
             allocate (z(nx, ny), x(nx), y(ny))
             r = 3
            do i=1,nx
                read(10, rec=r) x(i)
                r = r + 1
            end do
            do i=1,ny
                read(10, rec=r) y(i)
                r = r + 1
            end do
            do j=1,ny
                do i=1,nx
                    read(10, rec=r) z(i, j)
                    r = r + 1
                end do
            end do
         return
         end subroutine get_twod_grid

         subroutine get_twod_header(fn, x, y, z)
             implicit None
             character (len=128):: fn
             integer r, i, j
             integer*4 nx, ny
             real, dimension(:), allocatable :: x, y
             real, dimension(:, :), allocatable :: z
             open(10, access='DIRECT', recl=4, file=fn)
             read(10, rec=1) ny
             read(10, rec=2) nx
          return
         end subroutine get_twod_header

      SUBROUTINE write_twod_grid(fn, x, y, z, nx, ny)
          implicit none
          integer*4 :: nx, ny
          integer :: i, j, r
          real :: x(nx), y(ny), z(nx, ny)
          character(len=128) :: fn

          open(10, access='DIRECT', recl=4, file=fn)
          write(10, rec=1) nx
          write(10, rec=2) ny
             r = 3
            do i=1,nx
                write(10, rec=r) x(i)
                r = r + 1
            end do
            do i=1,ny
                write(10, rec=r) y(i)
                r = r + 1
            end do
            do j=1,ny
                do i=1,nx
                    write(10, rec=r) z(i, j)
                    r = r + 1
                end do
            end do
         return
      END SUBROUTINE write_twod_grid

        Function LinearInterp(dem,xx,yy,nx,ny,x,y) Result(InterP1)
        USE TYPES
        implicit none

        REAL(KIND=dp) :: dem(nx,ny),xx(nx),yy(ny)
        REAL(KIND=dp) :: Dx,Dy,DxDy
        Real(kind=dp) :: x,y,x_1,y_1,dist,B(4)
        Real(kind=dp) :: InterP1
        integer :: nx,ny,i,j
        integer :: nx_1,ny_1
        logical :: found
        
        IF ((nx <= 1) .OR. (ny <= 1)) THEN
            stop 'rrLI'
        END IF
        Dx=(xx(nx)-xx(1))/(nx-1)
        Dy=(yy(ny)-yy(1))/(ny-1)
        DxDy=Dx*Dy
        IF ((Dx <= 1) .OR. (Dy <= 1)) THEN
            stop 'rrLI'
        END IF
        ! lower left point in DEM
        nx_1=floor((x-xx(1))/Dx) + 1
        ny_1=floor((y-yy(1))/Dy) + 1
        nx_1=min(nx_1,nx-1)
        ny_1=min(ny_1,ny-1)
        nx_1=max(nx_1,1)
        ny_1=max(ny_1,1)

        x_1=xx(nx_1)
        y_1=yy(ny_1)


        ! DEM Value in surroundings points
        !       4 ----- 3
        !       |       |
        !       1 ----- 2
        B(1)=dem(nx_1,ny_1)
        B(2)=dem(nx_1+1,ny_1)
        B(3)=dem(nx_1+1,ny_1+1)
        B(4)=dem(nx_1,ny_1+1)

        found = .false.
        InterP1=0
        if (minval(B)/=-2e+9) then
            ! Linear Interpolation at Point x,y
            InterP1=(x-x_1)*(y-y_1)*(B(3)+B(1)-B(2)-B(4))/DxDy
            InterP1=InterP1+(x-x_1)*(B(2)-B(1))/Dx+(y-y_1)*(B(4)-B(1))/Dy+B(1)
            found = .true.
        else
            do i=0,1
                do j=0,1
                    dist = max( dabs(x-xx(nx_1+i)),dabs(y-yy(ny_1+j)) )
                    if (dist<=0.5*dx .and. dem(nx_1+i,ny_1+j)/=-2e+9) then
                        InterP1 = dem(nx_1+i,ny_1+j)
                        found = .true.
                    endif
                enddo
            enddo
        endif
        Return
        End

        Function LinearInterp_real(dem,xx,yy,nx,ny,x,y) Result(InterP1)
        USE TYPES
        implicit none

        REAL :: dem(nx,ny),xx(nx),yy(ny)
        REAL :: Dx,Dy,DxDy
        Real(kind=dp) :: x,y
        REAL :: x_1,y_1,dist,B(4)
        Real(kind=dp) :: InterP1
        integer :: nx,ny,i,j
        integer :: nx_1,ny_1
        logical :: found

        IF ((nx <= 1) .OR. (ny <= 1)) THEN
            stop 'rrLI'
        END IF
        Dx=(xx(nx)-xx(1))/(nx-1)
        Dy=(yy(ny)-yy(1))/(ny-1)
        DxDy=Dx*Dy
        IF ((Dx <= 1) .OR. (Dy <= 1)) THEN
            stop 'rrLI'
        END IF
        ! lower left point in DEM
        nx_1=floor((x-xx(1))/Dx) + 1
        ny_1=floor((y-yy(1))/Dy) + 1
        nx_1=min(nx_1,nx-1)
        ny_1=min(ny_1,ny-1)
        nx_1=max(nx_1,1)
        ny_1=max(ny_1,1)

        x_1=xx(nx_1)
        y_1=yy(ny_1)


        ! DEM Value in surroundings points
        !       4 ----- 3
        !       |       |
        !       1 ----- 2
        B(1)=dem(nx_1,ny_1)
        B(2)=dem(nx_1+1,ny_1)
        B(3)=dem(nx_1+1,ny_1+1)
        B(4)=dem(nx_1,ny_1+1)

        found = .false.
        InterP1=0
        if (minval(B)/=-2e+9) then
            ! Linear Interpolation at Point x,y
            InterP1=(x-x_1)*(y-y_1)*(B(3)+B(1)-B(2)-B(4))/DxDy
            InterP1=InterP1+(x-x_1)*(B(2)-B(1))/Dx+(y-y_1)*(B(4)-B(1))/Dy+B(1)
            found = .true.
        else
            do i=0,1
                do j=0,1
                    dist = max( dabs(x-xx(nx_1+i)),dabs(y-yy(ny_1+j)) )
                    if (dist<=0.5*dx .and. dem(nx_1+i,ny_1+j)/=-2e+9) then
                        InterP1 = dem(nx_1+i,ny_1+j)
                        found = .true.
                    endif
                enddo
            enddo
        endif
        Return
        End

        Function LID(dem,xx,yy,nx,ny,x,y,def) Result(InterP1)
        USE TYPES
        implicit none
        REAL :: dem(nx,ny),xx(nx),yy(ny)
        REAL :: Dx,Dy,DxDy
        Real(kind=dp) :: x,y
        REAL :: x_1,y_1,dist,B(4)
        Real :: def
        Real(kind=dp) :: InterP1
        integer :: nx,ny,i,j
        integer :: nx_1,ny_1
        logical :: found

        Dx=(xx(nx)-xx(1))/(nx-1)
        Dy=(yy(ny)-yy(1))/(ny-1)
        DxDy=Dx*Dy

        ! lower left point in DEM
        nx_1=floor((x-xx(1))/Dx) + 1
        ny_1=floor((y-yy(1))/Dy) + 1
        nx_1=min(nx_1,nx-1)
        ny_1=min(ny_1,ny-1)
        nx_1=max(nx_1,1)
        ny_1=max(ny_1,1)

        x_1=xx(nx_1)
        y_1=yy(ny_1)


        ! DEM Value in surroundings points
        !       4 ----- 3
        !       |       |
        !       1 ----- 2
        B(1)=dem(nx_1,ny_1)
        B(2)=dem(nx_1+1,ny_1)
        B(3)=dem(nx_1+1,ny_1+1)
        B(4)=dem(nx_1,ny_1+1)

        found = .false.
        InterP1 = def
        if (minval(B)/=-2e+9) then
            ! Linear Interpolation at Point x,y
            InterP1=(x-x_1)*(y-y_1)*(B(3)+B(1)-B(2)-B(4))/DxDy
            InterP1=InterP1+(x-x_1)*(B(2)-B(1))/Dx+(y-y_1)*(B(4)-B(1))/Dy+B(1)
            found = .true.
        else
            do i=0,1
                do j=0,1
                    dist = max( dabs(x-xx(nx_1+i)),dabs(y-yy(ny_1+j)) )
                    if (dist<=0.5*dx .and. dem(nx_1+i,ny_1+j)/=-2e+9) then
                        InterP1 = dem(nx_1+i,ny_1+j)
                        found = .true.
                    endif
                enddo
            enddo
        endif
        Return
        End

       end module read_routines

