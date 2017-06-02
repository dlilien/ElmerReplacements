!
! read_routines.f90
! Copyright (C) 2017 dlilien <dlilien@berens>
!
! Distributed under terms of the MIT license.
!

      module read_routines
          contains
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
            r = 6
            allocate (z(nx, ny, nz), x(nx), y(ny))
            do i=1,nx
                read(10, rec=r) x(i)
                r = r + 1
            end do
            do i=1,ny
                read(10, rec=6 + nx + i) y(i)
                r = r + 1
            end do
            do k=1,nz
                do i=1,nx
                    do j=1,ny
                        read(10, rec=r) z(i, j, k)
                        r = r + 1
                    end do
                end do
            end do
         return
         end subroutine get_threed_grid

         subroutine get_twod_grid(fn, x, y, z)
             implicit None
             character (len=128):: fn
             integer nx, ny, r, i, j
             real, dimension(:), allocatable :: x, y
             real, dimension(:, :), allocatable :: z
             open(10, access='DIRECT', recl=4, file=fn)
             read(10, rec=1) ny
             read(10, rec=3) nx
             allocate (z(nx, ny), x(nx), y(ny))
             r = 4
            do i=1,nx
                read(10, rec=r) x(i)
                r = r + 1
            end do
            do i=1,ny
                read(10, rec=6 + nx + i) y(i)
                r = r + 1
            end do
            do i=1,nx
                do j=1,ny
                    read(10, rec=r) z(i, j)
                    r = r + 1
                end do
            end do
         return
         end subroutine get_twod_grid

       end module read_routines


