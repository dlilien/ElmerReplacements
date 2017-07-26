!
! visc2damage.f90
! Copyright (C) 2017 davidl <davidl@ryder>
!
! Distributed under terms of the MIT license.
!

      PROGRAM visc2damage
          use read_routines
          implicit none
          real, dimension(:), allocatable :: x,y
          real, dimension(:, :), allocatable :: visc
          real, dimension(:, :), allocatable :: damage
          integer*4 nx, ny, i, j
          character (len=128) :: fnin
          character (len=128) :: dout='damage.xyz'
          character (len=128) :: dvout='damaged_visc.xyz'

		  CALL getarg(1, fnin)
          CALL get_twod_grid(fnin, x, y, visc)

          nx = SIZE(x)
          ny = SIZE(y)

          allocate (damage(nx, ny))

          do i=1,nx
              do j=1,ny
                  ! Calculate from Borstad equation 1
                  ! This is really easy because of how I defined the
                  ! viscosity
                  if (visc(i, j) ** 2. .LT. 1.0) then
                    damage(i, j) = 1.0 - visc(i, j) ** 2.
                    visc(i, j) = 1.0
                  else
                    damage(i, j) = 0.0
                  end if
              end do
          end do

          CALL write_twod_grid(dout, x, y, damage, nx, ny)
          CALL write_twod_grid(dvout, x, y, visc, nx, ny)

      END PROGRAM visc2damage
