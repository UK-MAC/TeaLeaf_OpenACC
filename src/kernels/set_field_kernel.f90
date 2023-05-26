!Crown Copyright 2014 AWE.
!
! This file is part of TeaLeaf.
!
! TeaLeaf is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! TeaLeaf is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! TeaLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Fortran set field kernel.
!>  @author David Beckingsale, Wayne Gaudin
!>  @author Douglas Shanks (OpenACC)
!>  @details Copies all of the final start of step filed data to the end of
!>  step data.

MODULE set_field_kernel_module

CONTAINS

SUBROUTINE set_field_kernel(x_min,x_max,y_min,y_max,halo_exchange_depth,    &
                            energy0,                    &
                            energy1)

  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) &
                          :: energy0, energy1

  INTEGER :: j,k

!$ACC DATA &
!$ACC PRESENT(energy0,energy1)

!$ACC KERNELS
!$ACC LOOP COLLAPSE(2) INDEPENDENT

  DO k=y_min,y_max
     DO j=x_min,x_max
        energy1(j,k)=energy0(j,k)
     ENDDO
  ENDDO

!$ACC END KERNELS

!$ACC END DATA

END SUBROUTINE set_field_kernel

END MODULE set_field_kernel_module
