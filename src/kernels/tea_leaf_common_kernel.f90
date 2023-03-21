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

!>  @brief Controls error reporting
!>  @author David Beckingsale, Wayne Gaudin
!>  @author Douglas Shanks (OpenACC)
!>  @details Common routines.



MODULE tea_leaf_common_kernel_module

  IMPLICIT NONE

  ! These need to be kept consistent with the data module to avoid use statement
  INTEGER,private,PARAMETER :: CHUNK_LEFT   =1    &
                            ,CHUNK_RIGHT  =2    &
                            ,CHUNK_BOTTOM =3    &
                            ,CHUNK_TOP    =4    &
                            ,EXTERNAL_FACE=-1

   INTEGER,PARAMETER        ::   TL_PREC_NONE       = 1 &
                                ,TL_PREC_JAC_DIAG   = 2 &
                                ,TL_PREC_JAC_BLOCK  = 3

   INTEGER,PRIVATE         ::    CONDUCTIVITY        = 1 &
                                ,RECIP_CONDUCTIVITY  = 2

  integer, private, parameter:: jac_block_size = 4

  INTEGER(KIND=4), parameter :: block_size=1
  INTEGER(KIND=4), parameter :: kstep = block_size*jac_block_size

CONTAINS

SUBROUTINE tea_leaf_common_init_kernel(x_min,  &
                           x_max,                  &
                           y_min,                  &
                           y_max,                  &
                           halo_exchange_depth,                  &
                           zero_boundary,       &
                           reflective_boundary,    &
                           density,                &
                           energy,                 &
                           u,                      &
                           u0,                     &
                           r,                      &
                           w,                      &
                           Kx,                     &
                           Ky,                     &
                           Di,                     &
                           cp,                     &
                           bfp,                    &
                           Mi,                     &
                           rx,                     &
                           ry,                     &
                           preconditioner_type,      &
                           coef)

  IMPLICIT NONE

  LOGICAL :: reflective_boundary
  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  LOGICAL, DIMENSION(4) :: zero_boundary
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) &
                          :: density, energy, u, r, w, Kx, Ky, Di, Mi, u0
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp

  INTEGER(KIND=4) :: coef
  INTEGER(KIND=4) :: j,k

  REAL(KIND=8) ::  rx, ry

!$ACC DATA &
!$ACC PRESENT(density, energy, u, r, w, Kx, Ky, Di, Mi, u0, cp, bfp)

!$ACC KERNELS
!$ACC LOOP COLLAPSE(2) INDEPENDENT
  DO k=y_min, y_max
    DO j=x_min, x_max
      u(j,k) = energy(j,k)*density(j,k)
      u0(j,k) = energy(j,k)*density(j,k)
    ENDDO
  ENDDO


  IF (coef .EQ. RECIP_CONDUCTIVITY) THEN
  
!$ACC LOOP COLLAPSE(2) INDEPENDENT
    ! use w as temp val
    DO k=y_min-halo_exchange_depth,y_max+halo_exchange_depth
      DO j=x_min-halo_exchange_depth,x_max+halo_exchange_depth
         w(j  ,k  )=1.0_8/density(j  ,k  )
      ENDDO
    ENDDO

  ELSE IF (coef .EQ. CONDUCTIVITY) THEN
  
!$ACC LOOP COLLAPSE(2) INDEPENDENT
    DO k=y_min-halo_exchange_depth,y_max+halo_exchange_depth
      DO j=x_min-halo_exchange_depth,x_max+halo_exchange_depth
         w(j  ,k  )=density(j  ,k  )
      ENDDO
    ENDDO

  ENDIF

!$ACC LOOP COLLAPSE(2) INDEPENDENT
  DO k=y_min-halo_exchange_depth + 1,y_max+halo_exchange_depth
    DO j=x_min-halo_exchange_depth + 1,x_max+halo_exchange_depth
      Kx(j,k)=(w(j-1,k  ) + w(j,k))/(2.0_8*w(j-1,k  )*w(j,k))
      Ky(j,k)=(w(j  ,k-1) + w(j,k))/(2.0_8*w(j  ,k-1)*w(j,k))
    ENDDO
  ENDDO

! Whether to apply reflective boundary conditions to all external faces
  IF (reflective_boundary .EQV. .FALSE.) THEN
    IF (zero_boundary(CHUNK_LEFT).EQV..TRUE.) THEN
    
!$ACC LOOP COLLAPSE(2) INDEPENDENT
      DO k=y_min-halo_exchange_depth,y_max+halo_exchange_depth
        DO j=x_min-halo_exchange_depth,x_min
          Kx(j,k)=0.0_8
        ENDDO
      ENDDO
    ENDIF
    IF (zero_boundary(CHUNK_RIGHT).EQV..TRUE.) THEN
!$ACC LOOP COLLAPSE(2) INDEPENDENT
      DO k=y_min-halo_exchange_depth,y_max+halo_exchange_depth
        DO j=x_max + 1,x_max+halo_exchange_depth
          Kx(j,k)=0.0_8
        ENDDO
      ENDDO
    ENDIF
    IF (zero_boundary(CHUNK_BOTTOM).EQV..TRUE.) THEN
!$ACC LOOP COLLAPSE(2) INDEPENDENT
      DO k=y_min-halo_exchange_depth,y_min
        DO j=x_min-halo_exchange_depth,x_max+halo_exchange_depth
          Ky(j,k)=0.0_8
        ENDDO
      ENDDO
    ENDIF
    IF (zero_boundary(CHUNK_TOP).EQV..TRUE.) THEN
!$ACC LOOP COLLAPSE(2) INDEPENDENT
      DO k=y_max + 1,y_max+halo_exchange_depth
        DO j=x_min-halo_exchange_depth,x_max+halo_exchange_depth
          Ky(j,k)=0.0_8
        ENDDO
      ENDDO
    ENDIF
  ENDIF

!Setup storage for the diagonal entries
!$ACC LOOP COLLAPSE(2) INDEPENDENT
  DO k=y_min-halo_exchange_depth+1,y_max+halo_exchange_depth-1
    DO j=x_min-halo_exchange_depth+1,x_max+halo_exchange_depth-1
      Di(j,k)=(1.0_8                                              &
                + ry*(Ky(j, k+1) + Ky(j, k))                      &
                + rx*(Kx(j+1, k) + Kx(j, k)))
    ENDDO
  ENDDO
!$ACC END KERNELS
  
  IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
    CALL tea_block_init(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                           cp, bfp, Kx, Ky, Di, rx, ry)
  ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
    CALL tea_diag_init(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                           Mi, Kx, Ky, Di, rx, ry)
  ENDIF

!$ACC KERNELS
!$ACC LOOP COLLAPSE(2) INDEPENDENT
    DO k=y_min,y_max
        DO j=x_min,x_max
            w(j, k) = Di(j,k)*u(j, k)                             &
                - ry*(Ky(j, k+1)*u(j, k+1) + Ky(j, k)*u(j, k-1))  &
                - rx*(Kx(j+1, k)*u(j+1, k) + Kx(j, k)*u(j-1, k))

            r(j, k) = u(j, k) - w(j, k)
            !r(j, k) = u(j, k)! This is required to make a zero initial guess to match petsc errant behaviour
                              ! Only works one timestep is run
        ENDDO
    ENDDO
!$ACC END KERNELS
!$ACC END DATA

END SUBROUTINE tea_leaf_common_init_kernel

SUBROUTINE tea_leaf_kernel_finalise(x_min,    &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           halo_exchange_depth,             &
                           energy,            &
                           density,           &
                           u)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) &
                          :: u, energy, density

  INTEGER(KIND=4) :: j,k

!$ACC DATA &
!$ACC PRESENT(u, energy, density)

!$ACC KERNELS
!$ACC LOOP COLLAPSE(2) INDEPENDENT
  DO k=y_min, y_max
    DO j=x_min, x_max
      energy(j,k) = u(j,k) / density(j,k)
    ENDDO
  ENDDO
!$ACC END KERNELS
!$ACC END DATA

END SUBROUTINE tea_leaf_kernel_finalise

SUBROUTINE tea_leaf_calc_residual_kernel(x_min,       &
                                  x_max,       &
                                  y_min,       &
                                  y_max,       &
                                  halo_exchange_depth,       &
                                  u ,          &
                                  u0,          &
                                  r,           &
                                  Kx,          &
                                  Ky,          &
                                  Di,          &
                                  rx, ry       )

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) &
                          :: Kx, u, r, Ky, u0, Di

  REAL(KIND=8) :: smvp, rx, ry

  INTEGER(KIND=4) :: j,k
  
!$ACC DATA &
!$ACC PRESENT(Kx, u, r, Ky, u0, Di)

!$ACC KERNELS
!$ACC LOOP COLLAPSE(2) INDEPENDENT
    DO k=y_min, y_max
      DO j=x_min, x_max
        smvp = Di(j,k)*u(j, k)                                &
            - ry*(Ky(j, k+1)*u(j, k+1) + Ky(j, k)*u(j, k-1))  &
            - rx*(Kx(j+1, k)*u(j+1, k) + Kx(j, k)*u(j-1, k))
        r(j, k) = u0(j, k) - smvp
      ENDDO
    ENDDO
!$ACC END KERNELS
!$ACC END DATA

END SUBROUTINE tea_leaf_calc_residual_kernel

SUBROUTINE tea_leaf_calc_2norm_kernel(x_min, &
                          x_max,             &
                          y_min,             &
                          y_max,             &
                          halo_exchange_depth,             &
                          arr,               &
                          norm)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: arr
  REAL(KIND=8) :: norm
  integer :: j, k

  norm = 0.0_8


!$ACC DATA &
!$ACC PRESENT(arr)

!$ACC KERNELS
!$ACC LOOP COLLAPSE(2) INDEPENDENT REDUCTION(+:norm) 
    DO k=y_min,y_max
        DO j=x_min,x_max
            norm = norm + arr(j, k)*arr(j, k)
        ENDDO
    ENDDO
!$ACC END KERNELS
!$ACC END DATA

END SUBROUTINE tea_leaf_calc_2norm_kernel

SUBROUTINE tea_diag_init(x_min,             &
                         x_max,             &
                         y_min,             &
                         y_max,             &
                         halo_exchange_depth,             &
                         Mi,                &
                         Kx, Ky, Di, rx, ry)

  IMPLICIT NONE

  INTEGER(KIND=4):: j, k
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) &
                          :: Kx, Ky, Di, Mi
  REAL(KIND=8) :: rx, ry

  REAL(KIND=8), PARAMETER :: omega=1.0_8

!$ACC DATA &
!$ACC PRESENT(Kx, Ky, Di, Mi)

!$ACC KERNELS
!$ACC LOOP COLLAPSE(2) INDEPENDENT
    DO k=y_min-halo_exchange_depth+1,y_max+halo_exchange_depth-1
      DO j=x_min-halo_exchange_depth+1,x_max+halo_exchange_depth-1
        IF (Di(j, k) /= 0.0_8) THEN
          Mi(j, k) = omega/Di(j, k)

        ELSE
          Mi(j, k) = 0.0_8
        ENDIF
      ENDDO
    ENDDO
!$ACC END KERNELS
!$ACC END DATA

END SUBROUTINE

SUBROUTINE tea_diag_solve(x_min,              &
                         x_max,               &
                         y_min,               &
                         y_max,               &
                         halo_exchange_depth, &
                         depth,               &
                         r,                   &
                         z,                   &
                         Mi)

  IMPLICIT NONE

  INTEGER(KIND=4):: j, k
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth,depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) &
                          :: r, z, Mi

!$ACC DATA &
!$ACC PRESENT(r, z, Mi)

!$ACC KERNELS
!$ACC LOOP COLLAPSE(2) INDEPENDENT
    DO k=y_min-depth,y_max+depth
      DO j=x_min-depth,x_max+depth
        z(j, k) = Mi(j, k)*r(j, k)
      ENDDO
    ENDDO
!$ACC END KERNELS
!$ACC END DATA


END SUBROUTINE

SUBROUTINE tea_block_init(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           halo_exchange_depth,             &
                           cp,                     &
                           bfp,                     &
                           Kx, Ky, Di, rx, ry)

  IMPLICIT NONE

  INTEGER(KIND=4):: j, ko, k, bottom, top
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: Kx, Ky, Di
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  REAL(KIND=8) :: rx, ry

!$ACC DATA &
!$ACC PRESENT(Kx, Ky, Di, cp, bfp)

!$ACC KERNELS
!$ACC LOOP INDEPENDENT
    DO ko=y_min,y_max,jac_block_size

      bottom = ko
      top = MIN(ko + jac_block_size - 1, y_max)
!$ACC LOOP INDEPENDENT
      DO j=x_min, x_max
        k = bottom
        cp(j,k) = (-Ky(j, k+1)*ry)/Di(j, k)

        DO k=bottom+1,top
            bfp(j, k) = 1.0_8/(Di(j,k) - (-Ky(j, k)*ry)*cp(j, k-1))
            cp(j, k) = (-Ky(j, k+1)*ry)*bfp(j, k)
        ENDDO
      ENDDO
    ENDDO
!$ACC END KERNELS
!$ACC END DATA

END SUBROUTINE

SUBROUTINE tea_block_solve(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           halo_exchange_depth,             &
                           r,                 &
                           z,                 &
                           cp,                     &
                           bfp,                     &
                           Kx, Ky, Di, rx, ry)

  IMPLICIT NONE

  INTEGER(KIND=4):: j, ko, k, bottom, top, ki, upper_k, k_extra
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: Kx, Ky, Di, r, z
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  REAL(KIND=8) :: rx, ry
  REAL(KIND=8), dimension(0:jac_block_size-1) :: dp_l, z_l

  dp_l = 0.0_8
  z_l = 0.0_8

  k_extra = y_max - MOD(y_max, kstep)

!$ACC DATA &
!$ACC PRESENT(Kx, Ky, Di, r, z, cp, bfp)

!$ACC KERNELS
!$ACC LOOP INDEPENDENT
    DO ko=y_min, k_extra, kstep
      upper_k = ko+kstep - jac_block_size

      DO ki=ko,upper_k,jac_block_size
        bottom = ki
        top = ki+jac_block_size - 1
        
!$ACC LOOP INDEPENDENT PRIVATE(dp_l, z_l)
        DO j=x_min,x_max
          k = bottom
          dp_l(k-bottom) = r(j, k)/Di(j, k)

          DO k=bottom+1,top
            dp_l(k-bottom) = (r(j, k) - (-Ky(j, k)*ry)*dp_l(k-bottom-1))*bfp(j, k)
          ENDDO

          k=top
          z_l(k-bottom) = dp_l(k-bottom)

          DO k=top-1, bottom, -1
            z_l(k-bottom) = dp_l(k-bottom) - cp(j, k)*z_l(k-bottom+1)
          ENDDO

          DO k=bottom,top
            z(j, k) = z_l(k-bottom)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
!$ACC END KERNELS

!$ACC KERNELS
!$ACC LOOP INDEPENDENT
    DO ki=k_extra+1, y_max, jac_block_size
      bottom = MIN(ki, y_max)
      top = MIN(ki+jac_block_size-1, y_max)

!$ACC LOOP INDEPENDENT PRIVATE(dp_l, z_l)
      DO j=x_min,x_max
        k = bottom
        dp_l(k-bottom) = r(j, k)/Di(j, k)

        DO k=bottom+1,top
          dp_l(k-bottom) = (r(j, k) - (-Ky(j, k)*ry)*dp_l(k-bottom-1))*bfp(j, k)
        ENDDO

        k=top
        z_l(k-bottom) = dp_l(k-bottom)

        DO k=top-1, bottom, -1
          z_l(k-bottom) = dp_l(k-bottom) - cp(j, k)*z_l(k-bottom+1)
        ENDDO

        DO k=bottom,top
          z(j, k) = z_l(k-bottom)
        ENDDO
      ENDDO
    ENDDO
!$ACC END KERNELS
!$ACC END DATA

END SUBROUTINE

END MODULE tea_leaf_common_kernel_module
