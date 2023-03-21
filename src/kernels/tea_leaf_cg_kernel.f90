!cROWn Copyright 2014 AWE.
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

!>  @brief Fortran heat conduction kernel
!>  @author Michael Boulton, Wayne Gaudin
!>  @author Douglas Shanks (OpenACC)
!>  @details Implicitly calculates the change in temperature using CG method

MODULE tea_leaf_cg_kernel_module

  USE definitions_module, only: tl_ppcg_active
  USE tea_leaf_common_kernel_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_cg_init_kernel(x_min,  &
                           x_max,                  &
                           y_min,                  &
                           y_max,                  &
                           halo_exchange_depth,                  &
                           p,                      &
                           r,                      &
                           Mi,                     &
                           z,                      &
                           Kx,                     &
                           Ky,                     &
                           Di,                     &
                           cp,                     &
                           bfp,                    &
                           rx,                     &
                           ry,                     &
                           rro,                    &
                           preconditioner_type)

  IMPLICIT NONE
  
  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth,ppcg_inner_iters
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: r, Kx, Ky, Di, z, Mi, p
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp

  INTEGER(KIND=4) :: j,k

  REAL(kind=8) :: rro
  REAL(KIND=8) :: rx, ry

! step 1 is a CG step, whereas steps 2 and 3 are partial steps from the PPCG algorithm to allow a middle step for the PP application

  rro = 0.0_8

!$ACC DATA &
!$ACC PRESENT(r, Kx, Ky, Di, z, Mi, p, cp, bfp)

!$ACC KERNELS
!$ACC LOOP COLLAPSE(2) INDEPENDENT REDUCTION(+:rro)
  DO k=y_min,y_max
    DO j=x_min,x_max
      p(j, k) = 0.0_8
      z(j, k) = 0.0_8
    ENDDO
  ENDDO
!$ACC END KERNELS

  IF (preconditioner_type .NE. TL_PREC_NONE ) THEN

    IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
      CALL tea_block_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                             r, z, cp, bfp, Kx, Ky, Di, rx, ry)
    ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
      CALL tea_diag_solve(x_min, x_max, y_min, y_max, halo_exchange_depth, 0,           &
                             r, z, Mi)
    ENDIF
!$ACC KERNELS
!$ACC LOOP COLLAPSE(2) INDEPENDENT
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k) = z(j, k)
        ENDDO
    ENDDO
!$ACC END KERNELS
    
  ELSE
!$ACC KERNELS          
!$ACC LOOP COLLAPSE(2) INDEPENDENT
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k) = r(j, k)
        ENDDO
    ENDDO
!$ACC END KERNELS
    
  ENDIF
!$ACC KERNELS  
!$ACC LOOP COLLAPSE(2) INDEPENDENT
  DO k=y_min,y_max
    DO j=x_min,x_max
      rro = rro + r(j, k)*p(j, k)
    ENDDO
  ENDDO
!$ACC END KERNELS

!$ACC END DATA

END SUBROUTINE tea_leaf_cg_init_kernel

SUBROUTINE tea_leaf_cg_calc_w_kernel(x_min,             &
                                                   x_max,             &
                                                   y_min,             &
                                                   y_max,             &
                                                   halo_exchange_depth,             &
                                                   p,                 &
                                                   w,                 &
                                                   Kx,                &
                                                   Ky,                &
                                                   Di,                &
                                                   rx,                &
                                                   ry,                &
                                                   pw                 )

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: w, Kx, Ky, p, Di

    REAL(KIND=8) ::  rx, ry

    INTEGER(KIND=4) :: j,k
    REAL(kind=8) :: pw

  pw = 0.0_8

!$ACC DATA &
!$ACC PRESENT(w, Kx, Ky, p, Di)

!$ACC KERNELS
!$ACC LOOP COLLAPSE(2) INDEPENDENT REDUCTION(+:pw)
    DO k=y_min,y_max
        DO j=x_min,x_max
            w(j, k) = Di(j,k)*p(j, k)                             &
                - ry*(Ky(j, k+1)*p(j, k+1) + Ky(j, k)*p(j, k-1))  &
                - rx*(Kx(j+1, k)*p(j+1, k) + Kx(j, k)*p(j-1, k))
            pw = pw + w(j, k)*p(j, k)
        ENDDO
    ENDDO
!$ACC END KERNELS

!$ACC END DATA

END SUBROUTINE tea_leaf_cg_calc_w_kernel

SUBROUTINE tea_leaf_cg_calc_w_kernel_norxy(x_min,             &
                                                   x_max,             &
                                                   y_min,             &
                                                   y_max,             &
                                                   halo_exchange_depth,             &
                                                   p,                 &
                                                   w,                 &
                                                   Kx,                &
                                                   Ky,                &
                                                   Di,                &
                                                   pw                 )

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: w, Kx, Ky, p, Di

    INTEGER(KIND=4) :: j,k
    REAL(kind=8) :: pw

  pw = 0.0_8

!$ACC DATA &
!$ACC PRESENT(w, Kx, Ky, p, Di)

!$ACC KERNELS
!$ACC LOOP COLLAPSE(2) INDEPENDENT REDUCTION(+:pw)
    DO k=y_min,y_max
        DO j=x_min,x_max
            w(j, k) = Di(j,k)*p(j, k)                             &
                - (Ky(j, k+1)*p(j, k+1) + Ky(j, k)*p(j, k-1))  &
                - (Kx(j+1, k)*p(j+1, k) + Kx(j, k)*p(j-1, k))
            pw = pw + w(j, k)*p(j, k)
        ENDDO
    ENDDO
!$ACC END KERNELS
!$ACC END DATA

END SUBROUTINE tea_leaf_cg_calc_w_kernel_norxy

SUBROUTINE tea_leaf_cg_calc_ur_kernel(x_min,             &
                                                    x_max,             &
                                                    y_min,             &
                                                    y_max,             &
                                                    halo_exchange_depth,             &
                                                    u,                 &
                                                    p,                 &
                                                    r,                 &
                                                    Mi,                &
                                                    w,                 &
                                                    z,                 &
                                                    cp,                     &
                                                    bfp,                     &
                                                    Kx, &
                                                    Ky, &
                                                    Di, &
                                                    rx, &
                                                    ry, &
                                                    alpha,             &
                                                    rrn,               &
                                                    preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: u, r, Mi, w, z, Kx, Ky, p, Di
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  REAL(KIND=8) :: rx, ry

  INTEGER(KIND=4) :: j,k
  REAL(kind=8) :: alpha, rrn

  rrn = 0.0_8


!$ACC DATA &
!$ACC PRESENT(u, r, Mi, w, z, Kx, Ky, p, Di, cp, bfp)

  IF (preconditioner_type .NE. TL_PREC_NONE) THEN

    IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
!$ACC KERNELS
!$ACC LOOP COLLAPSE(2) INDEPENDENT
      DO k=y_min,y_max
        DO j=x_min,x_max
          u(j, k) =  u(j, k) + alpha   *p(j, k)
        ENDDO
      ENDDO
!$ACC LOOP COLLAPSE(2) INDEPENDENT REDUCTION(+:rrn)
      DO k=y_min,y_max
        DO j=x_min,x_max
          r(j, k) =  r(j, k) - alpha   *w(j, k)
          z(j, k) =            Mi(j, k)*r(j, k)
          rrn     = rrn      +  r(j, k)*z(j, k)
        ENDDO
      ENDDO
!$ACC END KERNELS

    ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
!$ACC KERNELS
!$ACC LOOP COLLAPSE(2) INDEPENDENT
      DO k=y_min,y_max
        DO j=x_min,x_max
          u(j, k) = u(j, k) + alpha*p(j, k)
        ENDDO
      ENDDO
!$ACC LOOP COLLAPSE(2) INDEPENDENT
      DO k=y_min,y_max
        DO j=x_min,x_max
          r(j, k) = r(j, k) - alpha*w(j, k)
        ENDDO
      ENDDO
!$ACC END KERNELS

      CALL tea_block_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                             r, z, cp, bfp, Kx, Ky, Di, rx, ry)
!$ACC KERNELS
!$ACC LOOP COLLAPSE(2) INDEPENDENT REDUCTION(+:rrn)
      DO k=y_min,y_max
        DO j=x_min,x_max
          rrn = rrn + r(j, k)*z(j, k)
        ENDDO
      ENDDO
!$ACC END KERNELS

    ENDIF
  ELSE
!$ACC KERNELS          
!$ACC LOOP COLLAPSE(2) INDEPENDENT
    DO k=y_min,y_max
      DO j=x_min,x_max
        u(j, k) = u(j, k) + alpha*p(j, k)
      ENDDO
    ENDDO
!$ACC LOOP COLLAPSE(2) INDEPENDENT
    DO k=y_min,y_max
      DO j=x_min,x_max
        r(j, k) = r(j, k) - alpha*w(j, k)
        rrn     = rrn + r(j, k)*r(j, k)
      ENDDO
    ENDDO
!$ACC END KERNELS
    
  ENDIF

!$ACC END DATA

END SUBROUTINE tea_leaf_cg_calc_ur_kernel

SUBROUTINE tea_leaf_cg_calc_p_kernel(x_min,             &
                                                   x_max,             &
                                                   y_min,             &
                                                   y_max,             &
                                                   halo_exchange_depth,             &
                                                   p,                 &
                                                   r,                 &
                                                   z,                 &
                                                   beta,              &
                                                   preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: z, r, p

  INTEGER(KIND=4) :: j,k
  REAL(kind=8) :: beta

!$ACC DATA &
!$ACC PRESENT(z, r, p)

  IF (preconditioner_type .NE. TL_PREC_NONE .or. tl_ppcg_active) THEN
!$ACC KERNELS          
!$ACC LOOP COLLAPSE(2) INDEPENDENT
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k) = z(j, k) + beta*p(j, k)
        ENDDO
    ENDDO
!$ACC END KERNELS    
  ELSE
!$ACC KERNELS          
!$ACC LOOP COLLAPSE(2) INDEPENDENT
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k) = r(j, k) + beta*p(j, k)
        ENDDO
    ENDDO
!$ACC END KERNELS    
  ENDIF

!$ACC END DATA

END SUBROUTINE tea_leaf_cg_calc_p_kernel

END MODULE

