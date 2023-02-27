c ---------------------------------------------------------------------
c
c Copyright (c) 2019 - 2019 by the IBAMR developers
c All rights reserved.
c
c This file is part of IBAMR.
c
c IBAMR is free software and is distributed under the 3-clause BSD
c license. The full text of the license can be found in the file
c COPYRIGHT at the top level directory of IBAMR.
c
c ---------------------------------------------------------------------

c
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/fortran/pdat_m4arrdim3d.i $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for dimensioning 3d arrays in FORTRAN routines.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes cell centered Oldroyd-B type Convective Operator
c     for the square root evolution equation
c
c     where u is vector valued face centered velocity
c     and tau is symmetric sqrt tensor valued cell centered
c     c_data is u.grad(tau)
c     computes grad(u) using centered differences
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sqrt_upper_convective_op3d
     &        (dx, u_data_0, u_data_1, u_data_2,
     &        u_gcw, s_data, s_gcw, rhs_data, rhs_gcw,
     &        c_data, c_gcw, r_data, r_gcw,
     &        ilower0, iupper0, ilower1,
     &        iupper1, ilower2, iupper2)

      implicit none
      integer ilower0, iupper0
      integer ilower1, iupper1
      integer ilower2, iupper2
      double precision dx(0:2)
c
c    Velocity Data
c
      integer u_gcw
      double precision u_data_0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw,
     &          ilower2-u_gcw:iupper2+u_gcw)
      double precision u_data_1(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+1+u_gcw,
     &          ilower2-u_gcw:iupper2+u_gcw)
      double precision u_data_2(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw,
     &          ilower2-u_gcw:iupper2+1+u_gcw)
c
c    Tensor Data
c
      integer s_gcw
      double precision s_data(ilower0-s_gcw:iupper0+s_gcw,
     &          ilower1-s_gcw:iupper1+s_gcw,
     &          ilower2-s_gcw:iupper2+s_gcw,0:5)
c
c    RHS Data
c
      integer rhs_gcw
      double precision rhs_data(ilower0-rhs_gcw:iupper0+rhs_gcw,
     &          ilower1-rhs_gcw:iupper1+rhs_gcw,
     &          ilower2-rhs_gcw:iupper2+rhs_gcw,0:5)
c
c    Convec Data
c
      integer c_gcw
      double precision c_data(ilower0-c_gcw:iupper0+c_gcw,
     &          ilower1-c_gcw:iupper1+c_gcw,
     &          ilower2-c_gcw:iupper2+c_gcw,0:5)
c
c    Return Data
c
      integer r_gcw
      double precision r_data(ilower0-r_gcw:iupper0+r_gcw,
     &          ilower1-r_gcw:iupper1+r_gcw,
     &          ilower2-r_gcw:iupper2+r_gcw,0:5)

      integer i0, i1, i2
      double precision du_dx, dv_dx, dw_dx
      double precision du_dy, dv_dy, dw_dy
      double precision du_dz, dv_dz, dw_dz
      double precision scale_ux, scale_uy, scale_uz
      double precision scale_vx, scale_vy, scale_vz
      double precision scale_wx, scale_wy, scale_wz
      double precision qxx_ij, qyy_ij, qzz_ij
      double precision qyz_ij, qxz_ij, qxy_ij
      double precision det
      double precision g12, g23, g13

      scale_ux = 1.d0/dx(0)
      scale_uy = 1.d0/(4.d0*dx(1))
      scale_uz = 1.d0/(4.d0*dx(2))
      scale_vx = 1.d0/(4.d0*dx(0))
      scale_vy = 1.d0/dx(1)
      scale_vz = 1.d0/(4.d0*dx(2))
      scale_wx = 1.d0/(4.d0*dx(0))
      scale_wy = 1.d0/(4.d0*dx(1))
      scale_wz = 1.d0/dx(2)

      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0
            du_dx = 
     &        scale_ux*(u_data_0(i0+1,i1,i2)-u_data_0(i0,i1,i2))
            du_dy = 
     &        scale_uy*(u_data_0(i0+1,i1+1,i2)+u_data_0(i0,i1+1,i2)
     &              -u_data_0(i0+1,i1-1,i2)-u_data_0(i0,i1-1,i2))
            du_dz = 
     &        scale_uz*(u_data_0(i0+1,i1,i2+1)+u_data_0(i0,i1,i2+1)
     &              -u_data_0(i0+1,i1,i2-1)-u_data_0(i0,i1,i2-1))
            dv_dy = 
     &        scale_vy*(u_data_1(i0,i1+1,i2)-u_data_1(i0,i1,i2))
            dv_dx = 
     &        scale_vx*(u_data_1(i0+1,i1+1,i2)+u_data_1(i0+1,i1,i2)
     &              -u_data_1(i0-1,i1+1,i2) - u_data_1(i0-1,i1,i2))
            dv_dz = 
     &        scale_vz*(u_data_1(i0,i1+1,i2+1)+u_data_1(i0,i1,i2+1)
     &              -u_data_1(i0,i1+1,i2-1) - u_data_1(i0,i1,i2-1))
            dw_dx = 
     &        scale_wx*(u_data_2(i0+1,i1,i2+1)+u_data_2(i0+1,i1,i2)
     &              -u_data_2(i0-1,i1,i2+1)-u_data_2(i0-1,i1,i2))
            dw_dy = 
     &        scale_wy*(u_data_2(i0,i1+1,i2+1)+u_data_2(i0,i1+1,i2)
     &              -u_data_2(i0,i1-1,i2+1)-u_data_2(i0,i1-1,i2))
            dw_dz = 
     &        scale_wz*(u_data_2(i0,i1,i2+1)-u_data_2(i0,i1,i2))

            qxx_ij = s_data(i0,i1,i2,0)
            qyy_ij = s_data(i0,i1,i2,1)
            qzz_ij = s_data(i0,i1,i2,2)
            qyz_ij = s_data(i0,i1,i2,3)
            qxz_ij = s_data(i0,i1,i2,4)
            qxy_ij = s_data(i0,i1,i2,5)

            det = qxx_ij*qyy_ij*qzz_ij+
     &       2.d0*qxy_ij*qxz_ij*qyz_ij-
     &       qxz_ij**2*qyy_ij-
     &       qyz_ij**2*qxx_ij-
     &       qxy_ij**2*qzz_ij

            call find_G_vals(qxx_ij, qyy_ij, qzz_ij,
     &                       qyz_ij, qxz_ij, qxy_ij,
     &                       du_dx, du_dy, du_dz,
     &                       dv_dx, dv_dy, dv_dz,
     &                       dw_dx, dw_dy, dw_dz,
     &                       g12, g13, g23)

            r_data(i0,i1,i2,0) = 
     &        c_data(i0,i1,i2,0)
     &        - qxx_ij*du_dx - qxy_ij*du_dy - qxz_ij*du_dz
     &        - g12*qxy_ij - g13*qxz_ij
     &        -0.5d0/det*(rhs_data(i0,i1,i2,0)*(
     &        qyy_ij*qzz_ij-qyz_ij**2)+rhs_data(i0,i1,i2,4)*(
     &        qxy_ij*qyz_ij-qxz_ij*qyy_ij)+
     &        rhs_data(i0,i1,i2,5)*(qxz_ij*qyz_ij-qxy_ij*qzz_ij))
            r_data(i0,i1,i2,1) = 
     &        c_data(i0,i1,i2,1)
     &        - qxy_ij*dv_dx - qyy_ij*dv_dy - qyz_ij*dv_dz
     &        + g12*qxy_ij - g23*qyz_ij
     &        -0.5d0/det*(rhs_data(i0,i1,i2,1)*(
     &        qxx_ij*qzz_ij-qxz_ij**2)+rhs_data(i0,i1,i2,3)*(
     &        qxy_ij*qxz_ij-qxx_ij*qyz_ij)+
     &        rhs_data(i0,i1,i2,5)*(qxz_ij*qyz_ij-qxy_ij*qzz_ij))
            r_data(i0,i1,i2,2) = 
     &        c_data(i0,i1,i2,2)
     &        - qxz_ij*dw_dx - qyz_ij*dw_dy - qzz_ij*dw_dz
     &        + g13*qxz_ij + g23*qyz_ij
     &        -0.5d0/det*(rhs_data(i0,i1,i2,2)*(
     &        qxx_ij*qyy_ij-qxy_ij**2)+rhs_data(i0,i1,i2,3)*(
     &        qxy_ij*qxz_ij-qxx_ij*qyz_ij)+
     &        rhs_data(i0,i1,i2,4)*(qxy_ij*qyz_ij-qxz_ij*qyy_ij))
            r_data(i0,i1,i2,3) = 
     &        c_data(i0,i1,i2,3)
     &        - qxy_ij*dw_dx - qyy_ij*dw_dy - qyz_ij*dw_dz
     &        + g12*qxz_ij - g23*qzz_ij
     &        -0.5d0/det*(rhs_data(i0,i1,i2,2)*(
     &        qxy_ij*qxz_ij-qxx_ij*qyz_ij)+rhs_data(i0,i1,i2,3)
     &        *(qxx_ij*qzz_ij-qxz_ij**2)+rhs_data(i0,i1,i2,4)*(
     &        qxz_ij*qyz_ij-qxy_ij*qzz_ij))
            r_data(i0,i1,i2,4) = 
     &        c_data(i0,i1,i2,4)
     &        - qxx_ij*dw_dx - qxy_ij*dw_dy - qxz_ij*dw_dz
     &        - g12*qyz_ij - g13*qzz_ij
     &        -0.5d0/det*(rhs_data(i0,i1,i2,2)*(
     &        qxy_ij*qyz_ij-qxz_ij*qyy_ij)+rhs_data(i0,i1,i2,3)
     &        *(qxz_ij*qyz_ij-qxy_ij*qzz_ij)+rhs_data(i0,i1,i2,4)
     &        *(qyy_ij*qzz_ij-qyz_ij**2))
            r_data(i0,i1,i2,5) = 
     &        c_data(i0,i1,i2,5)
     &        - qxx_ij*dv_dx - qxy_ij*dv_dy - qxz_ij*dv_dz
     &        - g12*qyy_ij - g13*qyz_ij
     &        -0.5d0/det*(rhs_data(i0,i1,i2,1)*(
     &        qxz_ij*qyz_ij-qxy_ij*qzz_ij)+rhs_data(i0,i1,i2,3)
     &        *(qxy_ij*qyz_ij-qxz_ij*qyy_ij)+rhs_data(i0,i1,i2,5)
     &        *(qyy_ij*qzz_ij-qyz_ij**2))
          enddo
        enddo
      enddo
      end subroutine

      subroutine find_G_vals(q11,q22,q33,
     &           q23, q13, q12,
     &           dudx, dudy, dudz,
     &           dvdx, dvdy, dvdz,
     &           dwdx, dwdy, dwdz,
     &           g12, g13, g23)
      implicit none
c     INPUTS
      double precision q11, q22, q33
      double precision q23, q13, q12
      double precision dudx, dudy, dudz
      double precision dvdx, dvdy, dvdz
      double precision dwdx, dwdy, dwdz
c     TO RETURN
      double precision g12, g13, g23
c     SCRATCH
      double precision A(3,3), b(3)
      integer ipiv(3), info

      A(1,1) = q11 + q22; A(1,2) = q23; A(1,3) = -q13
      A(2,1) = q23; A(2,2) = q11 + q33; A(2,3) = q12
      A(3,1) = -q13; A(3,2) = q12; A(3,3) = q22+q33

      b(1) = q12*dudx-q11*dvdx+q22*dudy-q12*dvdy+q23*dudz-q13*dvdz
      b(2) = q13*dudx-q11*dwdx+q33*dudz-q13*dwdz+q23*dudy-q12*dwdy
      b(3) = q13*dvdx-q12*dwdx+q23*dvdy-q22*dwdy+q33*dvdz-q23*dwdz

      call dgesv(3, 1, A, 3, ipiv, b, 3, info)

      if (info /= 0) then
        print *, 'ERROR IN DGESV'
      endif

      g12 = b(1)
      g13 = b(2)
      g23 = b(3)

      endsubroutine
