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
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/fortran/pdat_m4arrdim2d.i $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for dimensioning 2d arrays in FORTRAN routines.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes cell centered Oldroyd-B type Convective Operator
c
c     where u is vector valued face centered velocity
c     and tau is the symmetric stress tensor cell centered
c     c_data is convective term u.grad(tau)
c     returns s_data
c     computes grad(u) using centered differencing
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine upper_convective_op2d
     &        (dx, u_data_0, u_data_1,
     &        u_gcw, s_data, s_gcw, rhs_data, rhs_gcw,
     &        c_data, c_gcw, r_data, r_gcw,
     &        ilower0, iupper0, ilower1, iupper1)
      implicit none
      integer ilower0, iupper0
      integer ilower1, iupper1
      double precision dx(0:1)

c
c     Velocity Data
c
      integer u_gcw
      double precision u_data_0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
      double precision u_data_1(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+1+u_gcw)
c
c     Tensor Data
c
      integer s_gcw
      double precision s_data(ilower0-s_gcw:iupper0+s_gcw,
     &          ilower1-s_gcw:iupper1+s_gcw,0:2)
c
c     RHS Data
c
      integer rhs_gcw
      double precision rhs_data(ilower0-rhs_gcw:iupper0+rhs_gcw,
     &          ilower1-rhs_gcw:iupper1+rhs_gcw,0:2)
c
c     Return Data
c
      integer r_gcw
      double precision r_data(ilower0-r_gcw:iupper0+r_gcw,
     &          ilower1-r_gcw:iupper1+r_gcw,0:2)
c
c     Convective Data
c
      integer c_gcw
      double precision c_data(ilower0-c_gcw:iupper0+c_gcw,
     &          ilower1-c_gcw:iupper1+c_gcw,0:2)

      integer i0, i1
      double precision du_dx, dv_dx
      double precision du_dy, dv_dy
      double precision scale_ux, scale_uy
      double precision scale_vx, scale_vy
      double precision qxx, qyy, qxy

      scale_vx = 1.d0/(4.d0*dx(0))
      scale_vy = 1.d0/dx(1)
      scale_ux = 1.d0/dx(0)
      scale_uy = 1.d0/(4.d0*dx(1))

      do i1 = ilower1, iupper1
         do i0 = ilower0, iupper0
!           2nd order approximations to derivatives
           du_dx = scale_ux*(u_data_0(i0+1,i1)-u_data_0(i0,i1))
           du_dy = scale_uy*(u_data_0(i0+1,i1+1)+u_data_0(i0,i1+1)
     &              -u_data_0(i0+1,i1-1)-u_data_0(i0,i1-1))
           dv_dy = scale_vy*(u_data_1(i0,i1+1)-u_data_1(i0,i1))
           dv_dx = scale_vx*(u_data_1(i0+1,i1+1)+u_data_1(i0+1,i1)
     &              -u_data_1(i0-1,i1) - u_data_1(i0-1,i1+1))
            qxx = s_data(i0,i1,0)
            qyy = s_data(i0,i1,1)
            qxy = s_data(i0,i1,2)

            r_data(i0,i1,0) = c_data(i0,i1,0)
     &        - 2.d0*du_dx*qxx - 2.d0*du_dy*qxy
     &        - rhs_data(i0,i1,0)
            r_data(i0,i1,1) = c_data(i0,i1,1)
     &        - 2.d0*dv_dx*qxy - 2.d0*dv_dy*qyy
     &        - rhs_data(i0,i1,1)
            r_data(i0,i1,2) = c_data(i0,i1,2)
     &        - qyy*du_dy - qxx*dv_dx
     &        - rhs_data(i0,i1,2)
         enddo
      enddo
      end subroutine
