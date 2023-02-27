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
c     for the log-evolution equation
c
c     where u is vector valued face centered velocity
c     and tau is symmetric log tensor valued cell centered
c     c_data is convective terms u.grad(tau)
c     computes grad(u) using centered differences
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine log_upper_convective_op3d
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
      double precision qxx, qyy, qzz
      double precision qyz, qxz, qxy

      double precision L(3,3), vecs(3,3)
      double precision vals(3)
      double precision convec_vals(3,3)
      double precision sigma(0:5)
      double precision temp(15)

      integer info

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

            qxx = s_data(i0,i1,i2,0)
            qyy = s_data(i0,i1,i2,1)
            qzz = s_data(i0,i1,i2,2)
            qyz = s_data(i0,i1,i2,3)
            qxz = s_data(i0,i1,i2,4)
            qxy = s_data(i0,i1,i2,5)

            L(1,1) = du_dx; L(1,2) = du_dy; L(1,3) = du_dz
            L(2,1) = dv_dx; L(2,2) = dv_dy; L(2,3) = dv_dz
            L(3,1) = dw_dx; L(3,2) = dw_dy; L(3,3) = dw_dz

            vecs(1,1) = qxx; vecs(1,2) = qxy; vecs(1,3) = qxz
            vecs(2,1) = qxy; vecs(2,2) = qyy; vecs(2,3) = qyz
            vecs(3,1) = qxz; vecs(3,2) = qyz; vecs(3,3) = qzz

            call dsyev('V','U',3,vecs,3,vals,temp,15,info)
            if (info /= 0) then
              print *, "ERROR IN DSYEV!!"
            endif
            call log_sum(L,vecs,vals,3,convec_vals)

            sigma(0) = vecs(1,1)**2*exp(-vals(1))
     &          +vecs(1,2)**2*exp(-vals(2))
     &          +vecs(1,3)**2*exp(-vals(3))
            sigma(1) = vecs(2,1)**2*exp(-vals(1))
     &          +vecs(2,2)**2*exp(-vals(2))
     &          +vecs(2,3)**2*exp(-vals(3))
            sigma(2) = vecs(3,1)**2*exp(-vals(1))
     &          +vecs(3,2)**2*exp(-vals(2))
     &          +vecs(3,3)**2*exp(-vals(3))
            sigma(3) = vecs(2,1)*vecs(3,1)*exp(-vals(1))
     &          +vecs(2,2)*vecs(3,2)*exp(-vals(2))
     &          +vecs(2,3)*vecs(3,3)*exp(-vals(3))
            sigma(4) = vecs(1,1)*vecs(3,1)*exp(-vals(1))
     &          +vecs(1,2)*vecs(3,2)*exp(-vals(2))
     &          +vecs(1,3)*vecs(3,3)*exp(-vals(3))
            sigma(5) = vecs(1,1)*vecs(2,1)*exp(-vals(1))
     &          +vecs(1,2)*vecs(2,2)*exp(-vals(2))
     &          +vecs(1,3)*vecs(2,3)*exp(-vals(3))

            r_data(i0,i1,i2,0) = 
     &        c_data(i0,i1,i2,0)
     &        - convec_vals(1,1) - (
     &        sigma(0)*rhs_data(i0,i1,i2,0) +
     &        sigma(4)*rhs_data(i0,i1,i2,4) +
     &        sigma(5)*rhs_data(i0,i1,i2,5))

            r_data(i0,i1,i2,1) =
     &        c_data(i0,i1,i2,1)
     &        - convec_vals(2,2) - (
     &        sigma(1)*rhs_data(i0,i1,i2,1) +
     &        sigma(3)*rhs_data(i0,i1,i2,3) +
     &        sigma(5)*rhs_data(i0,i1,i2,5))

            r_data(i0,i1,i2,2) =
     &        c_data(i0,i1,i2,2)
     &        - convec_vals(3,3) - (
     &        sigma(2)*rhs_data(i0,i1,i2,2) +
     &        sigma(3)*rhs_data(i0,i1,i2,3) +
     &        sigma(4)*rhs_data(i0,i1,i2,4))

            r_data(i0,i1,i2,3) =
     &        c_data(i0,i1,i2,3)
     &        - convec_vals(2,3) - (
     &        sigma(1)*rhs_data(i0,i1,i2,3) +
     &        sigma(3)*rhs_data(i0,i1,i2,2) +
     &        sigma(5)*rhs_data(i0,i1,i2,4))

            r_data(i0,i1,i2,4) =
     &        c_data(i0,i1,i2,4)
     &        - convec_vals(1,3) - (
     &        sigma(0)*rhs_data(i0,i1,i2,4) +
     &        sigma(4)*rhs_data(i0,i1,i2,2) +
     &        sigma(5)*rhs_data(i0,i1,i2,3))

            r_data(i0,i1,i2,5) =
     &        c_data(i0,i1,i2,5)
     &        - convec_vals(1,2) - (
     &        sigma(0)*rhs_data(i0,i1,i2,5) +
     &        sigma(4)*rhs_data(i0,i1,i2,3) +
     &        sigma(5)*rhs_data(i0,i1,i2,1))
          enddo
        enddo
      enddo
      end subroutine

      subroutine log_sum(L,vecs,vals,d,to_ret)

      implicit none
      integer d
      double precision L(d,d)
      double precision vecs(d,d)
      double precision vals(d)
      double precision to_ret(d,d)
      double precision Lij, Lji

      double precision exp_vals(d)
      integer i, j, k, kk, ii, jj

      do i=1,d
        do j=1,d
          to_ret(i,j) = 0.d0
        enddo
        exp_vals(i) = exp(vals(i))
      enddo

      do k=1,d
        do kk=1,d
          do i=1,d
            Lij = 0
            do ii=1,d
              do jj=1,d
                Lij = Lij + vecs(jj,i)*L(jj,ii)*vecs(ii,i)
              enddo
            enddo
            to_ret(k,kk) = to_ret(k,kk) +
     &        2*Lij*vecs(k,i)*vecs(kk,i)
          enddo
        enddo
      enddo
      do k=1,d
        do kk=1,d
          do i=1,d; do j=1,d
            if (i /= j) then
              Lij = 0.d0; Lji = 0.d0
              do ii=1,d
                do jj=1,d
                  Lij = Lij + vecs(jj,i)*L(jj,ii)*vecs(ii,j)
                  Lji = Lji + vecs(jj,j)*L(jj,ii)*vecs(ii,i)
                enddo
              enddo
              if(abs(vals(k)-vals(kk))<1d-10) then
                  to_ret(k,kk) = to_ret(k,kk) +
     &              (Lji+Lij)*vecs(k,i)*vecs(kk,j)
              else
                  to_ret(k,kk) = to_ret(k,kk) +
     &              (vals(i)-vals(j))/(exp_vals(i)-exp_vals(j))
     &              *(Lij*exp_vals(j)+Lji*exp_vals(i))
     &              *vecs(k,i)*vecs(kk,j)
              endif
            endif
          enddo; enddo
        enddo
      enddo
      endsubroutine
