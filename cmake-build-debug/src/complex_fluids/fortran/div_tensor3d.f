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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	Computes d = alpha div s
c
c	where d is vector valued side centered
c       and s is symmetric tensor valued cell centered
c       using centered differencing
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine div_tensor_c_to_s_3d(dx, d_data_0, d_data_1, d_data_2,
     &        d_gcw, s_data, s_gcw, ilower0, iupper0,
     &        ilower1,  iupper1, ilower2, iupper2, alpha)
      implicit none
c     INPUTS
      integer ilower0,  iupper0
      integer iupper1,  ilower1
      integer iupper2, ilower2
      integer s_gcw,  d_gcw

      double precision alpha
c     RETURNS
      double precision d_data_0(ilower0-d_gcw:iupper0+1+d_gcw,
     &          ilower1-d_gcw:iupper1+d_gcw,
     &          ilower2-d_gcw:iupper2+d_gcw)
      double precision d_data_1(ilower0-d_gcw:iupper0+d_gcw,
     &          ilower1-d_gcw:iupper1+1+d_gcw,
     &          ilower2-d_gcw:iupper2+d_gcw)
      double precision d_data_2(ilower0-d_gcw:iupper0+d_gcw,
     &          ilower1-d_gcw:iupper1+d_gcw,
     &          ilower2-d_gcw:iupper2+1+d_gcw)
c     TAU DATA
      double precision s_data(ilower0-s_gcw:iupper0+s_gcw,
     &          ilower1-s_gcw:iupper1+s_gcw,
     &          ilower2-s_gcw:iupper2+s_gcw,0:5)
      double precision dx(0:2)

      integer i0, i1, i2
      double precision scale0_x, scale0_y, scale0_z
      double precision scale1_x, scale1_y, scale1_z
      double precision scale2_x, scale2_y, scale2_z

      scale0_x = alpha/dx(0)
      scale0_y = alpha/(dx(1)*4.d0)
      scale0_z = alpha/(dx(2)*4.d0)
      scale1_y = alpha/dx(1)
      scale1_x = alpha/(dx(0)*4.d0)
      scale1_z = alpha/(dx(2)*4.d0)
      scale2_z = alpha/dx(2)
      scale2_x = alpha/(dx(0)*4.d0)
      scale2_y = alpha/(dx(1)*4.d0)

      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
          do i0 = ilower0, (iupper0+1)
            d_data_0(i0,i1,i2) =
     &       scale0_x*(s_data(i0, i1, i2, 0)-s_data(i0-1, i1, i2, 0))+
     &       scale0_y*(s_data(i0-1, i1+1, i2,5)+s_data(i0, i1+1, i2,5)
     &              -s_data(i0-1, i1-1, i2,5)-s_data(i0, i1-1, i2,5))+
     &       scale0_z*(s_data(i0-1, i1, i2+1,4)+s_data(i0, i1, i2+1,4)
     &              -s_data(i0-1,i1,i2-1,4) - s_data(i0,i1,i2-1,4))
          enddo
        enddo
      enddo
      do i2 = ilower2, iupper2
        do i1 = ilower1, (iupper1+1)
          do i0 = ilower0, iupper0
            d_data_1(i0,i1,i2) =
     &       scale1_y*(s_data(i0, i1, i2,1)-s_data(i0, i1-1, i2,1)) + 
     &       scale1_x*(s_data(i0+1, i1, i2,5)+s_data(i0+1, i1-1, i2,5)
     &              -s_data(i0-1, i1-1, i2,5)-s_data(i0-1, i1, i2,5))+
     &       scale1_z*(s_data(i0, i1,i2+1,3)+s_data(i0,i1-1,i2+1,3)
     &              -s_data(i0,i1,i2-1,3)-s_data(i0,i1-1,i2-1,3))
          enddo
        enddo
      enddo
      do i2 = ilower2, (iupper2+1)
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0
            d_data_2(i0,i1,i2) = 
     &       scale2_x*(s_data(i0+1,i1,i2,4)+s_data(i0+1,i1,i2-1,4)
     &                 -s_data(i0-1,i1,i2,4)-s_data(i0-1,i1,i2-1,4))+
     &       scale2_y*(s_data(i0,i1+1,i2,3)+s_data(i0,i1+1,i2-1,3)
     &                 -s_data(i0,i1-1,i2,3)-s_data(i0,i1-1,i2-1,3))+
     &       scale2_z*(s_data(i0,i1,i2,2)-s_data(i0,i1,i2-1,2))
          enddo
        enddo
      enddo
      end subroutine
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Computes d = alpha div s
c
c       where d is vector valued cell centered
c       and s is symmetric tensor valued cell centered
c       using centered differencing
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine div_tensor_c_to_c_3d(dx, d_data,
     &        d_gcw, s_data, s_gcw, ilower0, iupper0,
     &        ilower1,  iupper1, ilower2, iupper2, alpha)

      implicit none
c     INPUTS
      integer ilower0,  iupper0
      integer iupper1,  ilower1
      integer iupper2, ilower2
      integer s_gcw,  d_gcw

      double precision alpha
c     RETURNS
      double precision d_data(ilower0-d_gcw:iupper0+d_gcw,
     &          ilower1-d_gcw:iupper1+d_gcw,
     &          ilower2-d_gcw:iupper2+d_gcw,0:2)
c     TAU DATA
      double precision s_data(ilower0-s_gcw:iupper0+s_gcw,
     &          ilower1-s_gcw:iupper1+s_gcw,
     &          ilower2-s_gcw:iupper2+s_gcw,0:5)
      double precision dx(0:2)

      integer i0, i1, i2
      double precision scale_x, scale_y, scale_z

      scale_x = alpha/(2.d0*dx(0))
      scale_y = alpha/(2.d0*dx(1))
      scale_z = alpha/(2.d0*dx(2))

      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0
            d_data(i0,i1,i2,0) =
     &        scale_x*(s_data(i0+1,i1,i2,0)-s_data(i0-1,i1,i2,0))+
     &        scale_y*(s_data(i0,i1+1,i2,5)-s_data(i0,i1-1,i2,5))+
     &        scale_z*(s_data(i0,i1,i2+1,4)-s_data(i0,i1,i2-1,4))
            d_data(i0,i1,i2,1) =
     &        scale_y*(s_data(i0,i1+1,i2,1)-s_data(i0,i1-1,i2,1))+
     &        scale_x*(s_data(i0+1,i1,i2,5)-s_data(i0-1,i1,i2,5))+
     &        scale_z*(s_data(i0,i1,i2+1,3)-s_data(i0,i1,i2-1,3))
            d_data(i0,i1,i2,2) =
     &        scale_x*(s_data(i0+1,i1,i2,4)-s_data(i0-1,i1,i2,4))+
     &        scale_y*(s_data(i0,i1+1,i2,3)-s_data(i0,i1-1,i2,3))+
     &        scale_z*(s_data(i0,i1,i2+1,2)-s_data(i0,i1,i2-1,2))
          enddo
        enddo
      enddo
      end subroutine
