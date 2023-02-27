c ---------------------------------------------------------------------
c
c Copyright (c) 2014 - 2019 by the IBAMR developers
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
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute staggered-grid divergence of the stochastic stress tensor.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_stochastic_stress_div3d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     ifirst2,ilast2,
     &     scale,
     &     n_W_cc_gc0,n_W_cc_gc1,n_W_cc_gc2,
     &     W_cc,
     &     n_W_ec_gc0,n_W_ec_gc1,n_W_ec_gc2,
     &     W_ec0,W_ec1,W_ec2,
     &     n_divW_sc_gc0,n_divW_sc_gc1,n_divW_sc_gc2,
     &     divW_sc0,divW_sc1,divW_sc2)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0
      integer ifirst1,ilast1
      integer ifirst2,ilast2

      integer n_W_cc_gc0,n_W_cc_gc1,n_W_cc_gc2
      integer n_W_ec_gc0,n_W_ec_gc1,n_W_ec_gc2
      integer n_divW_sc_gc0,n_divW_sc_gc1,n_divW_sc_gc2

      double precision dx(0:3-1)

      double precision scale

      double precision W_cc(
     &     ifirst0-n_W_cc_gc0:ilast0+n_W_cc_gc0,
     &          ifirst1-n_W_cc_gc1:ilast1+n_W_cc_gc1,
     &          ifirst2-n_W_cc_gc2:ilast2+n_W_cc_gc2,
     &     0:3-1
     &     )

      double precision W_ec0(
     &     ifirst0-n_W_ec_gc0:ilast0+n_W_ec_gc0,
     &          ifirst1-n_W_ec_gc1:ilast1+1+n_W_ec_gc1,
     &          ifirst2-n_W_ec_gc2:ilast2+1+n_W_ec_gc2,
     &     0:1
     &     )
      double precision W_ec1(
     &     ifirst0-n_W_ec_gc0:ilast0+1+n_W_ec_gc0,
     &          ifirst1-n_W_ec_gc1:ilast1+n_W_ec_gc1,
     &          ifirst2-n_W_ec_gc2:ilast2+1+n_W_ec_gc2,
     &     0:1
     &     )
      double precision W_ec2(
     &     ifirst0-n_W_ec_gc0:ilast0+1+n_W_ec_gc0,
     &          ifirst1-n_W_ec_gc1:ilast1+1+n_W_ec_gc1,
     &          ifirst2-n_W_ec_gc2:ilast2+n_W_ec_gc2,
     &     0:1
     &     )
c
c     Output.
c
      double precision divW_sc0(
     &     ifirst0-n_divW_sc_gc0:ilast0+1+n_divW_sc_gc0,
     &          ifirst1-n_divW_sc_gc1:ilast1+n_divW_sc_gc1,
     &          ifirst2-n_divW_sc_gc2:ilast2+n_divW_sc_gc2
     &     )
      double precision divW_sc1(
     &     ifirst0-n_divW_sc_gc0:ilast0+n_divW_sc_gc0,
     &          ifirst1-n_divW_sc_gc1:ilast1+1+n_divW_sc_gc1,
     &          ifirst2-n_divW_sc_gc2:ilast2+n_divW_sc_gc2
     &     )
      double precision divW_sc2(
     &     ifirst0-n_divW_sc_gc0:ilast0+n_divW_sc_gc0,
     &          ifirst1-n_divW_sc_gc1:ilast1+n_divW_sc_gc1,
     &          ifirst2-n_divW_sc_gc2:ilast2+1+n_divW_sc_gc2
     &     )
c
c     Local variables.
c
      integer i0,i1,i2
c
c     Compute divW_sc0, the x component of div W.
c
      do i2=ifirst2,ilast2
         do i1=ifirst1,ilast1
            do i0=ifirst0,ilast0+1
               divW_sc0(i0,i1,i2) = scale*(
     &              (W_cc (i0,i1  ,i2  ,0)-W_cc (i0-1,i1,i2,0))/dx(0) +
     &              (W_ec2(i0,i1+1,i2  ,0)-W_ec2(i0  ,i1,i2,0))/dx(1) +
     &              (W_ec1(i0,i1  ,i2+1,0)-W_ec1(i0  ,i1,i2,0))/dx(2))
            enddo
         enddo
      enddo
c
c     Compute divW_sc1, the y component of div W.
c
      do i2=ifirst2,ilast2
         do i1=ifirst1,ilast1+1
            do i0=ifirst0,ilast0
               divW_sc1(i0,i1,i2) = scale*(
     &              (W_ec2(i0+1,i1,i2  ,1)-W_ec2(i0,i1  ,i2,1))/dx(0) +
     &              (W_cc (i0  ,i1,i2  ,1)-W_cc (i0,i1-1,i2,1))/dx(1) +
     &              (W_ec0(i0  ,i1,i2+1,0)-W_ec0(i0,i1  ,i2,0))/dx(2))
            enddo
         enddo
      enddo
c
c     Compute divW_sc2, the z component of div W.
c
      do i2=ifirst2,ilast2+1
         do i1=ifirst1,ilast1
            do i0=ifirst0,ilast0
               divW_sc2(i0,i1,i2) = scale*(
     &              (W_ec1(i0+1,i1  ,i2,1)-W_ec1(i0,i1,i2  ,1))/dx(0) +
     &              (W_ec0(i0  ,i1+1,i2,1)-W_ec0(i0,i1,i2  ,1))/dx(1) +
     &              (W_cc (i0  ,i1  ,i2,2)-W_cc (i0,i1,i2-1,2))/dx(2))
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
