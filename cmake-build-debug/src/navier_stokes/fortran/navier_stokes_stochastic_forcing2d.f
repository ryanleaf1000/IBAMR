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
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/fortran/pdat_m4arrdim2d.i $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for dimensioning 2d arrays in FORTRAN routines.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute staggered-grid divergence of the stochastic stress tensor.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_stochastic_stress_div2d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     scale,
     &     n_W_cc_gc0,n_W_cc_gc1,
     &     W_cc,
     &     n_W_nc_gc0,n_W_nc_gc1,
     &     W_nc,
     &     n_divW_sc_gc0,n_divW_sc_gc1,
     &     divW_sc0,divW_sc1)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0
      integer ifirst1,ilast1

      integer n_W_cc_gc0,n_W_cc_gc1
      integer n_W_nc_gc0,n_W_nc_gc1
      integer n_divW_sc_gc0,n_divW_sc_gc1

      double precision dx(0:2-1)

      double precision scale

      double precision W_cc(
     &     ifirst0-n_W_cc_gc0:ilast0+n_W_cc_gc0,
     &          ifirst1-n_W_cc_gc1:ilast1+n_W_cc_gc1,
     &     0:2-1
     &     )

      double precision W_nc(
     &     ifirst0-n_W_nc_gc0:ilast0+1+n_W_nc_gc0,
     &          ifirst1-n_W_nc_gc1:ilast1+1+n_W_nc_gc1,
     &     0:1
     &     )
c
c     Output.
c
      double precision divW_sc0(
     &     ifirst0-n_divW_sc_gc0:ilast0+1+n_divW_sc_gc0,
     &          ifirst1-n_divW_sc_gc1:ilast1+n_divW_sc_gc1
     &     )
      double precision divW_sc1(
     &     ifirst0-n_divW_sc_gc0:ilast0+n_divW_sc_gc0,
     &          ifirst1-n_divW_sc_gc1:ilast1+1+n_divW_sc_gc1
     &     )
c
c     Local variables.
c
      integer i0,i1
c
c     Compute divW_sc0, the x component of div W.
c
      do i1=ifirst1,ilast1
         do i0=ifirst0,ilast0+1
            divW_sc0(i0,i1) = scale*(
     &           (W_cc(i0,i1  ,0)-W_cc(i0-1,i1,0))/dx(0) +
     &           (W_nc(i0,i1+1,0)-W_nc(i0  ,i1,0))/dx(1))
         enddo
      enddo
c
c     Compute divW_sc1, the y component of div W.
c
      do i1=ifirst1,ilast1+1
         do i0=ifirst0,ilast0
            divW_sc1(i0,i1) = scale*(
     &           (W_nc(i0+1,i1,1)-W_nc(i0,i1  ,1))/dx(0) +
     &           (W_cc(i0  ,i1,1)-W_cc(i0,i1-1,1))/dx(1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
