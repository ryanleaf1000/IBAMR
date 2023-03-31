c ---------------------------------------------------------------------
c
c Copyright (c) 2011 - 2020 by the IBAMR developers
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

c     this is an m4 include, not a Fortran include
c ---------------------------------------------------------------------
c
c Copyright (c) 2011 - 2020 by the IBAMR developers
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the piecewise linear
c     "hat" function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_piecewise_linear_delta(r)
c
      implicit none
      double precision lagrangian_piecewise_linear_delta,r
c
      if ( r.lt.0.d0 ) then
         r = -r
      endif

      if ( r.lt.1.d0 ) then
         lagrangian_piecewise_linear_delta = 1.d0-r
      else
         lagrangian_piecewise_linear_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the broadened 4-point
c     piecewise linear "hat" function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_wide4_piecewise_linear_delta(r)
c
      implicit none
      external lagrangian_piecewise_linear_delta
      double precision lagrangian_piecewise_linear_delta
      double precision lagrangian_wide4_piecewise_linear_delta,r
c
      lagrangian_wide4_piecewise_linear_delta =
     &     0.5d0*lagrangian_piecewise_linear_delta(0.5d0*r)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the piecewise cubic
c     function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_piecewise_cubic_delta(r)
c
      implicit none
      double precision lagrangian_piecewise_cubic_delta,r
c
      if ( r.lt.0.d0 ) then
         r = -r
      endif

      if ( r.lt.1.d0 ) then
         lagrangian_piecewise_cubic_delta =
     &        1.d0-0.5d0*r-r*r+0.5d0*r*r*r
      else if ( r .lt. 2.d0 ) then
         lagrangian_piecewise_cubic_delta =
     &        1.d0-(11.d0/6.d0)*r+r*r-(1.d0/6.d0)*r*r*r
      else
         lagrangian_piecewise_cubic_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the broadened 8-point
c     piecewise cubic function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_wide8_piecewise_cubic_delta(r)
c
      implicit none
      external lagrangian_piecewise_cubic_delta
      double precision lagrangian_piecewise_cubic_delta
      double precision lagrangian_wide8_piecewise_cubic_delta,r
c
      lagrangian_wide8_piecewise_cubic_delta =
     &     0.5d0*lagrangian_piecewise_cubic_delta(0.5d0*r)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the IB 3-point delta
c     function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_ib_3_delta(r)
c
      implicit none
      double precision lagrangian_ib_3_delta,r
      double precision sixth, third
      parameter (sixth=0.16666666666667d0)
      parameter (third=0.333333333333333d0)
c
      if ( r.lt.0.d0 ) then
         r = -r
      endif

      if ( r.lt.0.5d0 ) then
         lagrangian_ib_3_delta =
     &        third*(1.d0+sqrt(1.d0-3.d0*r*r))
      else if ( r.lt.1.5d0 ) then
         lagrangian_ib_3_delta =
     &        sixth*(5.d0-3.d0*r-sqrt(1.d0-3.d0*(1.d0-r)*(1.d0-r)))
      else
         lagrangian_ib_3_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the broadened 6-point
c     version of the IB 3-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_wide6_ib_3_delta(r)
c
      implicit none
      external lagrangian_ib_3_delta
      double precision lagrangian_ib_3_delta
      double precision lagrangian_wide6_ib_3_delta,r
c
      lagrangian_wide6_ib_3_delta = 0.5d0*lagrangian_ib_3_delta(0.5d0*r)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the IB 4-point delta
c     function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_ib_4_delta(r)
c
      implicit none
      double precision lagrangian_ib_4_delta,r,t2,t6
c
      if ( r.lt.0.d0 ) then
         r = -r
      endif

      if ( r.lt.1.d0 ) then
         t2 = r * r
         t6 = sqrt(-0.4D1 * t2 + 0.4D1 * r + 0.1D1)
         lagrangian_ib_4_delta = -r / 0.4D1 + 0.3D1 / 0.8D1 + t6 / 0.8D1
      else if ( r.lt.2.d0 ) then
         t2 = r * r
         t6 = sqrt(0.12D2 * r - 0.7D1 - 0.4D1 * t2)
         lagrangian_ib_4_delta = -r / 0.4D1 + 0.5D1 / 0.8D1 - t6 / 0.8D1
      else
         lagrangian_ib_4_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the broadened 8-point
c     version of the IB 4-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_wide8_ib_4_delta(r)
c
      implicit none
      external lagrangian_ib_4_delta
      double precision lagrangian_ib_4_delta
      double precision lagrangian_wide8_ib_4_delta,r
c
      lagrangian_wide8_ib_4_delta = 0.5d0*lagrangian_ib_4_delta(0.5d0*r)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the broadened 16-point
c     version of the IB 4-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_wide16_ib_4_delta(r)
c
      implicit none
      external lagrangian_ib_4_delta
      double precision lagrangian_ib_4_delta
      double precision lagrangian_wide16_ib_4_delta,r
c
      lagrangian_wide16_ib_4_delta =
     &     0.25d0*lagrangian_ib_4_delta(0.25d0*r)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the 3-point B-spline.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_bspline_3_delta(x)
c
      implicit none
      double precision lagrangian_bspline_3_delta,modx,r,r2,x
c
      modx = dabs(x)
      r = modx + 1.5d0
      r2 = r*r
      if (modx .le. 0.5d0) then
         lagrangian_bspline_3_delta = 0.5d0*(-2.d0*r2 + 6.d0*r - 3.d0)
      elseif (modx .le. 1.5d0) then
         lagrangian_bspline_3_delta = 0.5d0*(r2 - 6.d0*r + 9.d0)
      else
         lagrangian_bspline_3_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the 4-point B-spline.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_bspline_4_delta(x)
c
      implicit none
      double precision lagrangian_bspline_4_delta,modx,r,r2,r3,x
c
      modx = dabs(x)
      r = modx + 2.d0
      r2 = r*r
      r3 = r2*r
      if (modx .le. 1.d0) then
         lagrangian_bspline_4_delta = (1.d0/6.d0)*(3.d0*r3 - 24.d0*r2 +
     &        60.d0*r - 44.d0)
      elseif (modx .le. 2.d0) then
         lagrangian_bspline_4_delta = (1.d0/6.d0)*(-r3 + 12.d0*r2 -
     &        48.d0*r + 64.d0)
      else
         lagrangian_bspline_4_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the 5-point B-spline.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_bspline_5_delta(x)
c
      implicit none
      double precision lagrangian_bspline_5_delta,modx,r,r2,r3,r4,x
c
      modx = dabs(x)
      r = modx + 2.5d0;
      r2 = r*r
      r3 = r2*r
      r4 = r3*r
      if (modx .le. 0.5d0) then
         lagrangian_bspline_5_delta = (1.d0/24.d0)*(6.d0*r4 - 60.d0*r3 +
     &        210.d0*r2 - 300.d0*r + 155.d0)
      elseif (modx .le. 1.5d0) then
         lagrangian_bspline_5_delta = (1.d0/24.d0)*(-4.d0*r4 + 60.d0*r3
     &        - 330.d0*r2 + 780.d0*r - 655.d0)
      elseif (modx .le. 2.5d0) then
         lagrangian_bspline_5_delta = (1.d0/24.d0)*(r4 - 20.d0*r3 +
     &        150.d0*r2 - 500.d0*r + 625.d0)
      else
         lagrangian_bspline_5_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the 6-point B-spline.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_bspline_6_delta(x)
c
      implicit none
      double precision lagrangian_bspline_6_delta,modx,r,r2,r3,r4,r5,x
c
      modx = dabs(x)
      r = modx + 3.d0;
      r2 = r*r
      r3 = r2*r
      r4 = r3*r
      r5 = r4*r
      if (modx .le. 1.d0) then
         lagrangian_bspline_6_delta = (1.d0/60.d0)*(2193.d0 - 3465.d0*r
     &        + 2130.d0*r2 - 630.d0*r3 + 90.d0*r4 - 5.d0*r5)
      elseif (modx .le. 2.d0) then
         lagrangian_bspline_6_delta = (1.d0/120.d0)*(-10974.d0 +
     &        12270.d0*r - 5340.d0*r2 + 1140.d0*r3 - 120.d0*r4 + 5.d0
     &        *r5)
      elseif (modx .le. 3.d0) then
         lagrangian_bspline_6_delta = (1.d0/120.d0)*(7776.d0 - 6480.d0*r
     $        + 2160.d0*r2 - 360.d0*r3 + 30.d0*r4 - r5)
      else
         lagrangian_bspline_6_delta = 0.d0
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Initializes the lookup table for the interpolation weight phi(r)
c     for the IB 4-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_table_init(ib_4_table,NTABLE)
c
      implicit none
      integer NTABLE,k
      double precision lagrangian_ib_4_delta,ib_4_table(0:NTABLE),x
c
      do k = 0,NTABLE
         x = 2.d0*dble(k)/dble(NTABLE)
         ib_4_table(k) = lagrangian_ib_4_delta(x)
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weights f(x,1/2), f(x,3/2), f(x,5/2),
c     and f(x,7/2) for the one-sided IB four point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_one_sided_ib_4_delta(f,x)
c
      implicit none
      double precision f(0:3),N,x
c
      f(0) = 0.75d0 - 0.25d0*x
      N = -4.d0*x*x+16.d0*x-14.d0
      if ( N.gt.0.d0 ) then
         f(0) = f(0) - 0.125d0*sqrt(N)
      endif
      f(1) =  1.5d0 - 0.5d0*x - f(0)
      f(2) =  0.5d0           - f(0)
      f(3) = -1.0d0 + 0.5d0*x + f(0)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weights f(x,3/2), f(x,5/2), f(x,7/2),
c     and f(x,9/2) for the one-sided IB four point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_alt_one_sided_ib_4_delta(f,x)
c
      implicit none
      double precision f(0:3),N,x
c
      f(0) = 1.d0 - 0.25d0*x
      N = -4.d0*x*x+24.d0*x-34.d0
      if ( N.gt.0.d0 ) then
         f(0) = f(0) - 0.125d0*sqrt(N)
      endif
      f(1) =  2.0d0 - 0.5d0*x - f(0)
      f(2) =  0.5d0           - f(0)
      f(3) = -1.5d0 + 0.5d0*x + f(0)
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the interpolation weight phi(r) for the 5-point IB delta
c     with three continuous derivatives
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function lagrangian_ib_5_delta(t)
c
      implicit none
      double precision lagrangian_ib_5_delta,r,x,t,K,phi
c
      x = dabs(t)
      K = (38.0d0 - sqrt(69.0d0))/60.0d0
c
      if (x .le. 0.5d0) then
c
          r = x
c
          phi = (136.0d0 - 40.0d0*K - 40.0d0*r**2 + sqrt(2.0d0)
     &          *sqrt(3123.0d0 - 6840.0d0*K + 3600.0d0*(K**2)
     &          - 12440.0d0*(r**2) + 25680.0d0*K*(r**2)
     &          - 12600.0d0*(K**2)*(r**2) + 8080.0d0*(r**4)
     &          - 8400.0d0*K*(r**4) - 1400.0d0*(r**6)))
     &          /280.0d0
c
         lagrangian_ib_5_delta = phi
c
      elseif (x .le. 1.5d0) then
c
          r = x - 1.0d0
c
          phi = (136.0d0 - 40.0d0*K - 40.0d0*r**2 + sqrt(2.0d0)
     &          *sqrt(3123.0d0 - 6840.0d0*K + 3600.0d0*(K**2)
     &          - 12440.0d0*(r**2) + 25680.0d0*K*(r**2)
     &          - 12600.0d0*(K**2)*(r**2) + 8080.0d0*(r**4)
     &          - 8400.0d0*K*(r**4) - 1400.0d0*(r**6)))
     &          /280.0d0
c
         lagrangian_ib_5_delta = (4.0d0 - 4.0d0*phi - K
     &          - 4.0d0*r + 3.0d0*K*r - r**2 + r**3)/6.0d0
c
      elseif (x .le. 2.5d0) then
c
          r = x - 2.0d0
c
          phi = (136.0d0 - 40.0d0*K - 40.0d0*r**2 + sqrt(2.0d0)
     &          *sqrt(3123.0d0 - 6840.0d0*K + 3600.0d0*(K**2)
     &          - 12440.0d0*(r**2) + 25680.0d0*K*(r**2)
     &          - 12600.0d0*(K**2)*(r**2) + 8080.0d0*(r**4)
     &          - 8400.0d0*K*(r**4) - 1400.0d0*(r**6)))
     &          /280.0d0
c
         lagrangian_ib_5_delta = (-2.0d0 + 2.0d0*phi + 2*K +
     &          r - 3.0d0*K*r + 2.0d0*r**2 - r**3)/12.0d0
c
      else
         lagrangian_ib_5_delta = 0.d0
      endif
c
      return
      end


c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     piecewise constant delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_constant_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Input.
c
      integer depth
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1
      integer nindices

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ic0,ic1
      integer d,l,s
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Use the piecewise constant delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic0 = NINT((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)-0.5d0)+ilower0
         ic1 = NINT((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)-0.5d0)+ilower1
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = u(ic0,ic1,d)
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the
c     piecewise constant delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_constant_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Input.
c
      integer depth
      integer nindices
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ic0,ic1
      integer d,l,s
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Use the piecewise constant delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic0 = NINT((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)-0.5d0)+ilower0
         ic1 = NINT((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)-0.5d0)+ilower1
c
c     Spread V onto u.
c
         do d = 0,depth-1
            u(ic0,ic1,d) = u(ic0,ic1,d) + V(d,s)/(dx(0)*dx(1))
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     discontinuous linear delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_discontinuous_linear_interp2d(
     &     dx,x_lower,x_upper,depth,axis,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Input.
c
      integer depth,axis
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1
      integer nindices

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer ic_trimmed_lower(0:2-1),ic_trimmed_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),X_shifted(0:2-1),w(0:2-1,0:1)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use the discontinuous linear delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Shift the point by X_shift.
c
         do d = 0,2-1
            X_shifted(d) = X(d,s)+Xshift(d,l)
         enddo
c
c     Determine the interpolation stencils and weights.
c
         do d = 0,2-1
            ic_center(d) = ilower(d) +
     &           NINT((X_shifted(d)-x_lower(d))/dx(d)-0.5d0)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
            if ( d.eq.axis ) then
               if ( X_shifted(d).lt.X_cell(d) ) then
                  ic_lower(d) = ic_center(d)-1
                  ic_upper(d) = ic_center(d)
                  w(d,0) = (X_cell(d)-X_shifted(d))/dx(d)
                  w(d,1) = 1.d0 - w(d,0)
               else
                  ic_lower(d) = ic_center(d)
                  ic_upper(d) = ic_center(d)+1
                  w(d,0) = 1.d0 + (X_cell(d)-X_shifted(d))/dx(d)
                  w(d,1) = 1.d0 - w(d,0)
               endif
            else
               w(d,0) = 1.d0
               ic_lower(d) = ic_center(d)
               ic_upper(d) = ic_center(d)
            endif

            ic_trimmed_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_trimmed_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
         enddo
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = 0.d0
            do ic1 = ic_trimmed_lower(1),ic_trimmed_upper(1)
               do ic0 = ic_trimmed_lower(0),ic_trimmed_upper(0)
                  V(d,s) = V(d,s)
     &                 +w(0,ic0-ic_lower(0))
     &                 *w(1,ic1-ic_lower(1))
     &                 *u(ic0,ic1,d)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the
c     discontinuous linear delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_discontinuous_linear_spread2d(
     &     dx,x_lower,x_upper,depth,axis,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Input.
c
      integer depth,axis
      integer nindices
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer ic_trimmed_lower(0:2-1),ic_trimmed_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),X_shifted(0:2-1),w(0:2-1,0:1)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use the discontinuous linear delta function to interpolate u onto
c     V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Shift the point by X_shift.
c
         do d = 0,2-1
            X_shifted(d) = X(d,s)+Xshift(d,l)
         enddo
c
c     Determine the interpolation stencils and weights.
c
         do d = 0,2-1
            ic_center(d) = ilower(d) +
     &           NINT((X_shifted(d)-x_lower(d))/dx(d)-0.5d0)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
            if ( d.eq.axis ) then
               if ( X_shifted(d).lt.X_cell(d) ) then
                  ic_lower(d) = ic_center(d)-1
                  ic_upper(d) = ic_center(d)
                  w(d,0) = (X_cell(d)-X_shifted(d))/dx(d)
                  w(d,1) = 1.d0 - w(d,0)
               else
                  ic_lower(d) = ic_center(d)
                  ic_upper(d) = ic_center(d)+1
                  w(d,0) = 1.d0 + (X_cell(d)-X_shifted(d))/dx(d)
                  w(d,1) = 1.d0 - w(d,0)
               endif
            else
               w(d,0) = 1.d0
               ic_lower(d) = ic_center(d)
               ic_upper(d) = ic_center(d)
            endif

            ic_trimmed_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_trimmed_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
         enddo
c
c     Spread V onto u.
c
         do d = 0,depth-1
            do ic1 = ic_trimmed_lower(1),ic_trimmed_upper(1)
               do ic0 = ic_trimmed_lower(0),ic_trimmed_upper(0)
                  u(ic0,ic1,d) = u(ic0,ic1,d)+(
     &                 w(0,ic0-ic_lower(0))*
     &                 w(1,ic1-ic_lower(1))*
     &                 V(d,s)/(dx(0)*dx(1)))
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     piecewise linear delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_linear_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Input.
c
      integer depth
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1
      integer nindices

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer ic_trimmed_lower(0:2-1),ic_trimmed_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),X_shifted(0:2-1),w(0:2-1,0:1)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use the piecewise linear delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Shift the point by X_shift.
c
         do d = 0,2-1
            X_shifted(d) = X(d,s)+Xshift(d,l)
         enddo
c
c     Determine the interpolation stencils and weights.
c
         do d = 0,2-1
            ic_center(d) = ilower(d) +
     &           NINT((X_shifted(d)-x_lower(d))/dx(d)-0.5d0)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)

            if ( X_shifted(d).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)
               w(d,0) = (X_cell(d)-X_shifted(d))/dx(d)
               w(d,1) = 1.d0 - w(d,0)
            else
               ic_lower(d) = ic_center(d)
               ic_upper(d) = ic_center(d)+1
               w(d,0) = 1.d0 + (X_cell(d)-X_shifted(d))/dx(d)
               w(d,1) = 1.d0 - w(d,0)
            endif

            ic_trimmed_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_trimmed_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
         enddo
c
c     Interpolate u onto V.
c
         do d = 0,depth-1
            V(d,s) = 0.d0
            do ic1 = ic_trimmed_lower(1),ic_trimmed_upper(1)
               do ic0 = ic_trimmed_lower(0),ic_trimmed_upper(0)
                  V(d,s) = V(d,s)
     &                 +w(0,ic0-ic_lower(0))
     &                 *w(1,ic1-ic_lower(1))
     &                 *u(ic0,ic1,d)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the
c     piecewise linear delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_linear_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Input.
c
      integer depth
      integer nindices
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer ic_trimmed_lower(0:2-1),ic_trimmed_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),X_shifted(0:2-1),w(0:2-1,0:1)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use the piecewise linear delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Shift the point by X_shift.
c
         do d = 0,2-1
            X_shifted(d) = X(d,s)+Xshift(d,l)
         enddo
c
c     Determine the interpolation stencils and weights.
c
         do d = 0,2-1
            ic_center(d) = ilower(d) +
     &           NINT((X_shifted(d)-x_lower(d))/dx(d)-0.5d0)
            X_cell(d) = x_lower(d) +
     &           (dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)

            if ( X_shifted(d).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)
               w(d,0) = (X_cell(d)-X_shifted(d))/dx(d)
               w(d,1) = 1.d0 - w(d,0)
            else
               ic_lower(d) = ic_center(d)
               ic_upper(d) = ic_center(d)+1
               w(d,0) = 1.d0 + (X_cell(d)-X_shifted(d))/dx(d)
               w(d,1) = 1.d0 - w(d,0)
            endif

            ic_trimmed_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_trimmed_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
         enddo
c
c     Spread V onto u.
c
         do d = 0,depth-1
            do ic1 = ic_trimmed_lower(1),ic_trimmed_upper(1)
               do ic0 = ic_trimmed_lower(0),ic_trimmed_upper(0)
                  u(ic0,ic1,d) = u(ic0,ic1,d)+(
     &                 w(0,ic0-ic_lower(0))*
     &                 w(1,ic1-ic_lower(1))*
     &                 V(d,s)/(dx(0)*dx(1)))
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     piecewise cubic delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_cubic_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      double precision lagrangian_piecewise_cubic_delta
c
c     Input.
c
      integer depth
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1
      integer nindices

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic,ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),w(0:2-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use the piecewise cubic delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,2-1
c
c     Determine the Cartesian cell in which X(s) is located.
c
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d)
     &           +(dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+1
            else
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)+2
            endif

            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic-ilower(d))+0.5d0)*dx(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_piecewise_cubic_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Interpolate u onto V.
c
            if (ic_upper(1) - ic_lower(1) == (4 - 1) .and.
     &       ic_upper(0) - ic_lower(0) == (4 - 1)) then
            do d = 0,depth-1
               V(d, s) = 0.d0
               do ic1 = ic_lower(1),(ic_lower(1) + 4 - 1)
                  do ic0 = ic_lower(0),(ic_lower(0) + 4 - 1)
                     V(d,s) = V(d,s)
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *u(ic0,ic1,d)
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               V(d, s) = 0.d0
               do ic1 = ic_lower(1),ic_upper(1)
                  do ic0 = ic_lower(0),ic_upper(0)
                     V(d,s) = V(d,s)
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *u(ic0,ic1,d)
                  enddo
               enddo
            enddo
         endif
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the
c     piecewise cubic delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_piecewise_cubic_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      double precision lagrangian_piecewise_cubic_delta
c
c     Input.
c
      integer depth
      integer nindices
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic,ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),w(0:2-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use the piecewise cubic delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d = 0,2-1
c
c     Determine the Cartesian cell in which X(s) is located.
c
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d)
     &           +(dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+1
            else
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)+2
            endif
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic-ilower(d))+0.5d0)*dx(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_piecewise_cubic_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Spread V onto u.
c
            if (ic_upper(1) - ic_lower(1) == (4 - 1) .and.
     &       ic_upper(0) - ic_lower(0) == (4 - 1)) then
            do d = 0,depth-1
               do ic1 = ic_lower(1),(ic_lower(1) + 4 - 1)
                  do ic0 = ic_lower(0),(ic_lower(0) + 4 - 1)
                     u(ic0,ic1,d) = u(ic0,ic1,d) + (
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *V(d,s)/(dx(0)*dx(1)))
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               do ic1 = ic_lower(1),ic_upper(1)
                  do ic0 = ic_lower(0),ic_upper(0)
                     u(ic0,ic1,d) = u(ic0,ic1,d) + (
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *V(d,s)/(dx(0)*dx(1)))
                  enddo
               enddo
            enddo
         endif
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the IB
c     3-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_3_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      double precision lagrangian_ib_3_delta
c
c     Input.
c
      integer depth
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1
      integer nindices

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic,ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),w(0:2-1,0:2)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use the IB 3-point delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d = 0,2-1
c
c     Determine the Cartesian cell in which X(s) is located.
c
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d)
     &           +(dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic-ilower(d))+0.5d0)*dx(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_ib_3_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Interpolate u onto V.
c
            if (ic_upper(1) - ic_lower(1) == (3 - 1) .and.
     &       ic_upper(0) - ic_lower(0) == (3 - 1)) then
            do d = 0,depth-1
               V(d, s) = 0.d0
               do ic1 = ic_lower(1),(ic_lower(1) + 3 - 1)
                  do ic0 = ic_lower(0),(ic_lower(0) + 3 - 1)
                     V(d,s) = V(d,s)
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *u(ic0,ic1,d)
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               V(d, s) = 0.d0
               do ic1 = ic_lower(1),ic_upper(1)
                  do ic0 = ic_lower(0),ic_upper(0)
                     V(d,s) = V(d,s)
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *u(ic0,ic1,d)
                  enddo
               enddo
            enddo
         endif
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the IB
c     3-point delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_3_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      double precision lagrangian_ib_3_delta
c
c     Input.
c
      integer depth
      integer nindices
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic,ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),w(0:2-1,0:2)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use the IB 3-point delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d = 0,2-1
c
c     Determine the Cartesian cell in which X(s) is located.
c
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d)
     &           +(dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic-ilower(d))+0.5d0)*dx(d)
               w(d, ic-ic_lower(d)) =
     &              lagrangian_ib_3_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Spread V onto u.
c
            if (ic_upper(1) - ic_lower(1) == (3 - 1) .and.
     &       ic_upper(0) - ic_lower(0) == (3 - 1)) then
            do d = 0,depth-1
               do ic1 = ic_lower(1),(ic_lower(1) + 3 - 1)
                  do ic0 = ic_lower(0),(ic_lower(0) + 3 - 1)
                     u(ic0,ic1,d) = u(ic0,ic1,d) + (
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *V(d,s)/(dx(0)*dx(1)))
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               do ic1 = ic_lower(1),ic_upper(1)
                  do ic0 = ic_lower(0),ic_upper(0)
                     u(ic0,ic1,d) = u(ic0,ic1,d) + (
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *V(d,s)/(dx(0)*dx(1)))
                  enddo
               enddo
            enddo
         endif
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the IB
c     4-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Input.
c
      integer depth
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1
      integer nindices

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer i0,i1,ic0,ic1
      integer ig_lower(0:2-1),ig_upper(0:2-1)
      integer ic_lower(0:2-1),ic_upper(0:2-1)
      integer istart0,istop0,istart1,istop1
      integer d,l,s

      double precision X_o_dx,q0,q1,r0,r1
      double precision w0(0:3),w1(0:3)
      double precision w(0:3,0:3),wy
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ilower0-nugc0
      ig_lower(1) = ilower1-nugc1
      ig_upper(0) = iupper0+nugc0
      ig_upper(1) = iupper1+nugc1
c
c     Use the IB 4-point delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell and compute the interpolation weights.
c
         X_o_dx = (X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)
         ic_lower(0) = NINT(X_o_dx)+ilower0-2
         ic_upper(0) = ic_lower(0) + 3
         r0 = X_o_dx - ((ic_lower(0)+1-ilower0)+0.5d0)
         q0 = sqrt(1.d0+4.d0*r0*(1.d0-r0))
         w0(0) = 0.125d0*(3.d0-2.d0*r0-q0)
         w0(1) = 0.125d0*(3.d0-2.d0*r0+q0)
         w0(2) = 0.125d0*(1.d0+2.d0*r0+q0)
         w0(3) = 0.125d0*(1.d0+2.d0*r0-q0)

         X_o_dx = (X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)
         ic_lower(1) = NINT(X_o_dx)+ilower1-2
         ic_upper(1) = ic_lower(1) + 3
         r1 = X_o_dx - ((ic_lower(1)+1-ilower1)+0.5d0)
         q1 = sqrt(1.d0+4.d0*r1*(1.d0-r1))
         w1(0) = 0.125d0*(3.d0-2.d0*r1-q1)
         w1(1) = 0.125d0*(3.d0-2.d0*r1+q1)
         w1(2) = 0.125d0*(1.d0+2.d0*r1+q1)
         w1(3) = 0.125d0*(1.d0+2.d0*r1-q1)
c
c     Compute the tensor product of the interpolation weights.
c
         do i1 = 0,3
            wy = w1(i1)
            do i0 = 0,3
               w(i0,i1) = w0(i0)*wy
            enddo
         enddo
c
c     Interpolate u onto V.
c
         istart0 =   max(ig_lower(0)-ic_lower(0),0)
         istop0  = 3-max(ic_upper(0)-ig_upper(0),0)
         istart1 =   max(ig_lower(1)-ic_lower(1),0)
         istop1  = 3-max(ic_upper(1)-ig_upper(1),0)
         do d = 0,depth-1
            V(d,s) = 0.d0
            do i1 = istart1,istop1
               ic1 = ic_lower(1)+i1
               do i0 = istart0,istop0
                  ic0 = ic_lower(0)+i0
                  V(d,s) = V(d,s) + w(i0,i1)*u(ic0,ic1,d)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the IB
c     4-point delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Input.
c
      integer depth
      integer nindices
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer i0,i1,ic0,ic1
      integer ig_lower(0:2-1),ig_upper(0:2-1)
      integer ic_lower(0:2-1),ic_upper(0:2-1)
      integer istart0,istop0,istart1,istop1
      integer d,l,s

      double precision X_o_dx,q0,q1,r0,r1
      double precision w0(0:3),w1(0:3)
      double precision w(0:3,0:3),wy
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ilower0-nugc0
      ig_lower(1) = ilower1-nugc1
      ig_upper(0) = iupper0+nugc0
      ig_upper(1) = iupper1+nugc1
c
c     Use the IB 4-point delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell and compute the interpolation weights.
c
         X_o_dx = (X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)
         ic_lower(0) = NINT(X_o_dx)+ilower0-2
         ic_upper(0) = ic_lower(0) + 3
         r0 = X_o_dx - ((ic_lower(0)+1-ilower0)+0.5d0)
         q0 = sqrt(1.d0+4.d0*r0*(1.d0-r0))
         w0(0) = 0.125d0*(3.d0-2.d0*r0-q0)
         w0(1) = 0.125d0*(3.d0-2.d0*r0+q0)
         w0(2) = 0.125d0*(1.d0+2.d0*r0+q0)
         w0(3) = 0.125d0*(1.d0+2.d0*r0-q0)

         X_o_dx = (X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)
         ic_lower(1) = NINT(X_o_dx)+ilower1-2
         ic_upper(1) = ic_lower(1) + 3
         r1 = X_o_dx - ((ic_lower(1)+1-ilower1)+0.5d0)
         q1 = sqrt(1.d0+4.d0*r1*(1.d0-r1))
         w1(0) = 0.125d0*(3.d0-2.d0*r1-q1)
         w1(1) = 0.125d0*(3.d0-2.d0*r1+q1)
         w1(2) = 0.125d0*(1.d0+2.d0*r1+q1)
         w1(3) = 0.125d0*(1.d0+2.d0*r1-q1)
c
c     Compute the tensor product of the scaled interpolation weights.
c
         do i1 = 0,3
            wy = w1(i1)/(dx(0)*dx(1))
            do i0 = 0,3
               w(i0,i1) = w0(i0)*wy
            enddo
         enddo
c
c     Spread V onto u.
c
         istart0 =   max(ig_lower(0)-ic_lower(0),0)
         istop0  = 3-max(ic_upper(0)-ig_upper(0),0)
         istart1 =   max(ig_lower(1)-ic_lower(1),0)
         istop1  = 3-max(ic_upper(1)-ig_upper(1),0)
         do d = 0,depth-1
            do i1 = istart1,istop1
               ic1 = ic_lower(1)+i1
               do i0 = istart0,istop0
                  ic0 = ic_lower(0)+i0
                  u(ic0,ic1,d) = u(ic0,ic1,d) + w(i0,i1)*V(d,s)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the IB
c     4-point delta function that has been broadened to have a support
c     of 8 meshwidths.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_w8_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Input.
c
      integer depth
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1
      integer nindices

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer i0,i1,ic0,ic1
      integer ig_lower(0:2-1),ig_upper(0:2-1)
      integer ic_lower(0:2-1),ic_upper(0:2-1)
      integer istart0,istop0,istart1,istop1
      integer d,l,s

      double precision X_o_dx,q0,q1,r0,r1
      double precision w0(0:7),w1(0:7)
      double precision w(0:7,0:7),wy
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ilower0-nugc0
      ig_lower(1) = ilower1-nugc1
      ig_upper(0) = iupper0+nugc0
      ig_upper(1) = iupper1+nugc1
c
c     Use a broadened version of the IB 4-point delta function to
c     interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell and compute the interpolation weights.
c
         X_o_dx = (X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)
         ic_lower(0) = NINT(X_o_dx)+ilower0-4
         ic_upper(0) = ic_lower(0) + 7
         r0 = 0.5d0*(X_o_dx - ((ic_lower(0)+3-ilower0)+0.5d0))
         q0 = sqrt(1.d0+4.d0*r0*(1.d0-r0))
         w0(1) = 0.0625d0*(3.d0-2.d0*r0-q0)
         w0(3) = 0.0625d0*(3.d0-2.d0*r0+q0)
         w0(5) = 0.0625d0*(1.d0+2.d0*r0+q0)
         w0(7) = 0.0625d0*(1.d0+2.d0*r0-q0)
         r0 = r0+0.5d0
         q0 = sqrt(1.d0+4.d0*r0*(1.d0-r0))
         w0(0) = 0.0625d0*(3.d0-2.d0*r0-q0)
         w0(2) = 0.0625d0*(3.d0-2.d0*r0+q0)
         w0(4) = 0.0625d0*(1.d0+2.d0*r0+q0)
         w0(6) = 0.0625d0*(1.d0+2.d0*r0-q0)

         X_o_dx = (X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)
         ic_lower(1) = NINT(X_o_dx)+ilower1-4
         ic_upper(1) = ic_lower(1) + 7
         r1 = 0.5d0*(X_o_dx - ((ic_lower(1)+3-ilower1)+0.5d0))
         q1 = sqrt(1.d0+4.d0*r1*(1.d0-r1))
         w1(1) = 0.0625d0*(3.d0-2.d0*r1-q1)
         w1(3) = 0.0625d0*(3.d0-2.d0*r1+q1)
         w1(5) = 0.0625d0*(1.d0+2.d0*r1+q1)
         w1(7) = 0.0625d0*(1.d0+2.d0*r1-q1)
         r1 = r1+0.5d0
         q1 = sqrt(1.d0+4.d0*r1*(1.d0-r1))
         w1(0) = 0.0625d0*(3.d0-2.d0*r1-q1)
         w1(2) = 0.0625d0*(3.d0-2.d0*r1+q1)
         w1(4) = 0.0625d0*(1.d0+2.d0*r1+q1)
         w1(6) = 0.0625d0*(1.d0+2.d0*r1-q1)
c
c     Compute the tensor product of the interpolation weights.
c
         do i1 = 0,7
            wy = w1(i1)
            do i0 = 0,7
               w(i0,i1) = w0(i0)*wy
            enddo
         enddo
c
c     Interpolate u onto V.
c
         istart0 =   max(ig_lower(0)-ic_lower(0),0)
         istop0  = 7-max(ic_upper(0)-ig_upper(0),0)
         istart1 =   max(ig_lower(1)-ic_lower(1),0)
         istop1  = 7-max(ic_upper(1)-ig_upper(1),0)
         do d = 0,depth-1
            V(d,s) = 0.d0
            do i1 = istart1,istop1
               ic1 = ic_lower(1)+i1
               do i0 = istart0,istop0
                  ic0 = ic_lower(0)+i0
                  V(d,s) = V(d,s) + w(i0,i1)*u(ic0,ic1,d)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the IB
c     4-point delta function that has been broadened to have a support
c     of 8 meshwidths.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_4_w8_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Input.
c
      integer depth
      integer nindices
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer i0,i1,ic0,ic1
      integer ig_lower(0:2-1),ig_upper(0:2-1)
      integer ic_lower(0:2-1),ic_upper(0:2-1)
      integer istart0,istop0,istart1,istop1
      integer d,l,s

      double precision X_o_dx,q0,q1,r0,r1
      double precision w0(0:7),w1(0:7)
      double precision w(0:7,0:7),wy
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ilower0-nugc0
      ig_lower(1) = ilower1-nugc1
      ig_upper(0) = iupper0+nugc0
      ig_upper(1) = iupper1+nugc1
c
c     Use a broadened version of the IB 4-point delta function to spread
c     V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell and compute the interpolation weights.
c
         X_o_dx = (X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)
         ic_lower(0) = NINT(X_o_dx)+ilower0-4
         ic_upper(0) = ic_lower(0) + 7
         r0 = 0.5d0*(X_o_dx - ((ic_lower(0)+3-ilower0)+0.5d0))
         q0 = sqrt(1.d0+4.d0*r0*(1.d0-r0))
         w0(1) = 0.0625d0*(3.d0-2.d0*r0-q0)
         w0(3) = 0.0625d0*(3.d0-2.d0*r0+q0)
         w0(5) = 0.0625d0*(1.d0+2.d0*r0+q0)
         w0(7) = 0.0625d0*(1.d0+2.d0*r0-q0)
         r0 = r0+0.5d0
         q0 = sqrt(1.d0+4.d0*r0*(1.d0-r0))
         w0(0) = 0.0625d0*(3.d0-2.d0*r0-q0)
         w0(2) = 0.0625d0*(3.d0-2.d0*r0+q0)
         w0(4) = 0.0625d0*(1.d0+2.d0*r0+q0)
         w0(6) = 0.0625d0*(1.d0+2.d0*r0-q0)

         X_o_dx = (X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)
         ic_lower(1) = NINT(X_o_dx)+ilower1-4
         ic_upper(1) = ic_lower(1) + 7
         r1 = 0.5d0*(X_o_dx - ((ic_lower(1)+3-ilower1)+0.5d0))
         q1 = sqrt(1.d0+4.d0*r1*(1.d0-r1))
         w1(1) = 0.0625d0*(3.d0-2.d0*r1-q1)
         w1(3) = 0.0625d0*(3.d0-2.d0*r1+q1)
         w1(5) = 0.0625d0*(1.d0+2.d0*r1+q1)
         w1(7) = 0.0625d0*(1.d0+2.d0*r1-q1)
         r1 = r1+0.5d0
         q1 = sqrt(1.d0+4.d0*r1*(1.d0-r1))
         w1(0) = 0.0625d0*(3.d0-2.d0*r1-q1)
         w1(2) = 0.0625d0*(3.d0-2.d0*r1+q1)
         w1(4) = 0.0625d0*(1.d0+2.d0*r1+q1)
         w1(6) = 0.0625d0*(1.d0+2.d0*r1-q1)
c
c     Compute the tensor product of the scaled interpolation weights.
c
         do i1 = 0,7
            wy = w1(i1)/(dx(0)*dx(1))
            do i0 = 0,7
               w(i0,i1) = w0(i0)*wy
            enddo
         enddo
c
c     Spread V onto u.
c
         istart0 =   max(ig_lower(0)-ic_lower(0),0)
         istop0  = 7-max(ic_upper(0)-ig_upper(0),0)
         istart1 =   max(ig_lower(1)-ic_lower(1),0)
         istop1  = 7-max(ic_upper(1)-ig_upper(1),0)
         do d = 0,depth-1
            do i1 = istart1,istop1
               ic1 = ic_lower(1)+i1
               do i0 = istart0,istop0
                  ic0 = ic_lower(0)+i0
                  u(ic0,ic1,d) = u(ic0,ic1,d) + w(i0,i1)*V(d,s)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     5-point IB delta with three continuous derivatives
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_5_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Input.
c
      integer depth
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1
      integer nindices

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)

c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      double precision r
      double precision phi
      double precision K
      integer i0,i1,ic0,ic1
      integer ig_lower(0:2-1),ig_upper(0:2-1)
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer istart0,istop0,istart1,istop1
      integer d,l,s

      double precision X_cell(0:2-1),w0(0:4),w1(0:4)

      PARAMETER (K = (38.0d0 - sqrt(69.0d0))/60.0d0)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ilower0-nugc0
      ig_lower(1) = ilower1-nugc1
      ig_upper(0) = iupper0+nugc0
      ig_upper(1) = iupper1+nugc1
c
c     Use a 5-point IB delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic_center(0) =
     &        floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ilower0
         ic_center(1) =
     &        floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ilower1

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ilower0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ilower1)+0.5d0)*dx(1)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c
         do d = 0,2-1
            ic_lower(d) = ic_center(d)-2
            ic_upper(d) = ic_center(d)+2
         enddo
c
c     Compute the interpolation weights.
c
         ic0 = ic_center(0)
         X_cell(0) = x_lower(0)+(dble(ic0-ilower0)+0.5d0)*dx(0)
         r = (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0)
         phi = (136.0d0 - 40.0d0*K - 40.0d0*r**2 + sqrt(2.0d0)
     &          *sqrt(3123.0d0 - 6840.0d0*K + 3600.0d0*(K**2)
     &          - 12440.0d0*(r**2) + 25680.0d0*K*(r**2)
     &          - 12600.0d0*(K**2)*(r**2) + 8080.0d0*(r**4)
     &          - 8400.0d0*K*(r**4) - 1400.0d0*(r**6)))
     &          /280.0d0

         w0(0) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K +
     &                         r - 3.0d0*K*r + 2.0d0*r**2 - r**3)
         w0(1) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K -
     &                   4.0d0*r + 3.0d0*K*r -       r**2 + r**3)
         w0(2) = phi
         w0(3) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K +
     &                   4.0d0*r - 3.0d0*K*r -       r**2 - r**3)
         w0(4) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K -
     &                         r + 3.0d0*K*r + 2.0d0*r**2 + r**3)

         ic1 = ic_center(1)
         X_cell(1) = x_lower(1)+(dble(ic1-ilower1)+0.5d0)*dx(1)
         r = (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1)
         phi = (136.0d0 - 40.0d0*K - 40.0d0*r**2 + sqrt(2.0d0)
     &          *sqrt(3123.0d0 - 6840.0d0*K + 3600.0d0*(K**2)
     &          - 12440.0d0*(r**2) + 25680.0d0*K*(r**2)
     &          - 12600.0d0*(K**2)*(r**2) + 8080.0d0*(r**4)
     &          - 8400.0d0*K*(r**4) - 1400.0d0*(r**6)))
     &          /280.0d0

         w1(0) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K +
     &                         r - 3.0d0*K*r + 2.0d0*r**2 - r**3)
         w1(1) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K -
     &                   4.0d0*r + 3.0d0*K*r -       r**2 + r**3)
         w1(2) = phi
         w1(3) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K +
     &                   4.0d0*r - 3.0d0*K*r -       r**2 - r**3)
         w1(4) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K -
     &                         r + 3.0d0*K*r + 2.0d0*r**2 + r**3)

c
c     Interpolate u onto V.
c
         istart0 =   max(ig_lower(0)-ic_lower(0),0)
         istop0  = 4-max(ic_upper(0)-ig_upper(0),0)
         istart1 =   max(ig_lower(1)-ic_lower(1),0)
         istop1  = 4-max(ic_upper(1)-ig_upper(1),0)
         do d = 0,depth-1
            V(d,s) = 0.d0
            do i1 = istart1,istop1
               ic1 = ic_lower(1)+i1
               do i0 = istart0,istop0
                  ic0 = ic_lower(0)+i0
                  V(d,s) = V(d,s) + w0(i0) * w1(i1) * u(ic0,ic1,d)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the
c     5-point IB delta with three continuous derivatives
c     using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_5_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Input.
c
      integer depth
      integer nindices
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      double precision r
      double precision phi
      double precision K
      integer i0,i1,ic0,ic1
      integer ig_lower(0:2-1),ig_upper(0:2-1)
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer istart0,istop0,istart1,istop1
      integer d,l,s

      double precision X_cell(0:2-1),w0(0:4),w1(0:4)

      PARAMETER (K = (38.0d0 - sqrt(69.0d0))/60.0d0)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ilower0-nugc0
      ig_lower(1) = ilower1-nugc1
      ig_upper(0) = iupper0+nugc0
      ig_upper(1) = iupper1+nugc1
c
c     Use a 5-point IB delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the Cartesian cell in which X(s) is located.
c
         ic_center(0) =
     &        floor((X(0,s)+Xshift(0,l)-x_lower(0))/dx(0))
     &        + ilower0
         ic_center(1) =
     &        floor((X(1,s)+Xshift(1,l)-x_lower(1))/dx(1))
     &        + ilower1

         X_cell(0) = x_lower(0)+(dble(ic_center(0)-ilower0)+0.5d0)*dx(0)
         X_cell(1) = x_lower(1)+(dble(ic_center(1)-ilower1)+0.5d0)*dx(1)
c
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c
         do d = 0,2-1
            ic_lower(d) = ic_center(d)-2
            ic_upper(d) = ic_center(d)+2
         enddo
c
c     Compute the spreading weights.
c

         ic0 = ic_center(0)
         X_cell(0) = x_lower(0)+(dble(ic0-ilower0)+0.5d0)*dx(0)
         r = (X(0,s)+Xshift(0,l)-X_cell(0))/dx(0)
         phi = (136.0d0 - 40.0d0*K - 40.0d0*r**2 + sqrt(2.0d0)
     &          *sqrt(3123.0d0 - 6840.0d0*K + 3600.0d0*(K**2)
     &          - 12440.0d0*(r**2) + 25680.0d0*K*(r**2)
     &          - 12600.0d0*(K**2)*(r**2) + 8080.0d0*(r**4)
     &          - 8400.0d0*K*(r**4) - 1400.0d0*(r**6)))
     &          /280.0d0

         w0(0) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K +
     &                         r - 3.0d0*K*r + 2.0d0*r**2 - r**3)
         w0(1) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K -
     &                   4.0d0*r + 3.0d0*K*r -       r**2 + r**3)
         w0(2) = phi
         w0(3) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K +
     &                   4.0d0*r - 3.0d0*K*r -       r**2 - r**3)
         w0(4) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K -
     &                         r + 3.0d0*K*r + 2.0d0*r**2 + r**3)

         ic1 = ic_center(1)
         X_cell(1) = x_lower(1)+(dble(ic1-ilower1)+0.5d0)*dx(1)
         r = (X(1,s)+Xshift(1,l)-X_cell(1))/dx(1)
         phi = (136.0d0 - 40.0d0*K - 40.0d0*r**2 + sqrt(2.0d0)
     &          *sqrt(3123.0d0 - 6840.0d0*K + 3600.0d0*(K**2)
     &          - 12440.0d0*(r**2) + 25680.0d0*K*(r**2)
     &          - 12600.0d0*(K**2)*(r**2) + 8080.0d0*(r**4)
     &          - 8400.0d0*K*(r**4) - 1400.0d0*(r**6)))
     &          /280.0d0

         w1(0) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K +
     &                         r - 3.0d0*K*r + 2.0d0*r**2 - r**3)
         w1(1) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K -
     &                   4.0d0*r + 3.0d0*K*r -       r**2 + r**3)
         w1(2) = phi
         w1(3) = (1.0d0/ 6.0d0) * ( 4.0d0 - 4.0d0*phi -       K +
     &                   4.0d0*r - 3.0d0*K*r -       r**2 - r**3)
         w1(4) = (1.0d0/12.0d0) * (-2.0d0 + 2.0d0*phi + 2.0d0*K -
     &                         r + 3.0d0*K*r + 2.0d0*r**2 + r**3)

c
c     Spread V onto u.
c
         istart0 =   max(ig_lower(0)-ic_lower(0),0)
         istop0  = 4-max(ic_upper(0)-ig_upper(0),0)
         istart1 =   max(ig_lower(1)-ic_lower(1),0)
         istop1  = 4-max(ic_upper(1)-ig_upper(1),0)
         do d = 0,depth-1
            do i1 = istart1,istop1
               ic1 = ic_lower(1)+i1
               do i0 = istart0,istop0
                  ic0 = ic_lower(0)+i0
                  u(ic0,ic1,d) = u(ic0,ic1,d)+(
     &                 w0(i0) * w1(i1) * V(d,s)/(dx(0)*dx(1)))
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the IB
c     6-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_6_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Input.
c
      integer depth
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1
      integer nindices

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer i0,i1,ic0,ic1
      integer ig_lower(0:2-1),ig_upper(0:2-1)
      integer ic_lower(0:2-1),ic_upper(0:2-1)
      integer istart0,istop0,istart1,istop1
      integer d,l,s

      double precision X_o_dx,r,alpha,beta,gamma,discr,K
      double precision pm3,pm2,pm1,p,pp1,pp2
      double precision w0(0:5),w1(0:5)
      double precision w(0:5,0:5),wy

      PARAMETER (K = (59.d0/60.d0)*(1.d0-sqrt(1.d0-(3220.d0/3481.d0))))
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ilower0-nugc0
      ig_lower(1) = ilower1-nugc1
      ig_upper(0) = iupper0+nugc0
      ig_upper(1) = iupper1+nugc1
c
c     Use the IB 6-point delta function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell and compute the interpolation weights.
c
         X_o_dx = (X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)
         ic_lower(0) = NINT(X_o_dx)+ilower0-3
         ic_upper(0) = ic_lower(0) + 5
         r = 1.d0 - X_o_dx + ((ic_lower(0)+2-ilower0)+0.5d0)

         alpha = 28.d0
         beta = (9.d0/4.d0)-(3.d0/2.d0)*(K+r**2)+((22.d0/3.d0)-7.d0*K)*r
     $        -(7.d0/3.d0)*r**3
         gamma = (1.d0/4.d0)*( ((161.d0/36.d0)-(59.d0/6.d0)*K+5.d0*K**2)
     $        *(1.d0/2.d0)*r**2 + (-(109.d0/24.d0)+5.d0*K)*(1.d0/3.d0)*r
     $        **4+ (5.d0/18.d0)*r**6 )
         discr = beta**2-4.d0*alpha*gamma

         pm3 = (-beta+sign(1.d0,(3.d0/2.d0)-K)*sqrt(discr))/(2.d0*alpha)
         pm2 =  -3.d0*pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) + (1.d0
     $        /12.d0)*(3.d0*K-1.d0)*r + (1.d0/12.d0)*r**3
         pm1 =   2.d0*pm3 + (1.d0/4.d0) +  (1.d0/6.d0)*(4.d0-3.d0*K)*r -
     $        (1.d0/6.d0)*r**3
         p   =   2.d0*pm3 + (5.d0/8.d0)  - (1.d0/4.d0)*(K+r**2)
         pp1 =  -3.d0*pm3 + (1.d0/4.d0) -  (1.d0/6.d0)*(4.d0-3.d0*K)*r +
     $        (1.d0/6.d0)*r**3
         pp2 =        pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) - (1.d0
     $        /12.d0)*(3.d0*K-1.d0)*r - (1.d0/12.d0)*r**3

         w0(0) = pm3
         w0(1) = pm2
         w0(2) = pm1
         w0(3) = p
         w0(4) = pp1
         w0(5) = pp2

         X_o_dx = (X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)
         ic_lower(1) = NINT(X_o_dx)+ilower1-3
         ic_upper(1) = ic_lower(1) + 5
         r = 1.d0 - X_o_dx + ((ic_lower(1)+2-ilower1)+0.5d0)

         alpha = 28.d0
         beta = (9.d0/4.d0)-(3.d0/2.d0)*(K+r**2)+((22.d0/3.d0)-7.d0*K)*r
     $        -(7.d0/3.d0)*r**3
         gamma = (1.d0/4.d0)*( ((161.d0/36.d0)-(59.d0/6.d0)*K+5.d0*K**2)
     $        *(1.d0/2.d0)*r**2 + (-(109.d0/24.d0)+5.d0*K)*(1.d0/3.d0)*r
     $        **4+ (5.d0/18.d0)*r**6 )
         discr = beta**2-4.d0*alpha*gamma

         pm3 = (-beta+sign(1.d0,(3.d0/2.d0)-K)*sqrt(discr))/(2.d0*alpha)
         pm2 =  -3.d0*pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) + (1.d0
     $        /12.d0)*(3.d0*K-1.d0)*r + (1.d0/12.d0)*r**3
         pm1 =   2.d0*pm3 + (1.d0/4.d0) +  (1.d0/6.d0)*(4.d0-3.d0*K)*r -
     $        (1.d0/6.d0)*r**3
         p   =   2.d0*pm3 + (5.d0/8.d0)  - (1.d0/4.d0)*(K+r**2)
         pp1 =  -3.d0*pm3 + (1.d0/4.d0) -  (1.d0/6.d0)*(4.d0-3.d0*K)*r +
     $        (1.d0/6.d0)*r**3
         pp2 =        pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) - (1.d0
     $        /12.d0)*(3.d0*K-1.d0)*r - (1.d0/12.d0)*r**3

         w1(0) = pm3
         w1(1) = pm2
         w1(2) = pm1
         w1(3) = p
         w1(4) = pp1
         w1(5) = pp2
c
c     Compute the tensor product of the interpolation weights.
c
         do i1 = 0,5
            wy = w1(i1)
            do i0 = 0,5
               w(i0,i1) = w0(i0)*wy
            enddo
         enddo
c
c     Interpolate u onto V.
c
         istart0 =   max(ig_lower(0)-ic_lower(0),0)
         istop0  = 5-max(ic_upper(0)-ig_upper(0),0)
         istart1 =   max(ig_lower(1)-ic_lower(1),0)
         istop1  = 5-max(ic_upper(1)-ig_upper(1),0)
         do d = 0,depth-1
            V(d,s) = 0.d0
            do i1 = istart1,istop1
               ic1 = ic_lower(1)+i1
               do i0 = istart0,istop0
                  ic0 = ic_lower(0)+i0
                  V(d,s) = V(d,s) + w(i0,i1)*u(ic0,ic1,d)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the IB
c     6-point delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_ib_6_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Input.
c
      integer depth
      integer nindices
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer i0,i1,ic0,ic1
      integer ig_lower(0:2-1),ig_upper(0:2-1)
      integer ic_lower(0:2-1),ic_upper(0:2-1)
      integer istart0,istop0,istart1,istop1
      integer d,l,s

      double precision X_o_dx,r,alpha,beta,gamma,discr,K
      double precision pm3,pm2,pm1,p,pp1,pp2
      double precision w0(0:5),w1(0:5)
      double precision w(0:5,0:5),wy

      PARAMETER (K = (59.d0/60.d0)*(1.d0-sqrt(1.d0-(3220.d0/3481.d0))))
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Compute the extents of the ghost box.
c
      ig_lower(0) = ilower0-nugc0
      ig_lower(1) = ilower1-nugc1
      ig_upper(0) = iupper0+nugc0
      ig_upper(1) = iupper1+nugc1
c
c     Use the IB 6-point delta function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell and compute the interpolation weights.
c
         X_o_dx = (X(0,s)+Xshift(0,l)-x_lower(0))/dx(0)
         ic_lower(0) = NINT(X_o_dx)+ilower0-3
         ic_upper(0) = ic_lower(0) + 5
         r = 1.d0 - X_o_dx + ((ic_lower(0)+2-ilower0)+0.5d0)

         alpha = 28.d0
         beta = (9.d0/4.d0)-(3.d0/2.d0)*(K+r**2)+((22.d0/3.d0)-7.d0*K)*r
     $        -(7.d0/3.d0)*r**3
         gamma = (1.d0/4.d0)*( ((161.d0/36.d0)-(59.d0/6.d0)*K+5.d0*K**2)
     $        *(1.d0/2.d0)*r**2 + (-(109.d0/24.d0)+5.d0*K)*(1.d0/3.d0)*r
     $        **4+ (5.d0/18.d0)*r**6 )
         discr = beta**2-4.d0*alpha*gamma

         pm3 = (-beta+sign(1.d0,(3.d0/2.d0)-K)*sqrt(discr))/(2.d0*alpha)
         pm2 =  -3.d0*pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) + (1.d0
     $        /12.d0)*(3.d0*K-1.d0)*r + (1.d0/12.d0)*r**3
         pm1 =   2.d0*pm3 + (1.d0/4.d0) +  (1.d0/6.d0)*(4.d0-3.d0*K)*r -
     $        (1.d0/6.d0)*r**3
         p   =   2.d0*pm3 + (5.d0/8.d0)  - (1.d0/4.d0)*(K+r**2)
         pp1 =  -3.d0*pm3 + (1.d0/4.d0) -  (1.d0/6.d0)*(4.d0-3.d0*K)*r +
     $        (1.d0/6.d0)*r**3
         pp2 =        pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) - (1.d0
     $        /12.d0)*(3.d0*K-1.d0)*r - (1.d0/12.d0)*r**3

         w0(0) = pm3
         w0(1) = pm2
         w0(2) = pm1
         w0(3) = p
         w0(4) = pp1
         w0(5) = pp2

         X_o_dx = (X(1,s)+Xshift(1,l)-x_lower(1))/dx(1)
         ic_lower(1) = NINT(X_o_dx)+ilower1-3
         ic_upper(1) = ic_lower(1) + 5
         r = 1.d0 - X_o_dx + ((ic_lower(1)+2-ilower1)+0.5d0)

         alpha = 28.d0
         beta = (9.d0/4.d0)-(3.d0/2.d0)*(K+r**2)+((22.d0/3.d0)-7.d0*K)*r
     $        -(7.d0/3.d0)*r**3
         gamma = (1.d0/4.d0)*( ((161.d0/36.d0)-(59.d0/6.d0)*K+5.d0*K**2)
     $        *(1.d0/2.d0)*r**2 + (-(109.d0/24.d0)+5.d0*K)*(1.d0/3.d0)*r
     $        **4+ (5.d0/18.d0)*r**6 )
         discr = beta**2-4.d0*alpha*gamma

         pm3 = (-beta+sign(1.d0,(3.d0/2.d0)-K)*sqrt(discr))/(2.d0*alpha)
         pm2 =  -3.d0*pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) + (1.d0
     $        /12.d0)*(3.d0*K-1.d0)*r + (1.d0/12.d0)*r**3
         pm1 =   2.d0*pm3 + (1.d0/4.d0) +  (1.d0/6.d0)*(4.d0-3.d0*K)*r -
     $        (1.d0/6.d0)*r**3
         p   =   2.d0*pm3 + (5.d0/8.d0)  - (1.d0/4.d0)*(K+r**2)
         pp1 =  -3.d0*pm3 + (1.d0/4.d0) -  (1.d0/6.d0)*(4.d0-3.d0*K)*r +
     $        (1.d0/6.d0)*r**3
         pp2 =        pm3 - (1.d0/16.d0) + (1.d0/8.d0)*(K+r**2) - (1.d0
     $        /12.d0)*(3.d0*K-1.d0)*r - (1.d0/12.d0)*r**3

         w1(0) = pm3
         w1(1) = pm2
         w1(2) = pm1
         w1(3) = p
         w1(4) = pp1
         w1(5) = pp2
c
c     Compute the tensor product of the scaled interpolation weights.
c
         do i1 = 0,5
            wy = w1(i1)/(dx(0)*dx(1))
            do i0 = 0,5
               w(i0,i1) = w0(i0)*wy
            enddo
         enddo
c
c     Spread V onto u.
c
         istart0 =   max(ig_lower(0)-ic_lower(0),0)
         istop0  = 5-max(ic_upper(0)-ig_upper(0),0)
         istart1 =   max(ig_lower(1)-ic_lower(1),0)
         istop1  = 5-max(ic_upper(1)-ig_upper(1),0)
         do d = 0,depth-1
            do i1 = istart1,istop1
               ic1 = ic_lower(1)+i1
               do i0 = istart0,istop0
                  ic0 = ic_lower(0)+i0
                  u(ic0,ic1,d) = u(ic0,ic1,d) + w(i0,i1)*V(d,s)
               enddo
            enddo
         enddo
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     3-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_3_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      double precision lagrangian_bspline_3_delta
c
c     Input.
c
      integer depth
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1
      integer nindices

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic,ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),w(0:2-1,0:2)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a 3-point B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,2-1
c
c     Determine the Cartesian cell in which X(s) is located.
c
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d)
     &           +(dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic-ilower(d))+0.5d0)*dx(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_bspline_3_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Interpolate u onto V.
c
            if (ic_upper(1) - ic_lower(1) == (3 - 1) .and.
     &       ic_upper(0) - ic_lower(0) == (3 - 1)) then
            do d = 0,depth-1
               V(d, s) = 0.d0
               do ic1 = ic_lower(1),(ic_lower(1) + 3 - 1)
                  do ic0 = ic_lower(0),(ic_lower(0) + 3 - 1)
                     V(d,s) = V(d,s)
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *u(ic0,ic1,d)
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               V(d, s) = 0.d0
               do ic1 = ic_lower(1),ic_upper(1)
                  do ic0 = ic_lower(0),ic_upper(0)
                     V(d,s) = V(d,s)
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *u(ic0,ic1,d)
                  enddo
               enddo
            enddo
         endif
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the 3-point
c     B-spline delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_3_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      double precision lagrangian_bspline_3_delta
c
c     Input.
c
      integer depth
      integer nindices
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic,ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),w(0:2-1,0:2)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a 3-point B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,2-1
c
c     Determine the Cartesian cell in which X(s) is located.
c
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d)
     &           +(dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic-ilower(d))+0.5d0)*dx(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_bspline_3_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Spread V onto u.
c
            if (ic_upper(1) - ic_lower(1) == (3 - 1) .and.
     &       ic_upper(0) - ic_lower(0) == (3 - 1)) then
            do d = 0,depth-1
               do ic1 = ic_lower(1),(ic_lower(1) + 3 - 1)
                  do ic0 = ic_lower(0),(ic_lower(0) + 3 - 1)
                     u(ic0,ic1,d) = u(ic0,ic1,d) + (
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *V(d,s)/(dx(0)*dx(1)))
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               do ic1 = ic_lower(1),ic_upper(1)
                  do ic0 = ic_lower(0),ic_upper(0)
                     u(ic0,ic1,d) = u(ic0,ic1,d) + (
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *V(d,s)/(dx(0)*dx(1)))
                  enddo
               enddo
            enddo
         endif
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     4-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_4_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      double precision lagrangian_bspline_4_delta
c
c     Input.
c
      integer depth
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1
      integer nindices

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic,ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),w(0:2-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a 4-point B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,2-1
c
c     Determine the Cartesian cell in which X(s) is located.
c
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d)
     &           +(dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+1
            else
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)+2
            endif
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic-ilower(d))+0.5d0)*dx(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_bspline_4_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Interpolate u onto V.
c
            if (ic_upper(1) - ic_lower(1) == (4 - 1) .and.
     &       ic_upper(0) - ic_lower(0) == (4 - 1)) then
            do d = 0,depth-1
               V(d, s) = 0.d0
               do ic1 = ic_lower(1),(ic_lower(1) + 4 - 1)
                  do ic0 = ic_lower(0),(ic_lower(0) + 4 - 1)
                     V(d,s) = V(d,s)
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *u(ic0,ic1,d)
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               V(d, s) = 0.d0
               do ic1 = ic_lower(1),ic_upper(1)
                  do ic0 = ic_lower(0),ic_upper(0)
                     V(d,s) = V(d,s)
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *u(ic0,ic1,d)
                  enddo
               enddo
            enddo
         endif
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the 4-point
c     B-spline delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_4_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      double precision lagrangian_bspline_4_delta
c
c     Input.
c
      integer depth
      integer nindices
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic,ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),w(0:2-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a 4-point B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,2-1
c
c     Determine the Cartesian cell in which X(s) is located.
c
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d)
     &           +(dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+1
            else
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)+2
            endif
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic-ilower(d))+0.5d0)*dx(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_bspline_4_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Spread V onto u.
c
            if (ic_upper(1) - ic_lower(1) == (4 - 1) .and.
     &       ic_upper(0) - ic_lower(0) == (4 - 1)) then
            do d = 0,depth-1
               do ic1 = ic_lower(1),(ic_lower(1) + 4 - 1)
                  do ic0 = ic_lower(0),(ic_lower(0) + 4 - 1)
                     u(ic0,ic1,d) = u(ic0,ic1,d) + (
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *V(d,s)/(dx(0)*dx(1)))
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               do ic1 = ic_lower(1),ic_upper(1)
                  do ic0 = ic_lower(0),ic_upper(0)
                     u(ic0,ic1,d) = u(ic0,ic1,d) + (
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *V(d,s)/(dx(0)*dx(1)))
                  enddo
               enddo
            enddo
         endif
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     5-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_5_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      double precision lagrangian_bspline_5_delta
c
c     Input.
c
      integer depth
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1
      integer nindices

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic,ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),w(0:2-1,0:4)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a 5-point B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,2-1
c
c     Determine the Cartesian cell in which X(s) is located.
c
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d)
     &           +(dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c
            ic_lower(d) = ic_center(d)-2
            ic_upper(d) = ic_center(d)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic-ilower(d))+0.5d0)*dx(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_bspline_5_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Interpolate u onto V.
c
            if (ic_upper(1) - ic_lower(1) == (5 - 1) .and.
     &       ic_upper(0) - ic_lower(0) == (5 - 1)) then
            do d = 0,depth-1
               V(d, s) = 0.d0
               do ic1 = ic_lower(1),(ic_lower(1) + 5 - 1)
                  do ic0 = ic_lower(0),(ic_lower(0) + 5 - 1)
                     V(d,s) = V(d,s)
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *u(ic0,ic1,d)
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               V(d, s) = 0.d0
               do ic1 = ic_lower(1),ic_upper(1)
                  do ic0 = ic_lower(0),ic_upper(0)
                     V(d,s) = V(d,s)
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *u(ic0,ic1,d)
                  enddo
               enddo
            enddo
         endif
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the 5-point
c     B-spline delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_5_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      double precision lagrangian_bspline_5_delta
c
c     Input.
c
      integer depth
      integer nindices
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic,ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),w(0:2-1,0:4)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a 5-point B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d = 0,2-1
c
c     Determine the Cartesian cell in which X(s) is located.
c
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d)
     &           +(dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c
            ic_lower(d) = ic_center(d)-2
            ic_upper(d) = ic_center(d)+2
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic-ilower(d))+0.5d0)*dx(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_bspline_5_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Spread V onto u.
c
            if (ic_upper(1) - ic_lower(1) == (5 - 1) .and.
     &       ic_upper(0) - ic_lower(0) == (5 - 1)) then
            do d = 0,depth-1
               do ic1 = ic_lower(1),(ic_lower(1) + 5 - 1)
                  do ic0 = ic_lower(0),(ic_lower(0) + 5 - 1)
                     u(ic0,ic1,d) = u(ic0,ic1,d) + (
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *V(d,s)/(dx(0)*dx(1)))
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               do ic1 = ic_lower(1),ic_upper(1)
                  do ic0 = ic_lower(0),ic_upper(0)
                     u(ic0,ic1,d) = u(ic0,ic1,d) + (
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *V(d,s)/(dx(0)*dx(1)))
                  enddo
               enddo
            enddo
         endif
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     6-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_6_interp2d(
     &     dx,x_lower,x_upper,depth,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      double precision lagrangian_bspline_6_delta
c
c     Input.
c
      integer depth
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1
      integer nindices

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic,ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),w(0:2-1,0:5)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a 6-point B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,2-1
c
c     Determine the Cartesian cell in which X(s) is located.
c
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d)
     &           +(dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-3
               ic_upper(d) = ic_center(d)+2
            else
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+3
            endif
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic-ilower(d))+0.5d0)*dx(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_bspline_6_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Interpolate u onto V.
c
            if (ic_upper(1) - ic_lower(1) == (6 - 1) .and.
     &       ic_upper(0) - ic_lower(0) == (6 - 1)) then
            do d = 0,depth-1
               V(d, s) = 0.d0
               do ic1 = ic_lower(1),(ic_lower(1) + 6 - 1)
                  do ic0 = ic_lower(0),(ic_lower(0) + 6 - 1)
                     V(d,s) = V(d,s)
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *u(ic0,ic1,d)
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               V(d, s) = 0.d0
               do ic1 = ic_lower(1),ic_upper(1)
                  do ic0 = ic_lower(0),ic_upper(0)
                     V(d,s) = V(d,s)
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *u(ic0,ic1,d)
                  enddo
               enddo
            enddo
         endif
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the 6-point
c     B-spline delta function using standard (double) precision
c     accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_bspline_6_spread2d(
     &     dx,x_lower,x_upper,depth,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      double precision lagrangian_bspline_6_delta
c
c     Input.
c
      integer depth
      integer nindices
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic,ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),w(0:2-1,0:5)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a 6-point B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,2-1
c
c     Determine the Cartesian cell in which X(s) is located.
c
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d)
     &           +(dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-3
               ic_upper(d) = ic_center(d)+2
            else
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+3
            endif
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic-ilower(d))+0.5d0)*dx(d)
               w(d,ic-ic_lower(d)) =
     &              lagrangian_bspline_6_delta(
     &              (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
            enddo
         enddo
c
c     Spread V onto u.
c
            if (ic_upper(1) - ic_lower(1) == (6 - 1) .and.
     &       ic_upper(0) - ic_lower(0) == (6 - 1)) then
            do d = 0,depth-1
               do ic1 = ic_lower(1),(ic_lower(1) + 6 - 1)
                  do ic0 = ic_lower(0),(ic_lower(0) + 6 - 1)
                     u(ic0,ic1,d) = u(ic0,ic1,d) + (
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *V(d,s)/(dx(0)*dx(1)))
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               do ic1 = ic_lower(1),ic_upper(1)
                  do ic0 = ic_lower(0),ic_upper(0)
                     u(ic0,ic1,d) = u(ic0,ic1,d) + (
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *V(d,s)/(dx(0)*dx(1)))
                  enddo
               enddo
            enddo
         endif
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     composite 3-point/2-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_32_interp2d(
     &     dx,x_lower,x_upper,depth,axis,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      double precision lagrangian_piecewise_linear_delta
      double precision lagrangian_bspline_3_delta
c
c     Input.
c
      integer depth,axis
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1
      integer nindices

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic,ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),w(0:2-1,0:2)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,2-1
c
c     Determine the Cartesian cell in which X(s) is located.
c
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d)
     &           +(dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic-ilower(d))+0.5d0)*dx(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_3_delta(
     &                 (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_piecewise_linear_delta(
     &                 (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
               endif
            enddo
         enddo
c
c     Interpolate u onto V.
c
            if (ic_upper(1) - ic_lower(1) == (3 - 1) .and.
     &       ic_upper(0) - ic_lower(0) == (3 - 1)) then
            do d = 0,depth-1
               V(d, s) = 0.d0
               do ic1 = ic_lower(1),(ic_lower(1) + 3 - 1)
                  do ic0 = ic_lower(0),(ic_lower(0) + 3 - 1)
                     V(d,s) = V(d,s)
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *u(ic0,ic1,d)
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               V(d, s) = 0.d0
               do ic1 = ic_lower(1),ic_upper(1)
                  do ic0 = ic_lower(0),ic_upper(0)
                     V(d,s) = V(d,s)
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *u(ic0,ic1,d)
                  enddo
               enddo
            enddo
         endif
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the composite
c     3-point/2-point B-spline delta function using standard (double)
c     precision accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_32_spread2d(
     &     dx,x_lower,x_upper,depth,axis,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      double precision lagrangian_piecewise_linear_delta
      double precision lagrangian_bspline_3_delta
c
c     Input.
c
      integer depth,axis
      integer nindices
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic,ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),w(0:2-1,0:2)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,2-1
c
c     Determine the Cartesian cell in which X(s) is located.
c
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d)
     &           +(dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c
            ic_lower(d) = ic_center(d)-1
            ic_upper(d) = ic_center(d)+1
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic-ilower(d))+0.5d0)*dx(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_3_delta(
     &                 (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_piecewise_linear_delta(
     &                 (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
               endif
            enddo
         enddo
c
c     Spread V onto u.
c
            if (ic_upper(1) - ic_lower(1) == (3 - 1) .and.
     &       ic_upper(0) - ic_lower(0) == (3 - 1)) then
            do d = 0,depth-1
               do ic1 = ic_lower(1),(ic_lower(1) + 3 - 1)
                  do ic0 = ic_lower(0),(ic_lower(0) + 3 - 1)
                     u(ic0,ic1,d) = u(ic0,ic1,d) + (
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *V(d,s)/(dx(0)*dx(1)))
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               do ic1 = ic_lower(1),ic_upper(1)
                  do ic0 = ic_lower(0),ic_upper(0)
                     u(ic0,ic1,d) = u(ic0,ic1,d) + (
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *V(d,s)/(dx(0)*dx(1)))
                  enddo
               enddo
            enddo
         endif
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Interpolate u onto V at the positions specified by X using the
c     composite 4-point/3-point B-spline delta function.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_43_interp2d(
     &     dx,x_lower,x_upper,depth,axis,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u,
     &     indices,Xshift,nindices,
     &     X,V)
c
      implicit none
c
c     Functions.
c
      double precision lagrangian_bspline_3_delta
      double precision lagrangian_bspline_4_delta
c
c     Input.
c
      integer depth,axis
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1
      integer nindices

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic,ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),w(0:2-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to interpolate u onto V.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,2-1
c
c     Determine the Cartesian cell in which X(s) is located.
c
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d)
     &           +(dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c
c     Determine the interpolation stencil corresponding to the position
c     of X(s) within the cell.
c
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+1
            else
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)+2
            endif
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the interpolation weights.
c
            do ic = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic-ilower(d))+0.5d0)*dx(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_4_delta(
     &                 (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_3_delta(
     &                 (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
               endif
            enddo
         enddo
c
c     Interpolate u onto V.
c
            if (ic_upper(1) - ic_lower(1) == (4 - 1) .and.
     &       ic_upper(0) - ic_lower(0) == (4 - 1)) then
            do d = 0,depth-1
               V(d, s) = 0.d0
               do ic1 = ic_lower(1),(ic_lower(1) + 4 - 1)
                  do ic0 = ic_lower(0),(ic_lower(0) + 4 - 1)
                     V(d,s) = V(d,s)
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *u(ic0,ic1,d)
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               V(d, s) = 0.d0
               do ic1 = ic_lower(1),ic_upper(1)
                  do ic0 = ic_lower(0),ic_upper(0)
                     V(d,s) = V(d,s)
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *u(ic0,ic1,d)
                  enddo
               enddo
            enddo
         endif
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Spread V onto u at the positions specified by X using the composite
c     4-point/3-point B-spline delta function using standard (double)
c     precision accumulation on the Cartesian grid.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lagrangian_composite_bspline_43_spread2d(
     &     dx,x_lower,x_upper,depth,axis,
     &     indices,Xshift,nindices,
     &     X,V,
     &     ilower0,iupper0,ilower1,iupper1,
     &     nugc0,nugc1,
     &     u)
c
      implicit none
c
c     Functions.
c
      double precision lagrangian_bspline_3_delta
      double precision lagrangian_bspline_4_delta
c
c     Input.
c
      integer depth,axis
      integer nindices
      integer ilower0,iupper0,ilower1,iupper1
      integer nugc0,nugc1

      integer indices(0:nindices-1)

      double precision Xshift(0:2-1,0:nindices-1)

      double precision dx(0:2-1),x_lower(0:2-1),x_upper(0:2-1)
      double precision u(ilower0-nugc0:iupper0+nugc0,
     &          ilower1-nugc1:iupper1+nugc1,0:depth-1)
      double precision X(0:2-1,0:*)
c
c     Input/Output.
c
      double precision V(0:depth-1,0:*)
c
c     Local variables.
c
      integer ilower(0:2-1),iupper(0:2-1)
      integer ic,ic0,ic1
      integer ic_center(0:2-1),ic_lower(0:2-1),ic_upper(0:2-1)
      integer d,l,s,nugc(0:2-1)

      double precision X_cell(0:2-1),w(0:2-1,0:3)
c
c     Prevent compiler warning about unused variables.
c
      x_upper(0) = x_upper(0)
c
c     Setup convenience arrays.
c
      ilower(0) = ilower0
      ilower(1) = ilower1

      iupper(0) = iupper0
      iupper(1) = iupper1

      nugc(0) = nugc0
      nugc(1) = nugc1
c
c     Use a composite B-spline function to spread V onto u.
c
      do l = 0,nindices-1
         s = indices(l)
         do d=0,2-1
c
c     Determine the Cartesian cell in which X(s) is located.
c
            ic_center(d) =
     &           floor((X(d,s)+Xshift(d,l)-x_lower(d))/dx(d))
     &           + ilower(d)
            X_cell(d) = x_lower(d)
     &           +(dble(ic_center(d)-ilower(d))+0.5d0)*dx(d)
c
c     Determine the spreading stencil corresponding to the position of
c     X(s) within the cell.
c
            if ( X(d,s).lt.X_cell(d) ) then
               ic_lower(d) = ic_center(d)-2
               ic_upper(d) = ic_center(d)+1
            else
               ic_lower(d) = ic_center(d)-1
               ic_upper(d) = ic_center(d)+2
            endif
            ic_lower(d) = max(ic_lower(d),ilower(d)-nugc(d))
            ic_upper(d) = min(ic_upper(d),iupper(d)+nugc(d))
c
c     Compute the spreading weights.
c
            do ic = ic_lower(d),ic_upper(d)
               X_cell(d) = x_lower(d)+(dble(ic-ilower(d))+0.5d0)*dx(d)
               if (d .eq. axis) then
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_4_delta(
     &                 (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
               else
                  w(d,ic-ic_lower(d)) =
     &                 lagrangian_bspline_3_delta(
     &                 (X(d,s)+Xshift(d,l)-X_cell(d))/dx(d))
               endif
            enddo
         enddo
c
c     Spread V onto u.
c
            if (ic_upper(1) - ic_lower(1) == (4 - 1) .and.
     &       ic_upper(0) - ic_lower(0) == (4 - 1)) then
            do d = 0,depth-1
               do ic1 = ic_lower(1),(ic_lower(1) + 4 - 1)
                  do ic0 = ic_lower(0),(ic_lower(0) + 4 - 1)
                     u(ic0,ic1,d) = u(ic0,ic1,d) + (
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *V(d,s)/(dx(0)*dx(1)))
                  enddo
               enddo
            enddo
         else
            do d = 0,depth-1
               do ic1 = ic_lower(1),ic_upper(1)
                  do ic0 = ic_lower(0),ic_upper(0)
                     u(ic0,ic1,d) = u(ic0,ic1,d) + (
     &                    +w(0,ic0-ic_lower(0))
     &                    *w(1,ic1-ic_lower(1))
     &                    *V(d,s)/(dx(0)*dx(1)))
                  enddo
               enddo
            enddo
         endif
c
c     End loop over points.
c
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
