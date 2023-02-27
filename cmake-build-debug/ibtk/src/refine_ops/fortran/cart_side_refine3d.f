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
c  File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/fortran/pdat_m4arrdim3d.i $
c  Package:     SAMRAI patchdata
c  Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
c  Revision:    $LastChangedRevision: 1917 $
c  Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
c  Description: m4 include file for dimensioning 3d arrays in FORTRAN routines.
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
c     The minmod function of three arguments.
c
c     The minmod function is a function of two or more arguments that
c     takes the value of the argument with the smallest modulus if all
c     arguments have the same sign.  Otherwise it takes the value zero.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      double precision function minmod3(a,b,c)
c
      implicit none
c
c     Input.
c
      double precision a,b,c
c
c     minmod(a,b,c)
c
      if     ( (a.ge.0.0d0).and.(b.ge.0.0d0).and.(c.ge.0.0d0) ) then
         minmod3 = dmin1(a,b,c)
      elseif ( (a.le.0.0d0).and.(b.le.0.0d0).and.(c.le.0.0d0) ) then
         minmod3 = dmax1(a,b,c)
      else
         minmod3 = 0.d0
      endif
c
      return
      end


c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform side-centered refine operation based on RT0 interpolation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cart_side_rt0_refine3d(
     &     u0_f,u1_f,u2_f,u_f_gcw,
     &     flower0,fupper0,
     &     flower1,fupper1,
     &     flower2,fupper2,
     &     u0_c,u1_c,u2_c,u_c_gcw,
     &     clower0,cupper0,
     &     clower1,cupper1,
     &     clower2,cupper2,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     ratio)
c
      implicit none
c
c     Input.
c
      integer u_f_gcw
      integer flower0,fupper0
      integer flower1,fupper1
      integer flower2,fupper2
      integer u_c_gcw
      integer clower0,cupper0
      integer clower1,cupper1
      integer clower2,cupper2
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer ilower2,iupper2
      integer ratio(0:3-1)

      double precision u0_c(clower0-u_c_gcw:cupper0+1+u_c_gcw,
     &          clower1-u_c_gcw:cupper1+u_c_gcw,
     &          clower2-u_c_gcw:cupper2+u_c_gcw)
      double precision u1_c(clower0-u_c_gcw:cupper0+u_c_gcw,
     &          clower1-u_c_gcw:cupper1+1+u_c_gcw,
     &          clower2-u_c_gcw:cupper2+u_c_gcw)
      double precision u2_c(clower0-u_c_gcw:cupper0+u_c_gcw,
     &          clower1-u_c_gcw:cupper1+u_c_gcw,
     &          clower2-u_c_gcw:cupper2+1+u_c_gcw)
c
c     Output.
c
      double precision u0_f(flower0-u_f_gcw:fupper0+1+u_f_gcw,
     &          flower1-u_f_gcw:fupper1+u_f_gcw,
     &          flower2-u_f_gcw:fupper2+u_f_gcw)
      double precision u1_f(flower0-u_f_gcw:fupper0+u_f_gcw,
     &          flower1-u_f_gcw:fupper1+1+u_f_gcw,
     &          flower2-u_f_gcw:fupper2+u_f_gcw)
      double precision u2_f(flower0-u_f_gcw:fupper0+u_f_gcw,
     &          flower1-u_f_gcw:fupper1+u_f_gcw,
     &          flower2-u_f_gcw:fupper2+1+u_f_gcw)
c
c     Local variables.
c
      integer i0,i1,i2
      integer i_c0,i_c1,i_c2
      integer i_f0,i_f1,i_f2
      double precision w0,w1
c
c     Refine data.
c
      do i2=ilower2,iupper2
         if (i2.lt.0) then
            i_c2=(i2+1)/ratio(2)-1
         else
            i_c2=i2/ratio(2)
         endif
         i_f2=i_c2*ratio(2)

         do i1=ilower1,iupper1
            if (i1.lt.0) then
            i_c1=(i1+1)/ratio(1)-1
         else
            i_c1=i1/ratio(1)
         endif
         i_f1=i_c1*ratio(1)

            do i0=ilower0,iupper0+1
               if (i0.lt.0) then
            i_c0=(i0+1)/ratio(0)-1
         else
            i_c0=i0/ratio(0)
         endif
         i_f0=i_c0*ratio(0)

               if ( i0 .eq. i_f0 ) then
                  u0_f(i0,i1,i2) = u0_c(i_c0,i_c1,i_c2)
               else
                  w1 = dble(i0-i_f0)/dble(ratio(0))
                  w0 = 1.d0-w1
                  u0_f(i0,i1,i2) =
     &                 w0*u0_c(i_c0  ,i_c1,i_c2) +
     &                 w1*u0_c(i_c0+1,i_c1,i_c2)
               endif
            enddo
         enddo
      enddo

      do i2=ilower2,iupper2
         if (i2.lt.0) then
            i_c2=(i2+1)/ratio(2)-1
         else
            i_c2=i2/ratio(2)
         endif
         i_f2=i_c2*ratio(2)

         do i1=ilower1,iupper1+1
            if (i1.lt.0) then
            i_c1=(i1+1)/ratio(1)-1
         else
            i_c1=i1/ratio(1)
         endif
         i_f1=i_c1*ratio(1)

            do i0=ilower0,iupper0
               if (i0.lt.0) then
            i_c0=(i0+1)/ratio(0)-1
         else
            i_c0=i0/ratio(0)
         endif
         i_f0=i_c0*ratio(0)

               if ( i1 .eq. i_f1 ) then
                  u1_f(i0,i1,i2) = u1_c(i_c0,i_c1,i_c2)
               else
                  w1 = dble(i1-i_f1)/dble(ratio(1))
                  w0 = 1.d0-w1
                  u1_f(i0,i1,i2) =
     &                 w0*u1_c(i_c0,i_c1  ,i_c2) +
     &                 w1*u1_c(i_c0,i_c1+1,i_c2)
               endif
            enddo
         enddo
      enddo

      do i2=ilower2,iupper2+1
         if (i2.lt.0) then
            i_c2=(i2+1)/ratio(2)-1
         else
            i_c2=i2/ratio(2)
         endif
         i_f2=i_c2*ratio(2)

         do i1=ilower1,iupper1
            if (i1.lt.0) then
            i_c1=(i1+1)/ratio(1)-1
         else
            i_c1=i1/ratio(1)
         endif
         i_f1=i_c1*ratio(1)

            do i0=ilower0,iupper0
               if (i0.lt.0) then
            i_c0=(i0+1)/ratio(0)-1
         else
            i_c0=i0/ratio(0)
         endif
         i_f0=i_c0*ratio(0)

               if ( i2 .eq. i_f2 ) then
                  u2_f(i0,i1,i2) = u2_c(i_c0,i_c1,i_c2)
               else
                  w1 = dble(i2-i_f2)/dble(ratio(2))
                  w0 = 1.d0-w1
                  u2_f(i0,i1,i2) =
     &                 w0*u2_c(i_c0,i_c1,i_c2  ) +
     &                 w1*u2_c(i_c0,i_c1,i_c2+1)
               endif
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform specialized refine operation that employs linear
c     interpolation in the normal direction and constant interpolation
c     in the tangential direction.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cart_side_specialized_linear_refine3d(
     &     u0_f,u1_f,u2_f,u_fine_gcw,
     &     flower0,fupper0,
     &     flower1,fupper1,
     &     flower2,fupper2,
     &     u0_c,u1_c,u2_c,u_coarse_gcw,
     &     clower0,cupper0,
     &     clower1,cupper1,
     &     clower2,cupper2,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     ratio)
c
      implicit none
c
c     External functions.
c
      double precision minmod3
c
c     Input.
c
      integer u_fine_gcw
      integer flower0,fupper0
      integer flower1,fupper1
      integer flower2,fupper2
      integer u_coarse_gcw
      integer clower0,cupper0
      integer clower1,cupper1
      integer clower2,cupper2
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer ilower2,iupper2
      integer ratio(0:3-1)

      double precision u0_c(clower0-u_coarse_gcw:cupper0+1+u_coarse_gcw,
     &          clower1-u_coarse_gcw:cupper1+u_coarse_gcw,
     &          clower2-u_coarse_gcw:cupper2+u_coarse_gcw)
      double precision u1_c(clower0-u_coarse_gcw:cupper0+u_coarse_gcw,
     &          clower1-u_coarse_gcw:cupper1+1+u_coarse_gcw,
     &          clower2-u_coarse_gcw:cupper2+u_coarse_gcw)
      double precision u2_c(clower0-u_coarse_gcw:cupper0+u_coarse_gcw,
     &          clower1-u_coarse_gcw:cupper1+u_coarse_gcw,
     &          clower2-u_coarse_gcw:cupper2+1+u_coarse_gcw)
c
c     Input/Output.
c
      double precision u0_f(flower0-u_fine_gcw:fupper0+1+u_fine_gcw,
     &          flower1-u_fine_gcw:fupper1+u_fine_gcw,
     &          flower2-u_fine_gcw:fupper2+u_fine_gcw)
      double precision u1_f(flower0-u_fine_gcw:fupper0+u_fine_gcw,
     &          flower1-u_fine_gcw:fupper1+1+u_fine_gcw,
     &          flower2-u_fine_gcw:fupper2+u_fine_gcw)
      double precision u2_f(flower0-u_fine_gcw:fupper0+u_fine_gcw,
     &          flower1-u_fine_gcw:fupper1+u_fine_gcw,
     &          flower2-u_fine_gcw:fupper2+1+u_fine_gcw)
c
c     Local variables.
c
      integer i0,i1,i2
      integer i_c0,i_c1,i_c2
      integer i_f0,i_f1,i_f2
      double precision w0,w1,w2
      double precision dx0,dx1,dx2
      double precision du0_dx1,du0_dx2,du1_dx0,du1_dx2,du2_dx0,du2_dx1
c
c     Refine data.
c
      dx0 = dble(ratio(0))      ! effective grid spacings on the coarse level
      dx1 = dble(ratio(1))
      dx2 = dble(ratio(2))

      do i2=ilower2,iupper2
         if (i2.lt.0) then
            i_c2=(i2+1)/ratio(2)-1
         else
            i_c2=i2/ratio(2)
         endif
         i_f2=i_c2*ratio(2)

         do i1=ilower1,iupper1
            if (i1.lt.0) then
            i_c1=(i1+1)/ratio(1)-1
         else
            i_c1=i1/ratio(1)
         endif
         i_f1=i_c1*ratio(1)

            do i0=ilower0,iupper0+1
               if (i0.lt.0) then
            i_c0=(i0+1)/ratio(0)-1
         else
            i_c0=i0/ratio(0)
         endif
         i_f0=i_c0*ratio(0)


               w0 = 1.d0-dble(i0-i_f0)/dx0
               w1 = dble(i1-i_f1)+0.5d0-0.5d0*dx1
               w2 = dble(i2-i_f2)+0.5d0-0.5d0*dx2

               du0_dx1 = minmod3(
     & 0.5d0*(u0_c(i_c0,i_c1+1,i_c2)-u0_c(i_c0,i_c1-1,i_c2)),
     & 2.d0*(u0_c(i_c0,i_c1,i_c2)-u0_c(i_c0,i_c1-1,i_c2)),
     & 2.d0*(u0_c(i_c0,i_c1+1,i_c2)-u0_c(i_c0,i_c1,i_c2)))/dx1
               du0_dx2 = minmod3(
     & 0.5d0*(u0_c(i_c0,i_c1,i_c2+1)-u0_c(i_c0,i_c1,i_c2-1)),
     & 2.d0*(u0_c(i_c0,i_c1,i_c2)-u0_c(i_c0,i_c1,i_c2-1)),
     & 2.d0*(u0_c(i_c0,i_c1,i_c2+1)-u0_c(i_c0,i_c1,i_c2)))/dx2
               u0_f(i0,i1,i2) =
     &              w0*(u0_c(i_c0,i_c1,i_c2)+w1*du0_dx1+w2*du0_dx2)

               w0 = 1.d0-w0

               du0_dx1 = minmod3(
     & 0.5d0*(u0_c(i_c0+1,i_c1+1,i_c2)-u0_c(i_c0+1,i_c1-1,i_c2)),
     & 2.d0*(u0_c(i_c0+1,i_c1,i_c2)-u0_c(i_c0+1,i_c1-1,i_c2)),
     & 2.d0*(u0_c(i_c0+1,i_c1+1,i_c2)-u0_c(i_c0+1,i_c1,i_c2)))/dx1
               du0_dx2 = minmod3(
     & 0.5d0*(u0_c(i_c0+1,i_c1,i_c2+1)-u0_c(i_c0+1,i_c1,i_c2-1)),
     & 2.d0*(u0_c(i_c0+1,i_c1,i_c2)-u0_c(i_c0+1,i_c1,i_c2-1)),
     & 2.d0*(u0_c(i_c0+1,i_c1,i_c2+1)-u0_c(i_c0+1,i_c1,i_c2)))/dx2
               u0_f(i0,i1,i2) = u0_f(i0,i1,i2)+
     &              w0*(u0_c(i_c0+1,i_c1,i_c2)+w1*du0_dx1+w2*du0_dx2)
            enddo
         enddo
      enddo

      do i2=ilower2,iupper2
         if (i2.lt.0) then
            i_c2=(i2+1)/ratio(2)-1
         else
            i_c2=i2/ratio(2)
         endif
         i_f2=i_c2*ratio(2)

         do i1=ilower1,iupper1
            if (i1.lt.0) then
            i_c1=(i1+1)/ratio(1)-1
         else
            i_c1=i1/ratio(1)
         endif
         i_f1=i_c1*ratio(1)

            do i0=ilower0,iupper0+1
               if (i0.lt.0) then
            i_c0=(i0+1)/ratio(0)-1
         else
            i_c0=i0/ratio(0)
         endif
         i_f0=i_c0*ratio(0)


               i_f0 = i_c0*ratio(0)
               i_f1 = i_c1*ratio(1)
               i_f2 = i_c2*ratio(2)

               w0 = dble(i0-i_f0)+0.5d0-0.5d0*dx0
               w1 = 1.d0-dble(i1-i_f1)/dx1
               w2 = dble(i2-i_f2)+0.5d0-0.5d0*dx2

               du1_dx0 = minmod3(
     & 0.5d0*(u1_c(i_c0+1,i_c1,i_c2)-u1_c(i_c0-1,i_c1,i_c2)),
     & 2.d0*(u1_c(i_c0,i_c1,i_c2)-u1_c(i_c0-1,i_c1,i_c2)),
     & 2.d0*(u1_c(i_c0+1,i_c1,i_c2)-u1_c(i_c0,i_c1,i_c2)))/dx0
               du1_dx2 = minmod3(
     & 0.5d0*(u1_c(i_c0,i_c1,i_c2+1)-u1_c(i_c0,i_c1,i_c2-1)),
     & 2.d0*(u1_c(i_c0,i_c1,i_c2)-u1_c(i_c0,i_c1,i_c2-1)),
     & 2.d0*(u1_c(i_c0,i_c1,i_c2+1)-u1_c(i_c0,i_c1,i_c2)))/dx2
               u1_f(i0,i1,i2) =
     &              w1*(u1_c(i_c0,i_c1,i_c2)+w0*du1_dx0+w2*du1_dx2)

               w1 = 1.d0-w1

               du1_dx0 = minmod3(
     & 0.5d0*(u1_c(i_c0+1,i_c1+1,i_c2)-u1_c(i_c0-1,i_c1+1,i_c2)),
     & 2.d0*(u1_c(i_c0,i_c1+1,i_c2)-u1_c(i_c0-1,i_c1+1,i_c2)),
     & 2.d0*(u1_c(i_c0+1,i_c1+1,i_c2)-u1_c(i_c0,i_c1+1,i_c2)))/dx0
               du1_dx2 = minmod3(
     & 0.5d0*(u1_c(i_c0,i_c1+1,i_c2+1)-u1_c(i_c0,i_c1+1,i_c2-1)),
     & 2.d0*(u1_c(i_c0,i_c1+1,i_c2)-u1_c(i_c0,i_c1+1,i_c2-1)),
     & 2.d0*(u1_c(i_c0,i_c1+1,i_c2+1)-u1_c(i_c0,i_c1+1,i_c2)))/dx2
               u1_f(i0,i1,i2) = u1_f(i0,i1,i2)+
     &              w1*(u1_c(i_c0,i_c1+1,i_c2)+w0*du1_dx0+w2*du1_dx2)
            enddo
         enddo
      enddo

      do i2=ilower2,iupper2
         if (i2.lt.0) then
            i_c2=(i2+1)/ratio(2)-1
         else
            i_c2=i2/ratio(2)
         endif
         i_f2=i_c2*ratio(2)

         do i1=ilower1,iupper1
            if (i1.lt.0) then
            i_c1=(i1+1)/ratio(1)-1
         else
            i_c1=i1/ratio(1)
         endif
         i_f1=i_c1*ratio(1)

            do i0=ilower0,iupper0+1
               if (i0.lt.0) then
            i_c0=(i0+1)/ratio(0)-1
         else
            i_c0=i0/ratio(0)
         endif
         i_f0=i_c0*ratio(0)


               i_f0 = i_c0*ratio(0)
               i_f1 = i_c1*ratio(1)
               i_f2 = i_c2*ratio(2)

               w0 = dble(i0-i_f0)+0.5d0-0.5d0*dx0
               w1 = dble(i1-i_f1)+0.5d0-0.5d0*dx1
               w2 = 1.d0-dble(i2-i_f2)/dx2

               du2_dx0 = minmod3(
     & 0.5d0*(u2_c(i_c0+1,i_c1,i_c2)-u2_c(i_c0-1,i_c1,i_c2)),
     & 2.d0*(u2_c(i_c0,i_c1,i_c2)-u2_c(i_c0-1,i_c1,i_c2)),
     & 2.d0*(u2_c(i_c0+1,i_c1,i_c2)-u2_c(i_c0,i_c1,i_c2)))/dx0
               du2_dx1 = minmod3(
     & 0.5d0*(u2_c(i_c0,i_c1+1,i_c2)-u2_c(i_c0,i_c1-1,i_c2)),
     & 2.d0*(u2_c(i_c0,i_c1,i_c2)-u2_c(i_c0,i_c1-1,i_c2)),
     & 2.d0*(u2_c(i_c0,i_c1+1,i_c2)-u2_c(i_c0,i_c1,i_c2)))/dx1
               u2_f(i0,i1,i2) =
     &              w2*(u2_c(i_c0,i_c1,i_c2)+w0*du2_dx0+w1*du2_dx1)

               w2 = 1.d0-w2

               du2_dx0 = minmod3(
     & 0.5d0*(u2_c(i_c0+1,i_c1,i_c2+1)-u2_c(i_c0-1,i_c1,i_c2+1)),
     & 2.d0*(u2_c(i_c0,i_c1,i_c2+1)-u2_c(i_c0-1,i_c1,i_c2+1)),
     & 2.d0*(u2_c(i_c0+1,i_c1,i_c2+1)-u2_c(i_c0,i_c1,i_c2+1)))/dx0
               du2_dx1 = minmod3(
     & 0.5d0*(u2_c(i_c0,i_c1+1,i_c2+1)-u2_c(i_c0,i_c1-1,i_c2+1)),
     & 2.d0*(u2_c(i_c0,i_c1,i_c2+1)-u2_c(i_c0,i_c1-1,i_c2+1)),
     & 2.d0*(u2_c(i_c0,i_c1+1,i_c2+1)-u2_c(i_c0,i_c1,i_c2+1)))/dx1
               u2_f(i0,i1,i2) = u2_f(i0,i1,i2)+
     &              w2*(u2_c(i_c0,i_c1,i_c2+1)+w0*du2_dx0+w1*du2_dx1)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
