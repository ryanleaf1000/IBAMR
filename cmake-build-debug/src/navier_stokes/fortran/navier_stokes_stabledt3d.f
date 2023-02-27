c ---------------------------------------------------------------------
c
c Copyright (c) 2006 - 2019 by the IBAMR developers
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
c     Determine the timestep size corresponding to a unit CFL number.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_cc_stabledt3d(
     &     dx,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     ngc0,ngc1,ngc2,
     &     U0,U1,U2,
     &     stabdt)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      integer ngc0,ngc1,ngc2

      double precision dx(0:3-1)

      double precision U0(ifirst0-ngc0:ilast0+ngc0,
     &          ifirst1-ngc1:ilast1+ngc1,
     &          ifirst2-ngc2:ilast2+ngc2)
      double precision U1(ifirst0-ngc0:ilast0+ngc0,
     &          ifirst1-ngc1:ilast1+ngc1,
     &          ifirst2-ngc2:ilast2+ngc2)
      double precision U2(ifirst0-ngc0:ilast0+ngc0,
     &          ifirst1-ngc1:ilast1+ngc1,
     &          ifirst2-ngc2:ilast2+ngc2)
c
c     Input/Output.
c
      double precision stabdt
c
c     Local variables.
c
      integer i0,i1,i2,d
      double precision maxspeed(0:3-1)
c
c     Determine the unit CFL number on the patch.
c
      do d = 0,3-1
         maxspeed(d) = 1.0d-12  ! avoid division by zero
      enddo

      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0
               maxspeed(0) = dmax1(maxspeed(0), dabs(U0(i0,i1,i2)))
            enddo
         enddo
      enddo

      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0
               maxspeed(1) = dmax1(maxspeed(1), dabs(U1(i0,i1,i2)))
            enddo
         enddo
      enddo

      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0
               maxspeed(2) = dmax1(maxspeed(2), dabs(U2(i0,i1,i2)))
            enddo
         enddo
      enddo

      stabdt = dmin1((dx(0)/maxspeed(0)),(dx(1)/maxspeed(1)),
     &               (dx(2)/maxspeed(2)))
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Determine the timestep size corresponding to a unit CFL number.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_fc_stabledt3d(
     &     dx,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     ngc0,ngc1,ngc2,
     &     U0,U1,U2,
     &     stabdt)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      integer ngc0,ngc1,ngc2

      double precision dx(0:3-1)

      double precision U0(ifirst0-ngc0:ilast0+1+ngc0,
     &          ifirst1-ngc1:ilast1+ngc1,
     &          ifirst2-ngc2:ilast2+ngc2)
      double precision U1(ifirst1-ngc1:ilast1+1+ngc1,
     &          ifirst2-ngc2:ilast2+ngc2,
     &          ifirst0-ngc0:ilast0+ngc0)
      double precision U2(ifirst2-ngc2:ilast2+1+ngc2,
     &          ifirst0-ngc0:ilast0+ngc0,
     &          ifirst1-ngc1:ilast1+ngc1)
c
c     Input/Output.
c
      double precision stabdt
c
c     Local variables.
c
      integer i0,i1,i2,d
      double precision maxspeed(0:3-1)
c
c     Determine the unit CFL number on the patch.
c
      do d = 0,3-1
         maxspeed(d) = 1.0d-12  ! avoid division by zero
      enddo

      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0+1
               maxspeed(0) = dmax1(maxspeed(0), dabs(U0(i0,i1,i2)))
            enddo
         enddo
      enddo

      do i0 = ifirst0,ilast0
         do i2 = ifirst2,ilast2
            do i1 = ifirst1,ilast1+1
               maxspeed(1) = dmax1(maxspeed(1), dabs(U1(i1,i2,i0)))
            enddo
         enddo
      enddo

      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
            do i2 = ifirst2,ilast2+1
               maxspeed(2) = dmax1(maxspeed(2), dabs(U2(i2,i0,i1)))
            enddo
         enddo
      enddo

      stabdt = dmin1((dx(0)/maxspeed(0)),(dx(1)/maxspeed(1)),
     &               (dx(2)/maxspeed(2)))
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Determine the timestep size corresponding to a unit CFL number.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_sc_stabledt3d(
     &     dx,
     &     ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &     ngc0,ngc1,ngc2,
     &     U0,U1,U2,
     &     stabdt)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2

      integer ngc0,ngc1,ngc2

      double precision dx(0:3-1)

      double precision U0(ifirst0-ngc0:ilast0+1+ngc0,
     &          ifirst1-ngc1:ilast1+ngc1,
     &          ifirst2-ngc2:ilast2+ngc2)
      double precision U1(ifirst0-ngc0:ilast0+ngc0,
     &          ifirst1-ngc1:ilast1+1+ngc1,
     &          ifirst2-ngc2:ilast2+ngc2)
      double precision U2(ifirst0-ngc0:ilast0+ngc0,
     &          ifirst1-ngc1:ilast1+ngc1,
     &          ifirst2-ngc2:ilast2+1+ngc2)
c
c     Input/Output.
c
      double precision stabdt
c
c     Local variables.
c
      integer i0,i1,i2,d
      double precision maxspeed(0:3-1)
c
c     Determine the unit CFL number on the patch.
c
      do d = 0,3-1
         maxspeed(d) = 1.0d-12  ! avoid division by zero
      enddo

      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0+1
               maxspeed(0) = dmax1(maxspeed(0), dabs(U0(i0,i1,i2)))
            enddo
         enddo
      enddo

      do i2 = ifirst2,ilast2
         do i1 = ifirst1,ilast1+1
            do i0 = ifirst0,ilast0
               maxspeed(1) = dmax1(maxspeed(1), dabs(U1(i0,i1,i2)))
            enddo
         enddo
      enddo

      do i2 = ifirst2,ilast2+1
         do i1 = ifirst1,ilast1
            do i0 = ifirst0,ilast0
               maxspeed(2) = dmax1(maxspeed(2), dabs(U2(i0,i1,i2)))
            enddo
         enddo
      enddo

      stabdt = dmin1((dx(0)/maxspeed(0)),(dx(1)/maxspeed(1)),
     &               (dx(2)/maxspeed(2)))
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
