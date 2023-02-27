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
c     Update a quantity using flux differencing.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine adv_diff_consdiff2d(
     &     dx,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nfluxgc0,nfluxgc1,
     &     nqvalgc0,nqvalgc1,
     &     flux0,flux1,
     &     qval)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0,ifirst1,ilast1

      integer nfluxgc0,nfluxgc1
      integer nqvalgc0,nqvalgc1

      double precision dx(0:2-1),dt

      double precision flux0(ifirst0-nfluxgc0:ilast0+1+nfluxgc0,
     &          ifirst1-nfluxgc1:ilast1+nfluxgc1)
      double precision flux1(ifirst1-nfluxgc1:ilast1+1+nfluxgc1,
     &          ifirst0-nfluxgc0:ilast0+nfluxgc0)
c
c     Input/Output.
c
      double precision qval(ifirst0-nqvalgc0:ilast0+nqvalgc0,
     &          ifirst1-nqvalgc1:ilast1+nqvalgc1)
c
c     Local variables.
c
      integer ic0,ic1,d
      double precision dtdx(0:2-1)
c
c     Update a quantity using flux differencing.
c
      do d = 0,2-1
         dtdx(d) = dt*dx(d)
      enddo

      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            qval(ic0,ic1) =
     &           -(flux0(ic0+1,ic1)-flux0(ic0,ic1))/dtdx(0)
     &           -(flux1(ic1+1,ic0)-flux1(ic1,ic0))/dtdx(1)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Update a quantity using flux differencing but include the proper
c     source term to account for a non-discretely divergence free
c     advection velocity.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine adv_diff_consdiffwithdivsource2d(
     &     dx,dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nfluxgc0,nfluxgc1,
     &     nqfluxgc0,nqfluxgc1,
     &     nufluxgc0,nufluxgc1,
     &     nqvalgc0,nqvalgc1,
     &     flux0,flux1,
     &     qflux0,qflux1,
     &     uflux0,uflux1,
     &     qval)
c
      implicit none
      double precision fourth
      parameter (fourth=0.25d0)
c
c     Input.
c
      integer ifirst0,ilast0,ifirst1,ilast1

      integer nfluxgc0,nfluxgc1
      integer nqfluxgc0,nqfluxgc1
      integer nufluxgc0,nufluxgc1
      integer nqvalgc0,nqvalgc1

      double precision dx(0:2-1),dt

      double precision flux0(ifirst0-nfluxgc0:ilast0+1+nfluxgc0,
     &          ifirst1-nfluxgc1:ilast1+nfluxgc1)
      double precision flux1(ifirst1-nfluxgc1:ilast1+1+nfluxgc1,
     &          ifirst0-nfluxgc0:ilast0+nfluxgc0)

      double precision qflux0(ifirst0-nqfluxgc0:ilast0+1+nqfluxgc0,
     &          ifirst1-nqfluxgc1:ilast1+nqfluxgc1)
      double precision qflux1(ifirst1-nqfluxgc1:ilast1+1+nqfluxgc1,
     &          ifirst0-nqfluxgc0:ilast0+nqfluxgc0)

      double precision uflux0(ifirst0-nufluxgc0:ilast0+1+nufluxgc0,
     &          ifirst1-nufluxgc1:ilast1+nufluxgc1)
      double precision uflux1(ifirst1-nufluxgc1:ilast1+1+nufluxgc1,
     &          ifirst0-nufluxgc0:ilast0+nufluxgc0)
c
c     Input/Output.
c
      double precision qval(ifirst0-nqvalgc0:ilast0+nqvalgc0,
     &          ifirst1-nqvalgc1:ilast1+nqvalgc1)
c
c     Local variables.
c
      integer ic0,ic1,d
      double precision dtdx(0:2-1),divsource
c
c     Update a quantity using flux differencing.
c
      do d = 0,2-1
         dtdx(d) = dt*dx(d)
      enddo

      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            divsource = (fourth/(dt**2.d0))*
     &           ( qflux0(ic0+1,ic1) + qflux0(ic0,ic1)
     &           + qflux1(ic1+1,ic0) + qflux1(ic1,ic0) )*
     &           ( (uflux0(ic0+1,ic1)-uflux0(ic0,ic1))/dx(0)
     &           + (uflux1(ic1+1,ic0)-uflux1(ic1,ic0))/dx(1) )

            qval(ic0,ic1) = divsource
     &           -(flux0(ic0+1,ic1)-flux0(ic0,ic1))/dtdx(0)
     &           -(flux1(ic1+1,ic0)-flux1(ic1,ic0))/dtdx(1)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
