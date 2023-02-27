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
c     Compute the source term F = -U max(q,0) that must be added to the
c     momentum equation to account for momentum loss due to internal
c     fluid sources.
c
c     NOTE: This is the source term that corresponds to the advective
c     (i.e., nonconservative) form of the equation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_adv_source2d(
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nQgc0,nQgc1,
     &     nFgc0,nFgc1,
     &     u0,u1,
     &     Q,
     &     F)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0,ifirst1,ilast1

      integer nugc0,nugc1
      integer nQgc0,nQgc1
      integer nFgc0,nFgc1

      double precision u0(ifirst0-nugc0:ilast0+1+nugc0,
     &          ifirst1-nugc1:ilast1+nugc1)
      double precision u1(ifirst1-nugc1:ilast1+1+nugc1,
     &          ifirst0-nugc0:ilast0+nugc0)

      double precision Q(ifirst0-nQgc0:ilast0+nQgc0,
     &          ifirst1-nQgc1:ilast1+nQgc1)
c
c     Input/Output.
c
      double precision F(ifirst0-nFgc0:ilast0+nFgc0,
     &          ifirst1-nFgc1:ilast1+nFgc1,0:2-1)
c
c     Local variables.
c
      integer ic0,ic1
c
c     Compute the source term F = -U max(Q,0).
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            F(ic0,ic1,0) = -0.5d0*(u0(ic0+1,ic1)+u0(ic0,ic1))*
     &           dmax1(Q(ic0,ic1),0.d0)
         enddo
      enddo

      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            F(ic0,ic1,1) = -0.5d0*(u1(ic1+1,ic0)+u1(ic1,ic0))*
     &           dmax1(Q(ic0,ic1),0.d0)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the source term F = +U min(q,0) that must be added to the
c     momentum equation to account for momentum loss due to internal
c     fluid sinks.
c
c     NOTE: This is the source term that corresponds to the conservative
c     form of the equation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_cons_source2d(
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nQgc0,nQgc1,
     &     nFgc0,nFgc1,
     &     u0,u1,
     &     Q,
     &     F)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0,ifirst1,ilast1

      integer nugc0,nugc1
      integer nQgc0,nQgc1
      integer nFgc0,nFgc1

      double precision u0(ifirst0-nugc0:ilast0+1+nugc0,
     &          ifirst1-nugc1:ilast1+nugc1)
      double precision u1(ifirst1-nugc1:ilast1+1+nugc1,
     &          ifirst0-nugc0:ilast0+nugc0)

      double precision Q(ifirst0-nQgc0:ilast0+nQgc0,
     &          ifirst1-nQgc1:ilast1+nQgc1)
c
c     Input/Output.
c
      double precision F(ifirst0-nFgc0:ilast0+nFgc0,
     &          ifirst1-nFgc1:ilast1+nFgc1,0:2-1)
c
c     Local variables.
c
      integer ic0,ic1
c
c     Compute the source term F = +U min(Q,0).
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            F(ic0,ic1,0) = +0.5d0*(u0(ic0+1,ic1)+u0(ic0,ic1))*
     &           dmin1(Q(ic0,ic1),0.d0)
         enddo
      enddo

      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            F(ic0,ic1,1) = +0.5d0*(u1(ic1+1,ic0)+u1(ic1,ic0))*
     &           dmin1(Q(ic0,ic1),0.d0)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the source term F = -0.5 U abs(q) that must be added to
c     the momentum equation to account for momentum loss due to internal
c     fluid sinks.
c
c     NOTE: This is the source term that corresponds to the
c     skew-symmetric form of the equation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_skew_sym_source2d(
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nQgc0,nQgc1,
     &     nFgc0,nFgc1,
     &     u0,u1,
     &     Q,
     &     F)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0,ifirst1,ilast1

      integer nugc0,nugc1
      integer nQgc0,nQgc1
      integer nFgc0,nFgc1

      double precision u0(ifirst0-nugc0:ilast0+1+nugc0,
     &          ifirst1-nugc1:ilast1+nugc1)
      double precision u1(ifirst1-nugc1:ilast1+1+nugc1,
     &          ifirst0-nugc0:ilast0+nugc0)

      double precision Q(ifirst0-nQgc0:ilast0+nQgc0,
     &          ifirst1-nQgc1:ilast1+nQgc1)
c
c     Input/Output.
c
      double precision F(ifirst0-nFgc0:ilast0+nFgc0,
     &          ifirst1-nFgc1:ilast1+nFgc1,0:2-1)
c
c     Local variables.
c
      integer ic0,ic1
c
c     Compute the source term F = -0.5 U abs(Q).
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            F(ic0,ic1,0) = -0.25d0*(u0(ic0+1,ic1)+u0(ic0,ic1))*
     &           abs(Q(ic0,ic1))
         enddo
      enddo

      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            F(ic0,ic1,1) = -0.25d0*(u1(ic1+1,ic0)+u1(ic1,ic0))*
     &           abs(Q(ic0,ic1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the source term F = -U max(q,0) that must be added to the
c     momentum equation to account for momentum loss due to internal
c     fluid sources.
c
c     NOTE: This is the source term that corresponds to the advective
c     (i.e., nonconservative) form of the equation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_staggered_adv_source2d(
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nQgc0,nQgc1,
     &     nfgc0,nfgc1,
     &     u0,u1,
     &     Q,
     &     f0,f1)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0,ifirst1,ilast1

      integer nugc0,nugc1
      integer nQgc0,nQgc1
      integer nfgc0,nfgc1

      double precision u0(ifirst0-nugc0:ilast0+1+nugc0,
     &          ifirst1-nugc1:ilast1+nugc1)
      double precision u1(ifirst0-nugc0:ilast0+nugc0,
     &          ifirst1-nugc1:ilast1+1+nugc1)

      double precision Q(ifirst0-nQgc0:ilast0+nQgc0,
     &          ifirst1-nQgc1:ilast1+nQgc1)
c
c     Input/Output.
c
      double precision f0(ifirst0-nfgc0:ilast0+1+nfgc0,
     &          ifirst1-nfgc1:ilast1+nfgc1)
      double precision f1(ifirst0-nfgc0:ilast0+nfgc0,
     &          ifirst1-nfgc1:ilast1+1+nfgc1)
c
c     Local variables.
c
      integer ic0,ic1
c
c     Compute the source term F = -U max(Q,0).
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0+1
            f0(ic0,ic1) = -u0(ic0,ic1)*
     &           dmax1(0.5d0*(Q(ic0-1,ic1)+Q(ic0,ic1)),0.d0)
         enddo
      enddo

      do ic1 = ifirst1,ilast1+1
         do ic0 = ifirst0,ilast0
            f1(ic0,ic1) = -u1(ic0,ic1)*
     &           dmax1(0.5d0*(Q(ic0,ic1-1)+Q(ic0,ic1)),0.d0)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the source term F = +U min(q,0) that must be added to the
c     momentum equation to account for momentum loss due to internal
c     fluid sinks.
c
c     NOTE: This is the source term that corresponds to the conservative
c     form of the equation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_staggered_cons_source2d(
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nQgc0,nQgc1,
     &     nfgc0,nfgc1,
     &     u0,u1,
     &     Q,
     &     f0,f1)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0,ifirst1,ilast1

      integer nugc0,nugc1
      integer nQgc0,nQgc1
      integer nfgc0,nfgc1

      double precision u0(ifirst0-nugc0:ilast0+1+nugc0,
     &          ifirst1-nugc1:ilast1+nugc1)
      double precision u1(ifirst0-nugc0:ilast0+nugc0,
     &          ifirst1-nugc1:ilast1+1+nugc1)

      double precision Q(ifirst0-nQgc0:ilast0+nQgc0,
     &          ifirst1-nQgc1:ilast1+nQgc1)
c
c     Input/Output.
c
      double precision f0(ifirst0-nfgc0:ilast0+1+nfgc0,
     &          ifirst1-nfgc1:ilast1+nfgc1)
      double precision f1(ifirst0-nfgc0:ilast0+nfgc0,
     &          ifirst1-nfgc1:ilast1+1+nfgc1)
c
c     Local variables.
c
      integer ic0,ic1
c
c     Compute the source term F = +U min(Q,0).
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0+1
            f0(ic0,ic1) = +u0(ic0,ic1)*
     &           dmin1(0.5d0*(Q(ic0-1,ic1)+Q(ic0,ic1)),0.d0)
         enddo
      enddo

      do ic1 = ifirst1,ilast1+1
         do ic0 = ifirst0,ilast0
            f1(ic0,ic1) = +u1(ic0,ic1)*
     &           dmin1(0.5d0*(Q(ic0,ic1-1)+Q(ic0,ic1)),0.d0)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the source term F = -0.5 U abs(q) that must be added to
c     the momentum equation to account for momentum loss due to internal
c     fluid sinks.
c
c     NOTE: This is the source term that corresponds to the
c     skew-symmetric form of the equation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine navier_stokes_staggered_skew_sym_source2d(
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nQgc0,nQgc1,
     &     nfgc0,nfgc1,
     &     u0,u1,
     &     Q,
     &     f0,f1)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0,ifirst1,ilast1

      integer nugc0,nugc1
      integer nQgc0,nQgc1
      integer nfgc0,nfgc1

      double precision u0(ifirst0-nugc0:ilast0+1+nugc0,
     &          ifirst1-nugc1:ilast1+nugc1)
      double precision u1(ifirst0-nugc0:ilast0+nugc0,
     &          ifirst1-nugc1:ilast1+1+nugc1)

      double precision Q(ifirst0-nQgc0:ilast0+nQgc0,
     &          ifirst1-nQgc1:ilast1+nQgc1)
c
c     Input/Output.
c
      double precision f0(ifirst0-nfgc0:ilast0+1+nfgc0,
     &          ifirst1-nfgc1:ilast1+nfgc1)
      double precision f1(ifirst0-nfgc0:ilast0+nfgc0,
     &          ifirst1-nfgc1:ilast1+1+nfgc1)
c
c     Local variables.
c
      integer ic0,ic1
c
c     Compute the source term F = - 0.5 U abs(Q).
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0+1
            f0(ic0,ic1) = -0.5d0*u0(ic0,ic1)*
     &           abs(0.5d0*(Q(ic0-1,ic1)+Q(ic0,ic1)))
         enddo
      enddo

      do ic1 = ifirst1,ilast1+1
         do ic0 = ifirst0,ilast0
            f1(ic0,ic1) = -0.5d0*u1(ic0,ic1)*
     &           abs(0.5d0*(Q(ic0,ic1-1)+Q(ic0,ic1)))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
