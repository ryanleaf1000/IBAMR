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
c     Computes the advective flux corresponding to a face centered value
c     and a face centered advective velocity.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_flux2d(
     &     dt,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nugc0,nugc1,
     &     nqhalfgc0,nqhalfgc1,
     &     nfluxgc0,nfluxgc1,
     &     u0,u1,
     &     qhalf0,qhalf1,
     &     flux0,flux1)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0,ifirst1,ilast1

      integer nugc0,nugc1
      integer nqhalfgc0,nqhalfgc1
      integer nfluxgc0,nfluxgc1

      double precision dt

      double precision u0(ifirst0-nugc0:ilast0+1+nugc0,
     &          ifirst1-nugc1:ilast1+nugc1)
      double precision u1(ifirst1-nugc1:ilast1+1+nugc1,
     &          ifirst0-nugc0:ilast0+nugc0)

      double precision qhalf0(ifirst0-nqhalfgc0:ilast0+1+nqhalfgc0,
     &          ifirst1-nqhalfgc1:ilast1+nqhalfgc1)
      double precision qhalf1(ifirst1-nqhalfgc1:ilast1+1+nqhalfgc1,
     &          ifirst0-nqhalfgc0:ilast0+nqhalfgc0)
c
c     Input/Output.
c
      double precision flux0(ifirst0-nfluxgc0:ilast0+1+nfluxgc0,
     &          ifirst1-nfluxgc1:ilast1+nfluxgc1)
      double precision flux1(ifirst1-nfluxgc1:ilast1+1+nfluxgc1,
     &          ifirst0-nfluxgc0:ilast0+nfluxgc0)
c
c     Local variables.
c
      integer ic0,ic1
c
c     Compute the time integral of the advective flux.
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0-1,ilast0
            flux0(ic0+1,ic1) = dt*u0(ic0+1,ic1)*qhalf0(ic0+1,ic1)
         enddo
      enddo

      do ic0 = ifirst0,ilast0
         do ic1 = ifirst1-1,ilast1
            flux1(ic1+1,ic0) = dt*u1(ic1+1,ic0)*qhalf1(ic1+1,ic0)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Update a quantity using flux differencing.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_consdiff2d(
     &     dx,
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

      double precision dx(0:2-1)

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
      integer ic0,ic1
c
c     Update a quantity using flux differencing.
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            qval(ic0,ic1) = qval(ic0,ic1)
     &           -(flux0(ic0+1,ic1)-flux0(ic0,ic1))/dx(0)
     &           -(flux1(ic1+1,ic0)-flux1(ic1,ic0))/dx(1)
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
      subroutine advect_consdiffwithdivsource2d(
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
      integer ic0,ic1
      double precision divsource
c
c     Update a quantity using flux differencing.
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            divsource = (fourth/dt)*
     &           ( qflux0(ic0+1,ic1) + qflux0(ic0,ic1)
     &           + qflux1(ic1+1,ic0) + qflux1(ic1,ic0) )*
     &           ( (uflux0(ic0+1,ic1)-uflux0(ic0,ic1))/dx(0)
     &           + (uflux1(ic1+1,ic0)-uflux1(ic1,ic0))/dx(1) )

            qval(ic0,ic1) = qval(ic0,ic1) + divsource
     &           -(flux0(ic0+1,ic1)-flux0(ic0,ic1))/dx(0)
     &           -(flux1(ic1+1,ic0)-flux1(ic1,ic0))/dx(1)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes the advective derivative N = [u_ADV*grad(q)].
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_derivative2d(
     &     dx,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nuadvgc0,nuadvgc1,
     &     nqgc0,nqgc1,
     &     uadv0,uadv1,
     &     q0,q1,
     &     nNgc0,nNgc1,
     &     N)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0,ifirst1,ilast1

      integer nuadvgc0,nuadvgc1
      integer nqgc0,nqgc1
      integer nNgc0,nNgc1

      double precision dx(0:2-1)

      double precision uadv0(ifirst0-nuadvgc0:ilast0+1+nuadvgc0,
     &          ifirst1-nuadvgc1:ilast1+nuadvgc1)
      double precision uadv1(ifirst1-nuadvgc1:ilast1+1+nuadvgc1,
     &          ifirst0-nuadvgc0:ilast0+nuadvgc0)

      double precision q0(ifirst0-nqgc0:ilast0+1+nqgc0,
     &          ifirst1-nqgc1:ilast1+nqgc1)
      double precision q1(ifirst1-nqgc1:ilast1+1+nqgc1,
     &          ifirst0-nqgc0:ilast0+nqgc0)
c
c     Input/Output.
c
      double precision N(ifirst0-nNgc0:ilast0+nNgc0,
     &          ifirst1-nNgc1:ilast1+nNgc1)
c
c     Local variables.
c
      integer ic0,ic1
      double precision U,V
      double precision Qx0,Qx1
c
c     Compute (U,V)*grad(q).
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            U = 0.5d0*(uadv0(ic0+1,ic1)+uadv0(ic0,ic1))
            Qx0 = (q0(ic0+1,ic1)-q0(ic0,ic1))/dx(0)
            N(ic0,ic1) = U*Qx0
         enddo
      enddo

      do ic0 = ifirst0,ilast0
         do ic1 = ifirst1,ilast1
            V = 0.5d0*(uadv1(ic1+1,ic0)+uadv1(ic1,ic0))
            Qx1 = (q1(ic1+1,ic0)-q1(ic1,ic0))/dx(1)
            N(ic0,ic1) = N(ic0,ic1) + V*Qx1
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes the convective derivative N = div[q*u_ADV].
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine convect_derivative2d(
     &     dx,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nuadvgc0,nuadvgc1,
     &     nqgc0,nqgc1,
     &     uadv0,uadv1,
     &     q0,q1,
     &     nNgc0,nNgc1,
     &     N)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0,ifirst1,ilast1

      integer nuadvgc0,nuadvgc1
      integer nqgc0,nqgc1
      integer nNgc0,nNgc1

      double precision dx(0:2-1)

      double precision uadv0(ifirst0-nuadvgc0:ilast0+1+nuadvgc0,
     &          ifirst1-nuadvgc1:ilast1+nuadvgc1)
      double precision uadv1(ifirst1-nuadvgc1:ilast1+1+nuadvgc1,
     &          ifirst0-nuadvgc0:ilast0+nuadvgc0)

      double precision q0(ifirst0-nqgc0:ilast0+1+nqgc0,
     &          ifirst1-nqgc1:ilast1+nqgc1)
      double precision q1(ifirst1-nqgc1:ilast1+1+nqgc1,
     &          ifirst0-nqgc0:ilast0+nqgc0)
c
c     Input/Output.
c
      double precision N(ifirst0-nNgc0:ilast0+nNgc0,
     &          ifirst1-nNgc1:ilast1+nNgc1)
c
c     Local variables.
c
      integer ic0,ic1
      double precision QUx0,QVx1
c
c     Compute div[q*(U,V)].
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            QUx0 = (uadv0(ic0+1,ic1)*q0(ic0+1,ic1)-
     &           uadv0(ic0,ic1)*q0(ic0,ic1))/dx(0)
            N(ic0,ic1) = QUx0
         enddo
      enddo

      do ic0 = ifirst0,ilast0
         do ic1 = ifirst1,ilast1
            QVx1 = (uadv1(ic1+1,ic0)*q1(ic1+1,ic0)-
     &           uadv1(ic1,ic0)*q1(ic1,ic0))/dx(1)
            N(ic0,ic1) = N(ic0,ic1) + QVx1
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes the skew-symmetric derivative N = 0.5([u_ADV*grad(q)] +
c     div[q*u_ADV]).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine skew_sym_derivative2d(
     &     dx,
     &     ifirst0,ilast0,ifirst1,ilast1,
     &     nuadvgc0,nuadvgc1,
     &     nqgc0,nqgc1,
     &     uadv0,uadv1,
     &     q0,q1,
     &     nNgc0,nNgc1,
     &     N)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0,ifirst1,ilast1

      integer nuadvgc0,nuadvgc1
      integer nqgc0,nqgc1
      integer nNgc0,nNgc1

      double precision dx(0:2-1)

      double precision uadv0(ifirst0-nuadvgc0:ilast0+1+nuadvgc0,
     &          ifirst1-nuadvgc1:ilast1+nuadvgc1)
      double precision uadv1(ifirst1-nuadvgc1:ilast1+1+nuadvgc1,
     &          ifirst0-nuadvgc0:ilast0+nuadvgc0)

      double precision q0(ifirst0-nqgc0:ilast0+1+nqgc0,
     &          ifirst1-nqgc1:ilast1+nqgc1)
      double precision q1(ifirst1-nqgc1:ilast1+1+nqgc1,
     &          ifirst0-nqgc0:ilast0+nqgc0)
c
c     Input/Output.
c
      double precision N(ifirst0-nNgc0:ilast0+nNgc0,
     &          ifirst1-nNgc1:ilast1+nNgc1)
c
c     Local variables.
c
      integer ic0,ic1
      double precision U,V
      double precision Qx0,Qx1
      double precision QUx0,QVx1
c
c     Compute 0.5*((U,V)*grad(q) + div[q*(U,V)]).
c
      do ic1 = ifirst1,ilast1
         do ic0 = ifirst0,ilast0
            U = 0.5d0*(uadv0(ic0+1,ic1)+uadv0(ic0,ic1))
            Qx0 = (q0(ic0+1,ic1)-q0(ic0,ic1))/dx(0)
            QUx0 = (uadv0(ic0+1,ic1)*q0(ic0+1,ic1)-
     &           uadv0(ic0,ic1)*q0(ic0,ic1))/dx(0)
            N(ic0,ic1) = 0.5d0*(U*Qx0+QUx0)
         enddo
      enddo

      do ic0 = ifirst0,ilast0
         do ic1 = ifirst1,ilast1
            V = 0.5d0*(uadv1(ic1+1,ic0)+uadv1(ic1,ic0))
            Qx1 = (q1(ic1+1,ic0)-q1(ic1,ic0))/dx(1)
            QVx1 = (uadv1(ic1+1,ic0)*q1(ic1+1,ic0)-
     &           uadv1(ic1,ic0)*q1(ic1,ic0))/dx(1)
            N(ic0,ic1) = N(ic0,ic1) + 0.5d0*(V*Qx1+QVx1)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
