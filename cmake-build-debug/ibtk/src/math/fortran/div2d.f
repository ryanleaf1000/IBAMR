c ---------------------------------------------------------------------
c
c Copyright (c) 2011 - 2019 by the IBAMR developers
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
c     Computes D = alpha div U.
c
c     Uses centered differences to compute the cell centered divergence
c     of a cell centered variable U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocdiv2d(
     &     D,D_gcw,
     &     alpha,
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer D_gcw,U_gcw

      double precision alpha

      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &          ilower1-U_gcw:iupper1+U_gcw,0:2-1)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision D(ilower0-D_gcw:iupper0+D_gcw,
     &          ilower1-D_gcw:iupper1+D_gcw)
c
c     Local variables.
c
      integer i0,i1
      double precision    fac0,fac1
c
c     Compute the cell centered divergence of U.
c
      fac0 = alpha/(2.d0*dx(0))
      fac1 = alpha/(2.d0*dx(1))

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            D(i0,i1) =
     &           fac0*(U(i0+1,i1,0)-U(i0-1,i1,0)) +
     &           fac1*(U(i0,i1+1,1)-U(i0,i1-1,1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes D = alpha div U + beta V.
c
c     Uses centered differences to compute the cell centered divergence
c     of a cell centered variable U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctocdivadd2d(
     &     D,D_gcw,
     &     alpha,
     &     U,U_gcw,
     &     beta,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer D_gcw,U_gcw,V_gcw

      double precision alpha

      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &          ilower1-U_gcw:iupper1+U_gcw,0:2-1)

      double precision beta

      double precision V(ilower0-V_gcw:iupper0+V_gcw,
     &          ilower1-V_gcw:iupper1+V_gcw)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision D(ilower0-D_gcw:iupper0+D_gcw,
     &          ilower1-D_gcw:iupper1+D_gcw)
c
c     Local variables.
c
      integer i0,i1
      double precision    fac0,fac1
c
c     Compute the cell centered divergence of U.
c
      fac0 = alpha/(2.d0*dx(0))
      fac1 = alpha/(2.d0*dx(1))

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            D(i0,i1) =
     &           fac0*(U(i0+1,i1,0)-U(i0-1,i1,0)) +
     &           fac1*(U(i0,i1+1,1)-U(i0,i1-1,1)) +
     &           beta*V(i0,i1)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes D = alpha div u.
c
c     Uses centered differences to compute the cell centered divergence
c     of a face centered variable u=(u0,u1).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ftocdiv2d(
     &     D,D_gcw,
     &     alpha,
     &     u0,u1,u_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer D_gcw,u_gcw

      double precision alpha

      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
      double precision u1(ilower1-u_gcw:iupper1+1+u_gcw,
     &          ilower0-u_gcw:iupper0+u_gcw)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision D(ilower0-D_gcw:iupper0+D_gcw,
     &          ilower1-D_gcw:iupper1+D_gcw)
c
c     Local variables.
c
      integer i0,i1
      double precision    fac0,fac1
c
c     Compute the cell centered divergence of u.
c
      fac0 = alpha/dx(0)
      fac1 = alpha/dx(1)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            D(i0,i1) =
     &           fac0*(u0(i0+1,i1)-u0(i0,i1)) +
     &           fac1*(u1(i1+1,i0)-u1(i1,i0))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes D = alpha div u + beta V.
c
c     Uses centered differences to compute the cell centered divergence
c     of a face centered variable u=(u0,u1).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ftocdivadd2d(
     &     D,D_gcw,
     &     alpha,
     &     u0,u1,u_gcw,
     &     beta,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer D_gcw,u_gcw,V_gcw

      double precision alpha

      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
      double precision u1(ilower1-u_gcw:iupper1+1+u_gcw,
     &          ilower0-u_gcw:iupper0+u_gcw)

      double precision beta

      double precision V(ilower0-V_gcw:iupper0+V_gcw,
     &          ilower1-V_gcw:iupper1+V_gcw)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision D(ilower0-D_gcw:iupper0+D_gcw,
     &          ilower1-D_gcw:iupper1+D_gcw)
c
c     Local variables.
c
      integer i0,i1
      double precision    fac0,fac1
c
c     Compute the cell centered divergence of u.
c
      fac0 = alpha/dx(0)
      fac1 = alpha/dx(1)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            D(i0,i1) =
     &           fac0*(u0(i0+1,i1)-u0(i0,i1)) +
     &           fac1*(u1(i1+1,i0)-u1(i1,i0)) +
     &           beta*V(i0,i1)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes D = alpha div u.
c
c     Uses centered differences to compute the cell centered divergence
c     of a side centered variable u=(u0,u1).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stocdiv2d(
     &     D,D_gcw,
     &     alpha,
     &     u0,u1,u_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer D_gcw,u_gcw

      double precision alpha

      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
      double precision u1(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+1+u_gcw)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision D(ilower0-D_gcw:iupper0+D_gcw,
     &          ilower1-D_gcw:iupper1+D_gcw)
c
c     Local variables.
c
      integer i0,i1
      double precision    fac0,fac1
c
c     Compute the cell centered divergence of u.
c
      fac0 = alpha/dx(0)
      fac1 = alpha/dx(1)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            D(i0,i1) =
     &        fac0*(u0(i0+1,i1)-u0(i0,i1)) +
     &        fac1*(u1(i0,i1+1)-u1(i0,i1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes D = alpha div u + beta V.
c
c     Uses centered differences to compute the cell centered divergence
c     of a side centered variable u=(u0,u1).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stocdivadd2d(
     &     D,D_gcw,
     &     alpha,
     &     u0,u1,u_gcw,
     &     beta,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer D_gcw,u_gcw,V_gcw

      double precision alpha

      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
      double precision u1(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+1+u_gcw)

      double precision beta

      double precision V(ilower0-V_gcw:iupper0+V_gcw,
     &          ilower1-V_gcw:iupper1+V_gcw)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision D(ilower0-D_gcw:iupper0+D_gcw,
     &          ilower1-D_gcw:iupper1+D_gcw)
c
c     Local variables.
c
      integer i0,i1
      double precision    fac0,fac1
c
c     Compute the cell centered divergence of u.
c
      fac0 = alpha/dx(0)
      fac1 = alpha/dx(1)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            D(i0,i1) =
     &           fac0*(u0(i0+1,i1)-u0(i0,i1)) +
     &           fac1*(u1(i0,i1+1)-u1(i0,i1)) +
     &           beta*V(i0,i1)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
