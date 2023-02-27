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
c     Computes W = curl U = dU1/dx0 - dU0/dx1.
c
c     Uses centered differences to compute the cell centered curl of a
c     cell centered vector field U=(U0,U1).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctoccurl2d(
     &     W,W_gcw,
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
      integer W_gcw,U_gcw

      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &          ilower1-U_gcw:iupper1+U_gcw,0:2-1)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision W(ilower0-W_gcw:iupper0+W_gcw,
     &          ilower1-W_gcw:iupper1+W_gcw)
c
c     Local variables.
c
      integer i0,i1
      double precision    dU0_dx1,dU1_dx0,fac01,fac10
c
c     Compute the cell centered curl of U=(U0,U1).
c
      fac01 = 0.5d0/dx(1)
      fac10 = 0.5d0/dx(0)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            dU0_dx1 = fac01*(U(i0  ,i1+1,0)-U(i0  ,i1-1,0))
            dU1_dx0 = fac10*(U(i0+1,i1  ,1)-U(i0-1,i1  ,1))
            W(i0,i1) = dU1_dx0-dU0_dx1
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes W = curl u = du1/dx0 - du0/dx1.
c
c     Uses centered differences to compute the cell centered curl of a
c     face centered vector field u=(u0,u1).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ftoccurl2d(
     &     W,W_gcw,
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
      integer W_gcw,u_gcw

      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
      double precision u1(ilower1-u_gcw:iupper1+1+u_gcw,
     &          ilower0-u_gcw:iupper0+u_gcw)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision W(ilower0-W_gcw:iupper0+W_gcw,
     &          ilower1-W_gcw:iupper1+W_gcw)
c
c     Local variables.
c
      integer i0,i1
      double precision    du0_dx1,du1_dx0,fac01,fac10
c
c     Compute the cell centered curl of u=(u0,u1).
c
      fac01 = 0.25d0/dx(1)
      fac10 = 0.25d0/dx(0)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            du0_dx1 = fac01*(
     &           +u0(i0  ,i1+1)+u0(i0+1,i1+1)
     &           -u0(i0  ,i1-1)-u0(i0+1,i1-1) )
            du1_dx0 = fac10*(
     &           +u1(i1  ,i0+1)+u1(i1+1,i0+1)
     &           -u1(i1  ,i0-1)-u1(i1+1,i0-1) )
            W(i0,i1) = du1_dx0-du0_dx1
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes W = curl u = du1/dx0 - du0/dx1.
c
c     Uses centered differences to compute the cell centered curl of a
c     side centered vector field u=(u0,u1).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stoccurl2d(
     &     W,W_gcw,
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
      integer W_gcw,u_gcw

      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
      double precision u1(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+1+u_gcw)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision W(ilower0-W_gcw:iupper0+W_gcw,
     &          ilower1-W_gcw:iupper1+W_gcw)
c
c     Local variables.
c
      integer i0,i1
      double precision    du0_dx1,du1_dx0,fac01,fac10
c
c     Compute the cell centered curl of u=(u0,u1).
c
      fac01 = 0.25d0/dx(1)
      fac10 = 0.25d0/dx(0)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            du0_dx1 = fac01*(
     &           +u0(i0  ,i1+1)+u0(i0+1,i1+1)
     &           -u0(i0  ,i1-1)-u0(i0+1,i1-1) )
            du1_dx0 = fac10*(
     &           +u1(i0+1,i1  )+u1(i0+1,i1+1)
     &           -u1(i0-1,i1  )-u1(i0-1,i1+1) )
            W(i0,i1) = du1_dx0-du0_dx1
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes W = curl u = du1/dx0 - du0/dx1.
c
c     Uses centered differences to compute the node centered curl of a
c     side centered vector field u=(u0,u1).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stoncurl2d(
     &     W,W_gcw,
     &     U0,U1,U_gcw,
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
      integer W_gcw,U_gcw

      double precision U0(ilower0-U_gcw:iupper0+1+U_gcw,
     &          ilower1-U_gcw:iupper1+U_gcw)
      double precision U1(ilower0-U_gcw:iupper0+U_gcw,
     &          ilower1-U_gcw:iupper1+1+U_gcw)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision W(ilower0-W_gcw:iupper0+1+W_gcw,
     &          ilower1-W_gcw:iupper1+1+W_gcw)
c
c     Local variables.
c
      integer i0,i1
c
c     Compute the node centered curl of U.
c

      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0+1
            W(i0,i1) = (U1(i0,i1) - U1(i0-1,i1))/dx(0) -
     &                 (U0(i0,i1) - U0(i0,i1-1))/dx(1) 
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
