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
c     Computes W = rot U = [dU/dx1, -dU/dx0].
c
c     Uses centered differences to compute the side centered rot of a
c     node centered scalar field U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ntosrot2d(
     &     W0,W1,W_gcw,
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

      double precision U(ilower0-U_gcw:iupper0+1+U_gcw,
     &          ilower1-U_gcw:iupper1+1+U_gcw)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision W0(ilower0-W_gcw:iupper0+1+W_gcw,
     &          ilower1-W_gcw:iupper1+W_gcw)
      double precision W1(ilower0-W_gcw:iupper0+W_gcw,
     &          ilower1-W_gcw:iupper1+1+W_gcw)
c
c     Local variables.
c
      integer i0,i1
c
c     Compute the side centered rot of node centered U.
c

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            W0(i0,i1) = (U(i0,i1+1)-U(i0,i1))/dx(1)
         enddo
      enddo

      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            W1(i0,i1) = -(U(i0+1,i1)-U(i0,i1))/dx(0)
         enddo
      enddo

c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes W = rot U = [dU/dy, -dU/dx]
c
c     Uses centered differences to compute the side centered rot of a
c     cell centered scalar field U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctosrot2d(
     &     W0,W1,W_gcw,
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
     &          ilower1-U_gcw:iupper1+U_gcw)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision W0(ilower0-W_gcw:iupper0+1+W_gcw,
     &          ilower1-W_gcw:iupper1+W_gcw)
      double precision W1(ilower0-W_gcw:iupper0+W_gcw,
     &          ilower1-W_gcw:iupper1+1+W_gcw)
c
c     Local variables.
c
      integer i0,i1
c
c     Compute the side centered rot of cell centered U.
c

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            W0(i0,i1) = (U(i0-1,i1+1)-U(i0-1,i1-1) + 
     &                   U(i0,i1+1)-U(i0,i1-1))/(0.4d1*dx(1))
         enddo
      enddo

      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            W1(i0,i1) = -(U(i0+1,i1-1)-U(i0-1,i1-1) + 
     &                    U(i0+1,i1)-U(i0-1,i1))/(0.4d1*dx(0))
         enddo
      enddo

c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
