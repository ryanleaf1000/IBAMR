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
c     Computes W = rot U = [du2/dx1 -du1/dx2, du0/dx2 - du2/dx0, du1/dx0 - du0/dx1].
c
c     Uses centered differences to compute the side centered rot of a
c     edge centered vector field U.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine etosrot3d(
     &     W0,W1,W2,W_gcw,
     &     U0,U1,U2,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx)
c
      implicit none
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer ilower2,iupper2
      integer W_gcw,U_gcw

      double precision U0(ilower0-U_gcw:iupper0+U_gcw,
     &          ilower1-U_gcw:iupper1+1+U_gcw,
     &          ilower2-U_gcw:iupper2+1+U_gcw)
      double precision U1(ilower0-U_gcw:iupper0+1+U_gcw,
     &          ilower1-U_gcw:iupper1+U_gcw,
     &          ilower2-U_gcw:iupper2+1+U_gcw)
      double precision U2(ilower0-U_gcw:iupper0+1+U_gcw,
     &          ilower1-U_gcw:iupper1+1+U_gcw,
     &          ilower2-U_gcw:iupper2+U_gcw)
      double precision dx(0:3-1)
c
c     Input/Output.
c
      double precision W0(ilower0-W_gcw:iupper0+1+W_gcw,
     &          ilower1-W_gcw:iupper1+W_gcw,
     &          ilower2-W_gcw:iupper2+W_gcw)
      double precision W1(ilower0-W_gcw:iupper0+W_gcw,
     &          ilower1-W_gcw:iupper1+1+W_gcw,
     &          ilower2-W_gcw:iupper2+W_gcw)
      double precision W2(ilower0-W_gcw:iupper0+W_gcw,
     &          ilower1-W_gcw:iupper1+W_gcw,
     &          ilower2-W_gcw:iupper2+1+W_gcw)
c
c     Local variables.
c
      integer i0,i1,i2
c
c     Compute the side centered rot of edge centered U.
c
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               W0(i0,i1,i2) = (U2(i0,i1+1,i2)-U2(i0,i1,i2))/dx(1) -
     &                        (U1(i0,i1,i2+1)-U1(i0,i1,i2))/dx(2)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0
               W1(i0,i1,i2) = (U0(i0,i1,i2+1)-U0(i0,i1,i2))/dx(2) -
     &                        (U2(i0+1,i1,i2)-U2(i0,i1,i2))/dx(0)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               W2(i0,i1,i2) = (U1(i0+1,i1,i2)-U1(i0,i1,i2))/dx(0) -
     &                        (U0(i0,i1+1,i2)-U0(i0,i1,i2))/dx(1)
            enddo
         enddo
      enddo

c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

