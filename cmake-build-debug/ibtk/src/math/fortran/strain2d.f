c ---------------------------------------------------------------------
c
c Copyright (c) 2016 - 2019 by the IBAMR developers
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
c     Computes E_diag    = (du0/dx0, du1/dx1).
c     Computes E_offDiag = 0.5*(du0/dx1 + du1/dx0)
c
c     Uses centered differences to compute the cell centered strain of a
c     side centered vector field u=(u0,u1).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stocstrain2d(
     &     E_diag,E_diag_gcw,
     &     E_off,E_off_gcw,
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
      integer E_diag_gcw,E_off_gcw,u_gcw

      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
      double precision u1(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+1+u_gcw)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision E_diag(ilower0-E_diag_gcw:iupper0+E_diag_gcw,
     &          ilower1-E_diag_gcw:iupper1+E_diag_gcw,0:2-1)
      double precision E_off(ilower0-E_off_gcw:iupper0+E_off_gcw,
     &          ilower1-E_off_gcw:iupper1+E_off_gcw)
c
c     Local variables.
c
      integer i0,i1
      double precision    du0_dx0,du1_dx1,du0_dx1,du1_dx0
      double precision    fac00,fac11,fac01,fac10
c
c     Compute the cell centered stress of u=(u0,u1).
c
      fac00 = 1.0d0/dx(0)
      fac11 = 1.0d0/dx(1)
      fac01 = 0.25d0/dx(1)
      fac10 = 0.25d0/dx(0)

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            du0_dx0 = fac00*(u0(i0+1,i1) - u0(i0,i1))
            du1_dx1 = fac11*(u1(i0,i1+1) - u1(i0,i1))
            
            du0_dx1 = fac01*(
     &           +u0(i0  ,i1+1)+u0(i0+1,i1+1)
     &           -u0(i0  ,i1-1)-u0(i0+1,i1-1) )
            du1_dx0 = fac10*(
     &           +u1(i0+1,i1  )+u1(i0+1,i1+1)
     &           -u1(i0-1,i1  )-u1(i0-1,i1+1) )
     
            E_diag(i0,i1,0) = du0_dx0
            E_diag(i0,i1,1) = du1_dx1
            
            E_off(i0,i1) = 0.5d0*(du1_dx0+du0_dx1)
         enddo
      enddo
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
