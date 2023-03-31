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
c     Apply the divergence- and gradient-preserving correction to values
c     refined from the next coarser level of the patch hierarchy.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine div_preserving_correction2d(
     &     u0,u1,u_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     correction_box_ilower0,correction_box_iupper0,
     &     correction_box_ilower1,correction_box_iupper1,
     &     ratio,dx_fine)
c
      implicit none
c
c     Input.
c
      integer u_gcw
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer correction_box_ilower0,correction_box_iupper0
      integer correction_box_ilower1,correction_box_iupper1
      integer ratio(0:2-1)
      double precision    dx_fine(0:2-1)
c
c     Input/Output.
c
      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
      double precision u1(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+1+u_gcw)
c
c     Local variables.
c
      integer d,i0,i1,i,j
      double precision u(-1:1,-1:1),u_xx
      double precision v(-1:1,-1:1),v_yy
      double precision dx,dy
c
c     Apply the divergence- and curl-preserving corrections.
c
      do d = 0,2-1
         if ( .not.(ratio(d).eq.2) ) then
            print *,'error: invalid refinement ratio'
            call abort
         endif
      enddo

      dx = dx_fine(0)           ! NOTE: These values are not used in the
      dy = dx_fine(1)           ! 2D code, but are in the 3D routine.

      do i1=correction_box_ilower1,correction_box_iupper1,ratio(1)
         do i0=correction_box_ilower0,correction_box_iupper0,ratio(0)
            u(-1,-1) = u0(i0  ,i1  )
            u( 1,-1) = u0(i0+2,i1  )
            u(-1, 1) = u0(i0  ,i1+1)
            u( 1, 1) = u0(i0+2,i1+1)

            v(-1,-1) = u1(i0  ,i1  )
            v( 1,-1) = u1(i0+1,i1  )
            v(-1, 1) = u1(i0  ,i1+2)
            v( 1, 1) = u1(i0+1,i1+2)

            u_xx = 0.25d0*(v(-1,-1)-v(-1, 1)-v( 1,-1)+v(1,1))
            v_yy = 0.25d0*(u(-1,-1)-u( 1,-1)-u(-1, 1)+u(1,1))

            do j = -1,1,2
               u(0,j) = 0.5d0*(u(1,j)+u(-1,j))+u_xx
            enddo

            u0(i0+1,i1  ) = u(0,-1)
            u0(i0+1,i1+1) = u(0, 1)

            do i = -1,1,2
               v(i,0) = 0.5d0*(v(i,1)+v(i,-1))+v_yy
            enddo

            u1(i0  ,i1+1) = v(-1,0)
            u1(i0+1,i1+1) = v( 1,0)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
