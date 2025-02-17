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
c     Apply the divergence- and gradient-preserving correction to values
c     refined from the next coarser level of the patch hierarchy.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine div_preserving_correction3d(
     &     u0,u1,u2,u_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     correction_box_ilower0,correction_box_iupper0,
     &     correction_box_ilower1,correction_box_iupper1,
     &     correction_box_ilower2,correction_box_iupper2,
     &     ratio,dx_fine)
c
      implicit none
c
c     Input.
c
      integer u_gcw
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer ilower2,iupper2
      integer correction_box_ilower0,correction_box_iupper0
      integer correction_box_ilower1,correction_box_iupper1
      integer correction_box_ilower2,correction_box_iupper2
      integer ratio(0:3-1)
      double precision    dx_fine(0:3-1)
c
c     Input/Output.
c
      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw,
     &          ilower2-u_gcw:iupper2+u_gcw)
      double precision u1(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+1+u_gcw,
     &          ilower2-u_gcw:iupper2+u_gcw)
      double precision u2(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw,
     &          ilower2-u_gcw:iupper2+1+u_gcw)
c
c     Local variables.
c
      integer d,i0,i1,i2,i,j,k
      double precision u(-2:2,-1:1,-1:1),u_xx,u_xyz
      double precision v(-1:1,-2:2,-1:1),v_yy,v_xyz
      double precision w(-1:1,-1:1,-2:2),w_zz,w_xyz
      double precision dx,dy,dz
c
c     Apply the divergence- and curl-preserving corrections.
c
      do d = 0,3-1
         if ( .not.(ratio(d).eq.2) ) then
            print *,'error: invalid refinement ratio'
            call abort
         endif
      enddo

      dx = dx_fine(0)
      dy = dx_fine(1)
      dz = dx_fine(2)

      do i2=correction_box_ilower2,correction_box_iupper2,ratio(2)
         do i1=correction_box_ilower1,correction_box_iupper1,ratio(1)
            do i0=correction_box_ilower0,correction_box_iupper0,ratio(0)

               u(-2,-1,-1) = u0(i0  ,i1  ,i2  )
               u( 2,-1,-1) = u0(i0+2,i1  ,i2  )
               u(-2, 1,-1) = u0(i0  ,i1+1,i2  )
               u( 2, 1,-1) = u0(i0+2,i1+1,i2  )
               u(-2,-1, 1) = u0(i0  ,i1  ,i2+1)
               u( 2,-1, 1) = u0(i0+2,i1  ,i2+1)
               u(-2, 1, 1) = u0(i0  ,i1+1,i2+1)
               u( 2, 1, 1) = u0(i0+2,i1+1,i2+1)

               v(-1,-2,-1) = u1(i0  ,i1  ,i2  )
               v( 1,-2,-1) = u1(i0+1,i1  ,i2  )
               v(-1, 2,-1) = u1(i0  ,i1+2,i2  )
               v( 1, 2,-1) = u1(i0+1,i1+2,i2  )
               v(-1,-2, 1) = u1(i0  ,i1  ,i2+1)
               v( 1,-2, 1) = u1(i0+1,i1  ,i2+1)
               v(-1, 2, 1) = u1(i0  ,i1+2,i2+1)
               v( 1, 2, 1) = u1(i0+1,i1+2,i2+1)

               w(-1,-1,-2) = u2(i0  ,i1  ,i2  )
               w( 1,-1,-2) = u2(i0+1,i1  ,i2  )
               w(-1, 1,-2) = u2(i0  ,i1+1,i2  )
               w( 1, 1,-2) = u2(i0+1,i1+1,i2  )
               w(-1,-1, 2) = u2(i0  ,i1  ,i2+2)
               w( 1,-1, 2) = u2(i0+1,i1  ,i2+2)
               w(-1, 1, 2) = u2(i0  ,i1+1,i2+2)
               w( 1, 1, 2) = u2(i0+1,i1+1,i2+2)

               u_xx = 0.d0;
               do k = -1,1,2
                  do j = -1,1,2
                     do i = -1,1,2
                        u_xx = u_xx + 0.125d0*(
     &                       dble(i)*dble(j)*v(i,2*j,k) +
     &                       dble(i)*dble(k)*w(i,j,2*k) )
                     enddo
                  enddo
               enddo

               v_yy = 0.d0;
               do k = -1,1,2
                  do j = -1,1,2
                     do i = -1,1,2
                        v_yy = v_yy + 0.125d0*(
     &                       dble(i)*dble(j)*u(2*i,j,k) +
     &                       dble(j)*dble(k)*w(i,j,2*k) )
                     enddo
                  enddo
               enddo

               w_zz = 0.d0;
               do k = -1,1,2
                  do j = -1,1,2
                     do i = -1,1,2
                        w_zz = w_zz + 0.125d0*(
     &                       dble(i)*dble(k)*u(2*i,j,k) +
     &                       dble(j)*dble(k)*v(i,2*j,k) )
                     enddo
                  enddo
               enddo

               u_xyz = 0.d0
               do k = -1,1,2
                  do j = -1,1,2
                     do i = -1,1,2
                        u_xyz = u_xyz + 0.125d0*(
     &                       dble(i)*dble(j)*dble(k)*u(2*i,j,k)/
     &                       (dy**2.d0 + dz**2.d0))
                     enddo
                  enddo
               enddo

               v_xyz = 0.d0
               do k = -1,1,2
                  do j = -1,1,2
                     do i = -1,1,2
                        v_xyz = v_xyz + 0.125d0*(
     &                       dble(i)*dble(j)*dble(k)*v(i,2*j,k)/
     &                       (dx**2.d0 + dz**2.d0))
                     enddo
                  enddo
               enddo

               w_xyz = 0.d0
               do k = -1,1,2
                  do j = -1,1,2
                     do i = -1,1,2
                        w_xyz = w_xyz + 0.125d0*(
     &                       dble(i)*dble(j)*dble(k)*w(i,j,2*k)/
     &                       (dx**2.d0 + dy**2.d0))
                     enddo
                  enddo
               enddo

               do k = -1,1,2
                  do j = -1,1,2
                     u(0,j,k) = 0.5d0*(u(2,j,k)+u(-2,j,k)) + u_xx +
     &                    dble(k)*dz**2.d0*v_xyz +
     &                    dble(j)*dy**2.d0*w_xyz
                  enddo
               enddo

               do k = -1,1,2
                  do i = -1,1,2
                     v(i,0,k) = 0.5d0*(v(i,2,k)+v(i,-2,k)) + v_yy +
     &                    dble(i)*dx**2.d0*w_xyz +
     &                    dble(k)*dz**2.d0*u_xyz
                  enddo
               enddo

               do j = -1,1,2
                  do i = -1,1,2
                     w(i,j,0) = 0.5d0*(w(i,j,2)+w(i,j,-2)) + w_zz +
     &                    dble(j)*dy**2.d0*u_xyz +
     &                    dble(i)*dx**2.d0*v_xyz
                  enddo
               enddo

               u0(i0+1,i1  ,i2  ) = u(0,-1,-1)
               u0(i0+1,i1+1,i2  ) = u(0, 1,-1)
               u0(i0+1,i1  ,i2+1) = u(0,-1, 1)
               u0(i0+1,i1+1,i2+1) = u(0, 1, 1)

               u1(i0  ,i1+1,i2  ) = v(-1,0,-1)
               u1(i0+1,i1+1,i2  ) = v( 1,0,-1)
               u1(i0  ,i1+1,i2+1) = v(-1,0, 1)
               u1(i0+1,i1+1,i2+1) = v( 1,0, 1)

               u2(i0  ,i1  ,i2+1) = w(-1,-1,0)
               u2(i0+1,i1  ,i2+1) = w( 1,-1,0)
               u2(i0  ,i1+1,i2+1) = w(-1, 1,0)
               u2(i0+1,i1+1,i2+1) = w( 1, 1,0)

            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
