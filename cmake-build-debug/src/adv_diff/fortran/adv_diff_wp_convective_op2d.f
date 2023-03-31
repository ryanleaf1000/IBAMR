c ---------------------------------------------------------------------
c
c Copyright (c) 2018 - 2019 by the IBAMR developers
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes r = u.grad(q)
c
c     where u is vector valued face centered velocity
c     q is cell centered with depth d
c     returns r_data at cell centeres
c     computes grad(q) using weno + wave propagation
c     interpolation coefficients and weights must be provided
c     currently only works for interp orders 3 (k=2) and 5 (k=3)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine adv_diff_wp_convective_op2d(
     &            q_data, q_gcw,
     &            u_data_0, u_data_1, u_gcw,
     &            r_data, r_gcw, d,
     &            ilower0, ilower1,
     &            iupper0, iupper1,
     &            dx, k)

      implicit none
      integer ilower0, iupper0
      integer ilower1, iupper1
      integer d

      integer q_gcw
      double precision q_data(ilower0-q_gcw:iupper0+q_gcw,
     &          ilower1-q_gcw:iupper1+q_gcw,0:(d-1))

      double precision s_data_0(ilower0:iupper0+1,
     &          ilower1:iupper1,0:1)
      double precision s_data_1(ilower1:iupper1+1,
     &          ilower0:iupper0,0:1)

      integer u_gcw
      double precision u_data_0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
      double precision u_data_1(ilower1-u_gcw:iupper1+1+u_gcw,
     &          ilower0-u_gcw:iupper0+u_gcw)

      integer r_gcw
      double precision r_data(ilower0-r_gcw:iupper0+r_gcw,
     &          ilower1-r_gcw:iupper1+r_gcw,0:(d-1))

      double precision dx(0:1)

      integer k, j
      integer i0, i1


      do j=0,(d-1)
      call reconstruct_data_on_patch_2d(q_data(:,:,j), q_gcw,
     &             s_data_0, s_data_1, 0,
     &             ilower0, ilower1, iupper0, iupper1, k)

      do i1 = ilower1, iupper1
        do i0 = ilower0, iupper0
         r_data(i0,i1,j) =
     &     1.d0/dx(0)*(max(u_data_0(i0,i1),0.d0)*
     &     (s_data_0(i0,i1,1)-s_data_0(i0,i1,0))
     &     + min(u_data_0(i0+1,i1),0.d0)*
     &     (s_data_0(i0+1,i1,1)-s_data_0(i0+1,i1,0))
     &     + 0.5d0*(u_data_0(i0+1,i1)+u_data_0(i0,i1))*
     &     (s_data_0(i0+1,i1,0)-s_data_0(i0,i1,1)))

         r_data(i0,i1,j) = r_data(i0,i1,j) +
     &     1.d0/dx(1)*(max(u_data_1(i1,i0),0.d0)*
     &     (s_data_1(i1,i0,1)-s_data_1(i1,i0,0))
     &     + min(u_data_1(i1+1,i0),0.d0)*
     &     (s_data_1(i1+1,i0,1)-s_data_1(i1+1,i0,0))
     &     + 0.5d0*(u_data_1(i1+1,i0)+u_data_1(i1,i0))*
     &     (s_data_1(i1+1,i0,0)-s_data_1(i1,i0,1)))
        enddo
      enddo
      enddo
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Reconstructs data on patches using a weno scheme
c       the convex and interpolation weights must be supplied
c
c       q_data is cell centered with depth 1
c       r_data_* are face centered with depth 2
c         and return the values reconstructed from each side
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine reconstruct_data_on_patch_2d(q_data, q_gcw,
     &            r_data_0, r_data_1, r_gcw,
     &            ilower0, ilower1, iupper0, iupper1, k)

       integer ilower0, iupper0
       integer ilower1, iupper1

       integer q_gcw
       double precision q_data(ilower0-q_gcw:iupper0+q_gcw,
     &          ilower1-q_gcw:iupper1+q_gcw)

       integer r_gcw
       double precision r_data_0(ilower0-r_gcw:iupper0+1+r_gcw,
     &          ilower1-r_gcw:iupper1+r_gcw,0:1)
       double precision r_data_1(ilower1-r_gcw:iupper1+1+r_gcw,
     &          ilower0-r_gcw:iupper0+r_gcw,0:1)


       integer i0, i1

       double precision WENO5_interp
c
c     Prevent compiler warning about unused variables.
c
      k = k

c     X DIRECTION
      do i1 = ilower1, iupper1
        do i0=ilower0,iupper0+1
           r_data_0(i0,i1,1) = WENO5_interp(q_data(i0+2:i0-2:-1,i1))
           r_data_0(i0,i1,0) = WENO5_interp(q_data(i0-3:i0+1,i1))
         enddo
       enddo

c      Y DIRECTION
       do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
           r_data_1(i1,i0,1) = WENO5_interp(q_data(i0,i1+2:i1-2:-1))
           r_data_1(i1,i0,0) = WENO5_interp(q_data(i0,i1-3:i1+1))
         enddo
       enddo
       end subroutine
