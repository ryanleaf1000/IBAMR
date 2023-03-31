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
c     Compute the divergence form of the convection term corresponding
c     to the given velocity field and concentration.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_centered_side_div_derivative2d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     n_U_gc0,n_U_gc1,
     &     U0,U1,
     &     n_C_gc0,n_C_gc1,
     &     C,
     &     n_N_gc0,n_N_gc1,
     &     N)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0
      integer ifirst1,ilast1

      integer n_U_gc0,n_U_gc1
      integer n_C_gc0,n_C_gc1
      integer n_N_gc0,n_N_gc1

      double precision dx(0:2-1)

      double precision U0(
     &     ifirst0-n_U_gc0:ilast0+1+n_U_gc0,
     &          ifirst1-n_U_gc1:ilast1+n_U_gc1
     &     )
      double precision U1(
     &     ifirst0-n_U_gc0:ilast0+n_U_gc0,
     &          ifirst1-n_U_gc1:ilast1+1+n_U_gc1
     &     )

      double precision C(
     &     ifirst0-n_C_gc0:ilast0+n_C_gc0,
     &          ifirst1-n_C_gc1:ilast1+n_C_gc1
     &     )
c
c     Output.
c
      double precision N(
     &     ifirst0-n_N_gc0:ilast0+n_N_gc0,
     &          ifirst1-n_N_gc1:ilast1+n_N_gc1
     &     )
c
c     Local variables.
c
      integer i0,i1
      double precision D,D_x,D_y
c
c     Compute N = div(UC).
c
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
            D_x = (0.5d0/dx(0))*(
     &           u0(i0+1,i1)*((C(i0+1,i1)+C(i0  ,i1))) -
     &           u0(i0  ,i1)*((C(i0  ,i1)+C(i0-1,i1))))

            D_y = (0.5d0/dx(1))*(
     &           u1(i0,i1+1)*((C(i0,i1+1)+C(i0,i1  ))) -
     &           u1(i0,i1  )*((C(i0,i1  )+C(i0,i1-1))))

            D = D_x + D_y

            N(i0,i1) = D
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the advective form of the convection term corresponding to
c     the given velocity field and concentration.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_centered_side_adv_derivative2d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     n_U_gc0,n_U_gc1,
     &     U0,U1,
     &     n_C_gc0,n_C_gc1,
     &     C,
     &     n_N_gc0,n_N_gc1,
     &     N)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0
      integer ifirst1,ilast1

      integer n_U_gc0,n_U_gc1
      integer n_C_gc0,n_C_gc1
      integer n_N_gc0,n_N_gc1

      double precision dx(0:2-1)

      double precision U0(
     &     ifirst0-n_U_gc0:ilast0+1+n_U_gc0,
     &          ifirst1-n_U_gc1:ilast1+n_U_gc1
     &     )
      double precision U1(
     &     ifirst0-n_U_gc0:ilast0+n_U_gc0,
     &          ifirst1-n_U_gc1:ilast1+1+n_U_gc1
     &     )

      double precision C(
     &     ifirst0-n_C_gc0:ilast0+n_C_gc0,
     &          ifirst1-n_C_gc1:ilast1+n_C_gc1
     &     )
c
c     Output.
c
      double precision N(
     &     ifirst0-n_N_gc0:ilast0+n_N_gc0,
     &          ifirst1-n_N_gc1:ilast1+n_N_gc1
     &     )
c
c     Local variables.
c
      integer i0,i1
      double precision A,A_x,A_y
c
c     Compute N = (U*grad)C.
c
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
            A_x = (0.5d0/dx(0))*(
     &           u0(i0+1,i1)*((C(i0+1,i1)-C(i0  ,i1))) +
     &           u0(i0  ,i1)*((C(i0  ,i1)-C(i0-1,i1))))

            A_y = (0.5d0/dx(1))*(
     &           u1(i0,i1+1)*((C(i0,i1+1)-C(i0,i1  ))) +
     &           u1(i0,i1  )*((C(i0,i1  )-C(i0,i1-1))))

            A = A_x + A_y

            N(i0,i1) = A
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the skew-symmetric form of the convection term
c     corresponding to the given velocity field and concentration.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_centered_side_skew_sym_derivative2d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     n_U_gc0,n_U_gc1,
     &     U0,U1,
     &     n_C_gc0,n_C_gc1,
     &     C,
     &     n_N_gc0,n_N_gc1,
     &     N)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0
      integer ifirst1,ilast1

      integer n_U_gc0,n_U_gc1
      integer n_C_gc0,n_C_gc1
      integer n_N_gc0,n_N_gc1

      double precision dx(0:2-1)

      double precision U0(
     &     ifirst0-n_U_gc0:ilast0+1+n_U_gc0,
     &          ifirst1-n_U_gc1:ilast1+n_U_gc1
     &     )
      double precision U1(
     &     ifirst0-n_U_gc0:ilast0+n_U_gc0,
     &          ifirst1-n_U_gc1:ilast1+1+n_U_gc1
     &     )

      double precision C(
     &     ifirst0-n_C_gc0:ilast0+n_C_gc0,
     &          ifirst1-n_C_gc1:ilast1+n_C_gc1
     &     )
c
c     Output.
c
      double precision N(
     &     ifirst0-n_N_gc0:ilast0+n_N_gc0,
     &          ifirst1-n_N_gc1:ilast1+n_N_gc1
     &     )
c
c     Local variables.
c
      integer i0,i1
      double precision A,A_x,A_y
      double precision D,D_x,D_y
c
c     Compute N = 0.5*(div(UC) + (U*grad)C).
c
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
            D_x = (0.5d0/dx(0))*(
     &           u0(i0+1,i1)*((C(i0+1,i1)+C(i0  ,i1))) -
     &           u0(i0  ,i1)*((C(i0  ,i1)+C(i0-1,i1))))

            D_y = (0.5d0/dx(1))*(
     &           u1(i0,i1+1)*((C(i0,i1+1)+C(i0,i1  ))) -
     &           u1(i0,i1  )*((C(i0,i1  )+C(i0,i1-1))))

            D = D_x + D_y

            A_x = (0.5d0/dx(0))*(
     &           u0(i0+1,i1)*((C(i0+1,i1)-C(i0  ,i1))) +
     &           u0(i0  ,i1)*((C(i0  ,i1)-C(i0-1,i1))))

            A_y = (0.5d0/dx(1))*(
     &           u1(i0,i1+1)*((C(i0,i1+1)-C(i0,i1  ))) +
     &           u1(i0,i1  )*((C(i0,i1  )-C(i0,i1-1))))

            A = A_x + A_y

            N(i0,i1) = 0.5d0*(D+A)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the divergence form of the convection term corresponding
c     to the given velocity field and concentration.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_centered_face_div_derivative2d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     n_U_gc0,n_U_gc1,
     &     U0,U1,
     &     n_C_gc0,n_C_gc1,
     &     C,
     &     n_N_gc0,n_N_gc1,
     &     N)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0
      integer ifirst1,ilast1

      integer n_U_gc0,n_U_gc1
      integer n_C_gc0,n_C_gc1
      integer n_N_gc0,n_N_gc1

      double precision dx(0:2-1)

      double precision U0(
     &     ifirst0-n_U_gc0:ilast0+1+n_U_gc0,
     &          ifirst1-n_U_gc1:ilast1+n_U_gc1
     &     )
      double precision U1(
     &     ifirst1-n_U_gc1:ilast1+1+n_U_gc1,
     &          ifirst0-n_U_gc0:ilast0+n_U_gc0
     &     )

      double precision C(
     &     ifirst0-n_C_gc0:ilast0+n_C_gc0,
     &          ifirst1-n_C_gc1:ilast1+n_C_gc1
     &     )
c
c     Output.
c
      double precision N(
     &     ifirst0-n_N_gc0:ilast0+n_N_gc0,
     &          ifirst1-n_N_gc1:ilast1+n_N_gc1
     &     )
c
c     Local variables.
c
      integer i0,i1
      double precision D,D_x,D_y
c
c     Compute N = div(UC).
c
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
            D_x = (0.5d0/dx(0))*(
     &           u0(i0+1,i1)*((C(i0+1,i1)+C(i0  ,i1))) -
     &           u0(i0  ,i1)*((C(i0  ,i1)+C(i0-1,i1))))

            D_y = (0.5d0/dx(1))*(
     &           u1(i1+1,i0)*((C(i0,i1+1)+C(i0,i1  ))) -
     &           u1(i1  ,i0)*((C(i0,i1  )+C(i0,i1-1))))

            D = D_x + D_y

            N(i0,i1) = D
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the advective form of the convection term corresponding to
c     the given velocity field and concentration.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_centered_face_adv_derivative2d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     n_U_gc0,n_U_gc1,
     &     U0,U1,
     &     n_C_gc0,n_C_gc1,
     &     C,
     &     n_N_gc0,n_N_gc1,
     &     N)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0
      integer ifirst1,ilast1

      integer n_U_gc0,n_U_gc1
      integer n_C_gc0,n_C_gc1
      integer n_N_gc0,n_N_gc1

      double precision dx(0:2-1)

      double precision U0(
     &     ifirst0-n_U_gc0:ilast0+1+n_U_gc0,
     &          ifirst1-n_U_gc1:ilast1+n_U_gc1
     &     )
      double precision U1(
     &     ifirst1-n_U_gc1:ilast1+1+n_U_gc1,
     &          ifirst0-n_U_gc0:ilast0+n_U_gc0
     &     )

      double precision C(
     &     ifirst0-n_C_gc0:ilast0+n_C_gc0,
     &          ifirst1-n_C_gc1:ilast1+n_C_gc1
     &     )
c
c     Output.
c
      double precision N(
     &     ifirst0-n_N_gc0:ilast0+n_N_gc0,
     &          ifirst1-n_N_gc1:ilast1+n_N_gc1
     &     )
c
c     Local variables.
c
      integer i0,i1
      double precision A,A_x,A_y
c
c     Compute N = (U*grad)C.
c
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
            A_x = (0.5d0/dx(0))*(
     &           u0(i0+1,i1)*((C(i0+1,i1)-C(i0  ,i1))) +
     &           u0(i0  ,i1)*((C(i0  ,i1)-C(i0-1,i1))))

            A_y = (0.5d0/dx(1))*(
     &           u1(i1+1,i0)*((C(i0,i1+1)-C(i0,i1  ))) +
     &           u1(i1  ,i0)*((C(i0,i1  )-C(i0,i1-1))))

            A = A_x + A_y

            N(i0,i1) = A
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the skew-symmetric form of the convection term
c     corresponding to the given velocity field and concentration.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_centered_face_skew_sym_derivative2d(
     &     dx,
     &     ifirst0,ilast0,
     &     ifirst1,ilast1,
     &     n_U_gc0,n_U_gc1,
     &     U0,U1,
     &     n_C_gc0,n_C_gc1,
     &     C,
     &     n_N_gc0,n_N_gc1,
     &     N)
c
      implicit none
c
c     Input.
c
      integer ifirst0,ilast0
      integer ifirst1,ilast1

      integer n_U_gc0,n_U_gc1
      integer n_C_gc0,n_C_gc1
      integer n_N_gc0,n_N_gc1

      double precision dx(0:2-1)

      double precision U0(
     &     ifirst0-n_U_gc0:ilast0+1+n_U_gc0,
     &          ifirst1-n_U_gc1:ilast1+n_U_gc1
     &     )
      double precision U1(
     &     ifirst1-n_U_gc1:ilast1+1+n_U_gc1,
     &          ifirst0-n_U_gc0:ilast0+n_U_gc0
     &     )

      double precision C(
     &     ifirst0-n_C_gc0:ilast0+n_C_gc0,
     &          ifirst1-n_C_gc1:ilast1+n_C_gc1
     &     )
c
c     Output.
c
      double precision N(
     &     ifirst0-n_N_gc0:ilast0+n_N_gc0,
     &          ifirst1-n_N_gc1:ilast1+n_N_gc1
     &     )
c
c     Local variables.
c
      integer i0,i1
      double precision A,A_x,A_y
      double precision D,D_x,D_y
c
c     Compute N = 0.5*(div(UC) + (U*grad)C).
c
      do i1 = ifirst1,ilast1
         do i0 = ifirst0,ilast0
            D_x = (0.5d0/dx(0))*(
     &           u0(i0+1,i1)*((C(i0+1,i1)+C(i0  ,i1))) -
     &           u0(i0  ,i1)*((C(i0  ,i1)+C(i0-1,i1))))

            D_y = (0.5d0/dx(1))*(
     &           u1(i1+1,i0)*((C(i0,i1+1)+C(i0,i1  ))) -
     &           u1(i1  ,i0)*((C(i0,i1  )+C(i0,i1-1))))

            D = D_x + D_y

            A_x = (0.5d0/dx(0))*(
     &           u0(i0+1,i1)*((C(i0+1,i1)-C(i0  ,i1))) +
     &           u0(i0  ,i1)*((C(i0  ,i1)-C(i0-1,i1))))

            A_y = (0.5d0/dx(1))*(
     &           u1(i1+1,i0)*((C(i0,i1+1)-C(i0,i1  ))) +
     &           u1(i1  ,i0)*((C(i0,i1  )-C(i0,i1-1))))

            A = A_x + A_y

            N(i0,i1) = 0.5d0*(D+A)
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
