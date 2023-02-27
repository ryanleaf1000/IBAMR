c ---------------------------------------------------------------------
c
c Copyright (c) 2011 - 2022 by the IBAMR developers
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
c     Perform a single Gauss-Seidel sweep for F = D div grad U +
c     C U. Both D and C coefficients are constant.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine smooth_gs_const_dc_3d(
     &     U,U_gcw,
     &     D,C,
     &     F,F_gcw,
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
      integer U_gcw,F_gcw

      double precision D,C

      double precision F(ilower0-F_gcw:iupper0+F_gcw,
     &     ilower1-F_gcw:iupper1+F_gcw,
     &     ilower2-F_gcw:iupper2+F_gcw)

      double precision dx(0:3-1)
c
c     Input/Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &     ilower1-U_gcw:iupper1+U_gcw,
     &     ilower2-U_gcw:iupper2+U_gcw)
c
c     Local variables.
c
      integer i0,i1,i2
      double precision    fac0,fac1,fac2,fac
c
c     Perform a single Gauss-Seidel sweep.
c
      fac0 = D/(dx(0)*dx(0))
      fac1 = D/(dx(1)*dx(1))
      fac2 = D/(dx(2)*dx(2))
      fac = 0.5d0/(fac0+fac1+fac2-0.5d0*C)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               U(i0,i1,i2) = fac*(
     &              fac0*(U(i0-1,i1,i2)+U(i0+1,i1,i2)) +
     &              fac1*(U(i0,i1-1,i2)+U(i0,i1+1,i2)) +
     &              fac2*(U(i0,i1,i2-1)+U(i0,i1,i2+1)) -
     &              F(i0,i1,i2))
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single "red" or "black" Gauss-Seidel sweep for F = D
c     div grad U + C U. Both D and C coefficients
c     are constant.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine smooth_gs_rb_const_dc_3d(
     &     U,U_gcw,
     &     D,C,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     red_or_black)
c
      implicit none
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer ilower2,iupper2
      integer U_gcw,F_gcw
      integer red_or_black

      double precision D,C

      double precision F(ilower0-F_gcw:iupper0+F_gcw,
     &     ilower1-F_gcw:iupper1+F_gcw,
     &     ilower2-F_gcw:iupper2+F_gcw)

      double precision dx(0:3-1)
c
c     Input/Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &     ilower1-U_gcw:iupper1+U_gcw,
     &     ilower2-U_gcw:iupper2+U_gcw)
c
c     Local variables.
c
      integer i0,i1,i2
      double precision    fac0,fac1,fac2,fac
c
c     Perform a single "red" or "black" Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1

      fac0 = D/(dx(0)*dx(0))
      fac1 = D/(dx(1)*dx(1))
      fac2 = D/(dx(2)*dx(2))
      fac = 0.5d0/(fac0+fac1+fac2-0.5d0*C)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               if ( mod(i0+i1+i2,2) .eq. red_or_black ) then
                  U(i0,i1,i2) = fac*(
     &                 fac0*(U(i0-1,i1,i2)+U(i0+1,i1,i2)) +
     &                 fac1*(U(i0,i1-1,i2)+U(i0,i1+1,i2)) +
     &                 fac2*(U(i0,i1,i2-1)+U(i0,i1,i2+1)) -
     &                 F(i0,i1,i2))
               endif
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single Gauss-Seidel sweep for F = alpha div grad U +
c     beta U with masking of certain degrees of freedom.
c
c     NOTE: The solution U is unmodified at masked degrees of freedom.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine gssmoothmask3d(
     &     U,U_gcw,
     &     alpha,beta,
     &     F,F_gcw,
     &     mask,mask_gcw,
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
      integer U_gcw,F_gcw,mask_gcw

      double precision alpha,beta

      double precision F(ilower0-F_gcw:iupper0+F_gcw,
     &     ilower1-F_gcw:iupper1+F_gcw,
     &     ilower2-F_gcw:iupper2+F_gcw)

      integer mask(ilower0-mask_gcw:iupper0+mask_gcw,
     &     ilower1-mask_gcw:iupper1+mask_gcw,
     &     ilower2-mask_gcw:iupper2+mask_gcw)

      double precision dx(0:3-1)
c
c     Input/Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &     ilower1-U_gcw:iupper1+U_gcw,
     &     ilower2-U_gcw:iupper2+U_gcw)
c
c     Local variables.
c
      integer i0,i1,i2
      double precision    fac0,fac1,fac2,fac
c
c     Perform a single Gauss-Seidel sweep.
c
      fac0 = alpha/(dx(0)*dx(0))
      fac1 = alpha/(dx(1)*dx(1))
      fac2 = alpha/(dx(2)*dx(2))
      fac = 0.5d0/(fac0+fac1+fac2-0.5d0*beta)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               if (mask(i0,i1,i2) .eq. 0) then
                  U(i0,i1,i2) = fac*(
     &                 fac0*(U(i0-1,i1,i2)+U(i0+1,i1,i2)) +
     &                 fac1*(U(i0,i1-1,i2)+U(i0,i1+1,i2)) +
     &                 fac2*(U(i0,i1,i2-1)+U(i0,i1,i2+1)) -
     &                 F(i0,i1,i2))
               endif
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single "red" or "black" Gauss-Seidel sweep for F = alpha
c     div grad U + beta U with masking of certain degrees of freedom.
c
c     NOTE: The solution U is unmodified at masked degrees of freedom.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine rbgssmoothmask3d(
     &     U,U_gcw,
     &     alpha,beta,
     &     F,F_gcw,
     &     mask,mask_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     red_or_black)
c
      implicit none
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer ilower2,iupper2
      integer U_gcw,F_gcw,mask_gcw
      integer red_or_black

      double precision alpha,beta

      double precision F(ilower0-F_gcw:iupper0+F_gcw,
     &     ilower1-F_gcw:iupper1+F_gcw,
     &     ilower2-F_gcw:iupper2+F_gcw)

      integer mask(ilower0-mask_gcw:iupper0+mask_gcw,
     &     ilower1-mask_gcw:iupper1+mask_gcw,
     &     ilower2-mask_gcw:iupper2+mask_gcw)

      double precision dx(0:3-1)
c
c     Input/Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &     ilower1-U_gcw:iupper1+U_gcw,
     &     ilower2-U_gcw:iupper2+U_gcw)
c
c     Local variables.
c
      integer i0,i1,i2
      double precision    fac0,fac1,fac2,fac
c
c     Perform a single "red" or "black" Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1

      fac0 = alpha/(dx(0)*dx(0))
      fac1 = alpha/(dx(1)*dx(1))
      fac2 = alpha/(dx(2)*dx(2))
      fac = 0.5d0/(fac0+fac1+fac2-0.5d0*beta)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               if ( (mod(i0+i1+i2,2) .eq. red_or_black) .and.
     &              (mask(i0,i1,i2) .eq. 0) ) then
                  U(i0,i1,i2) = fac*(
     &                 fac0*(U(i0-1,i1,i2)+U(i0+1,i1,i2)) +
     &                 fac1*(U(i0,i1-1,i2)+U(i0,i1+1,i2)) +
     &                 fac2*(U(i0,i1,i2-1)+U(i0,i1,i2+1)) -
     &                 F(i0,i1,i2))
               endif
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c              Variable coefficient patch smoothers

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single Gauss-Seidel sweep for F = D div grad U +
c     C U. C is a cell-centered variable and
c     D is constant.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine smooth_gs_const_d_var_c_3d(
     &     U,U_gcw,
     &     D,
     &     C,C_gcw,
     &     F,F_gcw,
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
      integer U_gcw,F_gcw,C_gcw

      double precision D
      double precision C(ilower0-C_gcw:iupper0+C_gcw,
     &          ilower1-C_gcw:iupper1+C_gcw,
     &          ilower2-C_gcw:iupper2+C_gcw)

      double precision F(ilower0-F_gcw:iupper0+F_gcw,
     &     ilower1-F_gcw:iupper1+F_gcw,
     &     ilower2-F_gcw:iupper2+F_gcw)

      double precision dx(0:3-1)
c
c     Input/Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &     ilower1-U_gcw:iupper1+U_gcw,
     &     ilower2-U_gcw:iupper2+U_gcw)
c
c     Local variables.
c
      integer i0,i1,i2
      double precision    fac0,fac1,fac2,fac
c
c     Perform a single Gauss-Seidel sweep.
c
      fac0 = D/(dx(0)*dx(0))
      fac1 = D/(dx(1)*dx(1))
      fac2 = D/(dx(2)*dx(2))

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               fac = 0.5d0/(fac0+fac1+fac2-0.5d0*C(i0,i1,i2))
               U(i0,i1,i2) = fac*(
     &              fac0*(U(i0-1,i1,i2)+U(i0+1,i1,i2)) +
     &              fac1*(U(i0,i1-1,i2)+U(i0,i1+1,i2)) +
     &              fac2*(U(i0,i1,i2-1)+U(i0,i1,i2+1)) -
     &              F(i0,i1,i2))
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single "red" or "black" Gauss-Seidel sweep for F = D
c     div grad U + C U. C is cell-centered variable and
c     D is constant.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine smooth_gs_rb_const_d_var_c_3d(
     &     U,U_gcw,
     &     D,
     &     C, C_gcw,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     red_or_black)
c
      implicit none
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer ilower2,iupper2
      integer U_gcw,F_gcw,C_gcw
      integer red_or_black

      double precision D
      double precision C(ilower0-C_gcw:iupper0+C_gcw,
     &          ilower1-C_gcw:iupper1+C_gcw,
     &          ilower2-C_gcw:iupper2+C_gcw)

      double precision F(ilower0-F_gcw:iupper0+F_gcw,
     &     ilower1-F_gcw:iupper1+F_gcw,
     &     ilower2-F_gcw:iupper2+F_gcw)

      double precision dx(0:3-1)
c
c     Input/Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &     ilower1-U_gcw:iupper1+U_gcw,
     &     ilower2-U_gcw:iupper2+U_gcw)
c
c     Local variables.
c
      integer i0,i1,i2
      double precision    fac0,fac1,fac2,fac
c
c     Perform a single "red" or "black" Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1

      fac0 = D/(dx(0)*dx(0))
      fac1 = D/(dx(1)*dx(1))
      fac2 = D/(dx(2)*dx(2))

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               if ( mod(i0+i1+i2,2) .eq. red_or_black ) then
                  fac = 0.5d0/(fac0+fac1+fac2-0.5d0*C(i0,i1,i2))
                  U(i0,i1,i2) = fac*(
     &                 fac0*(U(i0-1,i1,i2)+U(i0+1,i1,i2)) +
     &                 fac1*(U(i0,i1-1,i2)+U(i0,i1+1,i2)) +
     &                 fac2*(U(i0,i1,i2-1)+U(i0,i1,i2+1)) -
     &                 F(i0,i1,i2))
               endif
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single Gauss-Seidel sweep for F = div D grad U +
c     C U.
c
c     The smoother is written for cell-centered U and side-centered
c     D = (D0,D1) with constant C coefficient.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine smooth_gs_var_d_const_c_3d(
     &     U,U_gcw,
     &     D0,D1,D2,D_gcw,
     &     C,
     &     F,F_gcw,
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
      integer U_gcw,F_gcw,D_gcw

      double precision C

      double precision F(ilower0-F_gcw:iupper0+F_gcw,
     &     ilower1-F_gcw:iupper1+F_gcw,
     &     ilower2-F_gcw:iupper2+F_gcw)

      double precision D0(ilower0-D_gcw:iupper0+1+D_gcw,
     &          ilower1-D_gcw:iupper1+D_gcw,
     &          ilower2-D_gcw:iupper2+D_gcw)
      double precision D1(ilower0-D_gcw:iupper0+D_gcw,
     &          ilower1-D_gcw:iupper1+1+D_gcw,
     &          ilower2-D_gcw:iupper2+D_gcw)
      double precision D2(ilower0-D_gcw:iupper0+D_gcw,
     &          ilower1-D_gcw:iupper1+D_gcw,
     &          ilower2-D_gcw:iupper2+1+D_gcw)

      double precision dx(0:3-1)
c
c     Input/Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &     ilower1-U_gcw:iupper1+U_gcw,
     &     ilower2-U_gcw:iupper2+U_gcw)
c
c     Local variables.
c
      integer i0,i1,i2
      double precision    hx,hy,hz
      double precision    facu0,facl0
      double precision    facu1,facl1
      double precision    facu2,facl2
      double precision    fac
c
c     Perform a single Gauss-Seidel sweep.
c
      hx = dx(0)
      hy = dx(1)
      hz = dx(2)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               facu0 = D0(i0+1,i1,i2)/(hx*hx)
               facl0 = D0(i0,i1,i2)/(hx*hx)
               facu1 = D1(i0,i1+1,i2)/(hy*hy)
               facl1 = D1(i0,i1,i2)/(hy*hy)
               facu2 = D2(i0,i1,i2+1)/(hz*hz)
               facl2 = D2(i0,i1,i2)/(hz*hz)
               fac   = 1.d0/(facu0+facl0+facu1+facl1+facu2+facl2-C)
               U(i0,i1,i2) = fac*(
     &             facu0*U(i0+1,i1,i2) +
     &             facl0*U(i0-1,i1,i2) +
     &             facu1*U(i0,i1+1,i2) +
     &             facl1*U(i0,i1-1,i2) +
     &             facu2*U(i0,i1,i2+1) +
     &             facl2*U(i0,i1,i2-1) -
     &             F(i0,i1,i2))
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single "red" or "black" Gauss-Seidel sweep for F =
c     div D grad U + C U.
c
c     The smoother is written for cell-centered U and side-centered
c     D = (D0,D1) with constant C coefficient.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine smooth_gs_rb_var_d_const_c_3d(
     &     U,U_gcw,
     &     D0,D1,D2,D_gcw,
     &     C,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     red_or_black)
c
      implicit none
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer ilower2,iupper2
      integer U_gcw,F_gcw,D_gcw
      integer red_or_black

      double precision C

      double precision F(ilower0-F_gcw:iupper0+F_gcw,
     &     ilower1-F_gcw:iupper1+F_gcw,
     &     ilower2-F_gcw:iupper2+F_gcw)

      double precision D0(ilower0-D_gcw:iupper0+1+D_gcw,
     &          ilower1-D_gcw:iupper1+D_gcw,
     &          ilower2-D_gcw:iupper2+D_gcw)
      double precision D1(ilower0-D_gcw:iupper0+D_gcw,
     &          ilower1-D_gcw:iupper1+1+D_gcw,
     &          ilower2-D_gcw:iupper2+D_gcw)
      double precision D2(ilower0-D_gcw:iupper0+D_gcw,
     &          ilower1-D_gcw:iupper1+D_gcw,
     &          ilower2-D_gcw:iupper2+1+D_gcw)

      double precision dx(0:3-1)
c
c     Input/Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &     ilower1-U_gcw:iupper1+U_gcw,
     &     ilower2-U_gcw:iupper2+U_gcw)
c
c     Local variables.
c
      integer i0,i1,i2
      double precision    hx,hy,hz
      double precision    facu0,facl0
      double precision    facu1,facl1
      double precision    facu2,facl2
      double precision    fac
c
c     Perform a single "red" or "black" Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1

      hx = dx(0)
      hy = dx(1)
      hz = dx(2)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               if ( mod(i0+i1+i2,2) .eq. red_or_black ) then
                  facu0 = D0(i0+1,i1,i2)/(hx*hx)
                  facl0 = D0(i0,i1,i2)/(hx*hx)
                  facu1 = D1(i0,i1+1,i2)/(hy*hy)
                  facl1 = D1(i0,i1,i2)/(hy*hy)
                  facu2 = D2(i0,i1,i2+1)/(hz*hz)
                  facl2 = D2(i0,i1,i2)/(hz*hz)
                  fac  = 1.d0/(facu0+facl0+facu1+facl1+facu2+facl2-C)
                  U(i0,i1,i2) = fac*(
     &                facu0*U(i0+1,i1,i2) +
     &                facl0*U(i0-1,i1,i2) +
     &                facu1*U(i0,i1+1,i2) +
     &                facl1*U(i0,i1-1,i2) +
     &                facu2*U(i0,i1,i2+1) +
     &                facl2*U(i0,i1,i2-1) -
     &                F(i0,i1,i2))
               endif
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single Gauss-Seidel sweep for F = div D grad U +
c     C U.
c
c     The smoother is written for cell-centered U, side-centered
c     D = (D0,D1) and cell-centered C coefficient.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine smooth_gs_var_dc_3d(
     &     U,U_gcw,
     &     D0,D1,D2,D_gcw,
     &     C,C_gcw,
     &     F,F_gcw,
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
      integer U_gcw,F_gcw,D_gcw,C_gcw

      double precision C(ilower0-C_gcw:iupper0+C_gcw,
     &          ilower1-C_gcw:iupper1+C_gcw,
     &          ilower2-C_gcw:iupper2+C_gcw)

      double precision F(ilower0-F_gcw:iupper0+F_gcw,
     &     ilower1-F_gcw:iupper1+F_gcw,
     &     ilower2-F_gcw:iupper2+F_gcw)

      double precision D0(ilower0-D_gcw:iupper0+1+D_gcw,
     &          ilower1-D_gcw:iupper1+D_gcw,
     &          ilower2-D_gcw:iupper2+D_gcw)
      double precision D1(ilower0-D_gcw:iupper0+D_gcw,
     &          ilower1-D_gcw:iupper1+1+D_gcw,
     &          ilower2-D_gcw:iupper2+D_gcw)
      double precision D2(ilower0-D_gcw:iupper0+D_gcw,
     &          ilower1-D_gcw:iupper1+D_gcw,
     &          ilower2-D_gcw:iupper2+1+D_gcw)

      double precision dx(0:3-1)
c
c     Input/Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &     ilower1-U_gcw:iupper1+U_gcw,
     &     ilower2-U_gcw:iupper2+U_gcw)
c
c     Local variables.
c
      integer i0,i1,i2
      double precision    hx,hy,hz
      double precision    facu0,facl0
      double precision    facu1,facl1
      double precision    facu2,facl2
      double precision    fac
c
c     Perform a single Gauss-Seidel sweep.
c
      hx = dx(0)
      hy = dx(1)
      hz = dx(2)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               facu0 = D0(i0+1,i1,i2)/(hx*hx)
               facl0 = D0(i0,i1,i2)/(hx*hx)
               facu1 = D1(i0,i1+1,i2)/(hy*hy)
               facl1 = D1(i0,i1,i2)/(hy*hy)
               facu2 = D2(i0,i1,i2+1)/(hz*hz)
               facl2 = D2(i0,i1,i2)/(hz*hz)
               fac   = 1.d0/(facu0+facl0+facu1+facl1+facu2+facl2
     &                  -C(i0,i1,i2))
               U(i0,i1,i2) = fac*(
     &             facu0*U(i0+1,i1,i2) +
     &             facl0*U(i0-1,i1,i2) +
     &             facu1*U(i0,i1+1,i2) +
     &             facl1*U(i0,i1-1,i2) +
     &             facu2*U(i0,i1,i2+1) +
     &             facl2*U(i0,i1,i2-1) -
     &             F(i0,i1,i2))
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Perform a single "red" or "black" Gauss-Seidel sweep for F =
c     div D grad U + C U.
c
c     The smoother is written for cell-centered U, side-centered
c     D = (D0,D1) and cell-centered C coefficient.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine smooth_gs_rb_var_dc_3d(
     &     U,U_gcw,
     &     D0,D1,D2,D_gcw,
     &     C,C_gcw,
     &     F,F_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     red_or_black)
c
      implicit none
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer ilower2,iupper2
      integer U_gcw,F_gcw,D_gcw,C_gcw
      integer red_or_black

      double precision C(ilower0-C_gcw:iupper0+C_gcw,
     &          ilower1-C_gcw:iupper1+C_gcw,
     &          ilower2-C_gcw:iupper2+C_gcw)

      double precision F(ilower0-F_gcw:iupper0+F_gcw,
     &     ilower1-F_gcw:iupper1+F_gcw,
     &     ilower2-F_gcw:iupper2+F_gcw)

      double precision D0(ilower0-D_gcw:iupper0+1+D_gcw,
     &          ilower1-D_gcw:iupper1+D_gcw,
     &          ilower2-D_gcw:iupper2+D_gcw)
      double precision D1(ilower0-D_gcw:iupper0+D_gcw,
     &          ilower1-D_gcw:iupper1+1+D_gcw,
     &          ilower2-D_gcw:iupper2+D_gcw)
      double precision D2(ilower0-D_gcw:iupper0+D_gcw,
     &          ilower1-D_gcw:iupper1+D_gcw,
     &          ilower2-D_gcw:iupper2+1+D_gcw)

      double precision dx(0:3-1)
c
c     Input/Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &     ilower1-U_gcw:iupper1+U_gcw,
     &     ilower2-U_gcw:iupper2+U_gcw)
c
c     Local variables.
c
      integer i0,i1,i2
      double precision    hx,hy,hz
      double precision    facu0,facl0
      double precision    facu1,facl1
      double precision    facu2,facl2
      double precision    fac
c
c     Perform a single "red" or "black" Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1

      hx = dx(0)
      hy = dx(1)
      hz = dx(2)

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               if ( mod(i0+i1+i2,2) .eq. red_or_black ) then
                  facu0 = D0(i0+1,i1,i2)/(hx*hx)
                  facl0 = D0(i0,i1,i2)/(hx*hx)
                  facu1 = D1(i0,i1+1,i2)/(hy*hy)
                  facl1 = D1(i0,i1,i2)/(hy*hy)
                  facu2 = D2(i0,i1,i2+1)/(hz*hz)
                  facl2 = D2(i0,i1,i2)/(hz*hz)
                  fac  = 1.d0/(facu0+facl0+facu1+facl1+facu2+facl2
     &                   -C(i0,i1,i2))
                  U(i0,i1,i2) = fac*(
     &                facu0*U(i0+1,i1,i2) +
     &                facl0*U(i0-1,i1,i2) +
     &                facu1*U(i0,i1+1,i2) +
     &                facl1*U(i0,i1-1,i2) +
     &                facu2*U(i0,i1,i2+1) +
     &                facl2*U(i0,i1,i2-1) -
     &                F(i0,i1,i2))
               endif
            enddo
         enddo
      enddo
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Perform a single Gauss-Seidel sweep for
c     (f0,f1,f2) = alpha div mu (grad (u0,u1,u2) + grad (u0, u1,u2)^T) + beta c (u0,u1,u2).
c
c  The smoother is written for side-centered vector fields (u0, u1, u2) and (f0, f1, f2)
c  with edge-centered coefficient mu and side-centered coefficient (c0,c1,c2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vcgssmooth3d(
     &     u0,u1,u2,u_gcw,
     &     f0,f1,f2,f_gcw,
     &     c0,c1,c2,c_gcw,
     &     mu0,mu1,mu2,mu_gcw,
     &     alpha,beta,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     var_c,
     &     use_harmonic_interp)
c
      implicit none
c
c     Functions.
c
      double precision a_avg12, h_avg12
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer ilower2,iupper2
      integer u_gcw,f_gcw,c_gcw,mu_gcw
      integer var_c,use_harmonic_interp

      double precision alpha,beta

      double precision mu0(ilower0-mu_gcw:iupper0+mu_gcw,
     &          ilower1-mu_gcw:iupper1+1+mu_gcw,
     &          ilower2-mu_gcw:iupper2+1+mu_gcw)
      double precision mu1(ilower0-mu_gcw:iupper0+1+mu_gcw,
     &          ilower1-mu_gcw:iupper1+mu_gcw,
     &          ilower2-mu_gcw:iupper2+1+mu_gcw)
      double precision mu2(ilower0-mu_gcw:iupper0+1+mu_gcw,
     &          ilower1-mu_gcw:iupper1+1+mu_gcw,
     &          ilower2-mu_gcw:iupper2+mu_gcw)

      double precision f0(ilower0-f_gcw:iupper0+1+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f1(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f2(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+1+f_gcw)

      double precision c0(ilower0-c_gcw:iupper0+1+c_gcw,
     &          ilower1-c_gcw:iupper1+c_gcw,
     &          ilower2-c_gcw:iupper2+c_gcw)
      double precision c1(ilower0-c_gcw:iupper0+c_gcw,
     &          ilower1-c_gcw:iupper1+1+c_gcw,
     &          ilower2-c_gcw:iupper2+c_gcw)
      double precision c2(ilower0-c_gcw:iupper0+c_gcw,
     &          ilower1-c_gcw:iupper1+c_gcw,
     &          ilower2-c_gcw:iupper2+1+c_gcw)

      double precision dx(0:3-1)

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
      integer i0,i1,i2
      double precision fac0,fac1,fac2,fac,nmr,dnr,mu_lower,mu_upper,c
c
c     Perform a single Gauss-Seidel sweep.
c
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))
      fac2 = 1.d0/(dx(2))

      fac = 2.d0*fac0**2.d0
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1

               c = beta
               if (var_c .eq. 1) then
                  c = c0(i0,i1,i2)*beta
               endif

               if (use_harmonic_interp .eq. 1) then
                  mu_upper = h_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),
     &                             mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                             mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                             mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                             mu2(i0+1,i1+1,i2))

                  mu_lower = h_avg12(mu0(i0-1,i1,i2),mu0(i0-1,i1+1,i2),
     &                            mu0(i0-1,i1,i2+1),mu0(i0-1,i1+1,i2+1),
     &                            mu1(i0-1,i1,i2),mu1(i0,i1,i2),
     &                            mu1(i0-1,i1,i2+1),mu1(i0,i1,i2+1),
     &                            mu2(i0-1,i1,i2),mu2(i0,i1,i2),
     &                            mu2(i0-1,i1+1,i2),mu2(i0,i1+1,i2))
               else
                  mu_upper = a_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),
     &                             mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                             mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                             mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                             mu2(i0+1,i1+1,i2))

                  mu_lower = a_avg12(mu0(i0-1,i1,i2),mu0(i0-1,i1+1,i2),
     &                            mu0(i0-1,i1,i2+1),mu0(i0-1,i1+1,i2+1),
     &                            mu1(i0-1,i1,i2),mu1(i0,i1,i2),
     &                            mu1(i0-1,i1,i2+1),mu1(i0,i1,i2+1),
     &                            mu2(i0-1,i1,i2),mu2(i0,i1,i2),
     &                            mu2(i0-1,i1+1,i2),mu2(i0,i1+1,i2))
               endif

               dnr =  alpha*(fac*(mu_upper + mu_lower) +
     &             fac1**2.d0*(mu2(i0,i1+1,i2) + mu2(i0,i1,i2))+
     &             fac2**2.d0*(mu1(i0,i1,i2+1) + mu1(i0,i1,i2))) - c

               nmr = -f0(i0,i1,i2) + alpha*(fac*(
     &           mu_upper*u0(i0+1,i1,i2) + mu_lower*u0(i0-1,i1,i2))+
     &           fac1**2.d0*(mu2(i0,i1+1,i2)*u0(i0,i1+1,i2)+
     &             mu2(i0,i1,i2)*u0(i0,i1-1,i2))+
     &           fac0*fac1*(mu2(i0,i1+1,i2)*(u1(i0,i1+1,i2)-
     &            u1(i0-1,i1+1,i2))-mu2(i0,i1,i2)*(u1(i0,i1,i2)-
     &           u1(i0-1,i1,i2)))+
     &           fac2**2.d0*(mu1(i0,i1,i2+1)*u0(i0,i1,i2+1)+
     &             mu1(i0,i1,i2)*u0(i0,i1,i2-1))+
     &           fac0*fac2*(mu1(i0,i1,i2+1)*(u2(i0,i1,i2+1)-
     &             u2(i0-1,i1,i2+1))-mu1(i0,i1,i2)*(u2(i0,i1,i2)-
     &           u2(i0-1,i1,i2))))

               u0(i0,i1,i2) = nmr/dnr
            enddo
         enddo
      enddo

      fac = 2.d0*fac1**2.d0
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0

               c = beta
               if (var_c .eq. 1) then
                  c = c1(i0,i1,i2)*beta
               endif

               if (use_harmonic_interp .eq. 1) then
                 mu_upper = h_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                              mu0(i0,i1,i2+1),
     &                              mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                              mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                              mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                              mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                              mu2(i0+1,i1+1,i2))

                 mu_lower = h_avg12(mu0(i0,i1-1,i2),mu0(i0,i1,i2),
     &                              mu0(i0,i1-1,i2+1),mu0(i0,i1,i2+1),
     &                              mu1(i0,i1-1,i2),mu1(i0+1,i1-1,i2),
     &                           mu1(i0,i1-1,i2+1),mu1(i0+1,i1-1,i2+1),
     &                              mu2(i0,i1-1,i2),mu2(i0+1,i1-1,i2),
     &                              mu2(i0,i1,i2),mu2(i0+1,i1,i2))
               else
                 mu_upper = a_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                              mu0(i0,i1,i2+1),
     &                              mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                              mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                              mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                              mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                              mu2(i0+1,i1+1,i2))

                 mu_lower = a_avg12(mu0(i0,i1-1,i2),mu0(i0,i1,i2),
     &                              mu0(i0,i1-1,i2+1),mu0(i0,i1,i2+1),
     &                              mu1(i0,i1-1,i2),mu1(i0+1,i1-1,i2),
     &                           mu1(i0,i1-1,i2+1),mu1(i0+1,i1-1,i2+1),
     &                              mu2(i0,i1-1,i2),mu2(i0+1,i1-1,i2),
     &                              mu2(i0,i1,i2),mu2(i0+1,i1,i2))
               endif

               dnr = alpha*(fac*(mu_upper + mu_lower)+
     &            fac0**2.d0*(mu2(i0+1,i1,i2) + mu2(i0,i1,i2))+
     &            fac2**2.d0*(mu0(i0,i1,i2+1) + mu0(i0,i1,i2))) - c

               nmr = -f1(i0,i1,i2) + alpha*(fac*(
     &          mu_upper*u1(i0,i1+1,i2) + mu_lower*u1(i0,i1-1,i2))+
     &          fac0**2.d0*(mu2(i0+1,i1,i2)*u1(i0+1,i1,i2) +
     &             mu2(i0,i1,i2)*u1(i0-1,i1,i2))+
     &          fac0*fac1*(mu2(i0+1,i1,i2)*(u0(i0+1,i1,i2)-
     &             u0(i0+1,i1-1,i2))-
     &             mu2(i0,i1,i2)*(u0(i0,i1,i2)-u0(i0,i1-1,i2)))+
     &          fac2**2.d0*(mu0(i0,i1,i2+1)*u1(i0,i1,i2+1)+
     &             mu0(i0,i1,i2)*u1(i0,i1,i2-1))+
     &           fac1*fac2*(mu0(i0,i1,i2+1)*(u2(i0,i1,i2+1)-
     &             u2(i0,i1-1,i2+1))-mu0(i0,i1,i2)*(u2(i0,i1,i2)-
     &           u2(i0,i1-1,i2))))

               u1(i0,i1,i2) = nmr/dnr

            enddo
         enddo
      enddo

      fac = 2.d0*fac2**2.d0
      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0

               c = beta
               if (var_c .eq. 1) then
                  c = c2(i0,i1,i2)*beta
               endif
               if (use_harmonic_interp .eq. 1) then
                  mu_upper = h_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),
     &                             mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                             mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                             mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                             mu2(i0+1,i1+1,i2))

                  mu_lower = h_avg12(mu0(i0,i1,i2-1),mu0(i0,i1+1,i2-1),
     &                             mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu1(i0,i1,i2-1),
     &                             mu1(i0+1,i1,i2-1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu2(i0,i1,i2-1),
     &                             mu2(i0+1,i1,i2-1),mu2(i0,i1+1,i2-1),
     &                             mu2(i0+1,i1+1,i2-1))
               else
                  mu_upper = a_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),
     &                             mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                             mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                             mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                             mu2(i0+1,i1+1,i2))

                  mu_lower = a_avg12(mu0(i0,i1,i2-1),mu0(i0,i1+1,i2-1),
     &                             mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu1(i0,i1,i2-1),
     &                             mu1(i0+1,i1,i2-1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu2(i0,i1,i2-1),
     &                             mu2(i0+1,i1,i2-1),mu2(i0,i1+1,i2-1),
     &                             mu2(i0+1,i1+1,i2-1))
                endif

               dnr = alpha*(fac*(mu_upper + mu_lower)+
     &          fac0**2.d0*(mu1(i0+1,i1,i2) + mu1(i0,i1,i2))+
     &          fac1**2.d0*(mu0(i0,i1+1,i2) + mu0(i0,i1,i2))) - c

               nmr = -f2(i0,i1,i2) + alpha*(fac*(
     &         mu_upper*u2(i0,i1,i2+1) + mu_lower*u2(i0,i1,i2-1))+
     &         fac1**2.d0*(mu0(i0,i1+1,i2)*u2(i0,i1+1,i2)+
     &           mu0(i0,i1,i2)*u2(i0,i1-1,i2))+
     &         fac1*fac2*(mu0(i0,i1+1,i2)*(u1(i0,i1+1,i2)-
     &          u1(i0,i1+1,i2-1))-mu0(i0,i1,i2)*(u1(i0,i1,i2)-
     &         u1(i0,i1,i2-1)))+
     &         fac0**2.d0*(mu1(i0+1,i1,i2)*u2(i0+1,i1,i2)+
     &           mu1(i0,i1,i2)*u2(i0-1,i1,i2))+
     &         fac0*fac2*(mu1(i0+1,i1,i2)*(u0(i0+1,i1,i2)-
     &           u0(i0+1,i1,i2-1))-mu1(i0,i1,i2)*(u0(i0,i1,i2)-
     &         u0(i0,i1,i2-1))))

               u2(i0,i1,i2) = nmr/dnr

            enddo
         enddo
      enddo

c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Perform a single "red" or "black" Gauss-Seidel sweep for
c     (f0,f1,f2) = alpha div mu (grad (u0,u1,u2) + grad (u0, u1,u2)^T) + beta c (u0,u1,u2).
c
c  The smoother is written for side-centered vector fields (u0, u1, u2) and (f0, f1, f2)
c  with edge-centered coefficient mu and side-centered coefficient (c0,c1,c2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vcrbgssmooth3d(
     &     u0,u1,u2,u_gcw,
     &     f0,f1,f2,f_gcw,
     &     c0,c1,c2,c_gcw,
     &     mu0,mu1,mu2,mu_gcw,
     &     alpha,beta,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     var_c,
     &     use_harmonic_interp,
     &     red_or_black)
c
      implicit none
c
c     Functions.
c
      double precision a_avg12, h_avg12
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer ilower2,iupper2
      integer u_gcw,f_gcw,c_gcw,mu_gcw
      integer var_c,use_harmonic_interp
      integer red_or_black

      double precision alpha,beta

      double precision mu0(ilower0-mu_gcw:iupper0+mu_gcw,
     &          ilower1-mu_gcw:iupper1+1+mu_gcw,
     &          ilower2-mu_gcw:iupper2+1+mu_gcw)
      double precision mu1(ilower0-mu_gcw:iupper0+1+mu_gcw,
     &          ilower1-mu_gcw:iupper1+mu_gcw,
     &          ilower2-mu_gcw:iupper2+1+mu_gcw)
      double precision mu2(ilower0-mu_gcw:iupper0+1+mu_gcw,
     &          ilower1-mu_gcw:iupper1+1+mu_gcw,
     &          ilower2-mu_gcw:iupper2+mu_gcw)

      double precision f0(ilower0-f_gcw:iupper0+1+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f1(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f2(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+1+f_gcw)

      double precision c0(ilower0-c_gcw:iupper0+1+c_gcw,
     &          ilower1-c_gcw:iupper1+c_gcw,
     &          ilower2-c_gcw:iupper2+c_gcw)
      double precision c1(ilower0-c_gcw:iupper0+c_gcw,
     &          ilower1-c_gcw:iupper1+1+c_gcw,
     &          ilower2-c_gcw:iupper2+c_gcw)
      double precision c2(ilower0-c_gcw:iupper0+c_gcw,
     &          ilower1-c_gcw:iupper1+c_gcw,
     &          ilower2-c_gcw:iupper2+1+c_gcw)

      double precision dx(0:3-1)

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
      integer i0,i1,i2
      double precision fac0,fac1,fac2,fac,nmr,dnr,mu_lower,mu_upper,c
c
c     Perform a single "red" or "black"  Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1

      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))
      fac2 = 1.d0/(dx(2))

      fac = 2.d0*fac0**2.d0
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               if ( mod(i0+i1+i2,2) .eq. red_or_black ) then

                  c = beta
                  if (var_c .eq. 1) then
                     c = c0(i0,i1,i2)*beta
                  endif

               if (use_harmonic_interp .eq. 1) then
                  mu_upper = h_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),
     &                             mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                             mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                             mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                             mu2(i0+1,i1+1,i2))

                  mu_lower = h_avg12(mu0(i0-1,i1,i2),mu0(i0-1,i1+1,i2),
     &                            mu0(i0-1,i1,i2+1),mu0(i0-1,i1+1,i2+1),
     &                            mu1(i0-1,i1,i2),mu1(i0,i1,i2),
     &                            mu1(i0-1,i1,i2+1),mu1(i0,i1,i2+1),
     &                            mu2(i0-1,i1,i2),mu2(i0,i1,i2),
     &                            mu2(i0-1,i1+1,i2),mu2(i0,i1+1,i2))
               else
                  mu_upper = a_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),
     &                             mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                             mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                             mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                             mu2(i0+1,i1+1,i2))

                  mu_lower = a_avg12(mu0(i0-1,i1,i2),mu0(i0-1,i1+1,i2),
     &                            mu0(i0-1,i1,i2+1),mu0(i0-1,i1+1,i2+1),
     &                            mu1(i0-1,i1,i2),mu1(i0,i1,i2),
     &                            mu1(i0-1,i1,i2+1),mu1(i0,i1,i2+1),
     &                            mu2(i0-1,i1,i2),mu2(i0,i1,i2),
     &                            mu2(i0-1,i1+1,i2),mu2(i0,i1+1,i2))
               endif

                  dnr =  alpha*(fac*(mu_upper + mu_lower) +
     &                fac1**2.d0*(mu2(i0,i1+1,i2) + mu2(i0,i1,i2))+
     &                fac2**2.d0*(mu1(i0,i1,i2+1) + mu1(i0,i1,i2))) - c

                  nmr = -f0(i0,i1,i2) + alpha*(fac*(
     &              mu_upper*u0(i0+1,i1,i2) + mu_lower*u0(i0-1,i1,i2))+
     &              fac1**2.d0*(mu2(i0,i1+1,i2)*u0(i0,i1+1,i2)+
     &                mu2(i0,i1,i2)*u0(i0,i1-1,i2))+
     &              fac0*fac1*(mu2(i0,i1+1,i2)*(u1(i0,i1+1,i2)-
     &               u1(i0-1,i1+1,i2))-mu2(i0,i1,i2)*(u1(i0,i1,i2)-
     &              u1(i0-1,i1,i2)))+
     &              fac2**2.d0*(mu1(i0,i1,i2+1)*u0(i0,i1,i2+1)+
     &                mu1(i0,i1,i2)*u0(i0,i1,i2-1))+
     &              fac0*fac2*(mu1(i0,i1,i2+1)*(u2(i0,i1,i2+1)-
     &                u2(i0-1,i1,i2+1))-mu1(i0,i1,i2)*(u2(i0,i1,i2)-
     &              u2(i0-1,i1,i2))))

                  u0(i0,i1,i2) = nmr/dnr
               endif
            enddo
         enddo
      enddo

      fac = 2.d0*fac1**2.d0
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0
               if ( mod(i0+i1+i2,2) .eq. red_or_black ) then

                  c = beta
                  if (var_c .eq. 1) then
                     c = c1(i0,i1,i2)*beta
                  endif

                  if (use_harmonic_interp .eq. 1) then
                    mu_upper = h_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),
     &                             mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                             mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                             mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                             mu2(i0+1,i1+1,i2))

                    mu_lower = h_avg12(mu0(i0,i1-1,i2),mu0(i0,i1,i2),
     &                             mu0(i0,i1-1,i2+1),mu0(i0,i1,i2+1),
     &                             mu1(i0,i1-1,i2),mu1(i0+1,i1-1,i2),
     &                          mu1(i0,i1-1,i2+1),mu1(i0+1,i1-1,i2+1),
     &                             mu2(i0,i1-1,i2),mu2(i0+1,i1-1,i2),
     &                             mu2(i0,i1,i2),mu2(i0+1,i1,i2))
                  else
                    mu_upper = a_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),
     &                             mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                             mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                             mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                             mu2(i0+1,i1+1,i2))

                    mu_lower = a_avg12(mu0(i0,i1-1,i2),mu0(i0,i1,i2),
     &                             mu0(i0,i1-1,i2+1),mu0(i0,i1,i2+1),
     &                             mu1(i0,i1-1,i2),mu1(i0+1,i1-1,i2),
     &                          mu1(i0,i1-1,i2+1),mu1(i0+1,i1-1,i2+1),
     &                             mu2(i0,i1-1,i2),mu2(i0+1,i1-1,i2),
     &                             mu2(i0,i1,i2),mu2(i0+1,i1,i2))
                  endif


                  dnr = alpha*(fac*(mu_upper + mu_lower)+
     &               fac0**2.d0*(mu2(i0+1,i1,i2) + mu2(i0,i1,i2))+
     &               fac2**2.d0*(mu0(i0,i1,i2+1) + mu0(i0,i1,i2))) - c

                  nmr = -f1(i0,i1,i2) + alpha*(fac*(
     &             mu_upper*u1(i0,i1+1,i2) + mu_lower*u1(i0,i1-1,i2))+
     &             fac0**2.d0*(mu2(i0+1,i1,i2)*u1(i0+1,i1,i2) +
     &                mu2(i0,i1,i2)*u1(i0-1,i1,i2))+
     &             fac0*fac1*(mu2(i0+1,i1,i2)*(u0(i0+1,i1,i2)-
     &                u0(i0+1,i1-1,i2))-
     &                mu2(i0,i1,i2)*(u0(i0,i1,i2)-u0(i0,i1-1,i2)))+
     &             fac2**2.d0*(mu0(i0,i1,i2+1)*u1(i0,i1,i2+1)+
     &                mu0(i0,i1,i2)*u1(i0,i1,i2-1))+
     &              fac1*fac2*(mu0(i0,i1,i2+1)*(u2(i0,i1,i2+1)-
     &                u2(i0,i1-1,i2+1))-mu0(i0,i1,i2)*(u2(i0,i1,i2)-
     &              u2(i0,i1-1,i2))))

                  u1(i0,i1,i2) = nmr/dnr
               endif
            enddo
         enddo
      enddo

      fac = 2.d0*fac2**2.d0
      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               if ( mod(i0+i1+i2,2) .eq. red_or_black ) then

                  c = beta
                  if (var_c .eq. 1) then
                     c = c2(i0,i1,i2)*beta
                  endif

                  if (use_harmonic_interp .eq. 1) then
                      mu_upper = h_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                                 mu0(i0,i1,i2+1),
     &                                 mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                                 mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                                 mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                                 mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                                 mu2(i0+1,i1+1,i2))

                      mu_lower = h_avg12(mu0(i0,i1,i2-1),
     &                                 mu0(i0,i1+1,i2-1),
     &                                 mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                                 mu1(i0,i1,i2-1),
     &                                 mu1(i0+1,i1,i2-1),mu1(i0,i1,i2),
     &                                 mu1(i0+1,i1,i2),mu2(i0,i1,i2-1),
     &                              mu2(i0+1,i1,i2-1),mu2(i0,i1+1,i2-1),
     &                                 mu2(i0+1,i1+1,i2-1))
               else
                      mu_upper = a_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                                 mu0(i0,i1,i2+1),
     &                                 mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                                 mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                                 mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                                 mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                                 mu2(i0+1,i1+1,i2))

                      mu_lower = a_avg12(mu0(i0,i1,i2-1),
     &                                 mu0(i0,i1+1,i2-1),
     &                                 mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                                 mu1(i0,i1,i2-1),
     &                                 mu1(i0+1,i1,i2-1),mu1(i0,i1,i2),
     &                                 mu1(i0+1,i1,i2),mu2(i0,i1,i2-1),
     &                              mu2(i0+1,i1,i2-1),mu2(i0,i1+1,i2-1),
     &                                 mu2(i0+1,i1+1,i2-1))
                endif

                  dnr = alpha*(fac*(mu_upper + mu_lower)+
     &             fac0**2.d0*(mu1(i0+1,i1,i2) + mu1(i0,i1,i2))+
     &             fac1**2.d0*(mu0(i0,i1+1,i2) + mu0(i0,i1,i2))) - c

                  nmr = -f2(i0,i1,i2) + alpha*(fac*(
     &            mu_upper*u2(i0,i1,i2+1) + mu_lower*u2(i0,i1,i2-1))+
     &            fac1**2.d0*(mu0(i0,i1+1,i2)*u2(i0,i1+1,i2)+
     &              mu0(i0,i1,i2)*u2(i0,i1-1,i2))+
     &            fac1*fac2*(mu0(i0,i1+1,i2)*(u1(i0,i1+1,i2)-
     &             u1(i0,i1+1,i2-1))-mu0(i0,i1,i2)*(u1(i0,i1,i2)-
     &            u1(i0,i1,i2-1)))+
     &            fac0**2.d0*(mu1(i0+1,i1,i2)*u2(i0+1,i1,i2)+
     &              mu1(i0,i1,i2)*u2(i0-1,i1,i2))+
     &            fac0*fac2*(mu1(i0+1,i1,i2)*(u0(i0+1,i1,i2)-
     &              u0(i0+1,i1,i2-1))-mu1(i0,i1,i2)*(u0(i0,i1,i2)-
     &            u0(i0,i1,i2-1))))

                  u2(i0,i1,i2) = nmr/dnr
               endif
            enddo
         enddo
      enddo

c
      return
      end
c

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Perform a single Gauss-Seidel sweep for
c     (f0,f1,f2) = alpha div mu (grad (u0,u1,u2) + grad (u0, u1,u2)^T) + beta c (u0,u1,u2),
c  with masking of certain degrees of freedom.
c
c     NOTE: The solution (u0,u1,u2) is unmodified at masked degrees of freedom.
c
c  The smoother is written for side-centered vector fields (u0, u1, u2) and (f0, f1, f2)
c  with edge-centered coefficient mu and side-centered coefficient (c0,c1,c2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vcgssmoothmask3d(
     &     u0,u1,u2,u_gcw,
     &     f0,f1,f2,f_gcw,
     &     mask0,mask1,mask2,mask_gcw,
     &     c0,c1,c2,c_gcw,
     &     mu0,mu1,mu2,mu_gcw,
     &     alpha,beta,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     var_c,
     &     use_harmonic_interp)
c
      implicit none
c
c     Functions.
c
      double precision a_avg12, h_avg12
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer ilower2,iupper2
      integer u_gcw,f_gcw,c_gcw,mu_gcw,mask_gcw
      integer var_c,use_harmonic_interp

      double precision alpha,beta

      double precision mu0(ilower0-mu_gcw:iupper0+mu_gcw,
     &          ilower1-mu_gcw:iupper1+1+mu_gcw,
     &          ilower2-mu_gcw:iupper2+1+mu_gcw)
      double precision mu1(ilower0-mu_gcw:iupper0+1+mu_gcw,
     &          ilower1-mu_gcw:iupper1+mu_gcw,
     &          ilower2-mu_gcw:iupper2+1+mu_gcw)
      double precision mu2(ilower0-mu_gcw:iupper0+1+mu_gcw,
     &          ilower1-mu_gcw:iupper1+1+mu_gcw,
     &          ilower2-mu_gcw:iupper2+mu_gcw)

      double precision f0(ilower0-f_gcw:iupper0+1+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f1(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f2(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+1+f_gcw)

      integer mask0(ilower0-mask_gcw:iupper0+1+mask_gcw,
     &          ilower1-mask_gcw:iupper1+mask_gcw,
     &          ilower2-mask_gcw:iupper2+mask_gcw)
      integer mask1(ilower0-mask_gcw:iupper0+mask_gcw,
     &          ilower1-mask_gcw:iupper1+1+mask_gcw,
     &          ilower2-mask_gcw:iupper2+mask_gcw)
      integer mask2(ilower0-mask_gcw:iupper0+mask_gcw,
     &          ilower1-mask_gcw:iupper1+mask_gcw,
     &          ilower2-mask_gcw:iupper2+1+mask_gcw)

      double precision c0(ilower0-c_gcw:iupper0+1+c_gcw,
     &          ilower1-c_gcw:iupper1+c_gcw,
     &          ilower2-c_gcw:iupper2+c_gcw)
      double precision c1(ilower0-c_gcw:iupper0+c_gcw,
     &          ilower1-c_gcw:iupper1+1+c_gcw,
     &          ilower2-c_gcw:iupper2+c_gcw)
      double precision c2(ilower0-c_gcw:iupper0+c_gcw,
     &          ilower1-c_gcw:iupper1+c_gcw,
     &          ilower2-c_gcw:iupper2+1+c_gcw)

      double precision dx(0:3-1)

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
      integer i0,i1,i2
      double precision fac0,fac1,fac2,fac,nmr,dnr,mu_lower,mu_upper,c
c
c     Perform a single Gauss-Seidel sweep.
c
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))
      fac2 = 1.d0/(dx(2))

      fac = 2.d0*fac0**2.d0
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               if (mask0(i0,i1,i2) .eq. 0) then

                  c = beta
                  if (var_c .eq. 1) then
                     c = c0(i0,i1,i2)*beta
                  endif

               if (use_harmonic_interp .eq. 1) then
                  mu_upper = h_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),
     &                             mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                             mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                             mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                             mu2(i0+1,i1+1,i2))

                  mu_lower = h_avg12(mu0(i0-1,i1,i2),mu0(i0-1,i1+1,i2),
     &                            mu0(i0-1,i1,i2+1),mu0(i0-1,i1+1,i2+1),
     &                            mu1(i0-1,i1,i2),mu1(i0,i1,i2),
     &                            mu1(i0-1,i1,i2+1),mu1(i0,i1,i2+1),
     &                            mu2(i0-1,i1,i2),mu2(i0,i1,i2),
     &                            mu2(i0-1,i1+1,i2),mu2(i0,i1+1,i2))
               else
                  mu_upper = a_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),
     &                             mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                             mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                             mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                             mu2(i0+1,i1+1,i2))

                  mu_lower = a_avg12(mu0(i0-1,i1,i2),mu0(i0-1,i1+1,i2),
     &                            mu0(i0-1,i1,i2+1),mu0(i0-1,i1+1,i2+1),
     &                            mu1(i0-1,i1,i2),mu1(i0,i1,i2),
     &                            mu1(i0-1,i1,i2+1),mu1(i0,i1,i2+1),
     &                            mu2(i0-1,i1,i2),mu2(i0,i1,i2),
     &                            mu2(i0-1,i1+1,i2),mu2(i0,i1+1,i2))
               endif

                  dnr =  alpha*(fac*(mu_upper + mu_lower) +
     &                fac1**2.d0*(mu2(i0,i1+1,i2) + mu2(i0,i1,i2))+
     &                fac2**2.d0*(mu1(i0,i1,i2+1) + mu1(i0,i1,i2))) - c

                  nmr = -f0(i0,i1,i2) + alpha*(fac*(
     &              mu_upper*u0(i0+1,i1,i2) + mu_lower*u0(i0-1,i1,i2))+
     &              fac1**2.d0*(mu2(i0,i1+1,i2)*u0(i0,i1+1,i2)+
     &                mu2(i0,i1,i2)*u0(i0,i1-1,i2))+
     &              fac0*fac1*(mu2(i0,i1+1,i2)*(u1(i0,i1+1,i2)-
     &               u1(i0-1,i1+1,i2))-mu2(i0,i1,i2)*(u1(i0,i1,i2)-
     &              u1(i0-1,i1,i2)))+
     &              fac2**2.d0*(mu1(i0,i1,i2+1)*u0(i0,i1,i2+1)+
     &                mu1(i0,i1,i2)*u0(i0,i1,i2-1))+
     &              fac0*fac2*(mu1(i0,i1,i2+1)*(u2(i0,i1,i2+1)-
     &                u2(i0-1,i1,i2+1))-mu1(i0,i1,i2)*(u2(i0,i1,i2)-
     &              u2(i0-1,i1,i2))))

                  u0(i0,i1,i2) = nmr/dnr
               endif
            enddo
         enddo
      enddo

      fac = 2.d0*fac1**2.d0
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0
               if (mask1(i0,i1,i2) .eq. 0) then

                  c = beta
                  if (var_c .eq. 1) then
                     c = c1(i0,i1,i2)*beta
                  endif

                  if (use_harmonic_interp .eq. 1) then
                    mu_upper = h_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),
     &                             mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                             mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                             mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                             mu2(i0+1,i1+1,i2))

                    mu_lower = h_avg12(mu0(i0,i1-1,i2),mu0(i0,i1,i2),
     &                             mu0(i0,i1-1,i2+1),mu0(i0,i1,i2+1),
     &                             mu1(i0,i1-1,i2),mu1(i0+1,i1-1,i2),
     &                          mu1(i0,i1-1,i2+1),mu1(i0+1,i1-1,i2+1),
     &                             mu2(i0,i1-1,i2),mu2(i0+1,i1-1,i2),
     &                             mu2(i0,i1,i2),mu2(i0+1,i1,i2))
                  else
                    mu_upper = a_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),
     &                             mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                             mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                             mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                             mu2(i0+1,i1+1,i2))

                    mu_lower = a_avg12(mu0(i0,i1-1,i2),mu0(i0,i1,i2),
     &                             mu0(i0,i1-1,i2+1),mu0(i0,i1,i2+1),
     &                             mu1(i0,i1-1,i2),mu1(i0+1,i1-1,i2),
     &                          mu1(i0,i1-1,i2+1),mu1(i0+1,i1-1,i2+1),
     &                             mu2(i0,i1-1,i2),mu2(i0+1,i1-1,i2),
     &                             mu2(i0,i1,i2),mu2(i0+1,i1,i2))
                  endif

                  dnr = alpha*(fac*(mu_upper + mu_lower)+
     &               fac0**2.d0*(mu2(i0+1,i1,i2) + mu2(i0,i1,i2))+
     &               fac2**2.d0*(mu0(i0,i1,i2+1) + mu0(i0,i1,i2))) - c

                  nmr = -f1(i0,i1,i2) + alpha*(fac*(
     &             mu_upper*u1(i0,i1+1,i2) + mu_lower*u1(i0,i1-1,i2))+
     &             fac0**2.d0*(mu2(i0+1,i1,i2)*u1(i0+1,i1,i2) +
     &                mu2(i0,i1,i2)*u1(i0-1,i1,i2))+
     &             fac0*fac1*(mu2(i0+1,i1,i2)*(u0(i0+1,i1,i2)-
     &                u0(i0+1,i1-1,i2))-
     &                mu2(i0,i1,i2)*(u0(i0,i1,i2)-u0(i0,i1-1,i2)))+
     &             fac2**2.d0*(mu0(i0,i1,i2+1)*u1(i0,i1,i2+1)+
     &                mu0(i0,i1,i2)*u1(i0,i1,i2-1))+
     &              fac1*fac2*(mu0(i0,i1,i2+1)*(u2(i0,i1,i2+1)-
     &                u2(i0,i1-1,i2+1))-mu0(i0,i1,i2)*(u2(i0,i1,i2)-
     &              u2(i0,i1-1,i2))))

                  u1(i0,i1,i2) = nmr/dnr
               endif
            enddo
         enddo
      enddo

      fac = 2.d0*fac2**2.d0
      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               if (mask2(i0,i1,i2) .eq. 0) then

                  c = beta
                  if (var_c .eq. 1) then
                     c = c2(i0,i1,i2)*beta
                  endif

                  if (use_harmonic_interp .eq. 1) then
                      mu_upper = h_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                                 mu0(i0,i1,i2+1),
     &                                 mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                                 mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                                 mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                                 mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                                 mu2(i0+1,i1+1,i2))

                      mu_lower = h_avg12(mu0(i0,i1,i2-1),
     &                                 mu0(i0,i1+1,i2-1),
     &                                 mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                                 mu1(i0,i1,i2-1),
     &                                 mu1(i0+1,i1,i2-1),mu1(i0,i1,i2),
     &                                 mu1(i0+1,i1,i2),mu2(i0,i1,i2-1),
     &                              mu2(i0+1,i1,i2-1),mu2(i0,i1+1,i2-1),
     &                                 mu2(i0+1,i1+1,i2-1))
                  else
                      mu_upper = a_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                                 mu0(i0,i1,i2+1),
     &                                 mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                                 mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                                 mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                                 mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                                 mu2(i0+1,i1+1,i2))

                      mu_lower = a_avg12(mu0(i0,i1,i2-1),
     &                                 mu0(i0,i1+1,i2-1),
     &                                 mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                                 mu1(i0,i1,i2-1),
     &                                 mu1(i0+1,i1,i2-1),mu1(i0,i1,i2),
     &                                 mu1(i0+1,i1,i2),mu2(i0,i1,i2-1),
     &                              mu2(i0+1,i1,i2-1),mu2(i0,i1+1,i2-1),
     &                                 mu2(i0+1,i1+1,i2-1))
                  endif

                  dnr = alpha*(fac*(mu_upper + mu_lower)+
     &             fac0**2.d0*(mu1(i0+1,i1,i2) + mu1(i0,i1,i2))+
     &             fac1**2.d0*(mu0(i0,i1+1,i2) + mu0(i0,i1,i2))) - c

                  nmr = -f2(i0,i1,i2) + alpha*(fac*(
     &            mu_upper*u2(i0,i1,i2+1) + mu_lower*u2(i0,i1,i2-1))+
     &            fac1**2.d0*(mu0(i0,i1+1,i2)*u2(i0,i1+1,i2)+
     &              mu0(i0,i1,i2)*u2(i0,i1-1,i2))+
     &            fac1*fac2*(mu0(i0,i1+1,i2)*(u1(i0,i1+1,i2)-
     &             u1(i0,i1+1,i2-1))-mu0(i0,i1,i2)*(u1(i0,i1,i2)-
     &            u1(i0,i1,i2-1)))+
     &            fac0**2.d0*(mu1(i0+1,i1,i2)*u2(i0+1,i1,i2)+
     &              mu1(i0,i1,i2)*u2(i0-1,i1,i2))+
     &            fac0*fac2*(mu1(i0+1,i1,i2)*(u0(i0+1,i1,i2)-
     &              u0(i0+1,i1,i2-1))-mu1(i0,i1,i2)*(u0(i0,i1,i2)-
     &            u0(i0,i1,i2-1))))

                  u2(i0,i1,i2) = nmr/dnr
               endif
            enddo
         enddo
      enddo

c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Perform a single "red" or "black" Gauss-Seidel sweep for
c     (f0,f1,f2) = alpha div mu (grad (u0,u1,u2) + grad (u0,u1,u2)^T) + beta c (u0,u1,u2),
c  with masking of certain degrees of freedom.
c
c     NOTE: The solution (u0,u1,u2) is unmodified at masked degrees of freedom.
c
c  The smoother is written for side-centered vector fields (u0, u1, u2) and (f0, f1, f2)
c  with edge-centered coefficient mu and side-centered coefficient (c0,c1,c2)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vcrbgssmoothmask3d(
     &     u0,u1,u2,u_gcw,
     &     f0,f1,f2,f_gcw,
     &     mask0,mask1,mask2,mask_gcw,
     &     c0,c1,c2,c_gcw,
     &     mu0,mu1,mu2,mu_gcw,
     &     alpha,beta,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     var_c,
     &     use_harmonic_interp,
     &     red_or_black)
c
      implicit none
c
c     Functions.
c
      double precision a_avg12, h_avg12
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer ilower2,iupper2
      integer u_gcw,f_gcw,c_gcw,mu_gcw,mask_gcw
      integer var_c,use_harmonic_interp
      integer red_or_black

      double precision alpha,beta

      double precision mu0(ilower0-mu_gcw:iupper0+mu_gcw,
     &          ilower1-mu_gcw:iupper1+1+mu_gcw,
     &          ilower2-mu_gcw:iupper2+1+mu_gcw)
      double precision mu1(ilower0-mu_gcw:iupper0+1+mu_gcw,
     &          ilower1-mu_gcw:iupper1+mu_gcw,
     &          ilower2-mu_gcw:iupper2+1+mu_gcw)
      double precision mu2(ilower0-mu_gcw:iupper0+1+mu_gcw,
     &          ilower1-mu_gcw:iupper1+1+mu_gcw,
     &          ilower2-mu_gcw:iupper2+mu_gcw)

      double precision f0(ilower0-f_gcw:iupper0+1+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f1(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f2(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+1+f_gcw)

      integer mask0(ilower0-mask_gcw:iupper0+1+mask_gcw,
     &          ilower1-mask_gcw:iupper1+mask_gcw,
     &          ilower2-mask_gcw:iupper2+mask_gcw)
      integer mask1(ilower0-mask_gcw:iupper0+mask_gcw,
     &          ilower1-mask_gcw:iupper1+1+mask_gcw,
     &          ilower2-mask_gcw:iupper2+mask_gcw)
      integer mask2(ilower0-mask_gcw:iupper0+mask_gcw,
     &          ilower1-mask_gcw:iupper1+mask_gcw,
     &          ilower2-mask_gcw:iupper2+1+mask_gcw)

      double precision c0(ilower0-c_gcw:iupper0+1+c_gcw,
     &          ilower1-c_gcw:iupper1+c_gcw,
     &          ilower2-c_gcw:iupper2+c_gcw)
      double precision c1(ilower0-c_gcw:iupper0+c_gcw,
     &          ilower1-c_gcw:iupper1+1+c_gcw,
     &          ilower2-c_gcw:iupper2+c_gcw)
      double precision c2(ilower0-c_gcw:iupper0+c_gcw,
     &          ilower1-c_gcw:iupper1+c_gcw,
     &          ilower2-c_gcw:iupper2+1+c_gcw)

      double precision dx(0:3-1)

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
      integer i0,i1,i2
      double precision fac0,fac1,fac2,fac,nmr,dnr,mu_lower,mu_upper,c
c
c     Perform a single"red" or "black"  Gauss-Seidel sweep.
c
      red_or_black = mod(red_or_black,2) ! "red" = 0, "black" = 1

      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))
      fac2 = 1.d0/(dx(2))

      fac = 2.d0*fac0**2.d0
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1
               if ( (mod(i0+i1+i2,2) .eq. red_or_black) .and.
     &              (mask0(i0,i1,i2) .eq. 0) ) then

                  c = beta
                  if (var_c .eq. 1) then
                     c = c0(i0,i1,i2)*beta
                  endif

               if (use_harmonic_interp .eq. 1) then
                  mu_upper = h_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),
     &                             mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                             mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                             mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                             mu2(i0+1,i1+1,i2))

                  mu_lower = h_avg12(mu0(i0-1,i1,i2),mu0(i0-1,i1+1,i2),
     &                            mu0(i0-1,i1,i2+1),mu0(i0-1,i1+1,i2+1),
     &                            mu1(i0-1,i1,i2),mu1(i0,i1,i2),
     &                            mu1(i0-1,i1,i2+1),mu1(i0,i1,i2+1),
     &                            mu2(i0-1,i1,i2),mu2(i0,i1,i2),
     &                            mu2(i0-1,i1+1,i2),mu2(i0,i1+1,i2))
               else
                  mu_upper = a_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),
     &                             mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                             mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                             mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                             mu2(i0+1,i1+1,i2))

                  mu_lower = a_avg12(mu0(i0-1,i1,i2),mu0(i0-1,i1+1,i2),
     &                            mu0(i0-1,i1,i2+1),mu0(i0-1,i1+1,i2+1),
     &                            mu1(i0-1,i1,i2),mu1(i0,i1,i2),
     &                            mu1(i0-1,i1,i2+1),mu1(i0,i1,i2+1),
     &                            mu2(i0-1,i1,i2),mu2(i0,i1,i2),
     &                            mu2(i0-1,i1+1,i2),mu2(i0,i1+1,i2))
               endif

                  dnr =  alpha*(fac*(mu_upper + mu_lower) +
     &                fac1**2.d0*(mu2(i0,i1+1,i2) + mu2(i0,i1,i2))+
     &                fac2**2.d0*(mu1(i0,i1,i2+1) + mu1(i0,i1,i2))) - c

                  nmr = -f0(i0,i1,i2) + alpha*(fac*(
     &              mu_upper*u0(i0+1,i1,i2) + mu_lower*u0(i0-1,i1,i2))+
     &              fac1**2.d0*(mu2(i0,i1+1,i2)*u0(i0,i1+1,i2)+
     &                mu2(i0,i1,i2)*u0(i0,i1-1,i2))+
     &              fac0*fac1*(mu2(i0,i1+1,i2)*(u1(i0,i1+1,i2)-
     &               u1(i0-1,i1+1,i2))-mu2(i0,i1,i2)*(u1(i0,i1,i2)-
     &              u1(i0-1,i1,i2)))+
     &              fac2**2.d0*(mu1(i0,i1,i2+1)*u0(i0,i1,i2+1)+
     &                mu1(i0,i1,i2)*u0(i0,i1,i2-1))+
     &              fac0*fac2*(mu1(i0,i1,i2+1)*(u2(i0,i1,i2+1)-
     &                u2(i0-1,i1,i2+1))-mu1(i0,i1,i2)*(u2(i0,i1,i2)-
     &              u2(i0-1,i1,i2))))

                  u0(i0,i1,i2) = nmr/dnr
               endif
            enddo
         enddo
      enddo

      fac = 2.d0*fac1**2.d0
      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0
               if ( (mod(i0+i1+i2,2) .eq. red_or_black) .and.
     &              (mask1(i0,i1,i2) .eq. 0) ) then

                  c = beta
                  if (var_c .eq. 1) then
                     c = c1(i0,i1,i2)*beta
                  endif

                  if (use_harmonic_interp .eq. 1) then
                    mu_upper = h_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),
     &                             mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                             mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                             mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                             mu2(i0+1,i1+1,i2))

                    mu_lower = h_avg12(mu0(i0,i1-1,i2),mu0(i0,i1,i2),
     &                             mu0(i0,i1-1,i2+1),mu0(i0,i1,i2+1),
     &                             mu1(i0,i1-1,i2),mu1(i0+1,i1-1,i2),
     &                          mu1(i0,i1-1,i2+1),mu1(i0+1,i1-1,i2+1),
     &                             mu2(i0,i1-1,i2),mu2(i0+1,i1-1,i2),
     &                             mu2(i0,i1,i2),mu2(i0+1,i1,i2))
                  else
                    mu_upper = a_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),
     &                             mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                             mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                             mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                             mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                             mu2(i0+1,i1+1,i2))

                    mu_lower = a_avg12(mu0(i0,i1-1,i2),mu0(i0,i1,i2),
     &                             mu0(i0,i1-1,i2+1),mu0(i0,i1,i2+1),
     &                             mu1(i0,i1-1,i2),mu1(i0+1,i1-1,i2),
     &                          mu1(i0,i1-1,i2+1),mu1(i0+1,i1-1,i2+1),
     &                             mu2(i0,i1-1,i2),mu2(i0+1,i1-1,i2),
     &                             mu2(i0,i1,i2),mu2(i0+1,i1,i2))
                  endif

                  dnr = alpha*(fac*(mu_upper + mu_lower)+
     &               fac0**2.d0*(mu2(i0+1,i1,i2) + mu2(i0,i1,i2))+
     &               fac2**2.d0*(mu0(i0,i1,i2+1) + mu0(i0,i1,i2))) - c

                  nmr = -f1(i0,i1,i2) + alpha*(fac*(
     &             mu_upper*u1(i0,i1+1,i2) + mu_lower*u1(i0,i1-1,i2))+
     &             fac0**2.d0*(mu2(i0+1,i1,i2)*u1(i0+1,i1,i2) +
     &                mu2(i0,i1,i2)*u1(i0-1,i1,i2))+
     &             fac0*fac1*(mu2(i0+1,i1,i2)*(u0(i0+1,i1,i2)-
     &                u0(i0+1,i1-1,i2))-
     &                mu2(i0,i1,i2)*(u0(i0,i1,i2)-u0(i0,i1-1,i2)))+
     &             fac2**2.d0*(mu0(i0,i1,i2+1)*u1(i0,i1,i2+1)+
     &                mu0(i0,i1,i2)*u1(i0,i1,i2-1))+
     &              fac1*fac2*(mu0(i0,i1,i2+1)*(u2(i0,i1,i2+1)-
     &                u2(i0,i1-1,i2+1))-mu0(i0,i1,i2)*(u2(i0,i1,i2)-
     &              u2(i0,i1-1,i2))))

                  u1(i0,i1,i2) = nmr/dnr
               endif
            enddo
         enddo
      enddo

      fac = 2.d0*fac2**2.d0
      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               if ( (mod(i0+i1+i2,2) .eq. red_or_black) .and.
     &              (mask2(i0,i1,i2) .eq. 0) ) then

                  c = beta
                  if (var_c .eq. 1) then
                     c = c2(i0,i1,i2)*beta
                  endif

                  if (use_harmonic_interp .eq. 1) then
                      mu_upper = h_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                                 mu0(i0,i1,i2+1),
     &                                 mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                                 mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                                 mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                                 mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                                 mu2(i0+1,i1+1,i2))

                      mu_lower = h_avg12(mu0(i0,i1,i2-1),
     &                                 mu0(i0,i1+1,i2-1),
     &                                 mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                                 mu1(i0,i1,i2-1),
     &                                 mu1(i0+1,i1,i2-1),mu1(i0,i1,i2),
     &                                 mu1(i0+1,i1,i2),mu2(i0,i1,i2-1),
     &                              mu2(i0+1,i1,i2-1),mu2(i0,i1+1,i2-1),
     &                                 mu2(i0+1,i1+1,i2-1))
               else
                      mu_upper = a_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                                 mu0(i0,i1,i2+1),
     &                                 mu0(i0,i1+1,i2+1),mu1(i0,i1,i2),
     &                                 mu1(i0+1,i1,i2),mu1(i0,i1,i2+1),
     &                                 mu1(i0+1,i1,i2+1),mu2(i0,i1,i2),
     &                                 mu2(i0+1,i1,i2),mu2(i0,i1+1,i2),
     &                                 mu2(i0+1,i1+1,i2))

                      mu_lower = a_avg12(mu0(i0,i1,i2-1),
     &                                 mu0(i0,i1+1,i2-1),
     &                                 mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                                 mu1(i0,i1,i2-1),
     &                                 mu1(i0+1,i1,i2-1),mu1(i0,i1,i2),
     &                                 mu1(i0+1,i1,i2),mu2(i0,i1,i2-1),
     &                              mu2(i0+1,i1,i2-1),mu2(i0,i1+1,i2-1),
     &                                 mu2(i0+1,i1+1,i2-1))
                endif

                  dnr = alpha*(fac*(mu_upper + mu_lower)+
     &             fac0**2.d0*(mu1(i0+1,i1,i2) + mu1(i0,i1,i2))+
     &             fac1**2.d0*(mu0(i0,i1+1,i2) + mu0(i0,i1,i2))) - c

                  nmr = -f2(i0,i1,i2) + alpha*(fac*(
     &            mu_upper*u2(i0,i1,i2+1) + mu_lower*u2(i0,i1,i2-1))+
     &            fac1**2.d0*(mu0(i0,i1+1,i2)*u2(i0,i1+1,i2)+
     &              mu0(i0,i1,i2)*u2(i0,i1-1,i2))+
     &            fac1*fac2*(mu0(i0,i1+1,i2)*(u1(i0,i1+1,i2)-
     &             u1(i0,i1+1,i2-1))-mu0(i0,i1,i2)*(u1(i0,i1,i2)-
     &            u1(i0,i1,i2-1)))+
     &            fac0**2.d0*(mu1(i0+1,i1,i2)*u2(i0+1,i1,i2)+
     &              mu1(i0,i1,i2)*u2(i0-1,i1,i2))+
     &            fac0*fac2*(mu1(i0+1,i1,i2)*(u0(i0+1,i1,i2)-
     &              u0(i0+1,i1,i2-1))-mu1(i0,i1,i2)*(u0(i0,i1,i2)-
     &            u0(i0,i1,i2-1))))

                  u2(i0,i1,i2) = nmr/dnr
               endif
            enddo
         enddo
      enddo

c
      return
      end
c

c
