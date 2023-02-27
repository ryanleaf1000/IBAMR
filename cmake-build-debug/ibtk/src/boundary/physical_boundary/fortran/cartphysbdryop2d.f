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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     For cell-centered values, we follow a similar approach as that
c     implemented in class SAMRAI::solv::CartesianRobinBcHelper.  Let
c     u_g denote the ghost cell value and let u_i denote the
c     mirror-image interior cell value, and let n be the number of cell
c     widths separating the ghost cell center and the interior cell
c     center.  We define
c
c         u_b = (u_g + u_i)/2
c         u_n = (u_g - u_i)/(n*h)
c
c     If
c
c         a*u_b + b*u_n = g
c
c     then
c
c         u_g = (-(a*n*h-2*b)/(a*n*h+2*b))*u_i + (2*n*h/(a*n*h+2*b))*g
c             = f_i*u_i + f_g*g
c
c     with
c
c         f_i = -(a*n*h-2*b)/(a*n*h+2*b)
c         f_g = 2*n*h/(a*n*h+2*b)
c
c     For side-centered values, we follow a similar approach.  In this
c     case, however, u_b can be a degree of freedom of the problem, so
c     that
c
c         u_g = u_i + (-a*n*h/b)*u_b + (n*h/b)*g
c             = f_i*u_i + f_b*u_b + f_g*g
c
c     with
c
c         f_i = 1
c         f_b = -a*n*h/b
c         f_g = n*h/b
c
c     For Dirichlet boundary conditions, b=0, and the foregoing
c     expressions are ill defined.  Consequently, in this case, we
c     eliminate u_b and simply set
c
c         u_b = g/a
c         u_g = 2*u_b - u_i
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set cell centered boundary values using the supplied Robin
c     coefficients along the codimension 1 upper/lower x boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ccrobinphysbdryop1x2d(
     &     U,U_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower1,bupper1,
     &     dx,
     &     adjoint_op)
c
      implicit none
c
c     Input.
c
      integer U_gcw

      integer location_index

      integer ilower0,iupper0
      integer ilower1,iupper1

      integer blower1,bupper1

      double precision acoef(blower1:bupper1)
      double precision bcoef(blower1:bupper1)
      double precision gcoef(blower1:bupper1)

      double precision dx(0:2-1)

      integer adjoint_op
c
c     Input/Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &          ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      integer i,i_g,i_i
      integer j
      integer sgn
      double precision    a,b,g,h
      double precision    f_g,f_i,n,u_g,u_i
c
c     Initialize temporary variables to yield errors.
c
      u_g = 2.d0**15.d0
      u_i = 2.d0**15.d0
c
c     Set values along the upper/lower x side of the patch.
c
      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then

         h = dx(location_index/2)

         if (location_index .eq. 0) then
            sgn = -1
            i_g = ilower0-1     ! ghost    index
            i_i = ilower0       ! interior index
         else
            sgn = +1
            i_g = iupper0+1     ! ghost    index
            i_i = iupper0       ! interior index
         endif

         do j = blower1,bupper1
            a = acoef(j)
            b = bcoef(j)
            g = gcoef(j)
            do i = 0,U_gcw-1
               n = 1.d0+2.d0*i
               f_i = -(a*n*h-2.d0*b)/(a*n*h+2.d0*b)
               f_g = 2.d0*n*h/(a*n*h+2.d0*b)
               if (adjoint_op .eq. 1) then
                  u_g = U(i_g+sgn*i,j)
                  U(i_i-sgn*i,j) = U(i_i-sgn*i,j) + f_i*u_g
               else
                  u_i = U(i_i-sgn*i,j)
                  U(i_g+sgn*i,j) = f_i*u_i + f_g*g
               endif
            enddo
         enddo

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set cell centered boundary values using the supplied Robin
c     coefficients along the codimension 1 upper/lower y boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ccrobinphysbdryop1y2d(
     &     U,U_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower0,bupper0,
     &     dx,
     &     adjoint_op)
c
      implicit none
c
c     Input.
c
      integer U_gcw

      integer location_index

      integer ilower0,iupper0
      integer ilower1,iupper1

      integer blower0,bupper0

      double precision acoef(blower0:bupper0)
      double precision bcoef(blower0:bupper0)
      double precision gcoef(blower0:bupper0)

      double precision dx(0:2-1)

      integer adjoint_op
c
c     Input/Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &          ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      integer i
      integer j,j_g,j_i
      integer sgn
      double precision    a,b,g,h
      double precision    f_g,f_i,n,u_g,u_i
c
c     Initialize temporary variables to yield errors.
c
      u_g = 2.d0**15.d0
      u_i = 2.d0**15.d0
c
c     Set values along the upper/lower y side of the patch.
c
      if ( (location_index .eq. 2) .or.
     &     (location_index .eq. 3) ) then

         h = dx(location_index/2)

         if (location_index .eq. 2) then
            sgn = -1
            j_g = ilower1-1     ! ghost    index
            j_i = ilower1       ! interior index
         else
            sgn = +1
            j_g = iupper1+1     ! ghost    index
            j_i = iupper1       ! interior index
         endif

         do i = blower0,bupper0
            a = acoef(i)
            b = bcoef(i)
            g = gcoef(i)
            do j = 0,U_gcw-1
               n = 1.d0+2.d0*j
               f_i = -(a*n*h-2.d0*b)/(a*n*h+2.d0*b)
               f_g = 2.d0*n*h/(a*n*h+2.d0*b)
               if (adjoint_op .eq. 1) then
                  u_g = U(i,j_g+sgn*j)
                  U(i,j_i-sgn*j) = U(i,j_i-sgn*j) + f_i*u_g
               else
                  u_i = U(i,j_i-sgn*j)
                  U(i,j_g+sgn*j) = f_i*u_i + f_g*g
               endif
            enddo
         enddo

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set cell centered boundary values along the codimension 2 boundary
c     by extrapolating values from the codimension 1 boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ccrobinphysbdryop22d(
     &     U,U_gcw,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower0,bupper0,
     &     blower1,bupper1,
     &     adjoint_op)
c
      implicit none
c
c     Input.
c
      integer U_gcw

      integer location_index

      integer ilower0,iupper0
      integer ilower1,iupper1

      integer blower0,bupper0
      integer blower1,bupper1

      integer adjoint_op
c
c     Input/Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &          ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      double precision    U_g
      integer i,i_bdry,i_mirror
      integer j,j_bdry,j_mirror
      integer sgn_x,sgn_y
c
c     Set the codimension 2 boundary values via linear extrapolation.
c
      if     (location_index .eq. 0) then

         i_bdry = ilower0       ! lower x, lower y
         j_bdry = ilower1

         sgn_x = -1
         sgn_y = -1

      elseif (location_index .eq. 1) then

         i_bdry = iupper0       ! upper x, lower y
         j_bdry = ilower1

         sgn_x = +1
         sgn_y = -1

      elseif (location_index .eq. 2) then

         i_bdry = ilower0       ! lower x, upper y
         j_bdry = iupper1

         sgn_x = -1
         sgn_y = +1

      else

         i_bdry = iupper0       ! upper x, upper y
         j_bdry = iupper1

         sgn_x = +1
         sgn_y = +1

      endif

      do j = blower1,bupper1
         j_mirror = j_bdry+(j_bdry-j+sgn_y)
         do i = blower0,bupper0
            i_mirror = i_bdry+(i_bdry-i+sgn_x)
            if (adjoint_op .eq. 1) then
               U_g = U(i,j)
               U(i_mirror,j_mirror) = U(i_mirror,j_mirror) + U_g
               U(i,j_bdry)          = U(i,j_bdry)          + U_g
               U(i_bdry,j)          = U(i_bdry,j)          + U_g
               U(i_mirror,j_bdry)   = U(i_mirror,j_bdry)   - U_g
               U(i_bdry,j_mirror)   = U(i_bdry,j_mirror)   - U_g
            else
               U(i,j) = U(i_mirror,j_mirror)
     &              + (U(i,j_bdry)-U(i_mirror,j_bdry))
     &              + (U(i_bdry,j)-U(i_bdry,j_mirror))
            endif
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set side centered boundary values using the supplied Robin
c     coefficients along the codimension 1 upper/lower x boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrobinphysbdryop1x2d(
     &     u0,u_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower1,bupper1,
     &     dx,
     &     adjoint_op)
c
      implicit none
c
c     Input.
c
      integer u_gcw

      integer location_index

      integer ilower0,iupper0
      integer ilower1,iupper1

      integer blower1,bupper1

      double precision acoef(blower1:bupper1)
      double precision bcoef(blower1:bupper1)
      double precision gcoef(blower1:bupper1)

      double precision dx(0:2-1)

      integer adjoint_op
c
c     Input/Output.
c
      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
c
c     Local variables.
c
      integer i,i_b,i_i
      integer j
      integer sgn
      double precision    a,b,g,h,f_b,f_g,f_i,n,u_b,u_g,u_i
c
c     Initialize temporary variables to yield errors.
c
      u_b = 2.d0**15.d0
      u_g = 2.d0**15.d0
      u_i = 2.d0**15.d0
c
c     Set values along the upper/lower x side of the patch.
c
      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then

         h = dx(location_index/2)

         if (location_index .eq. 0) then
            sgn = -1
            i_b = ilower0       ! boundary index
            i_i = ilower0+1     ! interior index
         else
            sgn = +1
            i_b = iupper0+1     ! boundary index
            i_i = iupper0       ! interior index
         endif

         do j = blower1,bupper1
            a = acoef(j)
            b = bcoef(j)
            g = gcoef(j)
            if (abs(b) .lt. 1.d-12) then
c     Dirichlet boundary conditions
               u_b = g/a
               u0(i_b,j) = u_b
               do i = 1,u_gcw
                  f_i = -1.d0
                  f_b = 2.d0
                  if (adjoint_op .eq. 1) then
                     u_g = u0(i_b+sgn*i,j)
                     u0(i_b-sgn*i,j) = u0(i_b-sgn*i,j) + f_i*u_g
                     u0(i_b      ,j) = u0(i_b      ,j) + f_b*u_g
                  else
                     u_i = u0(i_b-sgn*i,j)
                     u0(i_b+sgn*i,j) = f_i*u_i + f_b*u_b
                  endif
               enddo
            else
c     Robin boundary conditions
               u_b = u0(i_b,j)
               do i = 1,u_gcw
                  n = 2.d0*i
                  f_i = 1.d0
                  f_b = -a*n*h/b
                  f_g = n*h/b
                  if (adjoint_op .eq. 1) then
                     u_g = u0(i_b+sgn*i,j)
                     u0(i_b-sgn*i,j) = u0(i_b-sgn*i,j) + f_i*u_g
                     u0(i_b      ,j) = u0(i_b      ,j) + f_b*u_g
                  else
                     u_i = u0(i_b-sgn*i,j)
                     u0(i_b+sgn*i,j) = f_i*u_i + f_b*u_b + f_g*g
                  endif
               enddo
            endif
         enddo

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set side centered boundary values using the supplied Robin
c     coefficients along the codimension 1 upper/lower y boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrobinphysbdryop1y2d(
     &     u1,u_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower0,bupper0,
     &     dx,
     &     adjoint_op)
c
      implicit none
c
c     Input.
c
      integer u_gcw

      integer location_index

      integer ilower0,iupper0
      integer ilower1,iupper1

      integer blower0,bupper0

      double precision acoef(blower0:bupper0)
      double precision bcoef(blower0:bupper0)
      double precision gcoef(blower0:bupper0)

      double precision dx(0:2-1)

      integer adjoint_op
c
c     Input/Output.
c
      double precision u1(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+1+u_gcw)
c
c     Local variables.
c
      integer i
      integer j,j_b,j_i
      integer sgn
      double precision    a,b,g,h,f_b,f_g,f_i,n,u_b,u_g,u_i
c
c     Initialize temporary variables to yield errors.
c
      u_b = 2.d0**15.d0
      u_g = 2.d0**15.d0
      u_i = 2.d0**15.d0
c
c     Set values along the upper/lower y side of the patch.
c
      if ( (location_index .eq. 2) .or.
     &     (location_index .eq. 3) ) then

         h = dx(location_index/2)

         if (location_index .eq. 2) then
            sgn = -1
            j_b = ilower1       ! boundary index
            j_i = ilower1+1     ! interior index
         else
            sgn = +1
            j_b = iupper1+1     ! boundary index
            j_i = iupper1       ! interior index
         endif

         do i = blower0,bupper0
            a = acoef(i)
            b = bcoef(i)
            g = gcoef(i)
            if (abs(b) .lt. 1.d-12) then
c     Dirichlet boundary conditions
               u_b = g/a
               u1(i,j_b) = u_b
               do j = 1,u_gcw
                  f_i = -1.d0
                  f_b = 2.d0
                  if (adjoint_op .eq. 1) then
                     u_g = u1(i,j_b+sgn*j)
                     u1(i,j_b-sgn*j) = u1(i,j_b-sgn*j) + f_i*u_g
                     u1(i,j_b      ) = u1(i,j_b      ) + f_b*u_g
                  else
                     u_i = u1(i,j_b-sgn*j)
                     u1(i,j_b+sgn*j) = f_i*u_i + f_b*u_b
                  endif
               enddo
            else
c     Robin boundary conditions
               u_b = u1(i,j_b)
               do j = 1,u_gcw
                  n = 2.d0*j
                  f_i = 1.d0
                  f_b = -a*n*h/b
                  f_g = n*h/b
                  if (adjoint_op .eq. 1) then
                     u_g = u1(i,j_b+sgn*j)
                     u1(i,j_b-sgn*j) = u1(i,j_b-sgn*j) + f_i*u_g
                     u1(i,j_b      ) = u1(i,j_b      ) + f_b*u_g
                  else
                     u_i = u1(i,j_b-sgn*j)
                     u1(i,j_b+sgn*j) = f_i*u_i + f_b*u_b + f_g*g
                  endif
               enddo
            endif
         enddo

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set side centered boundary values along the codimension 2 boundary
c     by extrapolating values from the codimension 1 boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrobinphysbdryop22d(
     &     u0,u1,u_gcw,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower0,bupper0,
     &     blower1,bupper1,
     &     adjoint_op)
c
      implicit none
c
c     Input.
c
      integer u_gcw

      integer location_index

      integer ilower0,iupper0
      integer ilower1,iupper1

      integer blower0,bupper0
      integer blower1,bupper1

      integer adjoint_op
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
      integer i,i_bdry,i_shift
      integer j,j_bdry,j_shift
      double precision u_g,del
c
c     Initialize index variables to yield errors in most cases.
c
      i       = 2**15
      i_bdry  = 2**15
      i_shift = 2**15

      j       = 2**15
      j_bdry  = 2**15
      j_shift = 2**15
c
c     Set the codimension 2 boundary values via linear extrapolation.
c
      if     (location_index .eq. 0) then

         i_bdry = ilower0       ! lower x, lower y
         j_bdry = ilower1
         i_shift = +1
         j_shift = +1

      elseif (location_index .eq. 1) then

         i_bdry = iupper0       ! upper x, lower y
         j_bdry = ilower1
         i_shift = -1
         j_shift = +1

      elseif (location_index .eq. 2) then

         i_bdry = ilower0       ! lower x, upper y
         j_bdry = iupper1
         i_shift = +1
         j_shift = -1

      else

         i_bdry = iupper0       ! upper x, upper y
         j_bdry = iupper1
         i_shift = -1
         j_shift = -1

      endif

      do j = blower1,bupper1
         do i = blower0,bupper0+1
            if ( (i .lt. ilower0) .or. (i .gt. iupper0+1) ) then
               del = dble(abs(j-j_bdry))
               if (adjoint_op .eq. 1) then
                  u_g = u0(i,j)
                  u0(i,j_bdry) = u0(i,j_bdry) + (1.d0+del)*u_g
                  u0(i,j_bdry+j_shift) = u0(i,j_bdry+j_shift) - del*u_g
               else
                  u0(i,j) = (1.d0+del)*u0(i,j_bdry)
     &                 - del*u0(i,j_bdry+j_shift)
               endif
            endif
         enddo
      enddo

      do j = blower1,bupper1+1
         if ( (j .lt. ilower1) .or. (j .gt. iupper1+1) ) then
            do i = blower0,bupper0
               del = dble(abs(i-i_bdry))
               if (adjoint_op .eq. 1) then
                  u_g = u1(i,j)
                  u1(i_bdry,j) = u1(i_bdry,j) + (1.d0+del)*u_g
                  u1(i_bdry+i_shift,j) = u1(i_bdry+i_shift,j) - del*u_g
               else
                  u1(i,j) = (1.d0+del)*u1(i_bdry,j)
     &                 - del*u1(i_bdry+i_shift,j)
               endif
            enddo
         endif
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set cell centered boundary values using the supplied Robin
c     coefficients along the codimension 1 upper/lower x boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ccrobinphysbdryquadop1x2d(
     &     U,U_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower1,bupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      integer U_gcw

      integer location_index

      integer ilower0,iupper0
      integer ilower1,iupper1

      integer blower1,bupper1

      double precision acoef(blower1:bupper1)
      double precision bcoef(blower1:bupper1)
      double precision gcoef(blower1:bupper1)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &          ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      integer i,i_g,i_i
      integer j
      integer sgn
      double precision    a,b,g,h
      double precision    f_g,f_i,f_i2,n,u_g,u_i,u_i2
      double precision    den
c
c     Initialize temporary variables to yield errors.
c
      u_g = 2.d0**15.d0
      u_i = 2.d0**15.d0
c
c     Set values along the upper/lower x side of the patch.
c
      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then

         h = dx(location_index/2)

         if (location_index .eq. 0) then
            sgn = -1
            i_g = ilower0-1     ! ghost    index
            i_i = ilower0       ! interior index
         else
            sgn = +1
            i_g = iupper0+1     ! ghost    index
            i_i = iupper0       ! interior index
         endif

         do j = blower1,bupper1
            a = acoef(j)
            b = bcoef(j)
            g = gcoef(j)
            do i = 0,U_gcw-1
               n = 1.d0+2.d0*i
               den = 1.d0/(4.d0*b*(1.d0+n)+a*h*n*(2.d0+n))
               f_i = den*(1.d0+n)*(4.d0*b-a*h*n*(2.d0+n))
               f_g = den*(4.d0*h*n*(1.d0+n))
               f_i2 = den*(a*h*n*n*n)
               u_i = U(i_i-sgn*i,j)
               u_i2 = U(i_i-sgn*(i+1),j)
               U(i_g+sgn*i,j) = f_i*u_i + f_g*g + f_i2*u_i2
            enddo
         enddo

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set cell centered boundary values using the supplied Robin
c     coefficients along the codimension 1 upper/lower y boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ccrobinphysbdryquadop1y2d(
     &     U,U_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower0,bupper0,
     &     dx)
c
      implicit none
c
c     Input.
c
      integer U_gcw

      integer location_index

      integer ilower0,iupper0
      integer ilower1,iupper1

      integer blower0,bupper0

      double precision acoef(blower0:bupper0)
      double precision bcoef(blower0:bupper0)
      double precision gcoef(blower0:bupper0)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &          ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      integer i
      integer j,j_g,j_i
      integer sgn
      double precision    a,b,g,h
      double precision    f_g,f_i,f_i2,n,u_g,u_i,u_i2
      double precision    den
c
c     Initialize temporary variables to yield errors.
c
      u_g = 2.d0**15.d0
      u_i = 2.d0**15.d0
c
c     Set values along the upper/lower y side of the patch.
c
      if ( (location_index .eq. 2) .or.
     &     (location_index .eq. 3) ) then

         h = dx(location_index/2)

         if (location_index .eq. 2) then
            sgn = -1
            j_g = ilower1-1     ! ghost    index
            j_i = ilower1       ! interior index
         else
            sgn = +1
            j_g = iupper1+1     ! ghost    index
            j_i = iupper1       ! interior index
         endif

         do i = blower0,bupper0
            a = acoef(i)
            b = bcoef(i)
            g = gcoef(i)
            do j = 0,U_gcw-1
               n = 1.d0+2.d0*j
               den = 1.d0/(4.d0*b*(1.d0+n)+a*h*n*(2.d0+n))
               f_i = den*(1.d0+n)*(4.d0*b-a*h*n*(2.d0+n))
               f_g = den*4.d0*h*n*(1.d0+n)
               f_i2 = den*a*h*n*n*n
               u_i = U(i,j_i-sgn*j)
               u_i2 = U(i,j_i-sgn*(j+1))
               U(i,j_g+sgn*j) = f_i*u_i + f_g*g + f_i2*u_i2
            enddo
         enddo

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set side centered boundary values using the supplied Robin
c     coefficients along the codimension 1 upper/lower x boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrobinphysbdryquadop1x2d(
     &     u0,u_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower1,bupper1,
     &     dx)
c
      implicit none
c
c     Input.
c
      integer u_gcw

      integer location_index

      integer ilower0,iupper0
      integer ilower1,iupper1

      integer blower1,bupper1

      double precision acoef(blower1:bupper1)
      double precision bcoef(blower1:bupper1)
      double precision gcoef(blower1:bupper1)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
c
c     Local variables.
c
      integer i,i_b,i_i
      integer j
      integer sgn
      double precision    a,b,g,h,f_b,f_g,f_i,f_i2,n,u_b,u_g,u_i,u_i2
c
c     Initialize temporary variables to yield errors.
c
      u_b = 2.d0**15.d0
      u_g = 2.d0**15.d0
      u_i = 2.d0**15.d0
      u_i2 = 2.d0**15.d0
c
c     Set values along the upper/lower x side of the patch.
c
      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then

         h = dx(location_index/2)

         if (location_index .eq. 0) then
            sgn = -1
            i_b = ilower0       ! boundary index
            i_i = ilower0+1     ! interior index
         else
            sgn = +1
            i_b = iupper0+1     ! boundary index
            i_i = iupper0       ! interior index
         endif

         do j = blower1,bupper1
            a = acoef(j)
            b = bcoef(j)
            g = gcoef(j)
            if (abs(b) .lt. 1.d-12) then
c     Dirichlet boundary conditions
               u_b = g/a
               u0(i_b,j) = u_b
               do i = 1,u_gcw
                  f_i2 = 1.d0
                  f_i = -3.d0
                  f_b = 3.d0
                  u_i = u0(i_b-sgn*i,j)
                  u_i2 = u0(i_b-sgn*(i+1),j)
                  u0(i_b+sgn*i,j) = f_i2*u_i2 + f_i*u_i + f_b*u_b
               enddo
            else
c     Robin boundary conditions
               u_b = u0(i_b,j)
               do i = 1,u_gcw
                  n = 2.d0*i
                  f_i = 1.d0
                  f_b = -a*n*h/b
                  f_g = n*h/b
                  u_i = u0(i_b-sgn*i,j)
                  u0(i_b+sgn*i,j) = f_i*u_i + f_b*u_b + f_g*g
               enddo
            endif
         enddo

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Set side centered boundary values using the supplied Robin
c     coefficients along the codimension 1 upper/lower y boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrobinphysbdryquadop1y2d(
     &     u1,u_gcw,
     &     acoef,bcoef,gcoef,
     &     location_index,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     blower0,bupper0,
     &     dx)
c
      implicit none
c
c     Input.
c
      integer u_gcw

      integer location_index

      integer ilower0,iupper0
      integer ilower1,iupper1

      integer blower0,bupper0

      double precision acoef(blower0:bupper0)
      double precision bcoef(blower0:bupper0)
      double precision gcoef(blower0:bupper0)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision u1(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+1+u_gcw)
c
c     Local variables.
c
      integer i
      integer j,j_b,j_i
      integer sgn
      double precision    a,b,g,h,f_b,f_g,f_i,f_i2,n,u_b,u_g,u_i,u_i2
c
c     Initialize temporary variables to yield errors.
c
      u_b = 2.d0**15.d0
      u_g = 2.d0**15.d0
      u_i = 2.d0**15.d0
      u_i2 = 2.d0**15.d0
c
c     Set values along the upper/lower y side of the patch.
c
      if ( (location_index .eq. 2) .or.
     &     (location_index .eq. 3) ) then

         h = dx(location_index/2)

         if (location_index .eq. 2) then
            sgn = -1
            j_b = ilower1       ! boundary index
            j_i = ilower1+1     ! interior index
         else
            sgn = +1
            j_b = iupper1+1     ! boundary index
            j_i = iupper1       ! interior index
         endif

         do i = blower0,bupper0
            a = acoef(i)
            b = bcoef(i)
            g = gcoef(i)
            if (abs(b) .lt. 1.d-12) then
c     Dirichlet boundary conditions
               u_b = g/a
               u1(i,j_b) = u_b
               do j = 1,u_gcw
                  f_i2 = 1.d0
                  f_i = -3.d0
                  f_b = 3.d0
                  u_i = u1(i,j_b-sgn*j)
                  u_i2 = u1(i,j_b-sgn*(j+1))
                  u1(i,j_b+sgn*j) = f_i2*u_i2 + f_i*u_i + f_b*u_b
               enddo
            else
c     Robin boundary conditions
               u_b = u1(i,j_b)
               do j = 1,u_gcw
                  n = 2.d0*j
                  f_i = 1.d0
                  f_b = -a*n*h/b
                  f_g = n*h/b
                  u_i = u1(i,j_b-sgn*j)
                  u1(i,j_b+sgn*j) = f_i*u_i + f_b*u_b + f_g*g
               enddo
            endif
         enddo

      endif
c
      return
      end
c
