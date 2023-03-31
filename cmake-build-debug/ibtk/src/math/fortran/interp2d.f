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
c     Compute the face centered normal vector field (u0,u1) from the
c     cell centered vector field V using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctofinterp2nd2d(
     &     u0,u1,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      integer u_gcw,V_gcw

      integer ilower0,iupper0
      integer ilower1,iupper1

      double precision V(ilower0-V_gcw:iupper0+V_gcw,
     &          ilower1-V_gcw:iupper1+V_gcw,0:2-1)
c
c     Output.
c
      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
      double precision u1(ilower1-u_gcw:iupper1+1+u_gcw,
     &          ilower0-u_gcw:iupper0+u_gcw)
c
c     Local variables.
c
      integer i0,i1
c
c     Compute the face centered vector field (u0,u1) from the cell
c     centered vector field V.
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            u0(i0,i1) = 0.5d0*(V(i0-1,i1,0)+V(i0,i1,0))
         enddo
      enddo
      do i0 = ilower0,iupper0
         do i1 = ilower1,iupper1+1
            u1(i1,i0) = 0.5d0*(V(i0,i1-1,1)+V(i0,i1,1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the face centered normal vector field (u0,u1) from the
c     cell centered scalar field V using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctofcwiseinterp2nd2d(
     &     u0,u1,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      integer u_gcw,V_gcw

      integer ilower0,iupper0
      integer ilower1,iupper1

      double precision V(ilower0-V_gcw:iupper0+V_gcw,
     &          ilower1-V_gcw:iupper1+V_gcw)
c
c     Output.
c
      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
      double precision u1(ilower1-u_gcw:iupper1+1+u_gcw,
     &          ilower0-u_gcw:iupper0+u_gcw)
c
c     Local variables.
c
      integer i0,i1
c
c     Compute the face centered vector field (u0,u1) from the cell
c     centered scalar field V.
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            u0(i0,i1) = 0.5d0*(V(i0-1,i1)+V(i0,i1))
         enddo
      enddo
      do i0 = ilower0,iupper0
         do i1 = ilower1,iupper1+1
            u1(i1,i0) = 0.5d0*(V(i0,i1-1)+V(i0,i1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the side centered normal vector field (u0,u1) from the
c     cell centered vector field V using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctosinterp2nd2d(
     &     u0,u1,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      integer u_gcw,V_gcw

      integer ilower0,iupper0
      integer ilower1,iupper1

      double precision V(ilower0-V_gcw:iupper0+V_gcw,
     &          ilower1-V_gcw:iupper1+V_gcw,0:2-1)
c
c     Output.
c
      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
      double precision u1(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+1+u_gcw)
c
c     Local variables.
c
      integer i0,i1
c
c     Compute the side centered vector field (u0,u1) from the cell
c     centered vector field V.
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            u0(i0,i1) = 0.5d0*(V(i0-1,i1,0)+V(i0,i1,0))
         enddo
      enddo
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            u1(i0,i1) = 0.5d0*(V(i0,i1-1,1)+V(i0,i1,1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the side centered normal vector field (u0,u1) from the
c     cell centered scalar field V using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctoscwiseinterp2nd2d(
     &     u0,u1,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      integer u_gcw,V_gcw

      integer ilower0,iupper0
      integer ilower1,iupper1

      double precision V(ilower0-V_gcw:iupper0+V_gcw,
     &          ilower1-V_gcw:iupper1+V_gcw)
c
c     Output.
c
      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
      double precision u1(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+1+u_gcw)
c
c     Local variables.
c
      integer i0,i1
c
c     Compute the side centered vector field (u0,u1) from the cell
c     centered scalar field V.
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            u0(i0,i1) = 0.5d0*(V(i0-1,i1)+V(i0,i1))
         enddo
      enddo
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            u1(i0,i1) = 0.5d0*(V(i0,i1-1)+V(i0,i1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the cell centered vector field U from the face centered
c     normal vector field (v0,v1) using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ftocinterp2nd2d(
     &     U,U_gcw,
     &     v0,v1,v_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      integer U_gcw,v_gcw

      integer ilower0,iupper0
      integer ilower1,iupper1

      double precision v0(ilower0-v_gcw:iupper0+1+v_gcw,
     &          ilower1-v_gcw:iupper1+v_gcw)
      double precision v1(ilower1-v_gcw:iupper1+1+v_gcw,
     &          ilower0-v_gcw:iupper0+v_gcw)
c
c     Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &          ilower1-U_gcw:iupper1+U_gcw,0:2-1)
c
c     Local variables.
c
      integer i0,i1
c
c     Compute the cell centered vector field U from the face centered
c     vector field (v0,v1).
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            U(i0,i1,0) = 0.5d0*(v0(i0,i1)+v0(i0+1,i1))
         enddo
      enddo
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            U(i0,i1,1) = 0.5d0*(v1(i1,i0)+v1(i1+1,i0))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the cell centered vector field U from the side centered
c     normal vector field (v0,v1) using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stocinterp2nd2d(
     &     U,U_gcw,
     &     v0,v1,v_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      integer U_gcw,v_gcw

      integer ilower0,iupper0
      integer ilower1,iupper1

      double precision v0(ilower0-v_gcw:iupper0+1+v_gcw,
     &          ilower1-v_gcw:iupper1+v_gcw)
      double precision v1(ilower0-v_gcw:iupper0+v_gcw,
     &          ilower1-v_gcw:iupper1+1+v_gcw)
c
c     Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &          ilower1-U_gcw:iupper1+U_gcw,0:2-1)
c
c     Local variables.
c
      integer i0,i1
c
c     Compute the cell centered vector field U from the side centered
c     vector field (v0,v1).
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            U(i0,i1,0) = 0.5d0*(v0(i0,i1)+v0(i0+1,i1))
         enddo
      enddo
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            U(i0,i1,1) = 0.5d0*(v1(i0,i1)+v1(i0,i1+1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the node centered vector field U from the face centered
c     normal vector field (v0,v1) using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ftoninterp2nd2d(
     &     U,U_gcw,
     &     v0,v1,v_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      integer U_gcw,v_gcw

      integer ilower0,iupper0
      integer ilower1,iupper1

      double precision v0(ilower0-v_gcw:iupper0+1+v_gcw,
     &          ilower1-v_gcw:iupper1+v_gcw)
      double precision v1(ilower1-v_gcw:iupper1+1+v_gcw,
     &          ilower0-v_gcw:iupper0+v_gcw)
c
c     Output.
c
      double precision U(ilower0-U_gcw:iupper0+1+U_gcw,
     &          ilower1-U_gcw:iupper1+1+U_gcw,0:2-1)
c
c     Local variables.
c
      integer i0,i1
c
c     Compute the node centered vector field U from the face centered
c     vector field (v0,v1).
c
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0+1
            U(i0,i1,0) = 0.5d0*(v0(i0,i1-1)+v0(i0,i1))
         enddo
      enddo
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0+1
            U(i0,i1,1) = 0.5d0*(v1(i1,i0-1)+v1(i1,i0))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the cell centered vector field U from the side centered
c     normal vector field (v0,v1) using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stocinterp2ndspecial2d(
     &     direction,
     &     U,U_gcw,
     &     alpha,
     &     v0,v1,v_gcw,
     &     beta,
     &     W,W_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      integer direction

      integer U_gcw,v_gcw,W_gcw

      integer ilower0,iupper0
      integer ilower1,iupper1

      double precision alpha

      double precision v0(ilower0-v_gcw:iupper0+1+v_gcw,
     &          ilower1-v_gcw:iupper1+v_gcw)
      double precision v1(ilower0-v_gcw:iupper0+v_gcw,
     &          ilower1-v_gcw:iupper1+1+v_gcw)

      double precision beta

      double precision W(ilower0-W_gcw:iupper0+W_gcw,
     &          ilower1-W_gcw:iupper1+W_gcw)
c
c     Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &          ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      integer i0,i1
c
c     Compute the cell centered vector field U from the side centered
c     vector field (v0,v1).
c
      if ( direction.eq.0 ) then
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               U(i0,i1) = 0.5d0*alpha*(v0(i0,i1)+v0(i0+1,i1))
     &              + beta*W(i0,i1)
            enddo
         enddo
      else
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0
               U(i0,i1) = 0.5d0*alpha*(v1(i0,i1)+v1(i0,i1+1))
     &              + beta*W(i0,i1)
            enddo
         enddo
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the node centered vector field U from the side centered
c     normal vector field (v0,v1) using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stoninterp2nd2d(
     &     U,U_gcw,
     &     v0,v1,v_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      integer U_gcw,v_gcw

      integer ilower0,iupper0
      integer ilower1,iupper1

      double precision v0(ilower0-v_gcw:iupper0+1+v_gcw,
     &          ilower1-v_gcw:iupper1+v_gcw)
      double precision v1(ilower0-v_gcw:iupper0+v_gcw,
     &          ilower1-v_gcw:iupper1+1+v_gcw)
c
c     Output.
c
      double precision U(ilower0-U_gcw:iupper0+1+U_gcw,
     &          ilower1-U_gcw:iupper1+1+U_gcw,0:2-1)
c
c     Local variables.
c
      integer i0,i1
c
c     Compute the node centered vector field U from the side centered
c     vector field (v0,v1).
c
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0+1
            U(i0,i1,0) = 0.5d0*(v0(i0,i1-1)+v0(i0,i1))
         enddo
      enddo
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0+1
            U(i0,i1,1) = 0.5d0*(v1(i0-1,i1)+v1(i0,i1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the cell centered field U from the node centered
c     field V using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ntocinterp2d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Input.
c
      integer U_gcw,V_gcw

      integer ilower0,iupper0
      integer ilower1,iupper1

      double precision V(ilower0-v_gcw:iupper0+1+v_gcw,
     &          ilower1-v_gcw:iupper1+1+v_gcw)
c
c     Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &          ilower1-U_gcw:iupper1+U_gcw)
c
c     Local variables.
c
      integer i0,i1
c
c     Compute the cell centered field U from the node centered
c     field v0.
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0
            U(i0,i1) = 0.25d0*(V(i0,i1)+V(i0+1,i1)+
     &                         V(i0,i1+1)+V(i0+1,i1+1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the node centered field U from the cell centered
c     field V using simple averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctoninterp2d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     U_ghost_interp)
c
      implicit none
c
c     Input.
c
      integer U_gcw,V_gcw
      integer U_ghost_interp

      integer ilower0,iupper0
      integer ilower1,iupper1

      double precision V(ilower0-V_gcw:iupper0+V_gcw,
     &          ilower1-V_gcw:iupper1+V_gcw)

c
c     Output.
c
      double precision U(ilower0-U_gcw:iupper0+1+U_gcw,
     &          ilower1-U_gcw:iupper1+1+U_gcw)
c
c     Local variables.
c
      integer i0,i1
      integer gcw_shift
c
c     Compute the node centered scalar field U from the cell centered
c     scalar field V.
c
      gcw_shift = 0
      if (U_ghost_interp .eq. 1) then
         gcw_shift = U_gcw
      endif

      do i1 = ilower1-gcw_shift,iupper1+gcw_shift+1
         do i0 = ilower0-gcw_shift,iupper0+gcw_shift+1
            U(i0,i1) = 0.25d0*(V(i0,i1)+V(i0-1,i1)
     &            +V(i0,i1-1)+V(i0-1,i1-1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the side centered normal vector field (u0,u1) from the
c     cell centered vector field V using harmonic averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctosharmonicinterp2nd2d(
     &     u0,u1,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Function.
c
      double precision h_avg2
c
c     Input.
c
      integer u_gcw,V_gcw

      integer ilower0,iupper0
      integer ilower1,iupper1

      double precision V(ilower0-V_gcw:iupper0+V_gcw,
     &          ilower1-V_gcw:iupper1+V_gcw,0:2-1)
c
c     Output.
c
      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
      double precision u1(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+1+u_gcw)
c
c     Local variables.
c
      integer i0,i1
c
c     Compute the side centered vector field (u0,u1) from the cell
c     centered vector field V.
c
      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            u0(i0,i1) = h_avg2(V(i0-1,i1,0),V(i0,i1,0))
         enddo
      enddo
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            u1(i0,i1) = h_avg2(V(i0,i1-1,1),V(i0,i1,1))
         enddo
      enddo
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the node centered field U from the cell centered
c     field V using harmonic averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctonharmonicinterp2d(
     &     U,U_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     U_ghost_interp)
c
      implicit none
c
c     Function.
c
      double precision h_avg4
c
c     Input.
c
      integer U_gcw,V_gcw
      integer U_ghost_interp

      integer ilower0,iupper0
      integer ilower1,iupper1

      double precision V(ilower0-V_gcw:iupper0+V_gcw,
     &          ilower1-V_gcw:iupper1+V_gcw)

c
c     Output.
c
      double precision U(ilower0-U_gcw:iupper0+1+U_gcw,
     &          ilower1-U_gcw:iupper1+1+U_gcw)
c
c     Local variables.
c
      integer i0,i1
      integer gcw_shift
c
c     Compute the node centered scalar field U from the cell centered
c     scalar field V.
c
      gcw_shift = 0
      if (U_ghost_interp .eq. 1) then
         gcw_shift = U_gcw
      endif

      do i1 = ilower1-gcw_shift,iupper1+gcw_shift+1
         do i0 = ilower0-gcw_shift,iupper0+gcw_shift+1
            U(i0,i1) = h_avg4(V(i0,i1),V(i0-1,i1),V(i0,i1-1),
     &                  V(i0-1,i1-1))
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute the side centered normal vector field (u0,u1) from the
c     cell centered scalar field V using harmonic averaging.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine ctoscwiseharmonicinterp2nd2d(
     &     u0,u1,u_gcw,
     &     V,V_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1)
c
      implicit none
c
c     Function.
c
      double precision h_avg2
c
c     Input.
c
      integer u_gcw,V_gcw

      integer ilower0,iupper0
      integer ilower1,iupper1

      double precision V(ilower0-V_gcw:iupper0+V_gcw,
     &          ilower1-V_gcw:iupper1+V_gcw)
c
c     Output.
c
      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
      double precision u1(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+1+u_gcw)
c
c     Local variables.
c
      integer i0,i1

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1
            u0(i0, i1) = h_avg2(V(i0-1,i1),V(i0,i1))
         enddo
      enddo
      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0
            u1(i0,i1) = h_avg2(V(i0,i1-1),V(i0,i1))
         enddo
      enddo
c
      return
      end
c
