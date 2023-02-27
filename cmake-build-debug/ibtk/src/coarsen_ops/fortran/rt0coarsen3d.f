c ---------------------------------------------------------------------
c
c Copyright (c) 2015 - 2020 by the IBAMR developers
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
c     Coarsen side-centered data via the adjoint of RT0 interpolation.
c
c     NOTE: Values at physical boundaries and coarse-fine interfaces
c     will need to be corrected.  This routine uses the correct stencil
c     only *away* from such boundaries.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrt0coarsen3d(
     &     U_crse0,U_crse1,U_crse2,U_crse_gcw,
     &     U_fine0,U_fine1,U_fine2,U_fine_gcw,
     &     ilowerc0,iupperc0,
     &     ilowerc1,iupperc1,
     &     ilowerc2,iupperc2,
     &     ilowerf0,iupperf0,
     &     ilowerf1,iupperf1,
     &     ilowerf2,iupperf2,
     &     ratio_to_coarser,
     &     coarse_box_lower,coarse_box_upper)
c
      implicit none
c
c     Input.
c
      integer ilowerc0,iupperc0
      integer ilowerc1,iupperc1
      integer ilowerc2,iupperc2
      integer ilowerf0,iupperf0
      integer ilowerf1,iupperf1
      integer ilowerf2,iupperf2
      integer U_crse_gcw,U_fine_gcw

      integer ratio_to_coarser(0:3-1)

      integer coarse_box_lower(0:3-1), coarse_box_upper(0:3-1)

      double precision U_fine0(
     &     ilowerf0-U_fine_gcw:iupperf0+1+U_fine_gcw,
     &          ilowerf1-U_fine_gcw:iupperf1+U_fine_gcw,
     &          ilowerf2-U_fine_gcw:iupperf2+U_fine_gcw
     &     )
      double precision U_fine1(
     &     ilowerf0-U_fine_gcw:iupperf0+U_fine_gcw,
     &          ilowerf1-U_fine_gcw:iupperf1+1+U_fine_gcw,
     &          ilowerf2-U_fine_gcw:iupperf2+U_fine_gcw
     &     )
      double precision U_fine2(
     &     ilowerf0-U_fine_gcw:iupperf0+U_fine_gcw,
     &          ilowerf1-U_fine_gcw:iupperf1+U_fine_gcw,
     &          ilowerf2-U_fine_gcw:iupperf2+1+U_fine_gcw
     &     )
c
c     Input/Output.
c
      double precision U_crse0(
     &     ilowerc0-U_crse_gcw:iupperc0+1+U_crse_gcw,
     &          ilowerc1-U_crse_gcw:iupperc1+U_crse_gcw,
     &          ilowerc2-U_crse_gcw:iupperc2+U_crse_gcw
     &     )
      double precision U_crse1(
     &     ilowerc0-U_crse_gcw:iupperc0+U_crse_gcw,
     &          ilowerc1-U_crse_gcw:iupperc1+1+U_crse_gcw,
     &          ilowerc2-U_crse_gcw:iupperc2+U_crse_gcw
     &     )
      double precision U_crse2(
     &     ilowerc0-U_crse_gcw:iupperc0+U_crse_gcw,
     &          ilowerc1-U_crse_gcw:iupperc1+U_crse_gcw,
     &          ilowerc2-U_crse_gcw:iupperc2+1+U_crse_gcw
     &     )
c
c     Local variables.
c
      integer i,i_c,i_f,j,j_c,j_f,k,k_c,k_f
      integer axis,d,istencil_lower(0:3-1),istencil_upper(0:3-1)
      double precision w,w_fac

c
c     Treat x components:
c
      axis = 0
      do d = 0,3-1
         if (d .eq. axis) then
            istencil_lower(d) = -ratio_to_coarser(d)+1
            istencil_upper(d) = +ratio_to_coarser(d)-1
         else
            istencil_lower(d) = 0
            istencil_upper(d) = ratio_to_coarser(d)-1
         endif
      enddo

      w_fac = 0.d0
      do k = istencil_lower(2),istencil_upper(2)
         do j = istencil_lower(1),istencil_upper(1)
            do i = istencil_lower(0),istencil_upper(0)
               w = 1.d0 - dabs(dble(i))/dble(ratio_to_coarser(axis))
               w_fac = w_fac + w
            enddo
         enddo
      enddo
      w_fac = 1.d0/w_fac

      do k_c = coarse_box_lower(2),coarse_box_upper(2)
         do j_c = coarse_box_lower(1),coarse_box_upper(1)
            do i_c = coarse_box_lower(0),coarse_box_upper(0)+1
               i_f = i_c*ratio_to_coarser(0)
               j_f = j_c*ratio_to_coarser(1)
               k_f = k_c*ratio_to_coarser(2)
               U_crse0(i_c,j_c,k_c) = 0.d0
               do k = istencil_lower(2),istencil_upper(2)
                  do j = istencil_lower(1),istencil_upper(1)
                     do i = istencil_lower(0),istencil_upper(0)
                        w = 1.d0 - dabs(dble(i))
     &                       /dble(ratio_to_coarser(axis))
                        U_crse0(i_c,j_c,k_c) = U_crse0(i_c,j_c,k_c) +
     &                       w*w_fac*U_fine0(i_f+i,j_f+j,k_f+k)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
c     Treat y components:
c
      axis = 1
      do d = 0,3-1
         if (d .eq. axis) then
            istencil_lower(d) = -ratio_to_coarser(d)+1
            istencil_upper(d) = +ratio_to_coarser(d)-1
         else
            istencil_lower(d) = 0
            istencil_upper(d) = ratio_to_coarser(d)-1
         endif
      enddo

      w_fac = 0.d0
      do k = istencil_lower(2),istencil_upper(2)
         do j = istencil_lower(1),istencil_upper(1)
            do i = istencil_lower(0),istencil_upper(0)
               w = 1.d0 - dabs(dble(j))/dble(ratio_to_coarser(axis))
               w_fac = w_fac + w
            enddo
         enddo
      enddo
      w_fac = 1.d0/w_fac

      do k_c = coarse_box_lower(2),coarse_box_upper(2)
         do j_c = coarse_box_lower(1),coarse_box_upper(1)+1
            do i_c = coarse_box_lower(0),coarse_box_upper(0)
               i_f = i_c*ratio_to_coarser(0)
               j_f = j_c*ratio_to_coarser(1)
               k_f = k_c*ratio_to_coarser(2)
               U_crse1(i_c,j_c,k_c) = 0.d0
               do k = istencil_lower(2),istencil_upper(2)
                  do j = istencil_lower(1),istencil_upper(1)
                     do i = istencil_lower(0),istencil_upper(0)
                        w = 1.d0 - dabs(dble(j))
     &                       /dble(ratio_to_coarser(axis))
                        U_crse1(i_c,j_c,k_c) = U_crse1(i_c,j_c,k_c) +
     &                       w*w_fac*U_fine1(i_f+i,j_f+j,k_f+k)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
c     Treat z components:
c
      axis = 2
      do d = 0,3-1
         if (d .eq. axis) then
            istencil_lower(d) = -ratio_to_coarser(d)+1
            istencil_upper(d) = +ratio_to_coarser(d)-1
         else
            istencil_lower(d) = 0
            istencil_upper(d) = ratio_to_coarser(d)-1
         endif
      enddo

      w_fac = 0.d0
      do k = istencil_lower(2),istencil_upper(2)
         do j = istencil_lower(1),istencil_upper(1)
            do i = istencil_lower(0),istencil_upper(0)
               w = 1.d0 - dabs(dble(k))/dble(ratio_to_coarser(axis))
               w_fac = w_fac + w
            enddo
         enddo
      enddo
      w_fac = 1.d0/w_fac

      do k_c = coarse_box_lower(2),coarse_box_upper(2)+1
         do j_c = coarse_box_lower(1),coarse_box_upper(1)
            do i_c = coarse_box_lower(0),coarse_box_upper(0)
               i_f = i_c*ratio_to_coarser(0)
               j_f = j_c*ratio_to_coarser(1)
               k_f = k_c*ratio_to_coarser(2)
               U_crse2(i_c,j_c,k_c) = 0.d0
               do k = istencil_lower(2),istencil_upper(2)
                  do j = istencil_lower(1),istencil_upper(1)
                     do i = istencil_lower(0),istencil_upper(0)
                        w = 1.d0 - dabs(dble(k))
     &                       /dble(ratio_to_coarser(axis))
                        U_crse2(i_c,j_c,k_c) = U_crse2(i_c,j_c,k_c) +
     &                       w*w_fac*U_fine2(i_f+i,j_f+j,k_f+k)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Fix up RT0 coarsening along a boundary.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scrt0coarsenbdry3d(
     &     U_crse0,U_crse1,U_crse2,U_crse_gcw,
     &     U_fine0,U_fine1,U_fine2,U_fine_gcw,
     &     ilowerc0,iupperc0,
     &     ilowerc1,iupperc1,
     &     ilowerc2,iupperc2,
     &     ilowerf0,iupperf0,
     &     ilowerf1,iupperf1,
     &     ilowerf2,iupperf2,
     &     ratio_to_coarser,
     &     coarse_box_lower,coarse_box_upper,
     &     bbox_ilowerf0,bbox_iupperf0,
     &     bbox_ilowerf1,bbox_iupperf1,
     &     bbox_ilowerf2,bbox_iupperf2,
     &     axis,upperlower)
c
      implicit none
c
c     Input.
c
      integer ilowerc0,iupperc0
      integer ilowerc1,iupperc1
      integer ilowerc2,iupperc2
      integer ilowerf0,iupperf0
      integer ilowerf1,iupperf1
      integer ilowerf2,iupperf2
      integer U_crse_gcw,U_fine_gcw

      integer ratio_to_coarser(0:3-1)

      integer coarse_box_lower(0:3-1), coarse_box_upper(0:3-1)

      integer bbox_ilowerf0,bbox_iupperf0
      integer bbox_ilowerf1,bbox_iupperf1
      integer bbox_ilowerf2,bbox_iupperf2
      integer axis,upperlower

      double precision U_fine0(
     &     ilowerf0-U_fine_gcw:iupperf0+1+U_fine_gcw,
     &          ilowerf1-U_fine_gcw:iupperf1+U_fine_gcw,
     &          ilowerf2-U_fine_gcw:iupperf2+U_fine_gcw
     &     )
      double precision U_fine1(
     &     ilowerf0-U_fine_gcw:iupperf0+U_fine_gcw,
     &          ilowerf1-U_fine_gcw:iupperf1+1+U_fine_gcw,
     &          ilowerf2-U_fine_gcw:iupperf2+U_fine_gcw
     &     )
      double precision U_fine2(
     &     ilowerf0-U_fine_gcw:iupperf0+U_fine_gcw,
     &          ilowerf1-U_fine_gcw:iupperf1+U_fine_gcw,
     &          ilowerf2-U_fine_gcw:iupperf2+1+U_fine_gcw
     &     )
c
c     Input/Output.
c
      double precision U_crse0(
     &     ilowerc0-U_crse_gcw:iupperc0+1+U_crse_gcw,
     &          ilowerc1-U_crse_gcw:iupperc1+U_crse_gcw,
     &          ilowerc2-U_crse_gcw:iupperc2+U_crse_gcw
     &     )
      double precision U_crse1(
     &     ilowerc0-U_crse_gcw:iupperc0+U_crse_gcw,
     &          ilowerc1-U_crse_gcw:iupperc1+1+U_crse_gcw,
     &          ilowerc2-U_crse_gcw:iupperc2+U_crse_gcw
     &     )
      double precision U_crse2(
     &     ilowerc0-U_crse_gcw:iupperc0+U_crse_gcw,
     &          ilowerc1-U_crse_gcw:iupperc1+U_crse_gcw,
     &          ilowerc2-U_crse_gcw:iupperc2+1+U_crse_gcw
     &     )
c
c     Local variables.
c
      integer i,i_c,i_f,j,j_c,j_f,k,k_c,k_f
c     NOTE: We require the boundary box to be a "side centered" box.
      integer bbox_ilowerc(0:3-1),bbox_iupperc(0:3-1)
      integer d,istencil_lower(0:3-1),istencil_upper(0:3-1)
      integer ibdryc(0:3-1)
      double precision w,w_fac

c
c     Prevent compiler warning about unused variables.
c
      coarse_box_lower(0) = coarse_box_lower(0)
      coarse_box_upper(0) = coarse_box_upper(0)

c
c     Setup boundary box extents.
c
      if (bbox_ilowerf0.lt.0) then
            bbox_ilowerc(0)=(bbox_ilowerf0+1)/ratio_to_coarser(0)-1
         else
            bbox_ilowerc(0)=bbox_ilowerf0/ratio_to_coarser(0)
         endif

      if (bbox_ilowerf1.lt.0) then
            bbox_ilowerc(1)=(bbox_ilowerf1+1)/ratio_to_coarser(1)-1
         else
            bbox_ilowerc(1)=bbox_ilowerf1/ratio_to_coarser(1)
         endif

      if (bbox_ilowerf2.lt.0) then
            bbox_ilowerc(2)=(bbox_ilowerf2+1)/ratio_to_coarser(2)-1
         else
            bbox_ilowerc(2)=bbox_ilowerf2/ratio_to_coarser(2)
         endif

      if (bbox_iupperf0.lt.0) then
            bbox_iupperc(0)=(bbox_iupperf0+1)/ratio_to_coarser(0)-1
         else
            bbox_iupperc(0)=bbox_iupperf0/ratio_to_coarser(0)
         endif

      if (bbox_iupperf1.lt.0) then
            bbox_iupperc(1)=(bbox_iupperf1+1)/ratio_to_coarser(1)-1
         else
            bbox_iupperc(1)=bbox_iupperf1/ratio_to_coarser(1)
         endif

      if (bbox_iupperf2.lt.0) then
            bbox_iupperc(2)=(bbox_iupperf2+1)/ratio_to_coarser(2)-1
         else
            bbox_iupperc(2)=bbox_iupperf2/ratio_to_coarser(2)
         endif

c
c     Setup stencil indices.
c
      do d = 0,3-1
         if (d .eq. axis) then
            if (upperlower .eq. 0)  then
               istencil_lower(d) = 0
               istencil_upper(d) = ratio_to_coarser(d)-1
            else if (upperlower .eq. 1)  then
               istencil_lower(d) = -ratio_to_coarser(d)+1
               istencil_upper(d) = 0
            endif
         else
            istencil_lower(d) = 0
            istencil_upper(d) = ratio_to_coarser(d)-1
         endif
      enddo
      ibdryc(axis) = bbox_ilowerc(axis)
c
c     Set values along the boundary.
c
      if (axis .eq. 0) then
c
c     Treat x boundaries:
c
         w_fac = 0.d0
         do k = istencil_lower(2),istencil_upper(2)
            do j = istencil_lower(1),istencil_upper(1)
               do i = istencil_lower(0),istencil_upper(0)
                  w = 1.d0 - dabs(dble(i))/dble(ratio_to_coarser(axis))
                  w_fac = w_fac + w
               enddo
            enddo
         enddo
         w_fac = 1.d0/w_fac

         do k_c = bbox_ilowerc(2),bbox_iupperc(2)
            do j_c = bbox_ilowerc(1),bbox_iupperc(1)
               i_c = ibdryc(0)
               i_f = i_c*ratio_to_coarser(0)
               j_f = j_c*ratio_to_coarser(1)
               k_f = k_c*ratio_to_coarser(2)
               U_crse0(i_c,j_c,k_c) = 0.d0
               do k = istencil_lower(2),istencil_upper(2)
                  do j = istencil_lower(1),istencil_upper(1)
                     do i = istencil_lower(0),istencil_upper(0)
                        w = 1.d0 - dabs(dble(i))
     &                     /dble(ratio_to_coarser(axis))
                        U_crse0(i_c,j_c,k_c) = U_crse0(i_c,j_c,k_c) +
     &                     w*w_fac*U_fine0(i_f+i,j_f+j,k_f+k)
                     enddo
                  enddo
               enddo
            enddo
         enddo

      else if (axis .eq. 1) then
c
c     Treat y boundaries:
c
         w_fac = 0.d0
         do k = istencil_lower(2),istencil_upper(2)
            do j = istencil_lower(1),istencil_upper(1)
               do i = istencil_lower(0),istencil_upper(0)
                  w = 1.d0 - dabs(dble(j))/dble(ratio_to_coarser(axis))
                  w_fac = w_fac + w
               enddo
            enddo
         enddo
         w_fac = 1.d0/w_fac

         do k_c = bbox_ilowerc(2),bbox_iupperc(2)
            j_c = ibdryc(1)
            do i_c = bbox_ilowerc(0),bbox_iupperc(0)
               i_f = i_c*ratio_to_coarser(0)
               j_f = j_c*ratio_to_coarser(1)
               k_f = k_c*ratio_to_coarser(2)
               U_crse1(i_c,j_c,k_c) = 0.d0
               do k = istencil_lower(2),istencil_upper(2)
                  do j = istencil_lower(1),istencil_upper(1)
                     do i = istencil_lower(0),istencil_upper(0)
                        w = 1.d0 - dabs(dble(j))
     &                     /dble(ratio_to_coarser(axis))
                        U_crse1(i_c,j_c,k_c) = U_crse1(i_c,j_c,k_c) +
     &                     w*w_fac*U_fine1(i_f+i,j_f+j,k_f+k)
                     enddo
                  enddo
               enddo
            enddo
         enddo

      else if (axis .eq. 2) then
c
c     Treat z boundaries:
c
         w_fac = 0.d0
         do k = istencil_lower(2),istencil_upper(2)
            do j = istencil_lower(1),istencil_upper(1)
               do i = istencil_lower(0),istencil_upper(0)
                  w = 1.d0 - dabs(dble(k))/dble(ratio_to_coarser(axis))
                  w_fac = w_fac + w
               enddo
            enddo
         enddo
         w_fac = 1.d0/w_fac

         k_c = ibdryc(2)
         do j_c = bbox_ilowerc(1),bbox_iupperc(1)
            do i_c = bbox_ilowerc(0),bbox_iupperc(0)
               i_f = i_c*ratio_to_coarser(0)
               j_f = j_c*ratio_to_coarser(1)
               k_f = k_c*ratio_to_coarser(2)
               U_crse2(i_c,j_c,k_c) = 0.d0
               do k = istencil_lower(2),istencil_upper(2)
                  do j = istencil_lower(1),istencil_upper(1)
                     do i = istencil_lower(0),istencil_upper(0)
                        w = 1.d0 - dabs(dble(k))
     &                     /dble(ratio_to_coarser(axis))
                        U_crse2(i_c,j_c,k_c) = U_crse2(i_c,j_c,k_c) +
     &                     w*w_fac*U_fine2(i_f+i,j_f+j,k_f+k)
                     enddo
                  enddo
               enddo
            enddo
         enddo

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
