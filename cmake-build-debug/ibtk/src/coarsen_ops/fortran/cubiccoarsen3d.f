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
c     Coarsen cell-centered data via cubic interpolation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cccubiccoarsen3d(
     &     U_crse,U_crse_gcw,
     &     U_fine,U_fine_gcw,
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
c     Constants.
c
      double precision wgts(-2:1)
      DATA wgts/-6.25d-2,5.625d-1,5.625d-1,-6.25d-2/
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

      double precision U_fine(ilowerf0-U_fine_gcw:iupperf0+U_fine_gcw,
     &          ilowerf1-U_fine_gcw:iupperf1+U_fine_gcw,
     &          ilowerf2-U_fine_gcw:iupperf2+U_fine_gcw)
c
c     Input/Output.
c
      double precision U_crse(ilowerc0-U_crse_gcw:iupperc0+U_crse_gcw,
     &          ilowerc1-U_crse_gcw:iupperc1+U_crse_gcw,
     &          ilowerc2-U_crse_gcw:iupperc2+U_crse_gcw)
c
c     Local variables.
c
      integer i,i_f,i_c
      integer j,j_f,j_c
      integer k,k_f,k_c
c
c     Coarsen the fine data via cubic interpolation.
c
      do k_c = coarse_box_lower(2),coarse_box_upper(2)
         do j_c = coarse_box_lower(1),coarse_box_upper(1)
            do i_c = coarse_box_lower(0),coarse_box_upper(0)
               U_crse(i_c,j_c,k_c) = 0.d0
            enddo
         enddo
      enddo
      do k_c = coarse_box_lower(2),coarse_box_upper(2)
         k_f = k_c*ratio_to_coarser(2)
         do k = -2,1
            do j_c = coarse_box_lower(1),coarse_box_upper(1)
               j_f = j_c*ratio_to_coarser(1)
               do j = -2,1
                  do i_c = coarse_box_lower(0),coarse_box_upper(0)
                     i_f = i_c*ratio_to_coarser(0)
                     do i = -2,1
                        U_crse(i_c,j_c,k_c) = U_crse(i_c,j_c,k_c) +
     &                       wgts(i)*wgts(j)*wgts(k)*U_fine(
     &                       i_f+ratio_to_coarser(0)/2+i,
     &                       j_f+ratio_to_coarser(1)/2+j,
     &                       k_f+ratio_to_coarser(2)/2+k)
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
c     Coarsen side-centered data via cubic interpolation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sccubiccoarsen3d(
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
c     Constants.
c
      double precision wgts(-2:1)
      DATA wgts/-6.25d-2,5.625d-1,5.625d-1,-6.25d-2/
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
      integer i,j,k,i_f,i_c,j_f,j_c,k_f,k_c
c
c     Coarsen the fine data via cubic interpolation.
c
      do k_c = coarse_box_lower(2),coarse_box_upper(2)
         do j_c = coarse_box_lower(1),coarse_box_upper(1)
            do i_c = coarse_box_lower(0),coarse_box_upper(0)+1
               U_crse0(i_c,j_c,k_c) = 0.d0
            enddo
         enddo
      enddo
      do k_c = coarse_box_lower(2),coarse_box_upper(2)
         k_f = k_c*ratio_to_coarser(2)
         do k = -2,1
            do j_c = coarse_box_lower(1),coarse_box_upper(1)
               j_f = j_c*ratio_to_coarser(1)
               do j = -2,1
                  do i_c = coarse_box_lower(0),coarse_box_upper(0)+1
                     i_f = i_c*ratio_to_coarser(0)
                     U_crse0(i_c,j_c,k_c) = U_crse0(i_c,j_c,k_c) +
     &                    wgts(j)*wgts(k)*U_fine0(
     &                    i_f,
     &                    j+j_f+ratio_to_coarser(1)/2,
     &                    k+k_f+ratio_to_coarser(2)/2)
                  enddo
               enddo
            enddo
         enddo
      enddo
      do k_c = coarse_box_lower(2),coarse_box_upper(2)
         do j_c = coarse_box_lower(1),coarse_box_upper(1)+1
            do i_c = coarse_box_lower(0),coarse_box_upper(0)
               U_crse1(i_c,j_c,k_c) = 0.d0
            enddo
         enddo
      enddo
      do k_c = coarse_box_lower(2),coarse_box_upper(2)
         k_f = k_c*ratio_to_coarser(2)
         do k = -2,1
            do j_c = coarse_box_lower(1),coarse_box_upper(1)+1
               j_f = j_c*ratio_to_coarser(1)
               do i_c = coarse_box_lower(0),coarse_box_upper(0)
                  i_f = i_c*ratio_to_coarser(0)
                  do i = -2,1
                     U_crse1(i_c,j_c,k_c) = U_crse1(i_c,j_c,k_c) +
     &                    wgts(i)*wgts(k)*U_fine1(
     &                    i+i_f+ratio_to_coarser(0)/2,
     &                    j_f,
     &                    k+k_f+ratio_to_coarser(2)/2)
                  enddo
               enddo
            enddo
         enddo
      enddo
      do k_c = coarse_box_lower(2),coarse_box_upper(2)+1
         do j_c = coarse_box_lower(1),coarse_box_upper(1)
            do i_c = coarse_box_lower(0),coarse_box_upper(0)
               U_crse2(i_c,j_c,k_c) = 0.d0
               enddo
            enddo
         enddo
      do k_c = coarse_box_lower(2),coarse_box_upper(2)+1
         k_f = k_c*ratio_to_coarser(2)
         do j_c = coarse_box_lower(1),coarse_box_upper(1)
            j_f = j_c*ratio_to_coarser(1)
            do j = -2,1
               do i_c = coarse_box_lower(0),coarse_box_upper(0)
                  i_f = i_c*ratio_to_coarser(0)
                  do i = -2,1
                     U_crse2(i_c,j_c,k_c) = U_crse2(i_c,j_c,k_c) +
     &                    wgts(i)*wgts(j)*U_fine2(
     &                    i+i_f+ratio_to_coarser(0)/2,
     &                    j+j_f+ratio_to_coarser(1)/2,
     &                    k_f)
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
