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
c     Coarsen cell-centered data via cubic interpolation.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cccubiccoarsen2d(
     &     U_crse,U_crse_gcw,
     &     U_fine,U_fine_gcw,
     &     ilowerc0,iupperc0,
     &     ilowerc1,iupperc1,
     &     ilowerf0,iupperf0,
     &     ilowerf1,iupperf1,
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
      integer ilowerf0,iupperf0
      integer ilowerf1,iupperf1
      integer U_crse_gcw,U_fine_gcw

      integer ratio_to_coarser(0:2-1)

      integer coarse_box_lower(0:2-1), coarse_box_upper(0:2-1)

      double precision U_fine(ilowerf0-U_fine_gcw:iupperf0+U_fine_gcw,
     &          ilowerf1-U_fine_gcw:iupperf1+U_fine_gcw)
c
c     Input/Output.
c
      double precision U_crse(ilowerc0-U_crse_gcw:iupperc0+U_crse_gcw,
     &          ilowerc1-U_crse_gcw:iupperc1+U_crse_gcw)
c
c     Local variables.
c
      integer i,i_f,i_c
      integer j,j_f,j_c
c
c     Coarsen the fine data via cubic interpolation.
c
      do j_c = coarse_box_lower(1),coarse_box_upper(1)
         do i_c = coarse_box_lower(0),coarse_box_upper(0)
            U_crse(i_c,j_c) = 0.d0
         enddo
      enddo
      do j_c = coarse_box_lower(1),coarse_box_upper(1)
         j_f = j_c*ratio_to_coarser(1)
         do j = -2,1
            do i_c = coarse_box_lower(0),coarse_box_upper(0)
               i_f = i_c*ratio_to_coarser(0)
               do i = -2,1
                  U_crse(i_c,j_c) = U_crse(i_c,j_c) +
     &                 wgts(i)*wgts(j)*U_fine(
     &                 i_f+ratio_to_coarser(0)/2+i,
     &                 j_f+ratio_to_coarser(1)/2+j)
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
      subroutine sccubiccoarsen2d(
     &     U_crse0,U_crse1,U_crse_gcw,
     &     U_fine0,U_fine1,U_fine_gcw,
     &     ilowerc0,iupperc0,
     &     ilowerc1,iupperc1,
     &     ilowerf0,iupperf0,
     &     ilowerf1,iupperf1,
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
      integer ilowerf0,iupperf0
      integer ilowerf1,iupperf1
      integer U_crse_gcw,U_fine_gcw

      integer ratio_to_coarser(0:2-1)

      integer coarse_box_lower(0:2-1), coarse_box_upper(0:2-1)

      double precision U_fine0(
     &     ilowerf0-U_fine_gcw:iupperf0+1+U_fine_gcw,
     &          ilowerf1-U_fine_gcw:iupperf1+U_fine_gcw
     &     )
      double precision U_fine1(
     &     ilowerf0-U_fine_gcw:iupperf0+U_fine_gcw,
     &          ilowerf1-U_fine_gcw:iupperf1+1+U_fine_gcw
     &     )
c
c     Input/Output.
c
      double precision U_crse0(
     &     ilowerc0-U_crse_gcw:iupperc0+1+U_crse_gcw,
     &          ilowerc1-U_crse_gcw:iupperc1+U_crse_gcw
     &     )
      double precision U_crse1(
     &     ilowerc0-U_crse_gcw:iupperc0+U_crse_gcw,
     &          ilowerc1-U_crse_gcw:iupperc1+1+U_crse_gcw
     &     )
c
c     Local variables.
c
      integer i,j,i_f,i_c,j_f,j_c
c
c     Coarsen the fine data via cubic interpolation.
c
      do j_c = coarse_box_lower(1),coarse_box_upper(1)
         do i_c = coarse_box_lower(0),coarse_box_upper(0)+1
            U_crse0(i_c,j_c) = 0.d0
         enddo
      enddo
      do j_c = coarse_box_lower(1),coarse_box_upper(1)
         j_f = j_c*ratio_to_coarser(1)
         do j = -2,1
            do i_c = coarse_box_lower(0),coarse_box_upper(0)+1
               i_f = i_c*ratio_to_coarser(0)
               U_crse0(i_c,j_c) = U_crse0(i_c,j_c) +
     &              wgts(j)*U_fine0(i_f,j+j_f+ratio_to_coarser(1)/2)
            enddo
         enddo
      enddo
      do j_c = coarse_box_lower(1),coarse_box_upper(1)+1
         do i_c = coarse_box_lower(0),coarse_box_upper(0)
            U_crse1(i_c,j_c) = 0.d0
         enddo
      enddo
      do j_c = coarse_box_lower(1),coarse_box_upper(1)+1
         j_f = j_c*ratio_to_coarser(1)
         do i_c = coarse_box_lower(0),coarse_box_upper(0)
            i_f = i_c*ratio_to_coarser(0)
            do i = -2,1
               U_crse1(i_c,j_c) = U_crse1(i_c,j_c) +
     &              wgts(i)*U_fine1(i+i_f+ratio_to_coarser(0)/2,j_f)
            enddo
         enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
