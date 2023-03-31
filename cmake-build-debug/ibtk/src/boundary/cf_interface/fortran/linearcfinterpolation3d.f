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
c     Compute a cell-centered linear interpolation at coarse-fine
c     interfaces.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cclinearnormalinterpolation3d(
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     location_index,ratio,
     &     blower,bupper)
c
      implicit none
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer ilower2,iupper2
      integer U_gcw

      integer location_index,ratio

      integer blower(0:3-1), bupper(0:3-1)
c
c     Input/Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &          ilower1-U_gcw:iupper1+U_gcw,
     &          ilower2-U_gcw:iupper2+U_gcw)
c
c     Local variables.
c
      integer i,ibeg,iend,i_p,i_q
      integer j,jbeg,jend,j_p,j_q
      integer k,kbeg,kend,k_p,k_q
      double precision R
c
c     Set the values along the appropriate side.
c
      ibeg = max(blower(0),ilower0)
      iend = min(bupper(0),iupper0)

      jbeg = max(blower(1),ilower1)
      jend = min(bupper(1),iupper1)

      kbeg = max(blower(2),ilower2)
      kend = min(bupper(2),iupper2)

      R = dble(ratio)

      if ( (location_index .eq. 0) .or.
     &     (location_index .eq. 1) ) then
c
c     Set the values along the upper/lower x side of the patch.
c
         if (location_index .eq. 0) then
            do j = jbeg,jend,ratio
               do j_p = j,j+ratio-1
                  j_q = j+(ratio-1+j-j_p)
                  do k = kbeg,kend,ratio
                     do k_p = k,k+ratio-1
                        k_q = k+(ratio-1+k-k_p)
                        U(ilower0-1,j_p,k_p) = (
     &                       2.d0*U(ilower0-1,j_p,k_p)
     &                       +  R*U(ilower0  ,j_p,k_p)
     &                       -    U(ilower0  ,j_q,k_q))/(R+1.d0)
                     enddo
                  enddo
               enddo
            enddo
         else
            do j = jbeg,jend,ratio
               do j_p = j,j+ratio-1
                  j_q = j+(ratio-1+j-j_p)
                  do k = kbeg,kend,ratio
                     do k_p = k,k+ratio-1
                        k_q = k+(ratio-1+k-k_p)
                        U(iupper0+1,j_p,k_p) = (
     &                       2.d0*U(iupper0+1,j_p,k_p)
     &                       +  R*U(iupper0  ,j_p,k_p)
     &                       -    U(iupper0  ,j_q,k_q))/(R+1.d0)
                     enddo
                  enddo
               enddo
            enddo
         endif

      elseif ( (location_index .eq. 2) .or.
     &         (location_index .eq. 3) ) then
c
c     Set the values along the upper/lower y side of the patch.
c
         if (location_index .eq. 2) then
            do i = ibeg,iend,ratio
               do i_p = i,i+ratio-1
                  i_q = i+(ratio-1+i-i_p)
                  do k = kbeg,kend,ratio
                     do k_p = k,k+ratio-1
                        k_q = k+(ratio-1+k-k_p)
                        U(i_p,ilower1-1,k_p) = (
     &                       2.d0*U(i_p,ilower1-1,k_p)
     &                       +  R*U(i_p,ilower1  ,k_p)
     &                       -    U(i_q,ilower1  ,k_q))/(R+1.d0)
                     enddo
                  enddo
               enddo
            enddo
         else
            do i = ibeg,iend,ratio
               do i_p = i,i+ratio-1
                  i_q = i+(ratio-1+i-i_p)
                  do k = kbeg,kend,ratio
                     do k_p = k,k+ratio-1
                        k_q = k+(ratio-1+k-k_p)
                        U(i_p,iupper1+1,k_p) = (
     &                       2.d0*U(i_p,iupper1+1,k_p)
     &                       +  R*U(i_p,iupper1  ,k_p)
     &                       -    U(i_q,iupper1  ,k_q))/(R+1.d0)
                     enddo
                  enddo
               enddo
            enddo
         endif

      elseif ( (location_index .eq. 4) .or.
     &         (location_index .eq. 5) ) then
c
c     Set the values along the upper/lower z side of the patch.
c
         if (location_index .eq. 4) then
            do i = ibeg,iend,ratio
               do i_p = i,i+ratio-1
                  i_q = i+(ratio-1+i-i_p)
                  do j = jbeg,jend,ratio
                     do j_p = j,j+ratio-1
                        j_q = j+(ratio-1+j-j_p)
                        U(i_p,j_p,ilower2-1) = (
     &                       2.d0*U(i_p,j_p,ilower2-1)
     &                       +  R*U(i_p,j_p,ilower2  )
     &                       -    U(i_q,j_q,ilower2  ))/(R+1.d0)
                     enddo
                  enddo
               enddo
            enddo
         else
            do i = ibeg,iend,ratio
               do i_p = i,i+ratio-1
                  i_q = i+(ratio-1+i-i_p)
                  do j = jbeg,jend,ratio
                     do j_p = j,j+ratio-1
                        j_q = j+(ratio-1+j-j_p)
                        U(i_p,j_p,iupper2+1) = (
     &                       2.d0*U(i_p,j_p,iupper2+1)
     &                       +  R*U(i_p,j_p,iupper2  )
     &                       -    U(i_q,j_q,iupper2  ))/(R+1.d0)
                     enddo
                  enddo
               enddo
            enddo
         endif

      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
