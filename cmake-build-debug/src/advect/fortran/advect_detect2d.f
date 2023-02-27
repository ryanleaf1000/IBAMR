c ---------------------------------------------------------------------
c
c Copyright (c) 2006 - 2019 by the IBAMR developers
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
c     Detect sharp spatial gradients.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine advect_detectgrad2d(
     &  ifirst0,ilast0,ifirst1,ilast1,
     &  vghost0,tagghost0,ttagghost0,
     &  vghost1,tagghost1,ttagghost1,
     &  dx,
     &  gradtol,
     &  dotag,
     &  var,
     &  tags,temptags)
c
      implicit none
c
      integer
     &  ifirst0,ifirst1,ilast0,ilast1,
     &  dotag,
     &  vghost0,vghost1,
     &  tagghost0,tagghost1,
     &  ttagghost0,ttagghost1
      double precision
     &  dx(0:2-1),
     &  gradtol
      double precision
     &  var(ifirst0-vghost0:ilast0+vghost0,
     &          ifirst1-vghost1:ilast1+vghost1)
      integer
     &  tags(ifirst0-tagghost0:ilast0+tagghost0,
     &          ifirst1-tagghost1:ilast1+tagghost1),
     &  temptags(ifirst0-ttagghost0:ilast0+ttagghost0,
     &          ifirst1-ttagghost1:ilast1+ttagghost1)
c
      double precision tol
      double precision facejump, loctol
      double precision presm1,presp1,diag01
      logical tagcell
      integer ic0,ic1
c
      tol = gradtol
      diag01 = sqrt(dx(0)**2+dx(1)**2)

      do ic1=ifirst1,ilast1
        do ic0=ifirst0,ilast0

          if (tags(ic0,ic1) .ne. 0) then
            loctol = 0.125d0*tol
          else
            loctol = tol
          endif

          tagcell = .false.

          presm1 = var(ic0-1,ic1)
          presp1 = var(ic0+1,ic1)
          facejump = abs(var(ic0,ic1)-presm1)
          facejump = max(facejump,abs(var(ic0,ic1)-presp1))
          tagcell = ((facejump).gt.(loctol*dx(0)))
          if (.not.tagcell) then
            presm1 = var(ic0,ic1-1)
            presp1 = var(ic0,ic1+1)
            facejump = abs(var(ic0,ic1)-presm1)
            facejump = max(facejump,abs(var(ic0,ic1)-presp1))
            tagcell = ((facejump).gt.(loctol*dx(1)))
          endif

          if (.not.tagcell) then
            presm1 = var(ic0-1,ic1-1)
            presp1 = var(ic0+1,ic1+1)
            facejump = abs(var(ic0,ic1)-presm1)
            facejump = max(facejump,abs(var(ic0,ic1)-presp1))
            tagcell = ((facejump).gt.(loctol*diag01))
          endif
          if (.not.tagcell) then
            presm1 = var(ic0-1,ic1+1)
            presp1 = var(ic0+1,ic1-1)
            facejump = abs(var(ic0,ic1)-presm1)
            facejump = max(facejump,abs(var(ic0,ic1)-presp1))
            tagcell = ((facejump).gt.(loctol*diag01))
          endif

          if ( tagcell ) then
            temptags(ic0,ic1) = dotag
          endif
        enddo
      enddo
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
