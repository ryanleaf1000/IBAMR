c ---------------------------------------------------------------------
c
c Copyright (c) 2017 - 2019 by the IBAMR developers
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Mollify indicator function using IB_4 kernel
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine mollify_ib_4_3d(
     &     V,V_gcw,
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2)
c
      implicit none
c ---------------------------------------------------------------------
c
c Copyright (c) 2019 - 2019 by the IBAMR developers
c All rights reserved.
c
c This file is part of IBAMR.
c
c IBAMR is free software and is distributed under the 3-clause BSD
c license. The full text of the license can be found in the file
c COPYRIGHT at the top level directory of IBAMR.
c
c ---------------------------------------------------------------------

      double precision zero,eighth,sixth,fourth,third,half,twothird,
     &  threefourth,fourthird,rt75,one,onept5,two,three,pi,
     &  four,seven,smallr
      parameter (zero=0.d0)
      parameter (eighth=0.125d0)
      parameter (sixth=0.16666666666667d0)
      parameter (fourth=0.25d0)
      parameter (third=0.333333333333333d0)
      parameter (half=0.5d0)
      parameter (twothird=0.66666666666667d0)
      parameter (threefourth=0.75d0)
      parameter (fourthird=1.3333333333333d0)
      parameter (rt75=0.8660254037844d0)
      parameter (one=1.d0)
      parameter (onept5=1.5d0)
      parameter (two=2.d0)
      parameter (three=3.d0)
      parameter (pi=3.14159265358979323846d0)
      parameter (four=4.d0)
      parameter (seven=7.d0)
      parameter (smallr=1.0d-32)
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer ilower2,iupper2
      integer V_gcw,U_gcw
    
c
c     Input/Output.
c
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &          ilower1-U_gcw:iupper1+U_gcw,
     &          ilower2-U_gcw:iupper2+U_gcw)
      double precision V(ilower0-V_gcw:iupper0+V_gcw,
     &          ilower1-V_gcw:iupper1+V_gcw,
     &          ilower2-V_gcw:iupper2+V_gcw)
      double precision w0(-1:1),w1(-1:1),w2(-1:1),wy,wz
      double precision w(-1:1,-1:1,-1:1)

c
c     Local variables.
c
      integer k0,k1,k2
      integer i0,i1,i2
      

c     Compute 1D weights.
      w0(-1) = fourth; w1(-1) = fourth; w2(-1) = fourth
      w0(0) = half; w1(0) = half; w2(0) = half
      w0(1) = fourth; w1(1) = fourth; w2(1) = fourth

c     Compute the tensor product weight
      do k2 = -1,1
         wz = w2(k2)
         do k1 = -1,1
           wy = w1(k1)
           do k0 = -1,1
              w(k0,k1,k2) = w0(k0)*wy*wz
           enddo
        enddo  
      enddo
    
c     Mollify U to V.
      do i2 = ilower2,iupper2
        do i1 = ilower1,iupper1
          do i0 = ilower0,iupper0

              V(i0,i1,i2) = 0.d0
              do k2 = -1,1
                do k1 = -1,1
                  do k0 = -1,1
                    V(i0,i1,i2) = V(i0,i1,i2) + 
     &                   U(i0+k0,i1+k1,i2+k2)*w(k0,k1,k2)
                  enddo
                enddo
              enddo
 
          enddo
        enddo
      enddo
      

      return
      end

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute gradient of the indicator function to estimate 
c     interface normal
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sc_normal_3d(
     &     N00,N01,N02,
     &     N10,N11,N12,
     &     N20,N21,N22,
     &     N_gcw,
     &     U,U_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx)
c
      implicit none
c ---------------------------------------------------------------------
c
c Copyright (c) 2019 - 2019 by the IBAMR developers
c All rights reserved.
c
c This file is part of IBAMR.
c
c IBAMR is free software and is distributed under the 3-clause BSD
c license. The full text of the license can be found in the file
c COPYRIGHT at the top level directory of IBAMR.
c
c ---------------------------------------------------------------------

      double precision zero,eighth,sixth,fourth,third,half,twothird,
     &  threefourth,fourthird,rt75,one,onept5,two,three,pi,
     &  four,seven,smallr
      parameter (zero=0.d0)
      parameter (eighth=0.125d0)
      parameter (sixth=0.16666666666667d0)
      parameter (fourth=0.25d0)
      parameter (third=0.333333333333333d0)
      parameter (half=0.5d0)
      parameter (twothird=0.66666666666667d0)
      parameter (threefourth=0.75d0)
      parameter (fourthird=1.3333333333333d0)
      parameter (rt75=0.8660254037844d0)
      parameter (one=1.d0)
      parameter (onept5=1.5d0)
      parameter (two=2.d0)
      parameter (three=3.d0)
      parameter (pi=3.14159265358979323846d0)
      parameter (four=4.d0)
      parameter (seven=7.d0)
      parameter (smallr=1.0d-32)
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer ilower2,iupper2
      integer N_gcw,U_gcw
    
c
c     Input/Output.
c
      double precision N00(ilower0-N_gcw:iupper0+1+N_gcw,
     &          ilower1-N_gcw:iupper1+N_gcw,
     &          ilower2-N_gcw:iupper2+N_gcw)
      double precision N01(ilower0-N_gcw:iupper0+1+N_gcw,
     &          ilower1-N_gcw:iupper1+N_gcw,
     &          ilower2-N_gcw:iupper2+N_gcw)
      double precision N02(ilower0-N_gcw:iupper0+1+N_gcw,
     &          ilower1-N_gcw:iupper1+N_gcw,
     &          ilower2-N_gcw:iupper2+N_gcw)
      double precision N10(ilower0-N_gcw:iupper0+N_gcw,
     &          ilower1-N_gcw:iupper1+1+N_gcw,
     &          ilower2-N_gcw:iupper2+N_gcw)
      double precision N11(ilower0-N_gcw:iupper0+N_gcw,
     &          ilower1-N_gcw:iupper1+1+N_gcw,
     &          ilower2-N_gcw:iupper2+N_gcw)
      double precision N12(ilower0-N_gcw:iupper0+N_gcw,
     &          ilower1-N_gcw:iupper1+1+N_gcw,
     &          ilower2-N_gcw:iupper2+N_gcw)
      double precision N20(ilower0-N_gcw:iupper0+N_gcw,
     &          ilower1-N_gcw:iupper1+N_gcw,
     &          ilower2-N_gcw:iupper2+1+N_gcw)
      double precision N21(ilower0-N_gcw:iupper0+N_gcw,
     &          ilower1-N_gcw:iupper1+N_gcw,
     &          ilower2-N_gcw:iupper2+1+N_gcw)
      double precision N22(ilower0-N_gcw:iupper0+N_gcw,
     &          ilower1-N_gcw:iupper1+N_gcw,
     &          ilower2-N_gcw:iupper2+1+N_gcw)
      double precision U(ilower0-U_gcw:iupper0+U_gcw,
     &          ilower1-U_gcw:iupper1+U_gcw,
     &          ilower2-U_gcw:iupper2+U_gcw)
     
      double precision dx(0:3-1)

c
c     Local variables.
c
      integer i0,i1,i2
      double precision fac0,fac1,fac2
      
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))
      fac2 = 1.d0/(dx(2))      

c
c     Find face normal gradients first and then interpolate 
c     face tangential gradients

c     Do N00.
      do i2 = ilower2 - 2, iupper2 + 2
        do i1 = ilower1 - 2, iupper1 + 2
          do i0 = ilower0 - 1, iupper0 + 2
             
            N00(i0,i1,i2) = fac0*(U(i0,i1,i2) - U(i0-1,i1,i2))

          enddo
        enddo
      enddo

      
c     Do N11.
      do i2 = ilower2 - 2, iupper2 + 2
        do i1 = ilower1 - 1, iupper1 + 2
          do i0 = ilower0 - 2, iupper0 + 2
             
              N11(i0,i1,i2) = fac1*(U(i0,i1,i2) - U(i0,i1-1,i2))

          enddo
        enddo
      enddo

c     Do N22.
      do i2 = ilower2 - 1, iupper2 + 2
        do i1 = ilower1 - 2, iupper1 + 2
          do i0 = ilower0 - 2, iupper0 + 2
             
              N22(i0,i1,i2) = fac2*(U(i0,i1,i2) - U(i0,i1,i2-1))

          enddo
        enddo
      enddo

c     Interpolate N11 to N01
      do i2 = ilower2 - 1, iupper2 + 1
        do i1 = ilower1 - 1, iupper1 + 1
          do i0 = ilower0 - 1, iupper0 + 2
             
              N01(i0,i1,i2) = fourth*(N11(i0-1,i1,i2) + N11(i0,i1,i2) + 
     &                         N11(i0-1,i1+1,i2) + N11(i0,i1+1,i2)) 
          enddo
        enddo
      enddo

c     Interpolate N22 to N02
      do i2 = ilower2 - 1, iupper2 + 1
        do i1 = ilower1 - 1, iupper1 + 1
          do i0 = ilower0 - 1, iupper0 + 2
             
              N02(i0,i1,i2) = fourth*(N22(i0-1,i1,i2) + N22(i0,i1,i2) + 
     &                         N22(i0-1,i1,i2+1) + N22(i0,i1,i2+1)) 
          enddo
        enddo
      enddo

c     Interpolate N00 to N10
      do i2 = ilower2 - 1, iupper2 + 1
        do i1 = ilower1 - 1, iupper1 + 2
          do i0 = ilower0 - 1, iupper0 + 1

              N10(i0,i1,i2) = fourth*(N00(i0,i1,i2) + N00(i0+1,i1,i2) +
     &                         N00(i0,i1-1,i2) + N00(i0+1,i1-1,i2))

          enddo
        enddo
      enddo

c     Interpolate N22 to N12
      do i2 = ilower2 - 1, iupper2 + 1
        do i1 = ilower1 - 1, iupper1 + 2
          do i0 = ilower0 - 1, iupper0 + 1

              N12(i0,i1,i2) = fourth*(N22(i0,i1,i2) + N22(i0,i1,i2+1) +
     &                         N22(i0,i1-1,i2) + N22(i0,i1-1,i2+1))

          enddo
        enddo
      enddo

c     Interpolate N00 to N20
      do i2 = ilower2 - 1, iupper2 + 2
        do i1 = ilower1 - 1, iupper1 + 1
          do i0 = ilower0 - 1, iupper0 + 1

              N20(i0,i1,i2) = fourth*(N00(i0,i1,i2) + N00(i0+1,i1,i2) +
     &                         N00(i0,i1,i2-1) + N00(i0+1,i1,i2-1))

          enddo
        enddo
      enddo

c     Interpolate N11 to N21
      do i2 = ilower2 - 1, iupper2 + 2
        do i1 = ilower1 - 1, iupper1 + 1
          do i0 = ilower0 - 1, iupper0 + 1
             
              N21(i0,i1,i2) = fourth*(N11(i0,i1,i2-1) + N11(i0,i1,i2) + 
     &                         N11(i0,i1+1,i2-1) + N11(i0,i1+1,i2)) 
          enddo
        enddo
      enddo

      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute cell center curvature of the interface.
c
c           K = - div (n/|n|)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine cc_curvature_3d(
     &     K,K_gcw,
     &     N00,N01,N02,
     &     N10,N11,N12,
     &     N20,N21,N22,
     &     N_gcw,
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
      integer K_gcw,N_gcw
    
c
c     Input/Output.
c
      double precision N00(ilower0-N_gcw:iupper0+1+N_gcw,
     &          ilower1-N_gcw:iupper1+N_gcw,
     &          ilower2-N_gcw:iupper2+N_gcw)
      double precision N01(ilower0-N_gcw:iupper0+1+N_gcw,
     &          ilower1-N_gcw:iupper1+N_gcw,
     &          ilower2-N_gcw:iupper2+N_gcw)
      double precision N02(ilower0-N_gcw:iupper0+1+N_gcw,
     &          ilower1-N_gcw:iupper1+N_gcw,
     &          ilower2-N_gcw:iupper2+N_gcw)
      double precision N10(ilower0-N_gcw:iupper0+N_gcw,
     &          ilower1-N_gcw:iupper1+1+N_gcw,
     &          ilower2-N_gcw:iupper2+N_gcw)
      double precision N11(ilower0-N_gcw:iupper0+N_gcw,
     &          ilower1-N_gcw:iupper1+1+N_gcw,
     &          ilower2-N_gcw:iupper2+N_gcw)
      double precision N12(ilower0-N_gcw:iupper0+N_gcw,
     &          ilower1-N_gcw:iupper1+1+N_gcw,
     &          ilower2-N_gcw:iupper2+N_gcw)
      double precision N20(ilower0-N_gcw:iupper0+N_gcw,
     &          ilower1-N_gcw:iupper1+N_gcw,
     &          ilower2-N_gcw:iupper2+1+N_gcw)
      double precision N21(ilower0-N_gcw:iupper0+N_gcw,
     &          ilower1-N_gcw:iupper1+N_gcw,
     &          ilower2-N_gcw:iupper2+1+N_gcw)
      double precision N22(ilower0-N_gcw:iupper0+N_gcw,
     &          ilower1-N_gcw:iupper1+N_gcw,
     &          ilower2-N_gcw:iupper2+1+N_gcw)
      double precision K(ilower0-K_gcw:iupper0+K_gcw,
     &          ilower1-K_gcw:iupper1+K_gcw,
     &          ilower2-K_gcw:iupper2+K_gcw)
     
      double precision dx(0:3-1)

c
c     Local variables.
c
      integer i0,i1,i2
      double precision fac0,fac1,fac2
      double precision norm_grad_upper,norm_grad_lower,eps
 
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))   
      fac2 = 1.d0/(dx(2))   
      eps = 1.d-10

c
c     Compute curvature K = -div (n/|n|) 
c     
      do i2 = ilower2 - 1, iupper2 + 1
        do i1 = ilower1 - 1, iupper1 + 1
          do i0 = ilower0 - 1, iupper0 + 1
            
c           compute -d/dx term.

            norm_grad_upper = sqrt(N00(i0+1,i1,i2)**2 +
     &                             N01(i0+1,i1,i2)**2 +
     &                             N02(i0+1,i1,i2)**2)
            norm_grad_lower = sqrt(N00(i0,i1,i2)**2 +
     &                             N01(i0,i1,i2)**2 +
     &                             N02(i0,i1,i2)**2)
            
            if (norm_grad_upper .gt. eps) then
                norm_grad_upper = 1.d0/norm_grad_upper
            else
                norm_grad_upper = 0.d0
            endif

            if (norm_grad_lower .gt. eps) then
                norm_grad_lower = 1.d0/norm_grad_lower
            else
                norm_grad_lower = 0.d0
            endif

c           Compute -div
            K(i0,i1,i2) = fac0*(N00(i0,i1,i2)*norm_grad_lower - 
     &              N00(i0+1,i1,i2)*norm_grad_upper)
            

c           compute -d/dy term.

            norm_grad_upper = sqrt(N10(i0,i1+1,i2)**2 +
     &                             N11(i0,i1+1,i2)**2 +
     &                             N12(i0,i1+1,i2)**2)
            norm_grad_lower = sqrt(N10(i0,i1,i2)**2 + 
     &                             N11(i0,i1,i2)**2 +
     &                             N12(i0,i1,i2)**2)

            if (norm_grad_upper .gt. eps) then
                norm_grad_upper = 1.d0/norm_grad_upper
            else
                norm_grad_upper = 0.d0
            endif
            
            if (norm_grad_lower .gt. eps) then
                norm_grad_lower = 1.d0/norm_grad_lower
            else
                norm_grad_lower = 0.d0
            endif    
              
c           Compute -div
            K(i0,i1,i2) = K(i0,i1,i2) + fac1*(N11(i0,i1,i2)
     &         * norm_grad_lower - N11(i0,i1+1,i2)*norm_grad_upper)

c           compute -d/dz term.

            norm_grad_upper = sqrt(N20(i0,i1,i2+1)**2 +
     &                             N21(i0,i1,i2+1)**2 +
     &                             N22(i0,i1,i2+1)**2)
            norm_grad_lower = sqrt(N20(i0,i1,i2)**2 + 
     &                             N21(i0,i1,i2)**2 +
     &                             N22(i0,i1,i2)**2)

            if (norm_grad_upper .gt. eps) then
                norm_grad_upper = 1.d0/norm_grad_upper
            else
                norm_grad_upper = 0.d0
            endif
            
            if (norm_grad_lower .gt. eps) then
                norm_grad_lower = 1.d0/norm_grad_lower
            else
                norm_grad_lower = 0.d0
            endif    
              
c           Compute -div
            K(i0,i1,i2) = K(i0,i1,i2) + fac2*(N22(i0,i1,i2)
     &         * norm_grad_lower - N22(i0,i1,i2+1)*norm_grad_upper)           

          enddo
        enddo
      enddo
      
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Compute surface tension forcing
c
c           F = sigma * K * grad(f)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine sc_surface_tension_force_3d(
     &     F0,F1,F2,F_gcw,
     &     K,K_gcw,
     &     N00,N11,N22,N_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     sigma)
c
      implicit none
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer ilower2,iupper2
      integer F_gcw,K_gcw,N_gcw
    
c
c     Input/Output.
c
      double precision F0(ilower0-F_gcw:iupper0+1+F_gcw,
     &          ilower1-F_gcw:iupper1+F_gcw,
     &          ilower2-F_gcw:iupper2+F_gcw)
      double precision F1(ilower0-F_gcw:iupper0+F_gcw,
     &          ilower1-F_gcw:iupper1+1+F_gcw,
     &          ilower2-F_gcw:iupper2+F_gcw)
      double precision F2(ilower0-F_gcw:iupper0+F_gcw,
     &          ilower1-F_gcw:iupper1+F_gcw,
     &          ilower2-F_gcw:iupper2+1+F_gcw)

      double precision N00(ilower0-N_gcw:iupper0+1+N_gcw,
     &          ilower1-N_gcw:iupper1+N_gcw,
     &          ilower2-N_gcw:iupper2+N_gcw)
      double precision N11(ilower0-N_gcw:iupper0+N_gcw,
     &          ilower1-N_gcw:iupper1+1+N_gcw,
     &          ilower2-N_gcw:iupper2+N_gcw)
      double precision N22(ilower0-N_gcw:iupper0+N_gcw,
     &          ilower1-N_gcw:iupper1+N_gcw,
     &          ilower2-N_gcw:iupper2+1+N_gcw)
      double precision K(ilower0-K_gcw:iupper0+K_gcw,
     &          ilower1-K_gcw:iupper1+K_gcw,
     &          ilower2-K_gcw:iupper2+K_gcw)
     
      double precision sigma

c
c     Local variables.
c
      integer i0,i1,i2
      double precision kappa
  
c
c     Compute F0  = sigma * K_x * N_x 
c     
      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0 + 1
            
              kappa = 0.5d0*(K(i0-1,i1,i2)+K(i0,i1,i2))
              F0(i0,i1,i2) = sigma*kappa*N00(i0,i1,i2)
            
          enddo
        enddo
      enddo

c
c     Compute F1  = sigma * K_y * N_y
c
      do i2 = ilower2, iupper2
        do i1 = ilower1, iupper1 + 1
          do i0 = ilower0, iupper0

              kappa = 0.5d0*(K(i0,i1-1,i2)+K(i0,i1,i2))
              F1(i0,i1,i2) = sigma*kappa*N11(i0,i1,i2)
                    
          enddo
        enddo
      enddo

c
c     Compute F2  = sigma * K_z * N_z
c
      do i2 = ilower2, iupper2 + 1
        do i1 = ilower1, iupper1
          do i0 = ilower0, iupper0

              kappa = 0.5d0*(K(i0,i1,i2-1)+K(i0,i1,i2))
              F2(i0,i1,i2) = sigma*kappa*N22(i0,i1,i2)
                    
          enddo
        enddo
      enddo
      
      return
      end




