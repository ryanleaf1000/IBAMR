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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Computes (f0,f1,f2) = alpha div mu (grad (u0,u1,u1) + grad(u0,u1,u2)^T) 
c     + beta (u0,u1,f2) + gamma (v0,v1,v2)
c
c     Computes the side-centered variable coefficient generalized
c     Laplacian, with edge-centered coefficient mu and side-centered
c     vector fields (u0,u1,u2) and (v0,v1,v2).
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stosvclaplace3d(
     &     f0,f1,f2,f_gcw,
     &     alpha,beta,
     &     mu0,mu1,mu2,mu_gcw,
     &     rho0,rho1,rho2,rho_gcw,
     &     u0,u1,u2,u_gcw,
     &     gamma,
     &     v0,v1,v2,v_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     ilower2,iupper2,
     &     dx,
     &     var_rho,
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
      integer f_gcw,mu_gcw,rho_gcw,u_gcw,v_gcw,var_rho
      integer use_harmonic_interp

      double precision alpha,beta,gamma

      double precision mu0(ilower0-mu_gcw:iupper0+mu_gcw,
     &          ilower1-mu_gcw:iupper1+1+mu_gcw,
     &          ilower2-mu_gcw:iupper2+1+mu_gcw)
      double precision mu1(ilower0-mu_gcw:iupper0+1+mu_gcw,
     &          ilower1-mu_gcw:iupper1+mu_gcw,
     &          ilower2-mu_gcw:iupper2+1+mu_gcw)
      double precision mu2(ilower0-mu_gcw:iupper0+1+mu_gcw,
     &          ilower1-mu_gcw:iupper1+1+mu_gcw,
     &          ilower2-mu_gcw:iupper2+mu_gcw)

      double precision rho0(ilower0-rho_gcw:iupper0+1+rho_gcw,
     &          ilower1-rho_gcw:iupper1+rho_gcw,
     &          ilower2-rho_gcw:iupper2+rho_gcw)
      double precision rho1(ilower0-rho_gcw:iupper0+rho_gcw,
     &          ilower1-rho_gcw:iupper1+1+rho_gcw,
     &          ilower2-rho_gcw:iupper2+rho_gcw)
      double precision rho2(ilower0-rho_gcw:iupper0+rho_gcw,
     &          ilower1-rho_gcw:iupper1+rho_gcw,
     &          ilower2-rho_gcw:iupper2+1+rho_gcw)

      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw,
     &          ilower2-u_gcw:iupper2+u_gcw)
      double precision u1(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+1+u_gcw,
     &          ilower2-u_gcw:iupper2+u_gcw)
      double precision u2(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw,
     &          ilower2-u_gcw:iupper2+1+u_gcw)

      double precision v0(ilower0-v_gcw:iupper0+1+v_gcw,
     &          ilower1-v_gcw:iupper1+v_gcw,
     &          ilower2-v_gcw:iupper2+v_gcw)
      double precision v1(ilower0-v_gcw:iupper0+v_gcw,
     &          ilower1-v_gcw:iupper1+1+v_gcw,
     &          ilower2-v_gcw:iupper2+v_gcw)
      double precision v2(ilower0-v_gcw:iupper0+v_gcw,
     &          ilower1-v_gcw:iupper1+v_gcw,
     &          ilower2-v_gcw:iupper2+1+v_gcw)
      
      double precision dx(0:3-1)
c
c     Input/Output.
c
      double precision f0(ilower0-f_gcw:iupper0+1+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f1(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+1+f_gcw,
     &          ilower2-f_gcw:iupper2+f_gcw)
      double precision f2(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw,
     &          ilower2-f_gcw:iupper2+1+f_gcw)
c
c     Local variables.
c
      integer i0,i1,i2
      double precision fac0,fac1,fac2,mu_lower,mu_upper,rho
c
c     Compute the discrete divergence of mu (grad (u0,u1,u2) + grad (u0,u1,u2)^T).
c
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))
      fac2 = 1.d0/(dx(2))

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0+1

               rho = beta
               if (var_rho .eq. 1) then
                  rho = rho0(i0,i1,i2)*beta
               endif

               if (use_harmonic_interp .eq. 1) then
                   mu_upper = h_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                              mu0(i0,i1,i2+1),mu0(i0,i1+1,i2+1),
     &                              mu1(i0,i1,i2),mu1(i0+1,i1,i2),
     &                              mu1(i0,i1,i2+1),mu1(i0+1,i1,i2+1),
     &                              mu2(i0,i1,i2),mu2(i0+1,i1,i2),
     &                              mu2(i0,i1+1,i2),mu2(i0+1,i1+1,i2))

                   mu_lower = h_avg12(mu0(i0-1,i1,i2),mu0(i0-1,i1+1,i2),
     &                            mu0(i0-1,i1,i2+1),mu0(i0-1,i1+1,i2+1),
     &                            mu1(i0-1,i1,i2),mu1(i0,i1,i2),
     &                            mu1(i0-1,i1,i2+1),mu1(i0,i1,i2+1),
     &                            mu2(i0-1,i1,i2),mu2(i0,i1,i2),
     &                            mu2(i0-1,i1+1,i2),mu2(i0,i1+1,i2))
               else
                    mu_upper = a_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                              mu0(i0,i1,i2+1),mu0(i0,i1+1,i2+1),
     &                              mu1(i0,i1,i2),mu1(i0+1,i1,i2),
     &                              mu1(i0,i1,i2+1),mu1(i0+1,i1,i2+1),
     &                              mu2(i0,i1,i2),mu2(i0+1,i1,i2),
     &                              mu2(i0,i1+1,i2),mu2(i0+1,i1+1,i2))

                    mu_lower = a_avg12(mu0(i0-1,i1,i2),mu0(i0-1,i1+1,i2)
     &                           ,mu0(i0-1,i1,i2+1),mu0(i0-1,i1+1,i2+1),
     &                            mu1(i0-1,i1,i2),mu1(i0,i1,i2),
     &                            mu1(i0-1,i1,i2+1),mu1(i0,i1,i2+1),
     &                            mu2(i0-1,i1,i2),mu2(i0,i1,i2),
     &                            mu2(i0-1,i1+1,i2),mu2(i0,i1+1,i2))
               endif

               f0(i0,i1,i2) = alpha*(
     &              2.d0*fac0**2.d0*(
     &              mu_upper*(u0(i0+1,i1,i2)-u0(i0,i1,i2))-
     &              mu_lower*(u0(i0,i1,i2)-u0(i0-1,i1,i2)))+
     &              fac1**2.d0*(
     &              mu2(i0,i1+1,i2)*(u0(i0,i1+1,i2)-u0(i0,i1,i2))-
     &              mu2(i0,i1,i2)*(u0(i0,i1,i2)-u0(i0,i1-1,i2)))+
     &              fac0*fac1*(
     &              mu2(i0,i1+1,i2)*(u1(i0,i1+1,i2)-u1(i0-1,i1+1,i2))-
     &              mu2(i0,i1,i2)*(u1(i0,i1,i2)-u1(i0-1,i1,i2)))+
     &              fac2**2.d0*(
     &              mu1(i0,i1,i2+1)*(u0(i0,i1,i2+1)-u0(i0,i1,i2))-
     &              mu1(i0,i1,i2)*(u0(i0,i1,i2)-u0(i0,i1,i2-1)))+
     &              fac0*fac2*(
     &              mu1(i0,i1,i2+1)*(u2(i0,i1,i2+1)-u2(i0-1,i1,i2+1))-
     &              mu1(i0,i1,i2)*(u2(i0,i1,i2)-u2(i0-1,i1,i2))))+
     &              rho*u0(i0,i1,i2) + gamma*v0(i0,i1,i2)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2
         do i1 = ilower1,iupper1+1
            do i0 = ilower0,iupper0

               rho = beta
               if (var_rho .eq. 1) then
                  rho = rho1(i0,i1,i2)*beta
               endif

               if (use_harmonic_interp .eq. 1) then
                  mu_upper = h_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),mu0(i0,i1+1,i2+1),
     &                             mu1(i0,i1,i2),mu1(i0+1,i1,i2),
     &                             mu1(i0,i1,i2+1),mu1(i0+1,i1,i2+1),
     &                             mu2(i0,i1,i2),mu2(i0+1,i1,i2),
     &                             mu2(i0,i1+1,i2),mu2(i0+1,i1+1,i2))
                  mu_lower = h_avg12(mu0(i0,i1-1,i2),mu0(i0,i1,i2),
     &                             mu0(i0,i1-1,i2+1),mu0(i0,i1,i2+1),
     &                             mu1(i0,i1-1,i2),mu1(i0+1,i1-1,i2),
     &                            mu1(i0,i1-1,i2+1),mu1(i0+1,i1-1,i2+1),
     &                            mu2(i0,i1-1,i2),mu2(i0+1,i1-1,i2),
     &                            mu2(i0,i1,i2),mu2(i0+1,i1,i2))
               else
                  mu_upper = a_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),mu0(i0,i1+1,i2+1),
     &                             mu1(i0,i1,i2),mu1(i0+1,i1,i2),
     &                             mu1(i0,i1,i2+1),mu1(i0+1,i1,i2+1),
     &                             mu2(i0,i1,i2),mu2(i0+1,i1,i2),
     &                             mu2(i0,i1+1,i2),mu2(i0+1,i1+1,i2))
                  mu_lower = a_avg12(mu0(i0,i1-1,i2),mu0(i0,i1,i2),
     &                             mu0(i0,i1-1,i2+1),mu0(i0,i1,i2+1),
     &                             mu1(i0,i1-1,i2),mu1(i0+1,i1-1,i2),
     &                            mu1(i0,i1-1,i2+1),mu1(i0+1,i1-1,i2+1),
     &                            mu2(i0,i1-1,i2),mu2(i0+1,i1-1,i2),
     &                            mu2(i0,i1,i2),mu2(i0+1,i1,i2))
                endif       

               f1(i0,i1,i2) = alpha*(
     &              2.d0*fac1**2.d0*(
     &              mu_upper*(u1(i0,i1+1,i2)-u1(i0,i1,i2))-
     &              mu_lower*(u1(i0,i1,i2)-u1(i0,i1-1,i2)))+
     &              fac0**2.d0*(
     &              mu2(i0+1,i1,i2)*(u1(i0+1,i1,i2)-u1(i0,i1,i2))-
     &              mu2(i0,i1,i2)*(u1(i0,i1,i2)-u1(i0-1,i1,i2)))+
     &              fac0*fac1*(
     &              mu2(i0+1,i1,i2)*(u0(i0+1,i1,i2)-u0(i0+1,i1-1,i2))-
     &              mu2(i0,i1,i2)*(u0(i0,i1,i2)-u0(i0,i1-1,i2)))+
     &              fac2**2.d0*(
     &              mu0(i0,i1,i2+1)*(u1(i0,i1,i2+1)-u1(i0,i1,i2))-
     &              mu0(i0,i1,i2)*(u1(i0,i1,i2)-u1(i0,i1,i2-1)))+
     &              fac1*fac2*(
     &              mu0(i0,i1,i2+1)*(u2(i0,i1,i2+1)-u2(i0,i1-1,i2+1))-
     &              mu0(i0,i1,i2)*(u2(i0,i1,i2)-u2(i0,i1-1,i2))))+
     &              rho*u1(i0,i1,i2) + gamma*v1(i0,i1,i2)
            enddo
         enddo
      enddo

      do i2 = ilower2,iupper2+1
         do i1 = ilower1,iupper1
            do i0 = ilower0,iupper0

               rho = beta
               if (var_rho .eq. 1) then
                  rho = rho2(i0,i1,i2)*beta
               endif

               if (use_harmonic_interp .eq. 1) then
                  mu_upper = h_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),mu0(i0,i1+1,i2+1),
     &                             mu1(i0,i1,i2),mu1(i0+1,i1,i2),
     &                             mu1(i0,i1,i2+1),mu1(i0+1,i1,i2+1),
     &                             mu2(i0,i1,i2),mu2(i0+1,i1,i2),
     &                             mu2(i0,i1+1,i2),mu2(i0+1,i1+1,i2))

                  mu_lower = h_avg12(mu0(i0,i1,i2-1),mu0(i0,i1+1,i2-1),
     &                             mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu1(i0,i1,i2-1),mu1(i0+1,i1,i2-1),
     &                             mu1(i0,i1,i2),mu1(i0+1,i1,i2),
     &                             mu2(i0,i1,i2-1),mu2(i0+1,i1,i2-1),
     &                            mu2(i0,i1+1,i2-1),mu2(i0+1,i1+1,i2-1))
               else
                  mu_upper = a_avg12(mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu0(i0,i1,i2+1),mu0(i0,i1+1,i2+1),
     &                             mu1(i0,i1,i2),mu1(i0+1,i1,i2),
     &                             mu1(i0,i1,i2+1),mu1(i0+1,i1,i2+1),
     &                             mu2(i0,i1,i2),mu2(i0+1,i1,i2),
     &                             mu2(i0,i1+1,i2),mu2(i0+1,i1+1,i2))

                  mu_lower = a_avg12(mu0(i0,i1,i2-1),mu0(i0,i1+1,i2-1),
     &                             mu0(i0,i1,i2),mu0(i0,i1+1,i2),
     &                             mu1(i0,i1,i2-1),mu1(i0+1,i1,i2-1),
     &                             mu1(i0,i1,i2),mu1(i0+1,i1,i2),
     &                             mu2(i0,i1,i2-1),mu2(i0+1,i1,i2-1),
     &                            mu2(i0,i1+1,i2-1),mu2(i0+1,i1+1,i2-1))
               endif

               f2(i0,i1,i2) = alpha*(
     &              2.d0*fac2**2.d0*(
     &              mu_upper*(u2(i0,i1,i2+1)-u2(i0,i1,i2))-
     &              mu_lower*(u2(i0,i1,i2)-u2(i0,i1,i2-1)))+
     &              fac1**2.d0*(
     &              mu0(i0,i1+1,i2)*(u2(i0,i1+1,i2)-u2(i0,i1,i2))-
     &              mu0(i0,i1,i2)*(u2(i0,i1,i2)-u2(i0,i1-1,i2)))+
     &              fac1*fac2*(
     &              mu0(i0,i1+1,i2)*(u1(i0,i1+1,i2)-u1(i0,i1+1,i2-1))-
     &              mu0(i0,i1,i2)*(u1(i0,i1,i2)-u1(i0,i1,i2-1)))+
     &              fac0**2.d0*(
     &              mu1(i0+1,i1,i2)*(u2(i0+1,i1,i2)-u2(i0,i1,i2))-
     &              mu1(i0,i1,i2)*(u2(i0,i1,i2)-u2(i0-1,i1,i2)))+
     &              fac0*fac2*(
     &              mu1(i0+1,i1,i2)*(u0(i0+1,i1,i2)-u0(i0+1,i1,i2-1))-
     &              mu1(i0,i1,i2)*(u0(i0,i1,i2)-u0(i0,i1,i2-1))))+
     &              rho*u2(i0,i1,i2) + gamma*v2(i0,i1,i2)
            enddo
         enddo
      enddo
c
      return
      end
c

