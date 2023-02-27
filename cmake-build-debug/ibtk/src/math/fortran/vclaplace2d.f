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
c     Computes (f0,f1) = alpha div mu (grad (u0,u1) + grad(u0,u1)^T) + beta (u0,u1) +
c     gamma (v0,v1).
c
c     Computes the side-centered variable coefficient generalized
c     Laplacian, with node-centered coefficient mu and side-centered
c     vector fields (u0,u1) and (v0,v1).
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine stosvclaplace2d(
     &     f0,f1,f_gcw,
     &     alpha,beta,
     &     mu,mu_gcw,
     &     rho0,rho1,rho_gcw,
     &     u0,u1,u_gcw,
     &     gamma,
     &     v0,v1,v_gcw,
     &     ilower0,iupper0,
     &     ilower1,iupper1,
     &     dx,
     &     var_rho,
     &     use_harmonic_interp)
c
      implicit none
c
c     Functions.
c
      double precision a_avg4, h_avg4
c
c     Input.
c
      integer ilower0,iupper0
      integer ilower1,iupper1
      integer f_gcw,mu_gcw,rho_gcw,u_gcw,v_gcw,var_rho
      integer use_harmonic_interp

      double precision alpha,beta,gamma

      double precision mu(ilower0-mu_gcw:iupper0+1+mu_gcw,
     &          ilower1-mu_gcw:iupper1+1+mu_gcw)

      double precision rho0(ilower0-rho_gcw:iupper0+1+rho_gcw,
     &          ilower1-rho_gcw:iupper1+rho_gcw)
      double precision rho1(ilower0-rho_gcw:iupper0+rho_gcw,
     &          ilower1-rho_gcw:iupper1+1+rho_gcw)
      
      double precision u0(ilower0-u_gcw:iupper0+1+u_gcw,
     &          ilower1-u_gcw:iupper1+u_gcw)
      double precision u1(ilower0-u_gcw:iupper0+u_gcw,
     &          ilower1-u_gcw:iupper1+1+u_gcw)

      double precision v0(ilower0-v_gcw:iupper0+1+v_gcw,
     &          ilower1-v_gcw:iupper1+v_gcw)
      double precision v1(ilower0-v_gcw:iupper0+v_gcw,
     &          ilower1-v_gcw:iupper1+1+v_gcw)

      double precision dx(0:2-1)
c
c     Input/Output.
c
      double precision f0(ilower0-f_gcw:iupper0+1+f_gcw,
     &          ilower1-f_gcw:iupper1+f_gcw)
      double precision f1(ilower0-f_gcw:iupper0+f_gcw,
     &          ilower1-f_gcw:iupper1+1+f_gcw)
c
c     Local variables.
c
      integer i0,i1
      double precision    fac0,fac1,rho
      double precision    mu_lower,mu_upper
c
c     Compute the discrete divergence of mu (grad (u0,u1) + grad (u0,u1)^T).
c
      fac0 = 1.d0/(dx(0))
      fac1 = 1.d0/(dx(1))

      do i1 = ilower1,iupper1
         do i0 = ilower0,iupper0+1

            rho = beta
            if (var_rho .eq. 1) then
               rho = rho0(i0,i1)*beta
            endif

            if (use_harmonic_interp .eq. 1) then
                mu_upper = h_avg4(mu(i0,i1),mu(i0+1,i1),
     &                          mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = h_avg4(mu(i0,i1),mu(i0-1,i1),
     &                       mu(i0,i1+1),mu(i0-1,i1+1))
            else
                mu_upper = a_avg4(mu(i0,i1),mu(i0+1,i1),
     &                          mu(i0,i1+1),mu(i0+1,i1+1))
                mu_lower = a_avg4(mu(i0,i1),mu(i0-1,i1),
     &                       mu(i0,i1+1),mu(i0-1,i1+1))
            endif

            f0(i0,i1) = alpha*(
     &           2.d0*fac0**2.d0*(
     &           mu_upper*
     &           (u0(i0+1,i1)-u0(i0,i1))-
     &           mu_lower*
     &           (u0(i0,i1)-u0(i0-1,i1)))+
     &           fac1**2.d0*(mu(i0,i1+1)*(u0(i0,i1+1)-u0(i0,i1))-
     &           mu(i0,i1)*(u0(i0,i1)-u0(i0,i1-1)))+
     &           fac0*fac1*(mu(i0,i1+1)*(u1(i0,i1+1)-u1(i0-1,i1+1))-
     &           mu(i0,i1)*(u1(i0,i1)-u1(i0-1,i1)))) +
     &           rho*u0(i0,i1) + gamma*v0(i0,i1)
         enddo
      enddo

      do i1 = ilower1,iupper1+1
         do i0 = ilower0,iupper0

            rho = beta
            if (var_rho .eq. 1) then
               rho = rho1(i0,i1)*beta
            endif

            if (use_harmonic_interp .eq. 1) then
               mu_upper = h_avg4(mu(i0,i1),mu(i0+1,i1),
     &                          mu(i0,i1+1),mu(i0+1,i1+1))

               mu_lower = h_avg4(mu(i0,i1),mu(i0+1,i1),
     &                          mu(i0,i1-1),mu(i0+1,i1-1))
            else
               mu_upper = a_avg4(mu(i0,i1),mu(i0+1,i1),
     &                          mu(i0,i1+1),mu(i0+1,i1+1))

               mu_lower = a_avg4(mu(i0,i1),mu(i0+1,i1),
     &                          mu(i0,i1-1),mu(i0+1,i1-1))
            endif


            f1(i0,i1) = alpha*(
     &           2.d0*fac1**2.d0*(
     &           mu_upper*
     &           (u1(i0,i1+1)-u1(i0,i1))-
     &           mu_lower*
     &           (u1(i0,i1)-u1(i0,i1-1)))+
     &           fac0**2.d0*(mu(i0+1,i1)*(u1(i0+1,i1)-u1(i0,i1))-
     &           mu(i0,i1)*(u1(i0,i1)-u1(i0-1,i1)))+
     &           fac0*fac1*(mu(i0+1,i1)*(u0(i0+1,i1)-u0(i0+1,i1-1))-
     &           mu(i0,i1)*(u0(i0,i1)-u0(i0,i1-1)))) +
     &           rho*u1(i0,i1) + gamma*v1(i0,i1)
         enddo
      enddo
c
      return
      end
c
