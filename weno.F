      module weno

      ! Provides weno style reconstructions at left cell face
      ! call as: 

#include "cppdefs.opt"

      use param
      use dimensions

      private

#  define WENO_Z

      public fluxes_uv

      contains

! ----------------------------------------------------------------------
      function weno5(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) ![
!
!     WENO5: Weighted Essentially Non-Oscillatory scheme
!            with 5th-order accuracy
!
! This function computes tracer reconstruction at the grid cell's left
! edge (u-point i-1/2 or v-point j-1/2). WENO5 uses a convex combination
! of the polynomials reconstructed on the three ENO3 stencils in order
! to achieve higher accuracy on smooth profiles. Both left and right
! combinations of mirror symetric stencils around i-1/2 (j-1/2) points
! are computed. For stability, the upstream stencil combination,
! identified by the sign of ua (u or v), is selected.
!
! The scheme includes improvements from Borges et al., 2008 (WENO_Z), i.e.,
! smoothness indicators of higher order (with new non-oscillatory weights)
! that provide a scheme with less dissipation, higher resolution and better
! monotonicity preservation than the classical WENO5. In the meantime, WENO_Z
! removes the need to tune parameter Eps for added dissipation (needed in the
! original scheme due to suboptimal performance on critical points).
!
! References:
!  Guang-Shan JIANG and Chi-Wang SHU, 1996: Efficient Implementation of
!  Weighted ENO Schemes. JOURNAL OF COMPUTATIONAL PHYSICS 126, 202–228.
!
!  Rong WANG and Raymond J. SPITERI, 2007: Linear instability of the
!  fifth-order WENO method, SIAM J. Numer. Anal., 45, 1871-1901.
!
!  Borges R., M. Carmona, B. Costa, W.S. Don, 2008: An improved weighted
!  essentially non-oscillatory scheme for hyperbolic conservation laws.
!  Journal of Computational Physics 227 (2008) 3191–3211
!
!  Implementation: P. Marchesiello and J. Demange 2013
!
      implicit none
      real    :: weno5
      real    :: q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua
      real    :: IS0, IS1, IS2
      real    :: d0, d1, d2
      real    :: a0, a1, a2
      real    :: w0, w1, w2
      real    :: p0, p1, p2
      real    :: Eps, cff1, cff2, T5


! Epsilon value for computing weigths
#  ifdef WENO_Z
      Eps = 1.e-40
#  else
!  --> from 1.e-7 (smoother) to 1.e-5 (sharper)
      Eps = 1.e-6
#  endif
! Weigths coefficients
      d0=1./10.
      d1=6./10.
      d2=3./10.
! other coefficients
      cff1=13./12.
      cff2=1./6.
!
      if (ua .ge. 0.) then        ! Take upstream stencils for stability
!
! === Reconstruction at u_i-1/2 using LEFT stencils ===
!
! Nonlinear Smoothness Indicators
        IS0 = cff1*(q_im3 - 2.*q_im2 +   q_im1)**2
     &      + 0.25*(q_im3 - 4.*q_im2 + 3*q_im1)**2
        IS1 = cff1*(q_im2 - 2.*q_im1 +   q_i  )**2
     &      + 0.25*(q_im2            -   q_i  )**2
        IS2 = cff1*(q_im1 - 2.*q_i   +   q_ip1)**2
     &   + 0.25*(3.*q_im1 - 4.*q_i   +   q_ip1)**2
#  ifdef WENO_Z
! Non-normalized stencil weights
! with parameter T5 of new smoothness indicator
        T5  = abs(IS2-IS0)
        a0  = d0*(1+T5/(Eps+IS0))
        a1  = d1*(1+T5/(Eps+IS1))
        a2  = d2*(1+T5/(Eps+IS2))
#  else
        a0  = d0/(Eps+IS0)**2
        a1  = d1/(Eps+IS1)**2
        a2  = d2/(Eps+IS2)**2
#  endif
! Normalized Weigths
        w0  = a0/(a0+a1+a2)
        w1  = a1/(a0+a1+a2)
        w2  = a2/(a0+a1+a2)
! Polynomials
        p0  = cff2*(2.*q_im3 - 7.*q_im2 + 11.*q_im1)
        p1  = cff2*(  -q_im2 + 5.*q_im1 +  2.*q_i  )
        p2  = cff2*(2.*q_im1 + 5.*q_i   -     q_ip1)

      else
!
! === Reconstruction at u_i-1/2 using RIGHT stencils ===
!
! Nonlinear Smoothness Indicators
        IS0 = cff1*(q_ip2 - 2.*q_ip1 +   q_i  )**2
     &      + 0.25*(q_ip2 - 4.*q_ip1 + 3*q_i  )**2
        IS1 = cff1*(q_ip1 - 2.*q_i   +   q_im1)**2
     &      + 0.25*(q_ip1            -   q_im1)**2
        IS2 = cff1*(q_i   - 2.*q_im1 +   q_im2)**2
     &   + 0.25*(3.*q_i   - 4.*q_im1 +   q_im2)**2
#  ifdef WENO_Z
! Non-normalized stencil weights
! with parameter T5 of new smoothness indicator
        T5  = abs(IS2-IS0)
        a0  = d0*(1+T5/(Eps+IS0))
        a1  = d1*(1+T5/(Eps+IS1))
        a2  = d2*(1+T5/(Eps+IS2))
#  else
        a0  = d0/(Eps+IS0)**2
        a1  = d1/(Eps+IS1)**2
        a2  = d2/(Eps+IS2)**2
#  endif
! Normalized Weigths
        w0  = a0/(a0+a1+a2)
        w1  = a1/(a0+a1+a2)
        w2  = a2/(a0+a1+a2)
! Polynomials
        p0  = cff2*(2.*q_ip2 - 7.*q_ip1 + 11.*q_i  )
        p1  = cff2*(  -q_ip1 + 5.*q_i   +  2.*q_im1)
        p2  = cff2*(2.*q_i   + 5.*q_im1 -     q_im2)

      endif

      ! Combination
      weno5 = w0*p0 + w1*p1 + w2*p2

      return
      end function weno5 !]
! ----------------------------------------------------------------------
      function weno3( q_im2, q_im1, q_i, q_ip1, ua)  ![
      ! This function computes tracer reconstruction at the grid cell's left
      ! edge (u-point i-1/2 or v-point j-1/2).
!
      implicit none
      real    :: weno3
      real :: f(-2:1) 
      real    :: q_im2, q_im1, q_i, q_ip1, ua
      real    :: IS0, IS1
      real    :: a0, a1
      real    :: w0, w1
      real    :: p0, p1
      real    :: Eps, d0,d1, T3
!
! Epsilon value for computing weigths
#  ifdef WENO_Z
      Eps = 1.e-40
#  else
!  --> from 1.e-7 (smoother) to 1.e-5 (sharper)
      Eps = 1.e-6
#  endif
      d0=1./3.
      d1=2./3.
!
!     if (ua==0)  weno3 = 0;return
!
!
      if (ua .ge. 0.) then        ! Take upstream stencils for stability

        IS0 = (q_im1-q_im2)**2
        IS1 = (q_im1-q_i)**2
!       IS0 = (f(-1)-f(-2))**2
!       IS1 = (f(-1)-f(0))**2
#  ifdef WENO_Z
        T3 = abs(IS1-IS0)
        a0  = d0*(1+T3/(Eps+IS0))
        a1  = d1*(1+T3/(Eps+IS1))
#  else
        a0  = 1./(3.*(Eps+IS0)**2)
        a1  = 2./(3.*(Eps+IS1)**2)
#  endif
        w0  = a0/(a0+a1)
        w1  = a1/(a0+a1)
        p0  = 1./2.*(3.*q_im1-q_im2)
        p1  = 1./2.*(q_im1+q_i)
!       p0  = 1./2.*(3.*f(-1)-q_im2)
!       p1  = 1./2.*(q_im1+q_i)

      else

        IS0 = (q_i-q_ip1)**2
        IS1 = (q_im1-q_i)**2
#  ifdef WENO_Z
        T3 = abs(IS1-IS0)
        a0  = d0*(1+T3/(Eps+IS0))
        a1  = d1*(1+T3/(Eps+IS1))
#  else
        a0  = 1./(3.*(Eps+IS0)**2)
        a1  = 2./(3.*(Eps+IS1)**2)
#  endif
        w0  = a0/(a0+a1)
        w1  = a1/(a0+a1)
        p0  = 1./2.*(3.*q_i-q_ip1)
        p1  = 1./2.*(q_im1+q_i)

      endif

      weno3 = w0*p0 + w1*p1

      end function weno3 !]
! ----------------------------------------------------------------------
      subroutine fluxes_uv(uu,uv,vu,vv,flxU,flxV,u,v) ![
!     fluxes of horizontal momentum at left/south cell faces
      implicit none

!     input/output
      real,intent(out),dimension(GLOBAL_2D_ARRAY) :: uu,uv,vu,vv
      real,intent(in) ,dimension(GLOBAL_2D_ARRAY) :: flxU,flxV
      real,intent(in) ,dimension(GLOBAL_2D_ARRAY) :: u,v

!     local
      integer :: i,j
      real    :: uval,vval

!    valu = weno5(u(i-3,j),u(i-2,j),u(i-1,j),u(i,j),u(i+1,j),u(i+2,j),uu(i,j))

      ! we need to reduce order near the physical boundaries and the masks

      do j=1,ny
        do i=1,nx+1
          uu(i,j) = 0.5*(FlxU(i,j)+FlxU(i-1,j))
          uv(i,j) = 0.5*(FlxU(i,j)+FlxU(i,j-1))
!         if (uu(i,j)>0) then
!          uval = weno3(u(i-2:i  ,j))
!         else
!          uval = weno3(u(i-1:i+1,j))
!         endif
          uval = weno3(u(i-2,j),u(i-1,j),u(i,j),u(i+1,j),uu(i,j))
          vval = weno3(v(i-2,j),v(i-1,j),v(i,j),v(i+1,j),uv(i,j))
          uu(i,j) = uu(i,j)*uval
          uv(i,j) = uv(i,j)*vval
        enddo
      enddo
      do j=1,ny+1
        do i=1,nx
          vu(i,j) = 0.5*(FlxV(i,j)+FlxV(i-1,j))
          vv(i,j) = 0.5*(FlxV(i,j)+FlxV(i,j-1))
          uval = weno3(u(i,j-2),u(i,j-1),u(i,j),u(i,j+1),vu(i,j))
          vval = weno3(v(i,j-2),v(i,j-1),v(i,j),v(i,j+1),vv(i,j))
          vu(i,j) = vu(i,j)*uval
          vv(i,j) = vv(i,j)*vval
        enddo
      enddo

      end subroutine fluxes_uv  !]
! ----------------------------------------------------------------------

      end module weno
