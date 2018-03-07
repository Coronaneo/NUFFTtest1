!#include "fintrf.h"cc Copyright (C) 2004-2009: Leslie Greengard and June-Yub Lee 
! Contact: greengard@cims.nyu.edu
! 
! This program is free software; you can redistribute it and/or modify 
! it under the terms of the GNU General Public License as published by 
! the Free Software Foundation; eithe version 2 of the License, or 
! (at your option) any later version.  This program is distributed in 
! the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
! even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
! PARTICULAR PURPOSE.  See the GNU General Public License for more 
! details. You should have received a copy of the GNU General Public 
! License along with this program; 
! if not, see <http://www.gnu.org/licenses/>.
!
!
!  NUFFT 1.2 release notes:
!
!  These codes are asymptotically fast (O(N log N)), but not optimizedc
!  1) We initialize the FFT on every call.
!
!  2) We do not precompute the exponentials involved in "fast Gaussian
!  gridding".
!
!  3) We do not block structure the code so that irregularly placed points
!  are interpolated (gridded) in a cache-aware fashion.
!
!  4) We use the Netlib FFT library (www.netlib.org) 
!     rather than the state of the art FFTW package (www.fftw.org).
!
!  Different applications have different needs, and we have chosen
!  to provide the simplest code as a reasonable efficient template.
!
       #include "fintrf.h"
 
        subroutine mexFunction(nlhs,plhs,nrhs,prhs)
        use mexprint
        implicit double precision (a-h,o-z)
        integer    :: nlhs,nrhs
        mwPointer  :: plhs(*),prhs(*)
         
        integer nj,iflag,ms,ier,nk,typ,m
        real*8  eps
        real*8,allocatable :: xj(:),sk(:)
        complex*16,allocatable  :: fk(:),cj(:)

        

        if (nrhs .eq. 9) then
          call mxCopyPtrToInteger4(mxGetPr(prhs(9)),typ,1*1)
          if (typ .eq. 1) then
           call mxCopyPtrToInteger4(mxGetPr(prhs(1)),nj,1*1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(4)),iflag,1*1)
           call mxCopyPtrToReal8(mxGetPr(prhs(5)),eps,1*1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(6)),ms,1*1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(8)),ier,1*1)
           m=mxGetM(prhs(2))
           
           allocate(xj(m))
           call mxCopyPtrToReal8(mxGetPr(prhs(2)),xj,m)
           m=mxGetM(prhs(2))
          
           allocate(cj(m))
           call mxCopyPtrToReal8(mxGetPr(prhs(3)),cj,m)
           m=mxGetM(prhs(2))
           
           allocate(fk(m))
           call mxCopyPtrToComplex16(mxGetPr(prhs(7)),fk,m)
           call  nufft1d1f90(nj,xj,cj,iflag,eps,ms,fk,ier)
           plhs(1)=mxCreateDoubleArray(m)
           call mxCopyComplex16ToPtr(fk,mxGetPr(plhs(1)),m)
           plhs(2)=mxCreateDoubleArray(1)
           call mxCopyInteger4ToPtr(ier,mxGetPr(plhs(2)),1)
          endif
          if (typ .eq. 2) then 
           call mxCopyPtrToInteger4(mxGetPr(prhs(1)),nj,1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(4)),iflag,1)
           call mxCopyPtrToReal8(mxGetPr(prhs(5)),eps,1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(6)),ms,1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(8)),ier,1)
           m=mxGetM(prhs(2))
        
           allocate(xj(m))
           call mxCopyPtrToReal8(mxGetPr(prhs(2)),xj,m)
           m=mxGetM(prhs(2))
        
           allocate(cj(m))
           call mxCopyPtrToReal8(mxGetPr(prhs(3)),cj,m)
           m=mxGetM(prhs(2))
        
           allocate(fk(m))
           call mxCopyPtrToComplex16(mxGetPr(prhs(7)),fk,m)
           call  nufft1d2f90(nj,xj,cj,iflag,eps,ms,fk,ier)
           plhs(1)=mxCreateDoubleArray(m)
           call mxCopyComplex16ToPtr(fk,mxGetPr(plhs(1)),m)
           plhs(2)=mxCreateDoubleArray(1)
           call mxCopyInteger4ToPtr(ier,mxGetPr(plhs(2)),1)
          endif
        endif
        if (nrhs .eq. 10) then
           call mxCopyPtrToInteger4(mxGetPr(prhs(1)),nj,1*1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(4)),iflag,1*1)
           call mxCopyPtrToReal8(mxGetPr(prhs(5)),eps,1*1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(6)),nk,1*1)
           
           call mxCopyPtrToInteger4(mxGetPr(prhs(9)),ier,1*1)
           allocate(sk(nk))
           call mxCopyPtrToReal8(mxGetPr(prhs(7)),sk,nk)
           m=mxGetM(prhs(2))
       
           allocate(xj(m))
           call mxCopyPtrToReal8(mxGetPr(prhs(2)),xj,m)
           m=mxGetM(prhs(2))
           
           allocate(cj(m))
           call mxCopyPtrToReal8(mxGetPr(prhs(3)),cj,m)
           m=mxGetM(prhs(2))
           
           allocate(fk(m))
           call mxCopyPtrToComplex16(mxGetPr(prhs(8)),fk,m)
           call  nufft1d3f90(nj,xj,cj,iflag,eps,nk,sk,fk,ier)
           plhs(1)=mxCreateDoubleArray(m)
           call mxCopyComplex16ToPtr(fk,mxGetPr(plhs(1)),m)
           plhs(2)=mxCreateDoubleArray(1)
           call mxCopyInteger4ToPtr(ier,mxGetPr(plhs(2)),1)
        endif
        deallocate(xj,cj,fk)
        return
        end subroutine


!*************************************************************************************************
      subroutine nufft1d1f90(nj,xj,ocj,iflag,eps,ms,fk,ier)
      implicit none
      integer ier,iflag,istart,iw1,iwtot,iwsav
      integer j,jb1,jb1u,jb1d,k1,ms,next235,nf1,nj,nspread
      real*8 ocross,ocross1,diff1,eps,hx,pi,rat,r2lamb,t1,tau
      real*8 xc(-47:47),xj(nj)
      parameter (pi=3.141592653589793d0)
      complex*16 ocj(nj),fk(-ms/2:(ms-1)/2),zz,occj
! ----------------------------------------------------------------------
      real*8, allocatable, save :: fw(:)
! ----------------------------------------------------------------------

      ier = 0
      if ((eps.lt.1d-13).or.(eps.gt.1d-1)) then
         ier = 1
         return
      endif
      if (eps.le.1d-11) then
         rat = 3.0d0
      else 
         rat = 2.0d0
      endif
      nspread = int(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0)
      nf1 = rat*ms
      if (2*nspread.gt.nf1) then
         nf1 = next235(2d0*nspread) 
      endif 

      r2lamb = rat*rat * nspread / (rat*(rat-.5d0))
      hx = 2*pi/nf1

      iw1 = 2*nf1
      iwsav = iw1+nspread+1
      iwtot = iwsav+4*nf1+15
      allocate ( fw(0:iwtot) )

      t1 = pi/r2lamb
      do k1 = 1, nspread
         fw(iw1+k1) = exp(-t1*k1**2)
      enddo
      call dcffti(nf1,fw(iwsav))

      do k1 = 0, 2*nf1-1
         fw(k1) = dcmplx(0d0,0d0)
      enddo

      do j = 1, nj
         occj = ocj(j)/dble(nj)

         jb1 = int((xj(j)+pi)/hx)
         diff1 = (xj(j)+pi)/hx - jb1
         jb1 = mod(jb1, nf1)
         if (jb1.lt.0) jb1=jb1+nf1

         xc(0) = exp(-t1*diff1**2)
         ocross = xc(0)
         ocross1 = exp(2d0*t1 * diff1)
         do k1 = 1, nspread
            ocross = ocross * ocross1
            xc(k1) = fw(iw1+k1)*ocross
         enddo
         ocross = xc(0)
         ocross1 = 1d0/ocross1
         do k1 = 1, nspread-1
            ocross = ocross * ocross1
            xc(-k1) = fw(iw1+k1)*ocross
         enddo

         jb1d = min(nspread-1, jb1)
         jb1u = min(nspread, nf1-jb1-1)
         do k1 = -nspread+1, -jb1d-1
	    istart = 2*(jb1+k1+nf1)
            zz=xc(k1)*occj
            fw(istart)=fw(istart)+dreal(zz)
            fw(istart+1)=fw(istart+1)+dimag(zz)
         enddo
         do k1 = -jb1d, jb1u
	    istart = 2*(jb1+k1)
            zz=xc(k1)*occj
            fw(istart)=fw(istart)+dreal(zz)
            fw(istart+1)=fw(istart+1)+dimag(zz)
         enddo
         do k1 = jb1u+1, nspread
	    istart = 2*(jb1+k1-nf1)
            zz=xc(k1)*occj
            fw(istart)=fw(istart)+dreal(zz)
            fw(istart+1)=fw(istart+1)+dimag(zz)
         enddo
      enddo

      if (iflag .ge. 0) then
         call dcfftb(nf1,fw(0),fw(iwsav))
      else
         call dcfftf(nf1,fw(0),fw(iwsav))
      endif

      tau = pi * r2lamb / dble(nf1)**2
      ocross1 = 1d0/sqrt(r2lamb)
      zz = dcmplx(fw(0),fw(1))
      fk(0) = ocross1*zz
      do k1 = 1, (ms-1)/2
         ocross1 = -ocross1
         ocross = ocross1*exp(tau*dble(k1)**2)
	 zz = dcmplx(fw(2*k1),fw(2*k1+1))
         fk(k1) = ocross*zz
	 zz = dcmplx(fw(2*(nf1-k1)),fw(2*(nf1-k1)+1))
         fk(-k1) = ocross*zz
      enddo
      if (ms/2*2.eq.ms) then
         ocross = -ocross1*exp(tau*dble(ms/2)**2)
         zz = dcmplx(fw(2*nf1-ms),fw(2*nf1-ms+1))
         fk(-ms/2) = ocross*zz
      endif
      deallocate(fw)
      return
      end

      subroutine nufft1d2f90(nj,xj,ocj, iflag,eps, ms,fk,ier)
      implicit none
      integer ier,iflag,iw1,iwsav,iwtot,j,jb1,jb1u,jb1d,k1
      integer ms,next235,nf1,nj,nspread,nw
      real*8 ocross,ocross1,diff1,eps,hx,pi,rat,r2lamb,t1
      real*8 xj(nj),xc(-47:47)
      parameter (pi=3.141592653589793d0)
      complex*16 ocj(nj), fk(-ms/2:(ms-1)/2)
      complex*16 zz
! ----------------------------------------------------------------------
      real*8, allocatable, save :: fw(:)
! ----------------------------------------------------------------------

      ier = 0
      if ((eps.lt.1d-13).or.(eps.gt.1d-1)) then
         ier = 1
         return
      endif
      if (eps.le.1d-11) then
         rat = 3.0d0
      else 
         rat = 2.0d0
      endif
      nspread = int(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0)
      nf1 = rat*ms
      if (2*nspread.gt.nf1) then
         nf1 = next235(2d0*nspread) 
      endif 

      ier = 0
      if ((eps.lt.1d-13).or.(eps.gt.1d-1)) then
         ier = 1
         return
      endif
      if (eps.le.1d-11) then
         rat = 3.0d0
      else 
         rat = 2.0d0
      endif

      nspread = int(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0)
      nf1 = rat*ms
      if (2*nspread.gt.nf1) then
         nf1 = next235(2d0*nspread) 
      endif 

      r2lamb = rat*rat * nspread / (rat*(rat-.5d0))
      hx = 2*pi/nf1

      iw1 = 2*nf1
      iwsav = iw1 + nspread+1
      iwtot = iwsav + 4*nf1 + 15
      allocate ( fw(0:iwtot))

      t1 = pi/r2lamb
      do k1 = 1, nspread
         fw(iw1+k1) = exp(-t1*k1**2)
      enddo
      call dcffti(nf1,fw(iwsav))

      t1 = pi * r2lamb / dble(nf1)**2
      ocross1 = 1d0/sqrt(r2lamb)
      zz = ocross1*fk(0)
      fw(0) = dreal(zz)
      fw(1) = dimag(zz)
      do k1 = 1, (ms-1)/2
         ocross1 = -ocross1
         ocross = ocross1*exp(t1*dble(k1)**2)
         zz = ocross*fk(k1)
         fw(2*k1) = dreal(zz)
         fw(2*k1+1) = dimag(zz)
         zz = ocross*fk(-k1)
         fw(2*(nf1-k1)) = dreal(zz)
         fw(2*(nf1-k1)+1) = dimag(zz)
      enddo
      ocross = -ocross1*exp(t1*dble(ms/2)**2)
      if (ms/2*2.eq.ms) then
	 zz = ocross*fk(-ms/2)
         fw(2*nf1-ms) = dreal(zz)
         fw(2*nf1-ms+1) = dimag(zz)
      endif
      do k1 = (ms+1)/2, nf1-ms/2-1
         fw(2*k1) = dcmplx(0d0, 0d0)
         fw(2*k1+1) = dcmplx(0d0, 0d0)
      enddo

      if (iflag .ge. 0) then
         call dcfftb(nf1,fw(0),fw(iwsav))
      else
         call dcfftf(nf1,fw(0),fw(iwsav))
      endif

      t1 = pi/r2lamb
      do j = 1, nj
         ocj(j) = dcmplx(0d0,0d0)
         jb1 = int((xj(j)+pi)/hx)
         diff1 = (xj(j)+pi)/hx - jb1
         jb1 = mod(jb1, nf1)
         if (jb1.lt.0) jb1=jb1+nf1
         xc(0) = exp(-t1*diff1**2)
         ocross = xc(0)
         ocross1 = exp(2d0*t1 * diff1)
         do k1 = 1, nspread
            ocross = ocross * ocross1
            xc(k1) = fw(iw1+k1)*ocross
         enddo
         ocross = xc(0)
         ocross1 = 1d0/ocross1
         do k1 = 1, nspread-1
            ocross = ocross * ocross1
            xc(-k1) = fw(iw1+k1)*ocross
         enddo

         jb1d = min(nspread-1, jb1)
         jb1u = min(nspread, nf1-jb1-1)
         do k1 = -nspread+1, -jb1d-1
	    zz = dcmplx(fw(2*(jb1+k1+nf1)),fw(2*(jb1+k1+nf1)+1))
            ocj(j) = ocj(j) + xc(k1)*zz
         enddo
         do k1 = -jb1d, jb1u
	    zz = dcmplx(fw(2*(jb1+k1)),fw(2*(jb1+k1)+1))
            ocj(j) = ocj(j) + xc(k1)*zz
         enddo
         do k1 = jb1u+1, nspread
	    zz = dcmplx(fw(2*(jb1+k1-nf1)),fw(2*(jb1+k1-nf1)+1))
            ocj(j) = ocj(j) + xc(k1)*zz
         enddo
      enddo
      deallocate(fw)
      return
      end

      subroutine nufft1d3f90(nj,xj,ocj, iflag,eps, nk,sk,fk,ier)
      implicit none
      integer ier,iw1,iwsave,iwtot,j,jb1,k1,kb1,kmax,nj,iflag,nk
      integer next235,nf1,nspread
      real*8 ang,ocross,ocross1,diff1,eps,hx,hs,rat,pi,r2lamb1
      real*8 sm,sb,t1,t2,xm,xb
      real*8 xc(-47:47), xj(nj), sk(nk)
      parameter (pi=3.141592653589793d0)
      complex*16 ocj(nj), fk(nk), zz, ocs

! ----------------------------------------------------------------------
      integer nw, istart
      real*8, allocatable, save :: fw(:)
! ----------------------------------------------------------------------

      ier = 0
      if ((eps.lt.1d-13).or.(eps.gt.1d-1)) then
         ier = 1
         return
      endif

      t1 = xj(1)
      t2 = xj(1)
      do j = 2, nj
         if (xj(j).gt.t2) then
             t2=xj(j)
         else if (xj(j).lt.t1) then
             t1=xj(j)
         endif
      enddo
      xb = (t1+t2) / 2d0
      xm = max(t2-xb,-t1+xb)  ! max(abs(t2-xb),abs(t1-xb))

      t1 = sk(1)
      t2 = sk(1)
      do k1 = 2, nk
         if (sk(k1).gt.t2) then
             t2=sk(k1)
         else if (sk(k1).lt.t1) then
             t1=sk(k1)
         endif
      enddo
      sb = (t1+t2) / 2d0
      sm = max(t2-sb,-t1+sb)

      if (eps.le.1d-11) then
         rat = sqrt(3.0d0)
      else 
         rat = sqrt(2.0d0)
      endif

      nspread = int(-log(eps)/(pi*(rat-1d0)/(rat-.5d0)) + .5d0)
      t1 = 2d0/pi * xm*sm
      nf1 = next235(rat*max(rat*t1+2*nspread,2*nspread/(rat-1)))
      rat = (sqrt(nf1*t1+nspread**2)-nspread)/t1

      r2lamb1 = rat*rat * nspread / (rat*(rat-.5d0))
      hx = pi/(rat*sm)
      hs = 2d0*pi/dble(nf1)/hx            ! hx hs = 2.pi/nf1

      kmax = int(nf1*(r2lamb1-nspread)/r2lamb1+.1d0)
      iw1 = 2*nf1
      iwsave = iw1 + nspread+1
      iwtot = iwsave + 16+4*nf1
      allocate ( fw(0:iwtot-1) )

      t1 = pi/r2lamb1
      do k1 = 1, nspread
         fw(iw1+k1) = exp(-t1*k1**2)
      enddo

      call dcffti(nf1,fw(iwsave))

      do k1 = 0, 2*nf1-1
         fw(k1) = dcmplx(0d0,0d0)
      enddo

      t1 = pi/r2lamb1
      if (iflag .lt. 0) sb = -sb
      do j = 1, nj
         jb1 = int(dble(nf1/2) + (xj(j)-xb)/hx)
         diff1 = dble(nf1/2) + (xj(j)-xb)/hx - jb1
         ang = sb*xj(j)
         ocs = dcmplx(cos(ang),sin(ang)) * ocj(j)

         xc(0) = exp(-t1*diff1**2)
         ocross = xc(0)
         ocross1 = exp(2d0*t1 * diff1)
         do k1 = 1, nspread
            ocross = ocross * ocross1
            xc(k1) = fw(iw1+k1)*ocross
         enddo
         ocross = xc(0)
         ocross1 = 1d0/ocross1
         do k1 = 1, nspread-1
            ocross = ocross * ocross1
            xc(-k1) = fw(iw1+k1)*ocross
         enddo

         do k1 = -nspread+1, nspread
	    istart = 2*(jb1+k1)
	    zz = xc(k1)*ocs
            fw(istart) = fw(istart) + dreal(zz)
            fw(istart+1) = fw(istart+1) + dimag(zz)
         enddo
      enddo
      if (iflag .lt. 0) sb = -sb

      t1 = pi * r2lamb1 / dble(nf1)**2
      ocross1 = (1d0-2d0*mod(nf1/2,2))/r2lamb1
      zz = dcmplx(fw(nf1),fw(nf1+1))
      zz = ocross1*zz
      fw(nf1) = dreal(zz)
      fw(nf1+1) = dimag(zz)
      do k1 = 1, kmax
         ocross1 = -ocross1
         ocross = ocross1*exp(t1*dble(k1)**2)
         zz = dcmplx(fw(nf1-2*k1),fw(nf1-2*k1+1))
         zz = ocross*zz
         fw(nf1-2*k1) = dreal(zz)
         fw(nf1-2*k1+1) = dimag(zz)
         zz = dcmplx(fw(nf1+2*k1),fw(nf1+2*k1+1))
         zz = ocross*zz
         fw(nf1+2*k1) = dreal(zz)
         fw(nf1+2*k1+1) = dimag(zz)
      enddo

      if (iflag .ge. 0) then
         call dcfftb(nf1,fw(0),fw(iwsave))
      else
         call dcfftf(nf1,fw(0),fw(iwsave))
      endif
      do k1 = 1, kmax+nspread, 2
         fw(nf1+2*k1) = -fw(nf1+2*k1)
         fw(nf1+2*k1+1) = -fw(nf1+2*k1+1)
         fw(nf1-2*k1) = -fw(nf1-2*k1)
         fw(nf1-2*k1+1) = -fw(nf1-2*k1+1)
      enddo

      t1 = pi/r2lamb1
      do j = 1, nk
         kb1 = int(dble(nf1/2) + (sk(j)-sb)/hs)
         diff1 = dble(nf1/2) + (sk(j)-sb)/hs - kb1

         ! exp(-t1*(diff1-k1)**2) = xc(k1)
         xc(0) = exp(-t1*diff1**2)
         ocross = xc(0)
         ocross1 = exp(2d0*t1 * diff1)
         do k1 = 1, nspread
            ocross = ocross * ocross1
            xc(k1) = fw(iw1+k1)*ocross
         enddo
         ocross = xc(0)
         ocross1 = 1d0/ocross1
         do k1 = 1, nspread-1
            ocross = ocross * ocross1
            xc(-k1) = fw(iw1+k1)*ocross
         enddo

         fk(j) = dcmplx(0d0,0d0)
         do k1 = -nspread+1, nspread
	    zz = dcmplx(fw(2*(kb1+k1)),fw(2*(kb1+k1)+1))
            fk(j) = fk(j) + xc(k1)*zz
         enddo
      enddo

      if (iflag .lt. 0) xb = -xb
      t1 = r2lamb1/(4d0*pi) * hx**2
      do j = 1, nk
         fk(j) = (exp(t1*(sk(j)-sb)**2))*fk(j)
         ang = (sk(j)-sb)*xb
         fk(j) = dcmplx(cos(ang),sin(ang)) * fk(j)
      enddo
      deallocate(fw)
      return
      end