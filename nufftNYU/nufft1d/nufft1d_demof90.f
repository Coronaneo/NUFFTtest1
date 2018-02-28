cc Copyright (C) 2004-2009: Leslie Greengard and June-Yub Lee 
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
c
      subroutine nufft1d1_demof90(nj,xj,cj,iflag,eps,ms,fk,ier)
      implicit none
      integer ier,iflag,q,num
      integer ms,nj
      real*8 eps,pi
      real*8 xj(nj)
      parameter (num=200)
      complex*16 cj(nj),fk(-ms/2:(ms-1)/2)


c     -----------------------
c     call 1D Type1 method
c     -----------------------
c
      do q=1,num
	   call nufft1d1f90(nj,xj,cj,iflag,eps, ms,fk,ier)
	enddo
	return
	end
c