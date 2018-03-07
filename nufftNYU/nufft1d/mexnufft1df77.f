       #include "fintrf.h"
 
        subroutine mexFunction(nlhs,plhs,nrhs,prhs)
        use mexprint
        implicit double precision (a-h,o-z)
        integer    :: nlhs,nrhs
        mwPointer  :: plhs(*),prhs(*)											  
         
        integer nj,iflag,ms,ier,nk,typ,m，lw,lused
        real*8  eps,fw(0:lw-1)
        real*8,allocatable :: xj(:),sk(:)
        complex*16,allocatable  :: fk(:),cj(:)

        

        if (nrhs .eq. 12) then
          call mxCopyPtrToInteger4(mxGetPr(prhs(12)),typ,1*1)
          if (typ .eq. 1) then
           call mxCopyPtrToInteger4(mxGetPr(prhs(1)),nj,1*1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(4)),iflag,1*1)
           call mxCopyPtrToReal8(mxGetPr(prhs(5)),eps,1*1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(6)),ms,1*1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(9)),lw,1*1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(10)),lused,1*1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(11)),ier,1*1)
           m=mxGetM(prhs(2))
           
           allocate(xj(m))
           call mxCopyPtrToReal8(mxGetPr(prhs(2)),xj,m)
           m=mxGetM(prhs(3))
          
           allocate(cj(m))
           call mxCopyPtrToComplex16(mxGetPr(prhs(3)),cj,m)
           m=mxGetM(prhs(7))
           
           allocate(fk(m))
           call mxCopyPtrToComplex16(mxGetPr(prhs(7)),fk,m)
           call mxCopyPtrToReal8(mxGetPr(prhs(8)),fw,lw)
           call  nufft1d1(nj,xj,cj,iflag,eps,ms,fk,fw,lw,lused,ier)
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
           call mxCopyPtrToInteger4(mxGetPr(prhs(9)),lw,1*1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(10)),lused,1*1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(11)),ier,1*1)
           m=mxGetM(prhs(2))
        
           allocate(xj(m))
           call mxCopyPtrToReal8(mxGetPr(prhs(2)),xj,m)
           m=mxGetM(prhs(3))
        
           allocate(cj(m))
           call mxCopyPtrToComplex16(mxGetPr(prhs(3)),cj,m)
           m=mxGetM(prhs(7))
        
           allocate(fk(m))
           call mxCopyPtrToComplex16(mxGetPr(prhs(7)),fk,m)
           call mxCopyPtrToReal8(mxGetPr(prhs(8)),fw,lw)
           call  nufft1d1(nj,xj,cj,iflag,eps,ms,fk,fw,lw,lused,ier)
           plhs(1)=mxCreateDoubleArray(m)
           call mxCopyComplex16ToPtr(fk,mxGetPr(plhs(1)),m)
           plhs(2)=mxCreateDoubleArray(1)
           call mxCopyInteger4ToPtr(ier,mxGetPr(plhs(2)),1)
          endif
        endif
        if (nrhs .eq. 13) then
           call mxCopyPtrToInteger4(mxGetPr(prhs(1)),nj,1*1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(4)),iflag,1*1)
           call mxCopyPtrToReal8(mxGetPr(prhs(5)),eps,1*1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(6)),nk,1*1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(10)),lw,1*1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(11)),lused,1*1)
           call mxCopyPtrToInteger4(mxGetPr(prhs(12)),ier,1*1)
           allocate(sk(nk))
           call mxCopyPtrToReal8(mxGetPr(prhs(7)),sk,nk)
           m=mxGetM(prhs(2))
       
           allocate(xj(m))
           call mxCopyPtrToReal8(mxGetPr(prhs(2)),xj,m)
           m=mxGetM(prhs(3))
           
           allocate(cj(m))
           call mxCopyPtrToComplex16(mxGetPr(prhs(3)),cj,m)
           m=mxGetM(prhs(8))
           
           allocate(fk(m))
           call mxCopyPtrToComplex16(mxGetPr(prhs(8)),fk,m)
           call mxCopyPtrToReal8(mxGetPr(prhs(9)),fk,lw)
           call  nufft1d3(nj,xj,cj,iflag,eps,nk,sk,fk,fw,lw,lused,ier)
           plhs(1)=mxCreateDoubleArray(m)
           call mxCopyComplex16ToPtr(fk,mxGetPr(plhs(1)),m)
           plhs(2)=mxCreateDoubleArray(1)
           call mxCopyInteger4ToPtr(ier,mxGetPr(plhs(2)),1)
        endif
        deallocate(xj,cj,fk)
        return
        end subroutine
