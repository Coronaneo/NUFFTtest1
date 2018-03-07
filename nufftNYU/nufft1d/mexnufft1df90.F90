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
