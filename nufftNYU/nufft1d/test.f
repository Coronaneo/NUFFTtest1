      #include "fintrf.h"  

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

      implicit none


      integer nlhs, nrhs
      mwPointer plhs(*), prhs(*)

      mwPointer mxGetPr, mxGetM, mxGetN
      mwPointer mxCreateDoubleMatrix

      mwPointer x_pr, y_pr
      mwPointer m, n
      mwSize size
 
      real*8 x(100), y

      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      size = m * n
      if(size .gt. 100) then
         call mexErrMsgIdAndTxt ('MATLAB:InputTooBig',+'Input #1: number of elements exceeds buffer')
      endif

      x_pr = mxGetPr(prhs(1))


      call mxCopyPtrToReal8(x_pr,x,size)


      call myprod(x,y,size)

      plhs(1) = mxCreateDoubleMatrix(1,1,0)
      y_pr = mxGetPr(plhs(1))

      call mxCopyReal8ToPtr(y,y_pr,1)      

      return
      end

      subroutine myprod(x, y, n)
      integer i, n
      real*8 x(*), y

      y = 1.0
      do i=1, n
         y = y * x(i) 
      end do     
      return
      end
