
      
      function next235(base)
      implicit none
      integer next235, numdiv
      real*8 base

      next235 = 2 * int(base/2d0+.9999d0)
      if (next235.le.0) next235 = 2

100   numdiv = next235
      do while (numdiv/2*2 .eq. numdiv)
         numdiv = numdiv /2
      enddo
      do while (numdiv/3*3 .eq. numdiv)
         numdiv = numdiv /3
      enddo
      do while (numdiv/5*5 .eq. numdiv)
         numdiv = numdiv /5
      enddo
      if (numdiv .eq. 1) return
      next235 = next235 + 2
      goto 100
      end

