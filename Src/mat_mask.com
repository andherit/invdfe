      common/mat_mask/maskgen,maskref,maskup,maskdw,mcrref
      integer maskgen(VITMAX)
      real*4 maskref(2*VITMAX+NXCMAX)
      real*4 maskdw(PARMAX),maskup(PARMAX)
      real*4 mcrref(PARMAX)
