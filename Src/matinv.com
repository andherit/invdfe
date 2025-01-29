c234567
      common/matinv/gcdg,ndim,cm,diag,fn
      real*4 gcdg(1600,1600),cm(1600,1600),fn(1600,1600)
      real*4 diag(1600)
      integer ndim
