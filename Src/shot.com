      common/shot/nsho,locsho
         integer nsho
         real*4 locsho(NSHOMAX,3)
