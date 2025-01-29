      common/station/nsta,locsta
         integer nsta
         real*4 locsta(NSTAMAX,3)
