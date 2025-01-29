      common/pick/tabpic,dumpic,refpic,dumref
         real*4 tabpic(NSTAMAX,NSHOMAX,2),dumpic(NSTAMAX,NSHOMAX)
         real*4 refpic(NSTAMAX,NSHOMAX,2),dumref(NSTAMAX,NSHOMAX)

c  max size for shot arrays : NSHOMAX
c  max size for station arrays : NSTAMAX
c  tabpic : real array containing the picking t1-t2
c     for the first arrivals (t1=-1 means no picking
c     and t1=0. only a t2 picking)
c  refpic : real array containing the picking t1-t2
c     for the reflected phases (t1=-1 means no picking
c     and t1=0. only a t2 picking)
