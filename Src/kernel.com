c234567
      common/ckernel/idatakernel,tk,tk0,dm_res
      integer idatakernel
      real*4 tk(NSHOMAX*NSTAMAX*2),tk0(NSHOMAX*NSTAMAX*2)
      real*4 dm_res
