      common/mat_int/matgeninf,tabinf,impact
      real*4 matgeninf(VITMAX),tabinf(NXCMAX)
      real*4 impact(NSTAMAX,2)
