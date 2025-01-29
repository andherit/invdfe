c234567
c shot and station maximum number
      integer NSHOMAX,NSTAMAX
      parameter(NSHOMAX=1000,NSTAMAX=1000)
c maximum number of horizontal cells in the velocity model
      integer NXCMAX
      parameter(NXCMAX=5000)
c maximum size of the velocity array [(nx+1)*(nz+1)]
      integer VITMAX
      parameter(VITMAX=2565070)
c maximum number of inversion parameters
      integer PARMAX
      parameter(PARMAX=1600)
