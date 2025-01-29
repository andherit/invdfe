      function amotry(psum,fac)
      implicit none  

      include 'prec_verb.com'

      real*4 psum(1600),fac
      real*8 amotry

      include 'geom_space.com'
      include 'simplex.com'
      include 'prec_mode.com'

      integer ndim,j
      real*4 fac1,fac2,sim_try(750)
      real*8 mis_try,comp_misfit

      ndim=nxt(1)*nzt(1)
      if (nm.eq.2) ndim=ndim+nxt(2)*nzt(2)+ni
      fac1=(1.-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        sim_try(j)=psum(j)*fac1-mat_sim(ihi,j)*fac2
11    continue
      mis_try=comp_misfit(sim_try,0)
      if (mis_try.lt.mat_mis(ihi)) then
         mat_mis(ihi)=mis_try
         do 12 j=1,ndim
            psum(j)=psum(j)-mat_sim(ihi,j)+sim_try(j)
            mat_sim(ihi,j)=sim_try(j)
12       continue
      endif
      amotry=mis_try
      return
      end

c###############################################################################
