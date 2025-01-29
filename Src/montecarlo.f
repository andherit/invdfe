c###############################################################################
      subroutine montecarlo(psum,fakir)
      implicit none

      real*4 psum(1600),fakir

      include 'geom_space.com'
      include 'mat_mask.com'
      include 'simplex.com'
      include 'prec_mode.com'
  
      integer ndim,i,j,nit
      real*4 sim_try(750),rdom,vmax,vmin,vmed
      real*8 mis_try,comp_misfit
      real*4 fac

      ndim=nxt(1)*nzt(1)
      if (nm.eq.2) ndim=ndim+nxt(2)*nzt(2)+ni
      fac=fakir
      nit=1
5     do 10 i=1,ndim
         vmed=mat_sim(ilo,i)
         vmax=maskup(i)
         vmin=maskdw(i)
         if (rdom(0).lt..5) then
            sim_try(i)=(vmin-vmed)*fac*rdom(0)+vmed
         else
            sim_try(i)=(vmax-vmed)*fac*rdom(0)+vmed
         endif
10    continue
      mis_try=comp_misfit(sim_try,0)
      if (mis_try.le.mat_mis(ihi)) then
         mat_mis(ihi)=mis_try
         do 30 j=1,ndim
            psum(j)=psum(j)-mat_sim(ihi,j)+sim_try(j)
            mat_sim(ihi,j)=sim_try(j)
30       continue
         return
      else
         nit=nit+1
         iter=iter+1
         if (nit.gt.10*ndim) then
            fac=fac*0.9
         endif
         goto 5
      endif
      return
      end

c###############################################################################
