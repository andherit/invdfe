c###############################################################################
c234567
      subroutine searchbysimplex
      implicit none

      include 'geom_space.com'
      include 'prec_stop.com'
      include 'prec_mode.com'
      include 'prec_verb.com'
      include 'simplex.com'

      integer ndim,n,m,i,j,imov,nbfree
      integer nbcarl,tymo
      real*4 sum,psum(1600),rtol
      real*8 mis_save,mis_try,amotry
      real*8 comp_misfit,lastcarl

      include 'searchbysimplex.verb'
      ndim=nxt(1)*nzt(1)
      if (nm.eq.2) ndim=ndim+nxt(2)*nzt(2)+ni
      nbcarl=iter+5*ndim
      tymo=0
      lastcarl=mat_mis(1)
      if (freeze.gt.0) then
         nbfree=freeze
3        if (nbfree.lt.iter) then
            nbfree=nbfree+freeze
            goto 3
         endif
      endif
1     do 12 n=1,ndim
         sum=0.
         do 11 m=1,ndim+1
            sum=sum+mat_sim(m,n)
11       continue
         psum(n)=sum
12    continue
2     ilo=1
      if (mat_mis(1).gt.mat_mis(2)) then
         ihi=1
         inhi=2
      else
         ihi=2
         inhi=1
      endif
      do 13 i=1,ndim+1
         if (mat_mis(i).le.mat_mis(ilo)) ilo=i
         if (mat_mis(i).gt.mat_mis(ihi)) then
            inhi=ihi
            ihi=i
         else
            if (mat_mis(i).gt.mat_mis(inhi)) then
               if (i.ne.ihi) inhi=i
            endif
         endif
13    continue
      include 'search_best.verb'
      if (stop_type.eq.2) then
         rtol=0.
         do 100 m=1,ndim+1
            sum=0.
            do 110 n=1,ndim
               sum=sum+abs(mat_sim(m,n)-psum(n)/float(ndim+1))
110         continue
            rtol=amax1(rtol,sum)
100      continue
         rtol=rtol/float(ndim)
         include 'writertol.verb'
         if (rtol.lt.stop_par) return
      endif
      if (iter.gt.nbfree.and.freeze.gt.0) then
         nbfree=nbfree+freeze
         include 'freeze.verb'
         call freeze_all(.true.)
      endif
      if (stop_type.eq.1.and.iter.ge.stop_par) return
      if (iter.gt.nbcarl) then
         if ((lastcarl-mat_mis(ilo))/lastcarl.lt.0.001) then
            call montecarlo(psum,.1)
            imov=imov+1
            tymo=4
            lastcarl=mat_mis(ilo)
            nbcarl=iter+5*ndim
            goto 2
         endif
         lastcarl=mat_mis(ilo)
         nbcarl=iter+ndim
      endif
      iter=iter+2
      mis_try=amotry(psum,-1.0)
      tymo=1
      imov=imov+2
      if (mis_try.le.mat_mis(ilo)) then
         mis_try=amotry(psum,2.0)
         tymo=2
         imov=imov+1
      else
         if (mis_try.ge.mat_mis(inhi)) then
            mis_save=mat_mis(ihi)
            mis_try=amotry(psum,.5)
            if (mis_try.ge.mis_save) then
               call montecarlo(psum,1.)
               tymo=4
            else
               tymo=3
            endif
         else
            iter=iter-1
            imov=imov-1
         endif
      endif
      goto 2
      end

c###############################################################################
