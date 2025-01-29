      common/simplex/mat_sim,mat_mis,iter,ilo,ihi,inhi
         real*4 mat_sim(PARMAX+1,PARMAX)
         real*8 mat_mis(PARMAX)
         integer iter,ilo,ihi,inhi
