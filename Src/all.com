      common/direct/dir_vit
         real*4 dir_vit(1002501)
      common/geom_space/dn,nxc,nzc,nxt,nzt,ni,nm
         real*4 dn
         integer nxc,nzc,nxt(2),nzt(2),ni,nm
      common/geom_topo/topo,topo_type,topo_alt,topo_zero,topo_file
         integer topo_type
         real*4 topo_alt,topo_zero,topo(3000)
         character*20 topo_file

      common/geom_water/water,water_type,water_vit,water_file
         integer water_type
         real*4 water_vit,water(3000)
         character*20 water_file
c234567
      common/ckernel/idatakernel,tk,tk0,dm_res
      integer idatakernel
      real*4 tk(10000),tk0(10000)
      real*4 dm_res
      common/mat_cstr/tabcstr,ntabcstr
         real*4 tabcstr(3,100)
         integer ntabcstr
      common/mat_int/matgeninf,tabinf,impact
      real*4 matgeninf(1002501),tabinf(3000)
      real*4 impact(300,2)
c234567
      common/matinv/gcdg,ndim,cm,diag,fn
      real*4 gcdg(1600,1600),cm(1600,1600),fn(1600,1600)
      real*4 diag(1600)
      integer ndim
      common/mat_mask/maskgen,maskref,maskup,maskdw,mcrref
      integer maskgen(1002501)
      real*4 maskref(2007000)
      real*4 maskdw(1600),maskup(1600)
      real*4 mcrref(1600)
c234567
      common/mat_resol/mat_best_rel,mat_best_abs,perc
      real*4 perc,mat_best_rel(750),mat_best_abs(750)
c234567
      common/var_mig/migration
      logical migration
      common/pick/tabpic,dumpic,refpic,dumref
         real*4 tabpic(300,100,2),dumpic(300,100)
         real*4 refpic(300,100,2),dumref(300,100)
         logical flagdir(100),flagref(100)

c  max size for shot arrays : 100
c  max size for station arrays : 300
c  tabpic : real array containing the picking t1-t2
c     for the first arrivals (t1=-1 means no picking
c     and t1=0. only a t2 picking)
c  refpic : real array containing the picking t1-t2
c     for the reflected phases (t1=-1 means no picking
c     and t1=0. only a t2 picking)
c  flagdir : boolean array indicating if there is a
c     first arrival picking on a station concerning
c     a given shot.
c  flagref : boolean array indicating if there is a
c     reflected phase picking on a station concerning
c     a given shot.
      common/prec_apri/apri_type,apri_vit,apri_grd,apri_file,apri_bnd
         integer apri_type
         real*4 apri_vit,apri_grd,apri_bnd
         character*20 apri_file
      common/prec_cstr/cstr_type,cstr_min,cstr_max,cstr_file
         integer cstr_type
         real*4 cstr_min,cstr_max
         character*20 cstr_file
      common/prec_mode/mode,mode_par,freeze,mode_file,fres
         integer freeze,mode,fres
         real*4 mode_par
         character*20 mode_file
      common/prec_prob/prob_type,prob_par
         integer prob_type,prob_par
      common/prec_resol/dm_res
         real*4 dm_res
      common/prec_search/search_type,search_par
         integer search_type,search_par
      common/prec_stop/stop_par,stop_type
         integer stop_type
         real*8 stop_par
      common/prec_verb/verb_type,verb_file
         integer verb_type
         character*20 verb_file
c234567
      common/mat_resolv/vecdat,tabidx,ndata
      real*4 vecdat(60000)
      integer tabidx(300,100,2),ndata
      common/shot/nsho,locsho
         integer nsho
         real*4 locsho(100,3)
      common/simplex/mat_sim,mat_mis,iter,ilo,ihi,inhi
         real*4 mat_sim(1601,1600)
         real*8 mat_mis(1600)
         integer iter,ilo,ihi,inhi
      common/station/nsta,locsta
         integer nsta
         real*4 locsta(300,3)
