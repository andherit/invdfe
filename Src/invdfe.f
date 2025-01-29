c234567
c  programme invdfe (inversion par differences
c    finies de l'Eikonal (podvin et lecomte 1991))
c
c  Version 2.5.2 (20/07/09)
c    Les dimensions maximales sont définies à travers des
c    statement parameter se trouvant dans maxval.com
c    passage des dimensions maximales pour les stations et les
c    shots a 1000.
c
c  Version 2.5.1 (15/07/04)
c    correction dans getprec afin de permettre a l'utilisateur
c    de ne pas mettre de graine dans les cas mode=start
c    et mode=joint. donc mode=start et mode=start,## sont
c    deux syntaxes acceptees.
c    29/10/2006 correction bug sur les dimensions (750 -> 1600)
c  Version 2.5.0 (19/01/04)
c    ajout du mode resolution calculant la matrice des
c    derivees partielles autour d'un modele donne. Ce
c    mode s'utilise comme le mode start. Le modele de
c    reference est introduit a travers la ligne de commande
c    apriori du fichier precision.
c    le calcul de la resolution introduite a moitie dans
c    la version 2.2.0 est definitivement abandonnee.
c
c  Version 2.4.1 (03/12/03)
c    correction d'un bug sur le back raytracing.
c
c  Version 2.4.0 (24/02/03)
c    Ajout d'une option water afin d'ajouter un masque sous
c    la topo de vitesse constante et d'epaisseur variable.
c
c  Version 2.3.1 (23/11/00)
c    Corrections de bugs. Fonction misfit mise a l'infini
c    si l'interface n'est pas totalement comprise entre la
c    topographie et le fond du modele ou si le point d'impact
c    de la reflection se situe sur un des bords du modele.
c    Bug dans comp_matvit2 pour le cas conjoint a vitesse
c    non constante. Le nombre de cp du milieu dw etait egal
c    a ceux du milieu up.
c
c  Version 2.3.0 (01/06/00)
c    Ajout d'une option de migration de l'interface dans
c    un milieu connu. Il s'agit d'un cas particulier du
c    mode joint. Ce cas est active quand nxt+nzt=0.
c    le milieu de reference est introduit par l'option
c    apriori dans le fichier precision.
c    implique une modification dans : getdata,compmisfit,
c    comp_matvit2
c
c  Version 2.2.1 (20/03/00)
c    modification du format de sortie du fichier hwd.dat
c    le programme sort maintenant un fichier par tir
c    (hwd??.dat). En plus, le format interne des fichiers
c    hwd passe de offset,temps a offset,profondeur relative,
c    a la topographie, temps
c
c  Version 2.2.0 (21/10/99)
c    introduction de la resolution par activation d'un
c    switch dans le fichier precision.
c  
c  Version 2.1.1 (22/09/99)
c    correction bugs
c
c  Version 2.1.0 (20/09/99)
c    decorrelation entre la discretisation des points
c    de controle entre milieux de vitesse et interface
c    nxt et nzt deviennent tableau a deux elements. le
c    nombre de points de controle de l'interface est
c    donne par ni (geometry.dat). introduction d'une
c    variable nm=mode+1 pour simplifier les tests entre
c    mode lisse et mode a interface
c
c  Version 2.0.1 (11/08/99)
c
c    modification de comp_misfit afin d'accelerer le
c    calcul. Le code gagne un facteur 6.
c
c  Version 2.0.0 (10/07/99)
c    introduction de l'inversion conjointe avec les
c    reflechies. Le milieu dans ce cas est discretise
c    par deux milieux de vitesses (au-dessus ou en dessous
c    de l'interface) + l'interface. le parametre mode devient
c      mode=start(0) inversion milieu lisse
c      mode=joint(1) inversion a interface
c      mode=continue(2) pour continuer une inversion quelque
c                       soit son type
c      mode=direct(3) forward modele en milieu lisse
c      mode=dirref(4) forward modele en milieu a interface
c   modification drastique de tout le code!
c
      program invdfe
      implicit none

      include 'maxval.com'
      include 'geom_space.com'
      include 'geom_water.com'
      include 'prec_cstr.com'
      include 'prec_apri.com'
      include 'prec_search.com'
      include 'prec_verb.com'
      include 'prec_mode.com'
      include 'mat_mask.com'
      include 'station.com'
      include 'shot.com'
      include 'simplex.com'

      character*50 filegeom,filedata,fileprec

      iter=0
c reading of the file names if different from default names
      call parse(filegeom,filedata,fileprec)
c reading of the main input file precision
      call getprec(fileprec)
c reading of the second input file geometry if the mode isn't "continue"
      if (mode.ne.2) call getgeom(filegeom)
c direct code branch
      if (mode.eq.3.or.mode.eq.4) then
         call get_dir
         call getdata(filedata)
         call comp_topo
         if (water_type.ne.0) call comp_water
         call pos_stasho
         call comp_dir
         stop
      endif
c inverse code branch
c reading secondary input files if mode isn't "continue"
      if (mode.lt.2..or.mode.eq.5) then
         call getdata(filedata)
         if (mode.ne.5) then
            include 'main_init.verb'
            include 'main_param.verb'
         endif
         call comp_topo
         if (water_type.ne.0) call comp_water
         call get_apri
         call pos_stasho
         call comp_mask
      endif
c resolution computation
      if (mode.eq.5) then
         call skernel
         stop
      endif
c non linear search by simplex
      if (search_type.eq.1) then
         if (mode.lt.2) call init_simplex
         if (mode.eq.2) then
            include 'main_cont.verb'
            call freeze_all(.false.)
         endif
         call searchbysimplex
         call dumpsimplex
      endif
      stop
      end

c###############################################################################
      subroutine skernel
      implicit none

      
      include 'maxval.com'
      include 'geom_space.com'
      include 'prec_verb.com'
      include 'prec_apri.com'
      include 'kernel.com'

      integer ndim,imodel,jmodel,idata,jdata,ik
c the kernel dimension is reduced.... for memory problem
      real*4 macrovit(PARMAX),kernel(NSHOMAX*NSTAMAX/16,PARMAX),dm
      real*8 vdum,comp_misfit,vdum0
      integer*4 nmodel,ndata
      real*4 rmodel,rdata,bascule
      equivalence (rmodel,nmodel),(rdata,ndata)

      dm=apri_bnd
      ndim=nzt(1)*nxt(1)
c calcul de mf0
      do jmodel=1,ndim
         macrovit(jmodel)=0.
      enddo
      vdum0=comp_misfit(macrovit)
      write(*,*) 'misfit value : ',vdum0
c     open(34,file='tk.dat')
      do idata=1,idatakernel
         tk0(idata)=tk(idata)
c        write(34,*) tk(idata)
      enddo
c     close(34)
c initialisation de la matrice derivee partielle
      do idata=1,idatakernel
         do jmodel=1,ndim
            kernel(idata,jmodel)=0
         enddo
      enddo
c calcul de la matrice derivee partielle
      do imodel=1,ndim
         write(*,'(i3.3,a1,i3.3,$)') imodel,'/',ndim
         do jmodel=1,ndim
            macrovit(jmodel)=0.
         enddo
         macrovit(imodel)=dm
         vdum=comp_misfit(macrovit)
         do jdata=1,idatakernel
            kernel(jdata,imodel)=(tk(jdata)-tk0(jdata))/dm
         enddo
         if (imodel.ne.ndim) then
         write(*,'(7a,$)') char(8),char(8),char(8),char(8),
     $char(8),char(8),char(8)
         endif
      enddo
      write(*,*)
      nmodel=ndim
      ndata=idatakernel
      open(10,file='kernel.bin',form='unformatted',access='direct',
     &recl=4)
      write(10,rec=1) rdata
      write(10,rec=2) rmodel
      ik=3
      do idata=1,idatakernel
         do jmodel=1,ndim
           write(10,rec=ik) kernel(idata,jmodel)*1000.
           ik=ik+1
         enddo
      enddo
      close(10)
      write(*,*)
      write(*,*) 'nmodel = ',ndim
      write(*,*) 'ndata  = ',idatakernel
c creation of the checker board test
      bascule=-dm
      do jmodel=1,ndim
         if (mod(jmodel,nxt(1)).ne.1) bascule=-bascule
         macrovit(jmodel)=bascule
         write(*,*) jmodel,macrovit(jmodel),mod(jmodel,nxt(1))
      enddo
      return
      end
         

      
c###############################################################################
      subroutine get_dir
      implicit none

      include 'maxval.com'
      include 'geom_space.com'
      include 'mat_mask.com'
      include 'prec_mode.com'
      include 'prec_verb.com'

      real*4 ztop,v,infinity
      integer i,j,idx,nx,nz,lnblnk,n
  
      infinity=1.e+06
      nx=nxc+1
      nz=nzc+1
      if (mode.eq.3) n=nx*nz
      if (mode.eq.4) n=2*nx*nz+nx
      open(30,file=mode_file(1:lnblnk(mode_file)),form='unformatted',
     &access='direct',recl=4,status='old',err=1000)
      do 10 i=1,n
         read(30,rec=i) maskref(i)
         if (i.le.2*nx*nz) then
            if (maskref(i).eq.0.) then
               maskref(i)=infinity
            else
               maskref(i)=dn/maskref(i)
            endif
         endif
10    continue
      close(30)
      return
1000  call msgexit('getdir',6)
      end

c###############################################################################
      subroutine mask_water(mat_vit)
      implicit none

      include 'maxval.com'
      real*4 mat_vit(VITMAX)

      include 'geom_space.com'
      include 'geom_water.com'
      include 'geom_topo.com'

      integer i,j,k,idx,nw
      real*4 x,y,v

      do 10 i=1,nxc
         nw=int(amin1(water(i),water(i+1))/dn)
         if (nw.eq.0) goto 10
         j=1
30       idx=i+(j-1)*(nxc+1)
         if (mat_vit(idx).lt.1.e+06) then
            do 20 k=1,nw
               idx=i+(k+j-2)*(nxc+1)
               mat_vit(idx)=dn/water_vit
20          continue
         else
            j=j+1
            if (j.le.nzc) goto 30
         endif
10    continue
      return
      end

c###############################################################################
      subroutine comp_dir
      implicit none

      include 'maxval.com'
      include 'geom_space.com'
      include 'geom_water.com'
      include 'geom_topo.com'
      include 'mat_mask.com'
      include 'simplex.com'
      include 'shot.com'
      include 'station.com'
      include 'picking.com'
      include 'prec_mode.com'
      include 'mat_int.com'

      integer ndim,i,j,idx,nx,nz
      integer status,time_2d_
      real*4 vx,vy,vv,poshotx,poshoty,time_ref
      real*4 mat_vit(VITMAX),timeonref(NXCMAX,2)
      real*4 mat_vitr(VITMAX),dummy
      real*4 mat_time(VITMAX),vect_hdw(NSHOMAX)
      real*4 mat_timer(VITMAX)
      character*11 filet

      nx=nxc+1
      nz=nzc+1
      filet='ref    .dat'
      do 10 i=1,nx*nz
         mat_vit(i)=maskref(i)
         if (mode.eq.4) mat_vitr(i)=maskref(i+nx*nz)
10    continue
      if (mode.eq.4) then
         do 20 i=1,nx
            tabinf(i)=maskref(i+2*nx*nz)
20       continue
      endif
      call mask_water(mat_vit)
      call dump_vit(mat_vit)
      do 60 i=1,nsho
         call init_time(mat_time)
         poshotx=locsho(i,1)/dn+1.
         poshoty=locsho(i,2)/dn+1.
         status=time_2d_(mat_vit,mat_time,
     &nx,nz,poshotx,poshoty,0.,0)
         if (status.ne.0) call msgexit('time_2d problem',15)
         call estim_time(mat_time,mat_vit,vect_hdw)
         do 70 j=1,nsta
            dumpic(j,i)=vect_hdw(j)
70       continue
         call dump_ray(mat_time,i)
c         call dump_time(mat_time,i)
         if (mode.eq.4) then
            call init_time(mat_timer)
            poshotx=locsho(i,1)/dn+1.
            poshoty=locsho(i,2)/dn+1.
            status=time_2d_(mat_vitr,mat_timer,
     &nx,nz,poshotx,poshoty,0.,0)
            if (status.ne.0) call msgexit('time_2d',7)
            dummy=time_ref(mat_timer,timeonref,0)
            write(filet(4:7),'(i4.4)') i
            open(60,file=filet)
            do 80 j=1,nsta
               call init_time(mat_time)
               poshotx=locsta(j,1)/dn+1.
               poshoty=locsta(j,2)/dn+1.
               status=time_2d_(mat_vitr,mat_time,
     &nx,nz,poshotx,poshoty,0.,0)
               if (status.ne.0) call msgexit('time_2d',7)
               dumref(j,i)=time_ref(mat_time,timeonref,j)
               if (dumref(j,i).gt.0..and.refpic(j,i,1).ge.0.)
     &call dump_ray_r(mat_timer,mat_time,i,j)
80          continue
            close(60)
         endif
60    continue
      call dump_hdw
      return
      end

c###############################################################################
      subroutine dumpsimplex
      implicit none

      include 'maxval.com'
      include 'geom_space.com'
      include 'geom_topo.com'
      include 'mat_mask.com'
      include 'mat_int.com'
      include 'simplex.com'
      include 'shot.com'
      include 'station.com'
      include 'picking.com'
      include 'prec_mode.com'
      include 'migration.com'

      integer ndim,i,j,idx,nx,nz,k
      integer status,time_2d_
      real*4 vx,vy,vv,poshotx,poshoty
      real*4 v1,v2,time_ref
      real*4 mat_vit(VITMAX),timeonref(NXCMAX,2)
      real*4 mat_vitr(VITMAX),macrovit(PARMAX)
      real*4 mat_time(VITMAX),vect_hdw(NSHOMAX)
      real*4 mat_timer(VITMAX),dummy
      character*11 filet

      ndim=nxt(1)*nzt(1)
      if (mode.eq.1) ndim=ndim+nxt(2)*nzt(2)+ni
      nx=nxc+1
      nz=nzc+1
      open(10,file='bestmodel.dat')
      if (.not.migration) then
         do 26 k=1,nm
            do 10 i=1,nxt(k)
               vx=(float(i)-.5)*(nxc*dn/float(nxt(k)))
               do 20 j=1,nzt(k)
               vy=topo_zero-(float(j)-.5)*(nzc*dn/float(nzt(k)))
                  idx=i+(j-1)*nxt(k)
                  if (k.eq.2) idx=idx+nxt(1)*nzt(1)
                  v1=mat_sim(ilo,idx)+mcrref(idx)
                  v2=mat_sim(ilo,idx)
                  write(10,*) vx,vy,v1,v2
20             continue
10          continue
26       continue
      endif
      if (nm.eq.2) then
         do 25 i=1,ni
            vx=(float(i)-.5)*(nxc*dn/float(ni))
            idx=nxt(1)*nzt(1)+nxt(2)*nzt(2)+i
            v1=mat_sim(ilo,idx)+mcrref(idx)
            v2=mat_sim(ilo,idx)
            write(10,*) vx,topo_zero-v1,v2
25       continue
      endif
      close(10)
      do 30 i=1,ndim
         macrovit(i)=mat_sim(ilo,i)
30    continue
      if (nm.eq.2) call comp_int(macrovit)
      call comp_matvit2(macrovit,mat_vit,mat_vitr,0)
      open(10,file='matv.bin',form='unformatted',
     &access='direct',recl=4)
      do 50 i=1,nx*nz
         vv=dn/mat_vit(i)
         if (vv.lt.1.) vv=0.
         write(10,rec=i) vv
50    continue
      if (mode.eq.1) then
         do 55 i=1,nx*nz
            idx=i+nx*nz
            vv=dn/mat_vitr(i)
            if (vv.lt.1.) vv=0.
            write(10,rec=idx) vv
55       continue
         do 56 i=1,nx
            idx=i+2*nx*nz
            write(10,rec=idx) tabinf(i)
56       continue
      endif
      close(10)
      call comp_matvit2(macrovit,mat_vit,mat_vitr,1)
      call dump_vit(mat_vit)
      write(*,*) 'best misfit : ',mat_mis(ilo)     
      write(*,*) 'worst misfit : ',mat_mis(ihi)     
      do 60 i=1,nsho
         if (.not.migration) then
            call init_time(mat_time)
           poshotx=locsho(i,1)/dn+1.
           poshoty=locsho(i,2)/dn+1.
           status=time_2d_(mat_vit,mat_time,
     &nx,nz,poshotx,poshoty,0.,0)
           if (status.ne.0) call msgexit('time_2d problem',15)
           call estim_time(mat_time,mat_vit,vect_hdw)
           do 70 j=1,nsta
              dumpic(j,i)=vect_hdw(j)
70         continue
           call dump_ray(mat_time,i)
c          call dump_time(mat_time,i)
         endif
         if (mode.eq.1) then
            call init_time(mat_timer)
            poshotx=locsho(i,1)/dn+1.
            poshoty=locsho(i,2)/dn+1.
            status=time_2d_(mat_vitr,mat_timer,
     &nx,nz,poshotx,poshoty,0.,0)
            if (status.ne.0) call msgexit('time_2d',7)
            dummy=time_ref(mat_timer,timeonref,0)
            filet='ref    .dat'
            write(filet(4:7),'(i4.4)') i
            open(60,file=filet)
            do 80 j=1,nsta
               call init_time(mat_time)
               poshotx=locsta(j,1)/dn+1.
               poshoty=locsta(j,2)/dn+1.
               status=time_2d_(mat_vitr,mat_time,
     &nx,nz,poshotx,poshoty,0.,0)
               if (status.ne.0) call msgexit('time_2d',7)
               dumref(j,i)=time_ref(mat_time,timeonref,j)
               if (dumref(j,i).gt.0..and.refpic(j,i,1).ge.0.)
     & call dump_ray_r(mat_timer,mat_time,i,j)
80          continue
            close(60)
         endif
60    continue
      call dump_hdw
      return
      end

c###############################################################################
      subroutine dump_time(mat_time,ish)
      implicit none

      include 'maxval.com'
      real*4 mat_time(VITMAX)
      integer ish

      include 'geom_space.com'
      include 'geom_topo.com'

      integer i,j,nx,nz
      real*4 x,y,v
      character*12 filet

      nx=nxc+1
      nz=nzc+1
      filet='time    .xyz'
      write(filet(5:8),'(i4.4)') ish
      open(50,file=filet)
      do 10 i=1,nx
         x=(i-1)*dn
         do 20 j=1,nz
            y=(j-1)*dn
            v=mat_time(i+(j-1)*nx)
            write(50,*) x,topo_zero-y,v
20       continue
10    continue
      close(50)
      return
      end

c###############################################################################
      subroutine dump_hdw
      implicit none

      include 'maxval.com'
      include 'shot.com'
      include 'station.com'
      include 'picking.com'
      include 'prec_mode.com'

      integer i,j,npic,km,lkm
      integer tpic(300),vtp
      real*4 va
      logical vswap,b1,b2,b3
      character*11 filed

      if (mode.eq.0.or.mode.eq.3) lkm=1
      if (mode.eq.1.or.mode.eq.4) lkm=2
      do 100 km=1,lkm
      if (km.eq.1) then
         filed='hwd????.dat'
      else
         filed='rwd????.dat'
      endif
      do 10 i=1,nsho
         npic=0
         write(filed(4:7),'(i4.4)') i
         open(10,file=filed)
         do 20 j=1,nsta
            if (km.eq.1) then
               va=tabpic(j,i,1)
            else
               va=refpic(j,i,1)
            endif
            if (va.ne.-1.) then
               npic=npic+1
               tpic(npic)=j
            endif
20       continue
         if (npic.eq.0) goto 10
         if (npic.gt.1) then
40          vswap=.false.
            do 30 j=1,npic-1
               b1=(locsta(tpic(j+1),1).lt.locsta(tpic(j),1))
               b2=(locsta(tpic(j+1),1).eq.locsta(tpic(j),1))
               b3=(locsta(tpic(j+1),2).lt.locsta(tpic(j),2))
               if (b1.or.(b2.and.b3)) then
                  vtp=tpic(j+1)
                  tpic(j+1)=tpic(j)
                  tpic(j)=vtp
                  vswap=.true.
               endif
30          continue
            if (vswap) goto 40
         endif
c le test est introduit pour securite en cas ou le  
c modele ne permet pas l'estimation d'une reflechie
c la ou elle est observee
         do 50 j=1,npic
            if (km.eq.1) then
               va=dumpic(tpic(j),i)
            else
               va=dumref(tpic(j),i)
            endif
            if (va.ge.0.) then
      write(10,*) locsta(tpic(j),1),locsta(tpic(j),3),va
            endif
50       continue
         close(10)
10    continue
100   continue
      close(10)
      return
      end

c###############################################################################
      subroutine dump_vit(mat_vit)
      implicit none

      include 'maxval.com'
      real*4 mat_vit(VITMAX)

      include 'geom_space.com'
      include 'geom_topo.com'

      integer i,j,idx
      real*4 x,y,v

      open(50,file='matv.xyz')
      do 10 i=1,nxc
         x=(float(i)-.5)*dn
         do 20 j=1,nzc
            y=(float(j)-.5)*dn
            idx=i+(j-1)*(nxc+1)
            if (mat_vit(idx).gt.1.e+06) then
               v=0.
            else
               v=dn/mat_vit(idx)
            endif
            write(50,*) x,topo_zero-y,v
20       continue
10    continue
      close(50)
      return
      end

c###############################################################################
      subroutine dump_ray(mat_time,ish)
      implicit none

      include 'maxval.com'
      real*4 mat_time(VITMAX)
      integer ish

      include 'geom_space.com'
      include 'geom_topo.com'
      include 'shot.com'
      include 'station.com'
      include 'picking.com'

      integer ks,xcell,zcell
      integer idx1,idx2,idx3,idx4,idum,icrt
      real*4 pas,xcrt,zcrt,xgrd,zgrd,vmod,wx,wz
      real*4 vclose
      character*11 filet

      pas=dn/2.
      vclose=1.5*dn
      filet='ray    .dat'
      write(filet(4:7),'(i4.4)') ish
      open(10,file=filet)
      do 2000 ks=1,nsta
         if (tabpic(ks,ish,1).eq.-1.) goto 2000
         write(10,'(a2,i3)') '> ',ks
         xcrt=locsta(ks,1)
         zcrt=locsta(ks,2)
         idum=0
         icrt=0
10       write(10,*) xcrt,topo_zero-zcrt
         xcell=int(xcrt/dn)+1
         zcell=int(zcrt/dn)+1
         idx1=xcell+(zcell-1)*(nxc+1)
         idx2=xcell+(zcell-1)*(nxc+1)+1
         idx3=xcell+zcell*(nxc+1)
         idx4=xcell+zcell*(nxc+1)+1
         xgrd=mat_time(idx4)+mat_time(idx2)
         xgrd=xgrd-mat_time(idx1)-mat_time(idx3)
         zgrd=mat_time(idx4)+mat_time(idx3)
         zgrd=zgrd-mat_time(idx1)-mat_time(idx2)
         vmod=sqrt(xgrd**2.+zgrd**2.)
         xcrt=xcrt-pas*xgrd/vmod
         zcrt=zcrt-pas*zgrd/vmod
         if (abs(xcrt-locsho(ish,1)).lt.vclose.
     &and.abs(zcrt-locsho(ish,2)).lt.vclose) then
            write(10,*) locsho(ish,1),topo_zero-locsho(ish,2)
            goto 2000
         else
            goto 10
         endif
2000  continue
      close(10)
      return
      end


c###############################################################################
      subroutine pos_stasho
      implicit none

      include 'maxval.com'
      include 'geom_space.com'
      include 'geom_topo.com'
      include 'station.com'
      include 'shot.com'

      integer i,cari
      real*4 ztop

      do 10 i=1,nsta
         cari=int(locsta(i,1)/dn)+1
         ztop=(locsta(i,1)/dn-(cari-1))*
     &(topo(cari+1)-topo(cari))+topo(cari)
         locsta(i,2)=ztop+locsta(i,3)
10    continue
      do 20 i=1,nsho
         cari=int(locsho(i,1)/dn)+1
         ztop=(locsho(i,1)/dn-(cari-1))*
     &(topo(cari+1)-topo(cari))+topo(cari)
         locsho(i,2)=ztop+locsho(i,3)
20    continue
      return
      end
  

c###############################################################################
      function comp_misfit(macrovit)
      implicit none

      include 'maxval.com'
      real*4 macrovit(PARMAX)
      integer type

      include 'geom_space.com'
      include 'geom_topo.com'
      include 'shot.com'
      include 'station.com'
      include 'prec_prob.com'
      include 'prec_mode.com'
      include 'prec_verb.com'
      include 'mat_int.com'
      include 'migration.com'
      include 'kernel.com'

      real*4 mat_vit(VITMAX),mat_time(VITMAX)
      real*4 mat_vitr(VITMAX)
      integer nx,nz,i,status,n_fit,t_fit,time_2d_,j,carj
      real*4 poshotx,poshoty,vect_hdw(NSTAMAX),vect_ref(NSTAMAX)
      real*4 timeonref(NXCMAX,NSTAMAX),estim_ref,dummy,x,z,y
      real*4 t1,t3,dy,timeonshot(NXCMAX)
      real*8 v_fit,comp_misfit,vmin,vmax
      real*8 s_fit
      logical verif_modin,verif_intin

      vmin=0.
      vmax=1.e+06
      if (mode.eq.5) idatakernel=1
      if (.not.verif_modin(macrovit).and.mode.ne.5) then
         if (prob_type.eq.1) comp_misfit=vmax
         if (prob_type.gt.1) comp_misfit=vmin
         return
      endif
      if (mode.eq.1) then 
         call comp_int(macrovit)
         if (.not.verif_intin()) then
            if (prob_type.eq.1) comp_misfit=vmax
            if (prob_type.gt.1) comp_misfit=vmin
            return
         endif
      endif
      call comp_matvit2(macrovit,mat_vit,mat_vitr,1)
      call dump_vit(mat_vit)
      nx=nxc+1
      nz=nzc+1
      if (mode.eq.1) then
         do 10 j=1,nsta
            call init_time(mat_time)
            poshotx=locsta(j,1)/dn+1.
            poshoty=locsta(j,2)/dn+1.
            status=time_2d_(mat_vitr,mat_time,
     &nx,nz,poshotx,poshoty,0.,0)
            if (status.ne.0) call msgexit('time_2d problem',15)
            do 20 i=1,nx
               x=(i-1)*dn
               y=tabinf(i)
               carj=int(y/dn)+1
               t1=mat_time(i+carj*nx)
               t3=mat_time(i+(carj-1)*nx)
               dy=y-(carj-1)*dn
               timeonref(i,j)=dy/dn*(t1-t3)+t3
20          continue
10       continue
      endif
      t_fit=0
      s_fit=0.
      do 100 i=1,nsho
         if (.not.migration) then
            call init_time(mat_time)
            poshotx=locsho(i,1)/dn+1.
            poshoty=locsho(i,2)/dn+1.
            status=time_2d_(mat_vit,mat_time,
     &nx,nz,poshotx,poshoty,0.,0)
            if (status.ne.0) call msgexit('time_2d problem',15)
            call estim_time(mat_time,mat_vit,vect_hdw)
            call comp_delta(vect_hdw,i,v_fit,n_fit,0)
            s_fit=s_fit+v_fit
            t_fit=t_fit+n_fit
         endif
         if (mode.eq.1) then
            call init_time(mat_time)
            poshotx=locsho(i,1)/dn+1.
            poshoty=locsho(i,2)/dn+1.
            status=time_2d_(mat_vitr,mat_time,
     &nx,nz,poshotx,poshoty,0.,0)
            if (status.ne.0) call msgexit('time_2d problem',15)
            do 30 j=1,nx
               x=(j-1)*dn
               y=tabinf(j)
               carj=int(y/dn)+1
               t1=mat_time(j+carj*nx)
               t3=mat_time(j+(carj-1)*nx)
               dy=y-(carj-1)*dn
               timeonshot(j)=dy/dn*(t1-t3)+t3
30          continue
            do 40 j=1,nsta
               vect_ref(j)=estim_ref(timeonref,timeonshot,j)
               if (vect_ref(j).lt.0.) then
                  if (prob_type.eq.1) comp_misfit=vmax
                  if (prob_type.gt.1) comp_misfit=vmin
                  return
               endif
40          continue
            call comp_delta(vect_ref,i,v_fit,n_fit,1)
            s_fit=s_fit+v_fit
            t_fit=t_fit+n_fit
         endif
100   continue
      comp_misfit=s_fit/dble(t_fit)
      if (mode.eq.5) idatakernel=idatakernel-1
      return
      end
      
c###############################################################################
      subroutine getprec(filep)
      implicit none

c     passing variables
         character*50 filep
      include 'maxval.com'
      include 'prec_mode.com'
      include 'prec_search.com'
      include 'prec_stop.com'
      include 'prec_apri.com'
      include 'prec_cstr.com'
      include 'prec_prob.com'
      include 'prec_verb.com'
      include 'kernel.com'

c     local variables

      character*20 keyword
      character*100 line,line2
      logical syntax
      integer k,ll,npar,lk,lpar(5),idx,lnblnk

c default parameter

      mode=0
      mode_par=0.
      search_type=1
      search_par=1
      stop_type=1
      stop_par=3000.
      apri_type=0
      apri_bnd=10.
      cstr_type=0
      prob_type=1
      prob_par=2
      freeze=250
      verb_type=1
      verb_file='spy'

      open(10,file=filep(1:lnblnk(filep)),status='old',err=1000)
10    read(10,'(a100)',end = 100) line
      if (line(1:1).eq.'#') goto 10
      if (.not.syntax(line,keyword,npar,line2,lk,lpar)) goto 2000
c mode parameter
      if (keyword(1:lk).eq."mode") then
         if (npar.lt.1.or.npar.gt.3) goto 2000
         mode=20
         if (line2(1:lpar(1)).eq."start") mode=0
         if (line2(1:lpar(1)).eq."joint") mode=1
         if (line2(1:lpar(1)).eq."continue") mode=2
         if (line2(1:lpar(1)).eq."direct") mode=3
         if (line2(1:lpar(1)).eq."dirref") mode=4
         if (line2(1:lpar(1)).eq."resolution") mode=5
         if (mode.eq.20) goto 2000
         if (npar.gt.1.and.mode.eq.2) goto 2000
         if (npar.eq.1.and.mode.eq.3) goto 2000
         if (npar.eq.1.and.mode.eq.4) goto 2000
         if (npar.ne.1.and.mode.eq.5) goto 2000
         if (npar.ge.2) then
            mode_file(1:lpar(2))=line2(21:20+lpar(2))
         endif
         if (npar.eq.2.and.mode.lt.2) 
     &read(mode_file(1:lpar(2)),*,err=2000) mode_par
         goto 10
      endif
c search parameter
      if (keyword(1:lk).eq."search") then
         if (npar.ne.2) goto 2000
         if (line2(1:lpar(1)).eq."simplex") then
            search_type=1
         else
            goto 2000
         endif
         search_par=0
         if (line2(21:20+lpar(2)).eq."random") search_par=1
         if (line2(21:20+lpar(2)).eq."auto") search_par=2
         if (search_par.eq.0) goto 2000
         goto 10
      endif
c stop parameter
      if (keyword(1:lk).eq."stop") then
         if (npar.ne.2) goto 2000
         stop_type=0
         if (line2(1:lpar(1)).eq."iteration") stop_type=1
         if (line2(1:lpar(1)).eq."precision") stop_type=2
         if (stop_type.eq.0) goto 2000
         read(line2(21:20+lpar(2)),*,err=2000) stop_par
         goto 10
      endif
c a priori model
      if (keyword(1:lk).eq.'apriori') then
         apri_type=10
         if (line2(1:lpar(1)).eq.'absolute') apri_type=0
         if (line2(1:lpar(1)).eq.'file') apri_type=1
         if (apri_type.eq.10) goto 2000
         if (apri_type.eq.0.and.npar.ne.1) goto 2000
         if (apri_type.eq.1.and.npar.ne.3) goto 2000
         if (apri_type.eq.1) then
            apri_file(1:lpar(2))=line2(21:20+lpar(2))
            read(line2(41:40+lpar(3)),*,err=2000) apri_bnd
         endif
         goto 10
      endif
c cstr parameter
      if (keyword(1:lk).eq.'cstr') then
         cstr_type=10
         if (line2(1:lpar(1)).eq.'none') cstr_type=0
         if (line2(1:lpar(1)).eq.'band') cstr_type=1
         if (cstr_type.eq.10) goto 2000
         if (cstr_type.eq.0.and.npar.ne.1) goto 2000
         if (cstr_type.eq.1.and.npar.ne.3) goto 2000
         if (cstr_type.eq.1) then
            read(line2(21:20+lpar(2)),*,err=2000) cstr_min
            read(line2(41:40+lpar(3)),*,err=2000) cstr_max
         endif
         goto 10
      endif
c prob parameter
      if (keyword(1:lk).eq.'prob') then
         if (npar.ne.2) goto 2000
         prob_type=0
         if (line2(1:lpar(1)).eq.'distance') prob_type=1
         if (line2(1:lpar(1)).eq.'gaussian') prob_type=2
         if (line2(1:lpar(1)).eq.'chi2') prob_type=3
         if (prob_type.eq.0) goto 2000
         read(line2(21:20+lpar(2)),*,err=2000) prob_par
         goto 10
      endif
c freeze parameter
      if (keyword(1:lk).eq."freeze") then
         if (npar.ne.1) goto 2000
         read(line2(1:lpar(1)),*,err=2000) freeze
         goto 10
      endif
c verbose parameter
      if (keyword(1:lk).eq."verbose") then
         verb_type=10
         if (line2(1:lpar(1)).eq."none") verb_type=0
         if (line2(1:lpar(1)).eq."minimum") verb_type=1
         if (line2(1:lpar(1)).eq."debug") verb_type=2
         if (verb_type.eq.10) goto 2000
         if (verb_type.eq.0.and.npar.ne.1) goto 2000
         if (verb_type.eq.1.and.npar.ne.2) goto 2000
         if (verb_type.eq.2.and.npar.ne.2) goto 2000
         if (verb_type.ne.0) verb_file(1:lpar(2))=line2(21:20+lpar(2))
         goto 10
      endif
      call msgexit('precision file: unknown keyword '
     &//keyword(1:lk),32+lk) 
100   close(10)
      return
1000  call msgexit('unknown precision file',22)
2000  call msgexit('precision file: syntax error',28)
      end
c###############################################################################
c  routine de lecture de la table de picking
c###############################################################################
      subroutine getdata(filedata)
      implicit none 

c     passing variables
         character*50 filedata

      include 'maxval.com'
      include 'shot.com'
      include 'station.com'
      include 'picking.com'
      include 'prec_mode.com'
      include 'migration.com'

      character*(NSHOMAX*10) line
      character*3 sdum
      integer i,j,k,lnblnk,lf,km,lkm,dkm
      real*4 dum(2*NSHOMAX+2)

      lf=lnblnk(filedata)
      dkm=1
      if (mode.eq.0.or.mode.eq.3.or.mode.eq.5) lkm=1
      if (mode.eq.1.or.mode.eq.4) lkm=2
      if (migration) dkm=2
      do 100 km=dkm,lkm
         if (km.eq.1) then
            open(10,file=filedata(1:lf),status='old'
     &,err=1000)
         else
            sdum=filedata(lf-2:lf)
            filedata(lf-2:lf)='ref'
            open(10,file=filedata(1:lf),status='old'
     &,err=1000)
            filedata(lf-2:lf)=sdum(1:3)
         endif
5        read(10,'(a)',err=2000) line
         if (line(1:1).eq.'#') goto 5
         if (km.eq.1.or.migration) then
            read(line,*) nsta,nsho
            read(line,*) (dum(i),i=1,(nsho+1)*2)
            do 10 i=1,nsho
               locsho(i,1)=dum((i+1)*2-1)
               locsho(i,3)=dum((i+1)*2)
10          continue
         endif
         do 20 i=1,nsta
            read(10,*) (dum(k),k=1,(nsho+1)*2)
            if (km.eq.1.or.migration) then
               locsta(i,1)=dum(1)
               locsta(i,3)=dum(2)
            endif
            if (km.eq.1) then
               do 30 j=1,nsho
                  tabpic(i,j,1)=dum((j+1)*2-1)
                  tabpic(i,j,2)=dum((j+1)*2)
30             continue
            else
               do 40 j=1,nsho
                  refpic(i,j,1)=dum((j+1)*2-1)
                  refpic(i,j,2)=dum((j+1)*2)
40             continue
            endif
20       continue
         close(10)
100   continue
      return
1000  call msgexit('data unknown file',17)
2000  call msgexit('data syntax error',17)
      end

c###############################################################################
c  routine d'association des parametres et leur valeur pour la geometrie
c###############################################################################
      subroutine getgeom(filegeom)
      implicit none

c     passing variables
         character*50 filegeom

      include 'maxval.com'
      include 'geom_space.com'
      include 'geom_topo.com'
      include 'geom_water.com'
      include 'prec_mode.com'
      include 'migration.com'

      character*100 line,line2
      character*20 keyword
      integer npar,lk,lpar(5),lnblnk
      logical syntax

c default parameters

      dn=50.
      nxc=300
      nzc=200
      nxt(1)=1
      nzt(1)=1
      nxt(2)=0
      nzt(2)=0
      ni=1
      if (mode.eq.1.or.mode.eq.4) then
         nm=2
      else
         nm=1
      endif
      topo_type=0
      topo_alt=0
      topo_zero=0
      water_type=0
      water_vit=0.
 
      open(10,file=filegeom(1:lnblnk(filegeom)),status='old',err=1000)
10    read(10,'(a100)',end=100) line
      if (line(1:1).eq.'#') goto 10
      if (.not.syntax(line,keyword,npar,line2,lk,lpar)) goto 2000
c dnum parameter
      if (keyword(1:lk).eq.'dnum') then
         if (npar.ne.1) goto 2000
         read(line2(1:lpar(1)),*,err=2000) dn
         goto 10
      endif
      if (keyword(1:lk).eq.'nxt') then
         if (npar.ne.1.and.nm.eq.1) goto 2000
         if (npar.ne.2.and.nm.eq.2) goto 2000
         read(line2(1:lpar(1)),*,err=2000) nxt(1)
         if (nm.eq.2.or.nm.eq.5) read(line2(21:20+lpar(2)),
     &*,err=2000) nxt(2)
         goto 10
      endif
      if (keyword(1:lk).eq.'nzt') then
         if (npar.ne.1.and.nm.eq.1) goto 2000
         if (npar.ne.2.and.nm.eq.2) goto 2000
         read(line2(1:lpar(1)),*,err=2000) nzt(1)
         if (nm.eq.2.or.nm.eq.5) read(line2(21:20+lpar(2))
     &,*,err=2000) nzt(2)
         goto 10
      endif
      if (keyword(1:lk).eq.'ni') then
         if (mode.ne.1) goto 2000
         if (npar.ne.1) goto 2000
         read(line2(1:lpar(1)),*,err=2000) ni
         goto 10
      endif
      if (keyword(1:lk).eq.'nxc') then
         if (npar.ne.1) goto 2000
         read(line2(1:lpar(1)),*,err=2000) nxc
         goto 10
      endif
      if (keyword(1:lk).eq.'nzc') then
         if (npar.ne.1) goto 2000
         read(line2(1:lpar(1)),*,err=2000) nzc
         goto 10
      endif
      if (keyword(1:lk).eq.'topo') then
         if (npar.ne.3) goto 2000
         topo_type=10
         if (line2(1:lpar(1)).eq.'flat') topo_type=0
         if (line2(1:lpar(1)).eq.'file') topo_type=1
         if (topo_type.eq.10) goto 2000
         if (topo_type.eq.0)
     & read(line2(21:20+lpar(2)),*,err=2000) topo_alt
         if (topo_type.eq.1) topo_file(1:lpar(2))=line2(21:20+lpar(2))
         read(line2(41:40+lpar(3)),*,err=2000) topo_zero
         goto 10
      endif
      if (keyword(1:lk).eq.'water') then
         water_type=10
         if (line2(1:lpar(1)).eq.'none') water_type=0
         if (line2(1:lpar(1)).eq.'file') water_type=1
         if (npar.ne.3.and.water_type.eq.1) goto 2000
         if (npar.ne.1.and.water_type.eq.0) goto 2000
         if (water_type.eq.1) then
            water_file(1:lpar(2))=line2(21:20+lpar(2))
            read(line2(41:40+lpar(3)),*,err=2000) water_vit
         endif
         goto 10
      endif
      call msgexit('geometry file: unknown parameter',32)
100   close(10)
      migration=(nxt(1)+nxt(2)+nzt(1)+nzt(2).eq.0)
      return
1000  call msgexit('unknown geometry file',21)
2000  call msgexit('geometry file: syntax error',27)
      end
      
c###############################################################################
c    routine d'association des fichiers d'input
c    Trois fichiers doivent etre fournis. Chaque fichier a un nom
c    par defaut qui sera pris si l'option n'est pas donnee.
c    I   fichier de la geometrie du probleme (par defaut: geometry.dat)
c            geom=....
c    II  fichier des donnees a inverser (par defaut: data.dat)
c            data=....
c    III fichier de la precision requise (par defaut: precision.dat)
c            prec=.... 
c    ex:
c            hostname%> invdfe
c    ou
c            hostname%> invdfe geom=flou.dat data=rien.dat prec=absolue.dat
c###############################################################################
      subroutine parse(filegeom,filedata,fileprec)
      implicit none

      integer iargc,n,i,k,lnblnk
      character*50 filegeom,filedata,fileprec
      character*100 line

      n=iargc()
      filegeom='geometry.dat'
      filedata='data.dat'
      fileprec='precision.dat'
      if (n.eq.0) return
      do 10 i=1,n
         call getarg(i,line)
         k=1
20       if (line(k:k).ne.'=') then
            k=k+1
            if (k.gt.lnblnk(line))
     & call msgexit('parse syntax error',18)
            goto 20
         endif
         if (line(1:k-1).eq.'geom') then
            filegeom=line(k+1:lnblnk(line))
            goto 10
         endif
         if (line(1:k-1).eq.'data') then
            filedata=line(k+1:lnblnk(line))
            goto 10
         endif
         if (line(1:k-1).eq.'prec') then
            fileprec=line(k+1:lnblnk(line))
            goto 10
         endif
         call msgexit('parse unknown option',20)
10    continue
      return
      end

c###############################################################################
c   routine d'erreur avec sortie du programme
c###############################################################################
      subroutine msgexit(message,ind)
      implicit none

      character*50 message
      integer ind

      write(*,'(a)') message(1:ind)
      stop 'program error exit'
      end

c###############################################################################
      function syntax(line,keyword,np,line2,lk,lpar)
      implicit none

      character*100 line
      character*100 line2
      character*20 keyword
      integer np,idx,i
      logical syntax
      integer k,ll,lk,kv,lnblnk,lpar(5)

      do 5 i=1,100
         line2(i:i)=''
5     continue
      syntax=.true.
      ll=lnblnk(line)
      k=1
10    if (line(k:k).ne.'=') then
         k=k+1
         if (k.gt.ll) then
            syntax=.false.
            return
         endif
         goto 10
      endif
      keyword(1:k-1)=line(1:k-1)
      lk=k-1
      np=0
30    np=np+1
      idx=(np-1)*20+1
      kv=k+1
20    if (line(kv:kv).ne.',') then
         kv=kv+1
         if (kv.gt.ll) then
            line2(idx:idx+kv-k-2)=line(k+1:kv-1)
            lpar(np)=kv-k-1
            return
         endif
         goto 20
      endif
      line2(idx:idx+kv-k-1)=line(k+1:kv-1)
      lpar(np)=kv-k-1
      k=kv
      goto 30
      end
c###############################################################################
      subroutine comp_water
      implicit none
      
      include 'maxval.com'
      include 'geom_space.com'
      include 'geom_water.com'
      include 'prec_verb.com'

      real*4 xtemp(10),ztemp(10)
      real*4 xwat(NXCMAX),ywat(NXCMAX)
      real*4 vx,vy,dy,vmax
      integer i,nt,j,k,nse,nx
      integer indmin(10),pmax


c
      open(50,file='spy',access='append')
      write(50,*) water_file
      write(50,*) water_vit
      write(50,*) water_type
c
      open(10,file=water_file,status='old',err=1000)
      i=1
20    read(10,*,end=100) xwat(i),ywat(i)
      i=i+1
      goto 20
100   close(10)
      nt=i-1
c
      write(50,*) nt
      close(50)
c
      nx=nxc+1
c
c degre de la polynomiale +1
c prevu pour etre parametre reglage de niveau 2 dans version ulterieure
c
      nse=5
c
c cas quelconque pour nt superieur a nse
c
      if (nt.gt.nse) then
      do 30 i=1,nx
         vx=(i-1)*dn
         vmax=0.
         do 40 j=1,nse
            indmin(j)=j
            if (abs(vx-xwat(indmin(j))).gt.vmax) then
               vmax=abs(vx-xwat(indmin(j)))
               pmax=j
            endif
40       continue
         do 50 j=nse+1,nt
            if (abs(vx-xwat(j)).lt.abs(vx-xwat(indmin(pmax)))) then
               indmin(pmax)=j
               vmax=0.
               do 60 k=1,nse
                  if (abs(vx-xwat(indmin(k))).gt.vmax) then
                     vmax=abs(vx-xwat(indmin(k)))
                     pmax=k
                  endif
60             continue
            endif
50       continue
         do 70 j=1,nse
            xtemp(j)=xwat(indmin(j))
            ztemp(j)=ywat(indmin(j))
70       continue
         call polint(xtemp,ztemp,nse,vx,vy,dy)
         water(i)=vy
         if (water(i).lt.0.) water(i)=0.
30    continue
      else
c
c cas quelconque pour nt inferieur a nse
c
         do 80 i=1,nx
            vx=(i-1)*dn
            call polint(xwat,ywat,nt,vx,vy,dy)
            water(i)=vy
            if (water(i).lt.0.) water(i)=0.
80       continue
      endif
      open(24,file='water.asc')
      do 500 i=1,nx
         vx=(i-1)*dn
         vy=water(i)
         write(24,*) vx,vy
500   continue
      close(24)
      return
1000  call msgexit('input water',11)
      end
c###############################################################################
      subroutine comp_topo
      implicit none
      
      include 'maxval.com'
      include 'geom_space.com'
      include 'geom_topo.com'
      include 'prec_verb.com'

      real*4 xtemp(10),ztemp(10)
      real*4 xtopo(NXCMAX),ytopo(NXCMAX)
      real*4 vx,vy,dy,vmax
      integer i,nt,j,k,nse,nx
      integer indmin(10),pmax

c
c  cas plan
c
      nx=nxc+1
      include 'comp_topo.verb'
      if (topo_type.eq.0) then
         do 10 i=1,nx
            topo(i)=topo_zero-topo_alt
            if (topo(i).lt.0.) call msgexit('topography',10)
10       continue
         return
      endif
      open(10,file=topo_file,status='old',err=1000)
      i=1
20    read(10,*,end=100) xtopo(i),ytopo(i)
      i=i+1
      goto 20
100   close(10)
      nt=i-1
c
c degre de la polynomiale +1
c prevu pour etre parametre reglage de niveau 2 dans version ulterieure
c
      nse=5
c
c cas quelconque pour nt superieur a nse
c
      if (nt.gt.nse) then
      do 30 i=1,nx
         vx=(i-1)*dn
         vmax=0.
         do 40 j=1,nse
            indmin(j)=j
            if (abs(vx-xtopo(indmin(j))).gt.vmax) then
               vmax=abs(vx-xtopo(indmin(j)))
               pmax=j
            endif
40       continue
         do 50 j=nse+1,nt
            if (abs(vx-xtopo(j)).lt.abs(vx-xtopo(indmin(pmax)))) then
               indmin(pmax)=j
               vmax=0.
               do 60 k=1,nse
                  if (abs(vx-xtopo(indmin(k))).gt.vmax) then
                     vmax=abs(vx-xtopo(indmin(k)))
                     pmax=k
                  endif
60             continue
            endif
50       continue
         do 70 j=1,nse
            xtemp(j)=xtopo(indmin(j))
            ztemp(j)=ytopo(indmin(j))
70       continue
         call polint(xtemp,ztemp,nse,vx,vy,dy)
         topo(i)=topo_zero-vy
         if (topo(i).lt.0.) then
            write(*,*) 'erreur de polint'
            write(*,*) 'pour le point :',vx
            write(*,*) 'altitude :',vy
            write(*,*) 'point de l interpolation :'
            do 200 j=1,nse
               write(*,*) xtemp(j),ztemp(j)
200         continue
            call msgexit('topography',10)
         endif
30    continue
      else
c
c cas quelconque pour nt inferieur a nse
c
         do 80 i=1,nx
            vx=(i-1)*dn
            call polint(xtopo,ytopo,nt,vx,vy,dy)
            topo(i)=topo_zero-vy
            if (topo(i).lt.0.) call msgexit('topography',10)
80       continue
      endif
      open(24,file='toto')
      do 500 i=1,nx
         vx=(i-1)*dn
         vy=topo_zero-topo(i)
         write(24,*) vx,vy
500   continue
      close(24)
      return
1000  call msgexit('topography',10)
      end
c###############################################################################
               
      subroutine comp_mask
      implicit none

      include 'maxval.com'
      include 'geom_space.com'
      include 'geom_topo.com'
      include 'geom_water.com'
      include 'mat_mask.com'
      include 'prec_mode.com'
      include 'prec_verb.com'
      include 'migration.com'

      integer nx,nz,i,j,idi,idj,idx,idij,k
      real*4 vmax,vmin,val,ztop,zprof,zwat

      nx=nxc+1
      nz=nzc+1
c creation of the topographic mask maskgen
      do 100 i=1,nxc
         ztop=amin1(topo(i),topo(i+1))
         zwat=amin1(water(i),water(i+1))
         do 200 j=1,nzc
            idx=i+(j-1)*nx
            zprof=j*dn
            if (zprof.le.ztop) then
               maskgen(idx)=0
            else
               if (zprof.le.ztop+zwat) then
                  maskgen(idx)=2
               else
                  maskgen(idx)=1
               endif
            endif
200      continue
         idx=i+nzc*nx
         maskgen(idx)=0
100   continue
      do 300 j=1,nz
         idx=j*nx
         maskgen(idx)=0
300   continue
c creation of masks on the interface
      if (nm.eq.2) then
         do 310 i=1,ni
           idi=int(((i-.5)*(nxc*dn/float(ni)))/dn)
           idx=nxt(1)*nzt(1)+nzt(2)*nxt(2)+i
           mcrref(idx)=maskref(idi+2*nz*nx)
           val=mcrref(idx)
           call estim_vit(vmin,vmax,val,idi)
           maskup(idx)=vmax
           maskdw(idx)=vmin
310     continue
      endif
      if (migration) return
c creation of masks for the medium
      do 700 k=1,nm
         do 500 i=1,nxt(k)
            idi=int(((i-.5)*(nxc*dn/float(nxt(k))))/dn)
            do 600 j=1,nzt(k)
               idj=int(((j-.5)*(nzc*dn/float(nzt(k))))/dn)
               idx=i+(j-1)*nxt(k)
               if (k.eq.2) idx=idx+nxt(1)*nzt(1)
               idij=idi+(idj-1)*nx
               if (k.eq.2) idij=idij+nx*nz
               mcrref(idx)=maskref(idij)
               val=mcrref(idx)
               call estim_vit(vmin,vmax,val,0)
               maskup(idx)=vmax
               maskdw(idx)=vmin
600         continue
500      continue
700   continue
      return
      end

c###############################################################################
      subroutine estim_vit(vmin,vmax,val,sw)
      implicit none

      real*4 vmin,vmax,val
      integer sw

      include 'maxval.com'
      include 'geom_space.com'
      include 'geom_topo.com'
      include 'prec_cstr.com'
      include 'prec_apri.com'

      real*4 tmin,tmax
      integer i,sct

      if (sw.ne.0) then
         sct=1
         tmin=amax1(topo(sw),topo(sw+1))
         tmax=nzc*dn
      else
         sct=cstr_type
         if (sct.eq.1) then
            tmin=cstr_min
            tmax=cstr_max
         endif
      endif
      if (sct.eq.0) then
         if (apri_type.eq.0)
     &        call msgexit('mask: no velocity constraints',29)
         vmin=-apri_bnd
         vmax=apri_bnd
         return
      endif
      if (apri_type.eq.0) then
         vmin=tmin
         vmax=tmax
         return
      endif
      if ((val-apri_bnd).gt.tmax.or.(val+apri_bnd).lt.tmin) then
           call msgexit('mask: constraint conflit',24)
      endif
      if ((val-apri_bnd).lt.tmin) then
         vmin=tmin-val
      else
         vmin=-apri_bnd
      endif
      if ((val+apri_bnd).gt.tmax) then
         vmax=tmax-val
      else
         vmax=apri_bnd
      endif
      return
      end

c###############################################################################
c234567
      subroutine searchbysimplex
      implicit none

      include 'maxval.com'
      include 'geom_space.com'
      include 'prec_stop.com'
      include 'prec_mode.com'
      include 'prec_verb.com'
      include 'simplex.com'

      integer ndim,n,m,i,j,imov,nbfree
      integer nbcarl,tymo
      real*4 sum,psum(PARMAX),rtol
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
      subroutine montecarlo(psum,fakir)
      implicit none

      include 'maxval.com'
      real*4 psum(PARMAX),fakir

      include 'geom_space.com'
      include 'mat_mask.com'
      include 'simplex.com'
      include 'prec_mode.com'
  
      integer ndim,i,j,nit
      real*4 sim_try(PARMAX),rdom,vmax,vmin,vmed
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
      mis_try=comp_misfit(sim_try)
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

      function amotry(psum,fac)
      implicit none  

      include 'prec_verb.com'

      include 'maxval.com'
      real*4 psum(PARMAX),fac
      real*8 amotry

      include 'geom_space.com'
      include 'simplex.com'
      include 'prec_mode.com'

      integer ndim,j
      real*4 fac1,fac2,sim_try(PARMAX)
      real*8 mis_try,comp_misfit

      ndim=nxt(1)*nzt(1)
      if (nm.eq.2) ndim=ndim+nxt(2)*nzt(2)+ni
      fac1=(1.-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        sim_try(j)=psum(j)*fac1-mat_sim(ihi,j)*fac2
11    continue
      mis_try=comp_misfit(sim_try)
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
c234567
      subroutine comp_delta(vect_hdw,ishot,v_fit,n_fit,sw)
      implicit none

      include 'maxval.com'
      real*4 vect_hdw(NSTAMAX)
      real*8 v_fit
      integer ishot,n_fit,sw

      include 'picking.com'
      include 'prec_prob.com'
      include 'prec_mode.com'
      include 'station.com'
      include 'kernel.com'

      integer i
      real*8 f
      real*4 d1,d2,t1,t2,t

      if (prob_type.gt.1) then
         call msgexit('not implemented prob',20)
      endif
      if (prob_type.eq.1) then
         n_fit=0
         v_fit=0.
         do 10 i=1,nsta
            if (sw.eq.0) then
               t1=tabpic(i,ishot,1)
               t2=tabpic(i,ishot,2)
            else
               t1=refpic(i,ishot,1)
               t2=refpic(i,ishot,2)
            endif
            t=vect_hdw(i)
            if (t1.lt.0.) goto 10
            d1=abs(t1-t)
            if (t2.eq.0.) t2=t1
            d2=abs(t2-t)
            if ((t.gt.t1).and.(t.lt.t2)) then
               f=0.
            else
               f=(amin1(d1,d2))**prob_par
            endif
            v_fit=v_fit+f
            n_fit=n_fit+1
            if (mode.eq.5) then
               tk(idatakernel)=t
               idatakernel=idatakernel+1
            endif
10       continue
         return
      endif
      end
c###############################################################################
      subroutine init_simplex
      implicit none

      include 'maxval.com'
      include 'simplex.com'
      include 'geom_space.com'
      include 'prec_verb.com'
      include 'prec_apri.com'
      include 'prec_mode.com'

      real*4 macrovit(PARMAX),a,rdom
      integer i,j,ival,k,ndim
      real*8 comp_misfit
      character*30 line

      if (mode_par.eq.0.) then
         call system('date > '//'/tmp/toto')
         open(10,file='/tmp/toto')
         read(10,'(a)') line
         close(10)
         call system('/bin/rm /tmp/toto')
         read(line(18:19),*) ival
      else
         ival=int(mode_par)
      endif
      if (ival.eq.0) ival=1
      a=rdom(ival)
      ndim=nxt(1)*nzt(1)
      if (nm.eq.2) ndim=ndim+nxt(2)*nzt(2)+ni
      include 'init_simplex.verb'
      do 10 i=1,ndim+1
         if (i.eq.1.and.apri_type.eq.1)  then
            do 15 j=1,ndim
               macrovit(j)=0.
15          continue
         else
            call comp_vit_al(macrovit)
         endif
         do 20 j=1,ndim
            mat_sim(i,j)=macrovit(j)
20       continue
         mat_mis(i)=comp_misfit(macrovit)
10    continue
      return
      end

c###############################################################################
      function verif_intin()
      implicit none

      logical verif_intin

      include 'maxval.com'
      include 'geom_space.com'
      include 'geom_topo.com'
      include 'mat_int.com'

      integer nzdn,nx,i

      nx=nxc+1
      nzdn=nzc*dn
      verif_intin=.true.
      do 10 i=1,nx
         if (tabinf(i).lt.topo(i).or
     &.abs(tabinf(i)-nzdn).lt.dn) then
            verif_intin=.false.
            return
         endif
10    continue
      return
      end

c###############################################################################
      function verif_modin(macrovit)
      implicit none

      include 'maxval.com'
      logical verif_modin
      real*4 macrovit(PARMAX)

      include 'geom_space.com'
      include 'mat_mask.com'
      include 'prec_apri.com'
      include 'prec_mode.com'

      integer i,n
      
      verif_modin=.true.
      n=nxt(1)*nzt(1)
      if (nm.eq.2) n=n+nxt(2)*nzt(2)+ni
      do 10 i=1,n
         if (macrovit(i).gt.maskup(i).or.
     & macrovit(i).lt.maskdw(i)) then
            verif_modin=.false.
            return
         endif
10    continue
      return
      end
c###############################################################################
      subroutine init_time(mat_time)
      implicit none

      include 'maxval.com'
      real*4 mat_time(VITMAX)

      include 'geom_space.com'

      real*4 infini
      integer i

      infini=1.e+06
      do 10 i=1,(nxc+1)*(nzc+1)
         mat_time(i)=infini
10    continue
      return
      end
c###############################################################################
      subroutine estim_time(mat_time,mat_vit,vect_hdw)
      implicit none

      include 'maxval.com'
      real*4 mat_time(VITMAX),vect_hdw(300)
      real*4 mat_vit(VITMAX)

      include 'geom_space.com'
      include 'station.com'
   
      integer i,cari,carj,nx
      real*4 x,y,v,t1,t2,xt,yt
      real*4 dst,t3,t4,xtm1,ytm1

      nx=nxc+1
      do 10 i=1,nsta
         x=locsta(i,1)
         y=locsta(i,2)
         cari=int(x/dn)+1
         carj=int(y/dn)+1
         v=mat_vit(cari+(carj-1)*nx)
         t1=mat_time(cari+carj*nx)
         t2=mat_time(cari+1+carj*nx)
         t3=mat_time(cari+(carj-1)*nx)
         t4=mat_time(cari+1+(carj-1)*nx)
         xtm1=(cari-1)*dn
         xt=cari*dn
         ytm1=(carj-1)*dn
         yt=carj*dn
         dst=(xt-x)*(yt-y)*t3
         dst=dst+(xt-x)*(y-ytm1)*t1
         dst=dst+(x-xtm1)*(yt-y)*t4
         dst=dst+(x-xtm1)*(y-ytm1)*t2
         vect_hdw(i)=dst/dn/dn
10    continue
      return
      end
c###############################################################################
      subroutine comp_vit_al(macrovit)
      implicit none

      include 'maxval.com'
      real*4 macrovit(PARMAX)

      include 'geom_space.com'
      include 'mat_mask.com'

      integer i,n
      real*4 a,rdom
     
      n=nxt(1)*nzt(1)
      if (nm.eq.2) n=n+nxt(2)*nzt(2)+ni
      do 10 i=1,n
         if (maskup(i).eq.maskdw(i)) then
             macrovit(i)=maskup(i)
             goto 10
         endif
         a=rdom(0)
         macrovit(i)=(maskup(i)-maskdw(i))*a+maskdw(i)
10    continue
      return
      end

c###############################################################################

      subroutine freeze_all(bwr)
      implicit none

      logical bwr

      include 'maxval.com'
      include 'geom_space.com'
      include 'geom_water.com'
      include 'geom_topo.com'
      include 'mat_mask.com'
      include 'shot.com'
      include 'station.com'
      include 'picking.com'
      include 'prec_mode.com'
      include 'prec_prob.com'
      include 'prec_search.com'
      include 'prec_verb.com'
      include 'simplex.com'
      include 'migration.com'
      common/fran/jran
        integer jran
  
      integer i,j,ncp,ngd,idx,k
      real*4 v1,v2,v3,v4,vx,vy

      open(30,file='freeze/freeze.dat',form='unformatted')
c mode
      if (bwr) then
         write(30) mode
      else
         read(30) mode
      endif
c iteration
      if (bwr) then
         write(30) iter
      else
         read(30) iter
      endif
c geom_space
      if (bwr) then
         write(30) dn,nxc,nzc,ni,nm
         write(30) nxt(1),nzt(1),nxt(2),nzt(2)
      else
         read(30) dn,nxc,nzc,ni,nm
         read(30) nxt(1),nzt(1),nxt(2),nzt(2)
         migration=(nxt(1)+nzt(1)+nxt(2)+nzt(2).eq.0)
      endif
c geom_topo
      if (bwr) then
         write(30) topo_zero
      else
         read(30) topo_zero
      endif
c mat_mask
      if (mode.eq.0) then
         ncp=nxt(1)*nzt(1)
         ngd=(nxc+1)*(nzc+1)
      else
         ncp=nxt(1)*nzt(1)+nxt(2)*nzt(2)+ni
         ngd=2*(nxc+1)*(nzc+1)+(nxc+1)
      endif
      if (bwr) then
         do 10 i=1,ngd
            write(30) maskgen(i)
            write(30) maskref(i)
10       continue
         do 20 i=1,ncp
            write(30) maskup(i)
            write(30) maskdw(i)
            write(30) mcrref(i)
20       continue
      else
         do 11 i=1,ngd
            read(30) maskgen(i)
            read(30) maskref(i)
11       continue
         do 21 i=1,ncp
            read(30) maskup(i)
            read(30) maskdw(i)
            read(30) mcrref(i)
21       continue
      endif
c shot
      if (bwr) then
         write(30) nsho
         do 30 i=1,nsho
            write(30) locsho(i,1)
            write(30) locsho(i,2)
30       continue
      else
         read(30) nsho
         do 31 i=1,nsho
            read(30) locsho(i,1)
            read(30) locsho(i,2)
31       continue
      endif
c station
      if (bwr) then
         write(30) nsta
         do 40 i=1,nsta
            write(30) locsta(i,1)
            write(30) locsta(i,2)
40       continue
      else
         read(30) nsta
         do 41 i=1,nsta
            read(30) locsta(i,1)
            read(30) locsta(i,2)
41       continue
      endif
c picking
      if (bwr) then
         do 50 i=1,nsta
            do 60 j=1,nsho
               if (.not.migration) then
                  write(30) tabpic(i,j,1)
                  write(30) tabpic(i,j,2)
               endif
               if (mode.eq.1) then
                  write(30) refpic(i,j,1)
                  write(30) refpic(i,j,2)
               endif
60          continue
50       continue
      else
         do 51 i=1,nsta
            do 61 j=1,nsho
               if (.not.migration) then
                  read(30) tabpic(i,j,1)
                  read(30) tabpic(i,j,2)
               endif
               if (mode.eq.1) then
                  read(30) refpic(i,j,1)
                  read(30) refpic(i,j,2)
               endif
61          continue
51       continue
      endif
c prec_prob
      if (bwr) then
         write(30) prob_type,prob_par
      else
         read(30) prob_type,prob_par
      endif
c prec_search
      if (bwr) then
         write(30) search_type,search_par
      else
         read(30) search_type,search_par
      endif
c prec_verb
      if (bwr) then
         write(30) verb_type,verb_file
      else
         read(30) verb_type,verb_file
      endif
c simplex
      if (bwr) then
         do 70 i=1,ncp+1
            do 80 j=1,ncp
               write(30) mat_sim(i,j)
80          continue
70       continue
         do 90 i=1,ncp+1
            write(30) mat_mis(i)
90       continue
      else
         do 71 i=1,ncp+1
            do 81 j=1,ncp
               read(30) mat_sim(i,j)
81          continue
71       continue
         do 91 i=1,ncp+1
            read(30) mat_mis(i)
91       continue
      endif
c random variable
      if (bwr) then
         write(30) jran
      else
         read(30) jran
      endif
c water variable
      if (bwr) then
         write(30) water_vit
      else
         read(30) water_vit
      endif
      close(30)
      if (.not.bwr) return
c only for download mode
      open(30,file='bestmodel.dat')
      if (.not.migration) then
      do 130 k=1,nm
      do 110 i=1,nxt(k)
         vx=(float(i)-.5)*(nxc*dn/float(nxt(k)))
         do 120 j=1,nzt(k)
            vy=topo_zero-(float(j)-.5)*(nzc*dn/float(nzt(k)))
            idx=i+(j-1)*nxt(k)
            if (k.eq.2) idx=idx+nxt(1)*nzt(1)
            v1=mat_sim(ilo,idx)+mcrref(idx)
            v2=mat_sim(ilo,idx)
            write(30,*) vx,vy,v1,v2
120       continue
110    continue
130   continue
      endif
      if (mode.eq.1) then
         do 125 i=1,ni
            vx=(float(i)-.5)*(nxc*dn/float(ni))
            idx=nxt(1)*nzt(1)+nxt(2)*nzt(2)+i
            v1=mat_sim(ilo,idx)+mcrref(idx)
            v2=mat_sim(ilo,idx)
            write(30,*) vx,topo_zero-v1,v2
125       continue
      endif
      close(30)
      return
      end

c###############################################################################

      function rdom(idum)
      implicit none

      integer idum,ia,im,ic
      real*4 rdom

      common/fran/jran
        integer jran

      ia=106
      im=6075
      ic=1283

      if (idum.ne.0) jran=idum
      jran=mod(jran*ia+ic,im)
      rdom=float(jran)/float(im)
      return
      end

c###############################################################################
      subroutine get_apri
      implicit none

      include 'maxval.com'
      include 'geom_space.com'
      include 'geom_topo.com'
      include 'mat_mask.com'
      include 'prec_apri.com'
      include 'prec_verb.com'
      include 'prec_mode.com'

      integer i,j,idx,nx,nz,lnblnk,n
  

      include 'get_apri.verb'
      nx=nxc+1
      nz=nzc+1
      if (mode.eq.0.or.mode.eq.5) n=nx*nz
      if (mode.eq.1) n=2*nx*nz+nx
      if (apri_type.eq.1) then
         open(30,file=apri_file(1:lnblnk(apri_file)),form='unformatted',
     &access='direct',recl=4,status='old',err=1000)
         do 10 i=1,n
            read(30,rec=i) maskref(i)
10       continue
         close(30)
         return
      endif
      do 20 i=1,n
         maskref(i)=0.
20    continue
      return
1000  call msgexit('get_apri',8)
      end
            

c###############################################################################
      subroutine comp_matvit2(macrovit,mat_vit,mat_vitr,bend)
      implicit none

      include 'maxval.com'
      real*4 macrovit(PARMAX),mat_vit(VITMAX)
      real*4 mat_vitr(VITMAX)
      real*4 mat_inf(VITMAX)
      integer bend

      include 'geom_space.com'
      include 'geom_water.com'
      include 'mat_mask.com'
      include 'prec_mode.com'
      include 'mat_int.com'
      include 'migration.com'

      integer i,j,idx
      integer nx,nz,k,kk
      real*4 curvit,infinity
      real*4 sp_tab(PARMAX),sp_x(50),sp_z(50)
      real*4 sp_dev(PARMAX),vx,vy,vd
 
      infinity=1.e+06
      nx=nxc+1
      nz=nzc+1
      do 1000 kk=1,nm
c constant case
         if ((nxt(kk).eq.1.and.nzt(kk).eq.1).or.migration) then
            if (kk.eq.1) then
               if (migration) then
                  vd=0.
               else
                  vd=macrovit(1)
               endif
               do 30 i=1,nx*nz
                  mat_vit(i)=vd+maskref(i)
c modif pour l'eau
c                 if (bend.eq.1) mat_vit(i)=mat_vit(i)*maskgen(i)
                  if (bend.eq.1) then
                     if (maskgen(i).eq.2) then
                        mat_vit(i)=water_vit
                     else
                        mat_vit(i)=mat_vit(i)*maskgen(i)
                     endif
                  endif
c fin modif
30             continue
            else
               if (migration) then
                  vd=0.
               else
                  vd=macrovit(nxt(1)*nzt(1)+1)
               endif
               do 40 i=1,nx*nz
                  mat_inf(i)=vd+maskref(i+nx*nz)
40             continue
            endif
         else
c spline initialisation
            if (kk.eq.1) then
               do 50 i=1,nxt(1)*nzt(1)
                  sp_tab(i)=macrovit(i)
50             continue
            else
               do 60 i=1,nxt(2)*nzt(2)
                  sp_tab(i)=macrovit(i+nxt(1)*nzt(1))
60             continue
            endif
            do 70 i=1,nxt(kk)
               sp_x(i)=(i-.5)*(nxc*dn/float(nxt(kk)))
70          continue
            do 80 j=1,nzt(kk)
               sp_z(j)=(j-.5)*(nzc*dn/float(nzt(kk)))
80          continue
            if (nxt(kk).eq.1) then
               call spline(sp_z,sp_tab,nzt(kk),1.e+30,
     &1.e+30,sp_dev)
            else
               call splie2(sp_x,sp_z,sp_tab,nxt(kk),
     &nzt(kk),sp_dev)
            endif
c spline computation
            do 90 i=1,nxc
               vx=(float(i)-.5)*dn
               do 100 j=1,nzc
                  vy=(float(j)-.5)*dn
                  idx=i+(j-1)*(nxc+1)
                  if (nxt(kk).eq.1) then
                     call splint(sp_z,sp_tab,sp_dev,
     &nzt(kk),vy,curvit)
                  else
                     call splin2(sp_x,sp_z,sp_tab,sp_dev,
     &nxt(kk),nzt(kk),vx,vy,curvit)
                  endif
                  if (kk.eq.1) then
                      mat_vit(idx)=curvit+maskref(idx)
c modif pour l'eau
c                     if (bend.eq.1) mat_vit(idx)=
c    &mat_vit(idx)*maskgen(idx)
                      if (bend.eq.1) then
                         if (maskgen(idx).eq.2) then
                            mat_vit(idx)=water_vit
                         else
                     mat_vit(idx)=mat_vit(idx)*maskgen(idx)
                         endif
                      endif
c fin modif
                  else
                      mat_inf(idx)=curvit+maskref(idx+nx*nz)
                  endif
100            continue
90          continue
         endif
1000  continue
c velocity matrix combination
      if (nm.eq.2) then
         do 110 i=1,nx*nz
            if (bend.eq.0) then
               mat_vitr(i)=mat_inf(i)
            else
               mat_vitr(i)=(1-matgeninf(i))*mat_vit(i)
               mat_vit(i)=mat_vitr(i)+mat_inf(i)*matgeninf(i)
            endif
110      continue
      endif
c velocity to slowness transformation
      do 120 i=1,nx*nz
         if (mat_vit(i).le.0.) then
             mat_vit(i)=infinity
         else
             mat_vit(i)=dn/mat_vit(i)
         endif
120   continue
      if (nm.eq.2) then
         do 130 i=1,nx*nz
            if (mat_vitr(i).le.0.) then
                mat_vitr(i)=infinity
            else
                mat_vitr(i)=dn/mat_vitr(i)
            endif
130      continue
      endif
c two borders put to infinity
      do 140 i=1,nxc+1
         mat_vit(i+nzc*(nxc+1))=infinity
140   continue
      do 150 j=1,nzc+1
         mat_vit(j*(nxc+1))=infinity
150   continue
      if (nm.eq.2) then
         do 160 i=1,nxc+1
            mat_vitr(i+nzc*(nxc+1))=infinity
160      continue
         do 170 j=1,nzc+1
            mat_vitr(j*(nxc+1))=infinity
170      continue
      endif
      return
      end

c###############################################################################
      subroutine comp_int(macrovit)
      implicit none

      include 'maxval.com'
      real*4 macrovit(PARMAX)
      include 'geom_space.com'
      include 'mat_mask.com'
      include 'mat_int.com'

      integer i,n,idx,j,nx,nz
      real*4 xinf(200),zinf(200),v,zcrt
      real*4 z2(200),x,nzdn

c cas ni = 1
      nx=nxc+1
      nz=nzc+1
      nzdn=nz*dn
      if (ni.eq.1) then
         do 5 i=1,nx
            tabinf(i)=macrovit(nxt(1)*nzt(1)+nxt(2)*nzt(2)+1)
5        continue
      else
         do 10 i=1,ni
            xinf(i)=((i-.5)*(nxc*dn/float(ni)))
            zinf(i)=macrovit(i+nxt(1)*nzt(1)+nxt(2)*nzt(2))
10       continue
         n=ni
         call spline(xinf,zinf,n,1.e32,1.e32,z2)
         do 20 i=1,nx
            x=(i-1)*dn
            call splint(xinf,zinf,z2,n,x,v)
            tabinf(i)=v+maskref(i+2*nx*nz)
            if (tabinf(i).gt.nzdn) tabinf(i)=nzdn
20       continue
      endif
      do 30 i=1,nx
         do 40 j=1,nz
            idx=i+(j-1)*nx
            zcrt=(j-1)*dn
            if (zcrt.lt.amax1(tabinf(i),tabinf(i+1))) then
               matgeninf(idx)=0.
            else
               matgeninf(idx)=1.
            endif
40       continue
30    continue
      return
      end
c###############################################################################
      function estim_ref(timeonref,timeonshot,sw)
      implicit none

      include 'maxval.com'
      real*4 timeonshot(NXCMAX),timeonref(NXCMAX,NSTAMAX)
      real*4 estim_ref
      integer sw

      include 'geom_space.com'
      include 'mat_int.com'

      integer i,carj,nx,idx
      real*4 x,y,t1,t3,dy

      nx=nxc+1
      estim_ref=10000.
      do 10 i=1,nx
         if (timeonref(i,sw)+timeonshot(i).lt.estim_ref) then
            estim_ref=timeonref(i,sw)+timeonshot(i)
            idx=i
         endif
10    continue
      if (idx.eq.1.or.idx.eq.nx) then
         estim_ref=-1.
      endif
      return
      end

c######################################################################
      subroutine dump_ray_r(mat_timer,mat_time,ish,ist)
      implicit none

      include 'maxval.com'
      include 'geom_topo.com'
      include 'geom_space.com'
      include 'shot.com'
      include 'station.com'
      include 'mat_int.com'

      real*4 mat_time(VITMAX)
      real*4 mat_timer(VITMAX)
      integer ish,sw,ist

      integer xcell,zcell,km
      integer idx1,idx2,idx3,idx4,idum,icrt
      real*4 pas,xcrt,zcrt,xgrd,zgrd,vmod,wx,wz,vclose
      real*4 xt,zt

      pas=dn/2.
      vclose=1.5*dn
      do 100 km=1,2
      write(60,'(a1)') '>'
      xcrt=impact(ist,1)
      zcrt=impact(ist,2)
      if (km.eq.1) then
         xt=locsho(ish,1)
         zt=locsho(ish,2)
      else
         xt=locsta(ist,1)
         zt=locsta(ist,2)
      endif
      idum=0
      icrt=0
10    write(60,*) xcrt,topo_zero-zcrt
      xcell=int(xcrt/dn)+1
      zcell=int(zcrt/dn)+1
      idx1=xcell+(zcell-1)*(nxc+1)
      idx2=xcell+(zcell-1)*(nxc+1)+1
      idx3=xcell+zcell*(nxc+1)
      idx4=xcell+zcell*(nxc+1)+1
      if (km.eq.1) then
      xgrd=mat_timer(idx4)+mat_timer(idx2)
      xgrd=xgrd-mat_timer(idx1)-mat_timer(idx3)
      zgrd=mat_timer(idx4)+mat_timer(idx3)
      zgrd=zgrd-mat_timer(idx1)-mat_timer(idx2)
      else
      xgrd=mat_time(idx4)+mat_time(idx2)
      xgrd=xgrd-mat_time(idx1)-mat_time(idx3)
      zgrd=mat_time(idx4)+mat_time(idx3)
      zgrd=zgrd-mat_time(idx1)-mat_time(idx2)
      endif
      vmod=sqrt(xgrd**2.+zgrd**2.)
      xcrt=xcrt-pas*xgrd/vmod
      zcrt=zcrt-pas*zgrd/vmod
      if (abs(xcrt-xt).lt.vclose.and.abs(zcrt-zt).lt.vclose) then
         write(60,*) xt,topo_zero-zt
      else
         goto 10
      endif
100   continue
      return
      end
c###############################################################################
      function time_ref(mat_time,timeonref,sw)
      implicit none

      include 'maxval.com'
      real*4 mat_time(VITMAX),timeonref(NXCMAX,2)
      real*4 time_ref
      integer sw

      include 'geom_space.com'
      include 'mat_int.com'

      integer i,carj,nx,idx
      real*4 x,y,t1,t3,dy

      nx=nxc+1
      time_ref=10000.
      do 10 i=1,nx
         x=(i-1)*dn
         y=tabinf(i)
         carj=int(y/dn)+1
         t1=mat_time(i+carj*nx)
         t3=mat_time(i+(carj-1)*nx)
         dy=y-(carj-1)*dn
         if (sw.eq.0) then
            timeonref(i,1)=dy/dn*(t1-t3)+t3
         else
            timeonref(i,2)=dy/dn*(t1-t3)+t3
c manu
c SP
c           if (1.8*timeonref(i,1)+timeonref(i,2).lt.time_ref) then
c              time_ref=1.8*timeonref(i,1)+timeonref(i,2)
c PS
            if (timeonref(i,1)+1.8*timeonref(i,2).lt.time_ref) then
               time_ref=timeonref(i,1)+1.8*timeonref(i,2)
c PP
c           if (timeonref(i,1)+timeonref(i,2).lt.time_ref) then
c              time_ref=timeonref(i,1)+timeonref(i,2)
               idx=i
            endif
         endif
10    continue
      if (sw.eq.0) then
         time_ref=0.
         return
      endif
      if (idx.eq.1.or.idx.eq.nx) then
         time_ref=-1.
      else
         impact(sw,1)=(idx-1)*dn
         impact(sw,2)=tabinf(idx)
      endif
      return
      end


      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      PARAMETER (NMAX=10) 
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF(DEN.EQ.0.)PAUSE
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M)THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE
      RETURN
      END


      SUBROUTINE SPLIE2(X1A,X2A,YA,M,N,Y2A)
      PARAMETER (NN=100)
      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),Y2TMP(NN)
      DO 13 J=1,M
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
11      CONTINUE
        CALL SPLINE(X2A,YTMP,N,1.E30,1.E30,Y2TMP)
        DO 12 K=1,N
          Y2A(J,K)=Y2TMP(K)
12      CONTINUE
13    CONTINUE
      RETURN
      END


      SUBROUTINE SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y)
      PARAMETER (NN=100)
      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),Y2TMP(NN),YYTMP(
     *NN)
      DO 12 J=1,M
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
          Y2TMP(K)=Y2A(J,K)
11      CONTINUE
        CALL SPLINT(X2A,YTMP,Y2TMP,N,X2,YYTMP(J))
12    CONTINUE
      CALL SPLINE(X1A,YYTMP,M,1.E30,1.E30,Y2TMP)
      CALL SPLINT(X1A,YYTMP,Y2TMP,M,X1,Y)
      RETURN
      END


      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      PARAMETER (NMAX=100)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END


      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      DIMENSION XA(N),YA(N),Y2A(N)
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) PAUSE 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END
