c###############################################################################
      subroutine comp_water
      implicit none
      
      include 'geom_space.com'
      include 'geom_water.com'
      include 'prec_verb.com'

      real*4 xtemp(10),ztemp(10)
      real*4 xwat(500),ywat(500)
      real*4 vx,vy,dy,vmax
      integer i,nt,j,k,nse,nx
      integer indmin(10),pmax

      open(10,file=water_file,status='old',err=1000)
      i=1
20    read(10,*,end=100) xwat(i),ywat(i)
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
            if (topo(i).lt.0.) water(i)=0.
80       continue
      endif
      return
1000  call msgexit('input water',11)
      end
c###############################################################################
