c###############################################################################
      subroutine mask_water(mat_vit)
      implicit none

      real*4 mat_vit(1002501)

      include 'geom_space.com'
      include 'geom_water.com'
      include 'geom_topo.com'

      integer i,j,k,idx
      real*4 x,y,v

      do 10 i=1,nxc
         nw=int(amin1(water(i),water(i+1)/dn)
         if (nw.eq.0) goto 10
         j=1
30       idx=i+(j-1)*(nxc+1)
         if (mat_vit(idx).lt.1.e+06) then
            do 20 k=1,nw
               idx=idx+(k-1)*(nxc+1)
               mat_vit(idx)=dn/water_vit
20          continue
         else
            j=j+1
            goto 30
         endif
10    continue
      return
      end

c###############################################################################
