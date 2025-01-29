      common/geom_water/water,water_type,water_vit,water_file
         integer water_type
         real*4 water_vit,water(NXCMAX)
         character*20 water_file
