      common/geom_topo/topo,topo_type,topo_alt,topo_zero,topo_file
         integer topo_type
         real*4 topo_alt,topo_zero,topo(NXCMAX)
         character*20 topo_file

