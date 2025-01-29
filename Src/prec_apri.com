      common/prec_apri/apri_type,apri_vit,apri_grd,apri_file,apri_bnd
         integer apri_type
         real*4 apri_vit,apri_grd,apri_bnd
         character*20 apri_file
