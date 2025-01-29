      common/prec_cstr/cstr_type,cstr_min,cstr_max,cstr_file
         integer cstr_type
         real*4 cstr_min,cstr_max
         character*20 cstr_file
