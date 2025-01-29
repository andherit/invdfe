      common/prec_mode/mode,mode_par,freeze,mode_file,fres
         integer freeze,mode,fres
         real*4 mode_par
         character*20 mode_file
