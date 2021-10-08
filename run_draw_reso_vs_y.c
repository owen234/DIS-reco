
#include "draw_reso_vs_y.c"


   void run_draw_reso_vs_y( const char* input_root_file = "plots-v4d-f2b/dnn-output.root" ) {

      draw_reso_vs_y( "x", "", input_root_file, "athena", 800 ) ;
      draw_reso_vs_y( "x", "has_norad", input_root_file, "athena", 600 ) ;
      draw_reso_vs_y( "x", "has_isr", input_root_file, "athena", 600 ) ;
      draw_reso_vs_y( "x", "has_fsr", input_root_file, "athena", 600 ) ;

      draw_reso_vs_y( "y", "", input_root_file, "athena", 800 ) ;
      draw_reso_vs_y( "y", "has_norad", input_root_file, "athena", 600 ) ;
      draw_reso_vs_y( "y", "has_isr", input_root_file, "athena", 600 ) ;
      draw_reso_vs_y( "y", "has_fsr", input_root_file, "athena", 600 ) ;

      draw_reso_vs_y( "q2", "", input_root_file, "athena", 800 ) ;
      draw_reso_vs_y( "q2", "has_norad", input_root_file, "athena", 600 ) ;
      draw_reso_vs_y( "q2", "has_isr", input_root_file, "athena", 600 ) ;
      draw_reso_vs_y( "q2", "has_fsr", input_root_file, "athena", 600 ) ;


   }
