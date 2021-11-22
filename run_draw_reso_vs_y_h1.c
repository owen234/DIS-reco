
#include "draw_reso_vs_y.c"


   void run_draw_reso_vs_y_h1( const char* input_root_file = "plots-h1-v4c-f3a/dnn-output-h1.root" ) {

      draw_reso_vs_y( "x", "", input_root_file, "h1", 800 ) ;
      draw_reso_vs_y( "x", "has_norad" , input_root_file, "h1", 600) ;
      draw_reso_vs_y( "x", "has_isr" , input_root_file, "h1", 600) ;
      draw_reso_vs_y( "x", "has_fsr" , input_root_file, "h1", 600) ;

      draw_reso_vs_y( "y", "" , input_root_file, "h1", 800) ;
      draw_reso_vs_y( "y", "has_norad" , input_root_file, "h1", 600) ;
      draw_reso_vs_y( "y", "has_isr" , input_root_file, "h1", 600) ;
      draw_reso_vs_y( "y", "has_fsr" , input_root_file, "h1", 600) ;

      draw_reso_vs_y( "q2", "" , input_root_file, "h1", 800) ;
      draw_reso_vs_y( "q2", "has_norad" , input_root_file, "h1", 600) ;
      draw_reso_vs_y( "q2", "has_isr" , input_root_file, "h1", 600) ;
      draw_reso_vs_y( "q2", "has_fsr" , input_root_file, "h1", 600) ;


   }
