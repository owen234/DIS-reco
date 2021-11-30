


void compare_dnn_reso_vs_y( const char* var = "x", const char* cuts = "",
       const char* input_root_file1 = "plots-h1-v4c-f3a/dnn-output-h1.root",
       const char* input_root_file2 = "plots-h1-v4c-f3a-django/dnn-output-h1.root",
       const char* experiment = "h1",
       const char* dnn_name1 = "RAPGAP",
       const char* dnn_name2 = "DJANGOH" ) {


   //char dnn_name1[100] ;
   //char dnn_name2[100] ;
   //sprintf( dnn_name1, "RAPGAP" ) ;
   //sprintf( dnn_name2, "DJANGO" ) ;

   gSystem -> Exec( "mkdir -p reso-plots" ) ;

   gStyle -> SetOptStat(0) ;

   gStyle -> SetPadBottomMargin(0.15) ;
   gStyle -> SetPadLeftMargin(0.15) ;
   gStyle -> SetPadTopMargin(0.10) ;
   gStyle -> SetPadRightMargin(0.05) ;

      gStyle -> SetPadLeftMargin(0.18) ;
      gStyle -> SetPadBottomMargin(0.17) ;
      gStyle -> SetTitleOffset( 1.00, "y" ) ;
      gStyle -> SetTitleOffset( 0.90, "x" ) ;
         gStyle -> SetLabelOffset( 0.015, "y" ) ;
         gStyle -> SetLabelOffset( 0.015, "x" ) ;


   TChain ch1("dnnout") ;
   TChain ch2("dnnout") ;

   char var_text[100] ;
   if ( strcmp( var, "q2") == 0 ) {
      sprintf( var_text, "Q^{2}" ) ;
   } else {
      sprintf( var_text, "%s", var ) ;
   }

   TText* tt_title = new TText() ;
   ////////////tt_title -> SetTextSize( 0.06 ) ;
      tt_title -> SetTextSize( 0.08 ) ;



   char cut_label[100] ;

   TString ts( cuts ) ;
   ts.ReplaceAll("_","-") ;
   ts.ReplaceAll(" ","-") ;
   if ( strlen( cuts ) > 0 ) {
      sprintf( cut_label, "-%s", ts.Data() ) ;
   } else {
      sprintf( cut_label, "-allevts" ) ;
   }



   char plot_title[1000] ;
   sprintf( plot_title, "" ) ;
   if ( strlen( cuts ) == 0 ) { sprintf( plot_title, "All events" ) ; }
   if ( strcmp( cuts, "has_norad" ) == 0 ) { sprintf( plot_title, "No QED radiation only" ) ; }
   if ( strcmp( cuts, "has_isr" ) == 0 ) { sprintf( plot_title, "ISR events only" ) ; }
   if ( strcmp( cuts, "has_fsr" ) == 0 ) { sprintf( plot_title, "FSR events only" ) ; }







   ch1.Add( input_root_file1 ) ;
   ch2.Add( input_root_file2 ) ;

   int nbins=20 ;
   float xaxis_min = 0.0 ;
   float xaxis_max = 0.8 ;
   float yaxis_min = 0. ;
   float yaxis_max = 2.0 ;

   char hname[1000] ;
   char htitle[1000] ;
   char arg1[1000] ;
   char label[100] ;

   char save_fname[1000] ;


   TCanvas* can = (TCanvas*) gDirectory -> FindObject("can") ;
   if ( can == 0x0 ) can = new TCanvas("can","", 50, 50, 900, 600 ) ;



   sprintf( hname, "tp_dnn1_%s", var ) ;
   sprintf( htitle, "%s resolution vs y, DNN1", var ) ;

   TProfile* tp_dnn1 = new TProfile( hname, htitle, nbins, xaxis_min, xaxis_max, yaxis_min, yaxis_max, "s" ) ;

   sprintf( arg1, "dnn_%s/true_%s:true_y>>%s", var, var, hname ) ;

   ch1.Draw( arg1, cuts ) ;



   sprintf( hname, "tp_dnn2_%s", var ) ;
   sprintf( htitle, "%s resolution vs y, DNN2", var ) ;

   TProfile* tp_dnn2 = new TProfile( hname, htitle, nbins, xaxis_min, xaxis_max, yaxis_min, yaxis_max, "s" ) ;

   sprintf( arg1, "dnn_%s/true_%s:true_y>>%s", var, var, hname ) ;

   ch2.Draw( arg1, cuts ) ;





   double y_vals[nbins] ;

   double dnn1_rms[nbins] ;
   double dnn2_rms[nbins] ;

   double dnn1_mean[nbins] ;
   double dnn2_mean[nbins] ;

   for ( int i=0; i<nbins; i++ ) {

      int bi = i+1 ;

      y_vals[i] = tp_dnn1 -> GetBinCenter( bi ) ;

      dnn1_rms[i] = tp_dnn1 -> GetBinError( bi ) ;
      dnn2_rms[i] = tp_dnn2 -> GetBinError( bi ) ;

      dnn1_mean[i] = tp_dnn1 -> GetBinContent( bi ) ;
      dnn2_mean[i] = tp_dnn2 -> GetBinContent( bi ) ;

      printf("  %3d : y = %8.3f,  DNN1 = %7.4f, DNN2 = %7.4f\n", bi, y_vals[i], dnn1_rms[i], dnn2_rms[i] ) ;
   }

   TGraph* tg_rms_dnn1 = new TGraph( nbins, y_vals, dnn1_rms ) ;
   TGraph* tg_rms_dnn2 = new TGraph( nbins, y_vals, dnn2_rms ) ;

   TGraph* tg_mean_dnn1 = new TGraph( nbins, y_vals, dnn1_mean ) ;
   TGraph* tg_mean_dnn2 = new TGraph( nbins, y_vals, dnn2_mean ) ;


   float marker_size = 2.0 ;

   tg_rms_dnn1 -> SetMarkerStyle(23) ;
   tg_rms_dnn1 -> SetLineColor(2) ;
   tg_rms_dnn1 -> SetLineWidth(3) ;
   tg_rms_dnn1 -> SetMarkerColor(2) ;
   tg_rms_dnn1 -> SetMarkerSize(marker_size) ;

   tg_rms_dnn2 -> SetMarkerStyle(22) ;
   tg_rms_dnn2 -> SetLineColor(4) ;
   tg_rms_dnn2 -> SetLineWidth(3) ;
   tg_rms_dnn2 -> SetMarkerColor(4) ;
   tg_rms_dnn2 -> SetMarkerSize(marker_size) ;


   tg_mean_dnn1 -> SetMarkerStyle(23) ;
   tg_mean_dnn1 -> SetLineColor(2) ;
   tg_mean_dnn1 -> SetLineWidth(3) ;
   tg_mean_dnn1 -> SetMarkerColor(2) ;
   tg_mean_dnn1 -> SetMarkerSize(marker_size) ;

   tg_mean_dnn2 -> SetMarkerStyle(22) ;
   tg_mean_dnn2 -> SetLineColor(4) ;
   tg_mean_dnn2 -> SetLineWidth(3) ;
   tg_mean_dnn2 -> SetMarkerColor(4) ;
   tg_mean_dnn2 -> SetMarkerSize(marker_size) ;









   TH2F* hd_rms = new TH2F( "hd_rms", "", 100, xaxis_min, xaxis_max, 100, 0., 0.25 ) ;

   /////////////hd_rms -> SetTitleSize( 0.06, "x" ) ;
   /////////////hd_rms -> SetTitleSize( 0.06, "y" ) ;
   /////////////hd_rms -> SetLabelSize( 0.05, "x" ) ;
   /////////////hd_rms -> SetLabelSize( 0.05, "y" ) ;

      hd_rms -> SetTitleSize( 0.085, "x" ) ;
      hd_rms -> SetTitleSize( 0.085, "y" ) ;
      hd_rms -> SetLabelSize( 0.075, "x" ) ;
      hd_rms -> SetLabelSize( 0.075, "y" ) ;

        hd_rms -> SetNdivisions( 805, "x" ) ;
        hd_rms -> SetNdivisions( 805, "y" ) ;


   hd_rms -> SetXTitle( "Gen y" ) ;
   sprintf( label, "RMS, %s / %s_{gen}", var_text, var_text ) ;
   hd_rms -> SetYTitle( label ) ;


   TH2F* hd_mean = new TH2F( "hd_mean", "", 100, xaxis_min, xaxis_max, 100, 0.9, 1.1 ) ;

   ////////////hd_mean -> SetTitleSize( 0.06, "x" ) ;
   ////////////hd_mean -> SetTitleSize( 0.06, "y" ) ;
   ////////////hd_mean -> SetLabelSize( 0.05, "x" ) ;
   ////////////hd_mean -> SetLabelSize( 0.05, "y" ) ;

      hd_mean -> SetTitleSize( 0.085, "x" ) ;
      hd_mean -> SetTitleSize( 0.085, "y" ) ;
      hd_mean -> SetLabelSize( 0.075, "x" ) ;
      hd_mean -> SetLabelSize( 0.075, "y" ) ;

        hd_mean -> SetNdivisions( 805, "x" ) ;
        hd_mean -> SetNdivisions( 805, "y" ) ;


   hd_mean -> SetXTitle( "Gen y" ) ;
   sprintf( label, "Mean, %s / %s_{gen}", var_text, var_text ) ;
   hd_mean -> SetYTitle( label ) ;



   float lx  ;
   float ly  ;
   float lw  ;
   float lh  ;



   TCanvas* can_rms = (TCanvas*) gDirectory -> FindObject("can_rms") ;
   if ( can_rms == 0x0 ) can_rms = new TCanvas("can_rms","", 50, 50, 900, 600 ) ;

   can_rms -> cd() ;

   hd_rms -> Draw() ;
   tg_rms_dnn1->Draw("pl") ;
   tg_rms_dnn2->Draw("pl") ;

   ///////////////gPad -> SetGridy(1) ;

   tt_title -> DrawTextNDC( 0.15, 0.93, plot_title ) ;

   lx = 0.64 ;
   ly = 0.70 ;
   lw = 0.27 ;
   lh = 0.18 ;

      gStyle -> SetLegendTextSize( 0.060 ) ;

   TLegend* tl_rms = new TLegend( lx, ly, lx+lw, ly+lh ) ;

   tl_rms -> AddEntry( tg_rms_dnn1, dnn_name1 ) ;
   tl_rms -> AddEntry( tg_rms_dnn2, dnn_name2 ) ;
   tl_rms -> Draw() ;


   sprintf( save_fname, "reso-plots/%s-comparison-rms-%s-vs-y%s.pdf", experiment, var, cut_label ) ;
   can_rms -> SaveAs( save_fname ) ;

   sprintf( save_fname, "reso-plots/%s-comparison-rms-%s-vs-y%s.png", experiment, var, cut_label ) ;
   can_rms -> SaveAs( save_fname ) ;

   can -> Update() ;
   can -> Draw() ;
   gSystem -> ProcessEvents() ;





   TCanvas* can_mean = (TCanvas*) gDirectory -> FindObject("can_mean") ;
   if ( can_mean == 0x0 ) can_mean = new TCanvas("can_mean","", 950, 50, 900, 600 ) ;

   can_mean -> cd() ;

   hd_mean -> Draw() ;
   tg_mean_dnn1->Draw("pl") ;
   tg_mean_dnn2->Draw("pl") ;

   ///////////////gPad -> SetGridy(1) ;

   tt_title -> DrawTextNDC( 0.15, 0.93, plot_title ) ;

   lx = 0.64 ;
   ly = 0.70 ;
   lw = 0.27 ;
   lh = 0.18 ;


   TLegend* tl_mean = new TLegend( lx, ly, lx+lw, ly+lh ) ;

   tl_mean -> AddEntry( tg_mean_dnn1, dnn_name1 ) ;
   tl_mean -> AddEntry( tg_mean_dnn2, dnn_name2 ) ;
   tl_mean -> Draw() ;


   sprintf( save_fname, "reso-plots/%s-comparison-mean-%s-vs-y%s.pdf", experiment, var, cut_label ) ;
   can_mean -> SaveAs( save_fname ) ;

   sprintf( save_fname, "reso-plots/%s-comparison-mean-%s-vs-y%s.png", experiment, var, cut_label ) ;
   can_mean -> SaveAs( save_fname ) ;

   can -> Update() ;
   can -> Draw() ;
   gSystem -> ProcessEvents() ;



   gDirectory->ls() ;

}






