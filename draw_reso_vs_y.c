


void draw_reso_vs_y( const char* var = "x", const char* cuts = "", const char* input_root_file = "dnn-output1b.root", const char* experiment = "athena", int can_size_y=800 ) {


   gSystem -> Exec( "mkdir -p reso-plots" ) ;

   gStyle -> SetOptStat(0) ;
   gStyle -> SetPadBottomMargin(0.15) ;
   gStyle -> SetPadLeftMargin(0.18) ;
   gStyle -> SetPadTopMargin(0.10) ;
   gStyle -> SetPadRightMargin(0.05) ;

   if ( can_size_y >= 800 ) {
      gStyle -> SetTitleOffset( 1.30, "y" ) ;
   } else {
      gStyle -> SetTitleOffset( 1.00, "y" ) ;
   }

   TChain ch("dnnout") ;

   char var_text[100] ;
   if ( strcmp( var, "q2") == 0 ) {
      sprintf( var_text, "Q^{2}" ) ;
   } else {
      sprintf( var_text, "%s", var ) ;
   }

   TText* tt_title = new TText() ;
   tt_title -> SetTextSize( 0.06 ) ;



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







   ch.Add( input_root_file ) ;

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


   //int can_size_y = 800 ;

   //TCanvas* can = (TCanvas*) gDirectory -> FindObject("can") ;
   //if ( can == 0x0 ) can = new TCanvas("can","", 50, 50, 900, can_size_y ) ;

   TCanvas* can = new TCanvas("can","", 50, 50, 900, can_size_y ) ;



   sprintf( hname, "tp_dnn_%s", var ) ;
   sprintf( htitle, "%s resolution vs y, DNN", var ) ;

   TProfile* tp_dnn = new TProfile( hname, htitle, nbins, xaxis_min, xaxis_max, yaxis_min, yaxis_max, "s" ) ;

   sprintf( arg1, "dnn_%s/true_%s:true_y>>%s", var, var, hname ) ;

   ch.Draw( arg1, cuts ) ;





   sprintf( hname, "tp_e_%s", var ) ;
   sprintf( htitle, "%s resolution vs y, e", var ) ;

   TProfile* tp_e = new TProfile( hname, htitle, nbins, xaxis_min, xaxis_max, yaxis_min, yaxis_max, "s" ) ;

   sprintf( arg1, "rec_%s_e/true_%s:true_y>>%s", var, var, hname ) ;

   ch.Draw( arg1, cuts ) ;




   sprintf( hname, "tp_da_%s", var ) ;
   sprintf( htitle, "%s resolution vs y, DA", var ) ;

   TProfile* tp_da = new TProfile( hname, htitle, nbins, xaxis_min, xaxis_max, yaxis_min, yaxis_max, "s" ) ;

   sprintf( arg1, "rec_%s_da/true_%s:true_y>>%s", var, var, hname ) ;

   ch.Draw( arg1, cuts ) ;





   sprintf( hname, "tp_h_%s", var ) ;
   sprintf( htitle, "%s resolution vs y, hadron", var ) ;

   TProfile* tp_h = new TProfile( hname, htitle, nbins, xaxis_min, xaxis_max, yaxis_min, yaxis_max, "s" ) ;

   sprintf( arg1, "rec_%s_h/true_%s:true_y>>%s", var, var, hname ) ;

   ch.Draw( arg1, cuts ) ;



   sprintf( hname, "tp_is_%s", var ) ;
   sprintf( htitle, "%s resolution vs y, ISigma", var ) ;

   TProfile* tp_is = new TProfile( hname, htitle, nbins, xaxis_min, xaxis_max, yaxis_min, yaxis_max, "s" ) ;

   sprintf( arg1, "rec_%s_is/true_%s:true_y>>%s", var, var, hname ) ;

   ch.Draw( arg1, cuts ) ;





   double y_vals[nbins] ;
   double dnn_rms[nbins] ;
   double e_rms[nbins] ;
   double da_rms[nbins] ;
   double h_rms[nbins] ;
   double is_rms[nbins] ;

   double dnn_mean[nbins] ;
   double e_mean[nbins] ;
   double da_mean[nbins] ;
   double h_mean[nbins] ;
   double is_mean[nbins] ;

   for ( int i=0; i<nbins; i++ ) {

      int bi = i+1 ;

      y_vals[i] = tp_dnn -> GetBinCenter( bi ) ;

      dnn_rms[i] = tp_dnn -> GetBinError( bi ) ;
      e_rms[i] = tp_e -> GetBinError( bi ) ;
      da_rms[i] = tp_da -> GetBinError( bi ) ;
      h_rms[i] = tp_h -> GetBinError( bi ) ;
      is_rms[i] = tp_is -> GetBinError( bi ) ;

      dnn_mean[i] = tp_dnn -> GetBinContent( bi ) ;
      e_mean[i] = tp_e -> GetBinContent( bi ) ;
      da_mean[i] = tp_da -> GetBinContent( bi ) ;
      h_mean[i] = tp_h -> GetBinContent( bi ) ;
      is_mean[i] = tp_is -> GetBinContent( bi ) ;


      printf("  %3d : y = %8.3f,  DNN = %7.3f, e = %7.3f, DA = %7.3f\n", bi, y_vals[i], dnn_rms[i], e_rms[i], da_rms[i] ) ;
   }

   TGraph* tg_rms_dnn = new TGraph( nbins, y_vals, dnn_rms ) ;
   TGraph* tg_rms_e = new TGraph( nbins, y_vals, e_rms ) ;
   TGraph* tg_rms_da = new TGraph( nbins, y_vals, da_rms ) ;
   TGraph* tg_rms_h = new TGraph( nbins, y_vals, h_rms ) ;
   TGraph* tg_rms_is = new TGraph( nbins, y_vals, is_rms ) ;

   TGraph* tg_mean_dnn = new TGraph( nbins, y_vals, dnn_mean ) ;
   TGraph* tg_mean_e = new TGraph( nbins, y_vals, e_mean ) ;
   TGraph* tg_mean_da = new TGraph( nbins, y_vals, da_mean ) ;
   TGraph* tg_mean_h = new TGraph( nbins, y_vals, h_mean ) ;
   TGraph* tg_mean_is = new TGraph( nbins, y_vals, is_mean ) ;


   float marker_size = 1.5 ;

   tg_rms_dnn -> SetMarkerStyle(20) ;
   tg_rms_dnn -> SetLineColor(2) ;
   tg_rms_dnn -> SetLineWidth(3) ;
   tg_rms_dnn -> SetMarkerColor(2) ;
   tg_rms_dnn -> SetMarkerSize(marker_size) ;

   tg_rms_e -> SetMarkerStyle(21) ;
   tg_rms_e -> SetLineColor(4) ;
   tg_rms_e -> SetLineWidth(3) ;
   tg_rms_e -> SetMarkerColor(4) ;
   tg_rms_e -> SetMarkerSize(marker_size) ;

   tg_rms_da -> SetMarkerStyle(22) ;
   tg_rms_da -> SetLineColor(6) ;
   tg_rms_da -> SetLineWidth(3) ;
   tg_rms_da -> SetMarkerColor(6) ;
   tg_rms_da -> SetMarkerSize(marker_size) ;

   tg_rms_h -> SetMarkerStyle(23) ;
   tg_rms_h -> SetLineColor(7) ;
   tg_rms_h -> SetLineWidth(3) ;
   tg_rms_h -> SetMarkerColor(7) ;
   tg_rms_h -> SetMarkerSize(marker_size) ;

   tg_rms_is -> SetMarkerStyle(33) ;
   tg_rms_is -> SetLineColor(8) ;
   tg_rms_is -> SetLineWidth(3) ;
   tg_rms_is -> SetMarkerColor(8) ;
   tg_rms_is -> SetMarkerSize(marker_size*1.3) ;




   tg_mean_dnn -> SetMarkerStyle(20) ;
   tg_mean_dnn -> SetLineColor(2) ;
   tg_mean_dnn -> SetLineWidth(3) ;
   tg_mean_dnn -> SetMarkerColor(2) ;
   tg_mean_dnn -> SetMarkerSize(marker_size) ;

   tg_mean_e -> SetMarkerStyle(21) ;
   tg_mean_e -> SetLineColor(4) ;
   tg_mean_e -> SetLineWidth(3) ;
   tg_mean_e -> SetMarkerColor(4) ;
   tg_mean_e -> SetMarkerSize(marker_size) ;

   tg_mean_da -> SetMarkerStyle(22) ;
   tg_mean_da -> SetLineColor(6) ;
   tg_mean_da -> SetLineWidth(3) ;
   tg_mean_da -> SetMarkerColor(6) ;
   tg_mean_da -> SetMarkerSize(marker_size) ;

   tg_mean_h -> SetMarkerStyle(23) ;
   tg_mean_h -> SetLineColor(7) ;
   tg_mean_h -> SetLineWidth(3) ;
   tg_mean_h -> SetMarkerColor(7) ;
   tg_mean_h -> SetMarkerSize(marker_size) ;

   tg_mean_is -> SetMarkerStyle(33) ;
   tg_mean_is -> SetLineColor(8) ;
   tg_mean_is -> SetLineWidth(3) ;
   tg_mean_is -> SetMarkerColor(8) ;
   tg_mean_is -> SetMarkerSize(marker_size*1.3) ;











   TH2F* hd_rms = new TH2F( "hd_rms", "", 100, xaxis_min, xaxis_max, 100, 0., 0.5 ) ;
   hd_rms -> SetTitleSize( 0.06, "x" ) ;
   hd_rms -> SetTitleSize( 0.06, "y" ) ;
   hd_rms -> SetLabelSize( 0.05, "x" ) ;
   hd_rms -> SetLabelSize( 0.05, "y" ) ;
   hd_rms -> SetXTitle( "Gen y" ) ;
   sprintf( label, "RMS, %s / %s_{gen}", var_text, var_text ) ;
   hd_rms -> SetYTitle( label ) ;


   TH2F* hd_mean = new TH2F( "hd_mean", "", 100, xaxis_min, xaxis_max, 100, 0.8, 1.2 ) ;
   hd_mean -> SetTitleSize( 0.06, "x" ) ;
   hd_mean -> SetTitleSize( 0.06, "y" ) ;
   hd_mean -> SetLabelSize( 0.05, "x" ) ;
   hd_mean -> SetLabelSize( 0.05, "y" ) ;
   hd_mean -> SetXTitle( "Gen y" ) ;
   sprintf( label, "Mean, %s / %s_{gen}", var_text, var_text ) ;
   hd_mean -> SetYTitle( label ) ;



   float lx  ;
   float ly  ;
   float lw  ;
   float lh  ;



   //TCanvas* can_rms = (TCanvas*) gDirectory -> FindObject("can_rms") ;
   //if ( can_rms == 0x0 ) can_rms = new TCanvas("can_rms","", 50, 50, 900, can_size_y ) ;

   TCanvas* can_rms = new TCanvas("can_rms","", 50, 50, 900, can_size_y ) ;

   can_rms -> cd() ;

   hd_rms -> Draw() ;
   tg_rms_e->Draw("pl") ;
   tg_rms_da->Draw("pl") ;
   tg_rms_h->Draw("pl") ;
   tg_rms_is->Draw("pl") ;
   tg_rms_dnn->Draw("pl") ;

   gPad -> SetGridy(1) ;

   tt_title -> DrawTextNDC( 0.15, 0.93, plot_title ) ;

   lx = 0.35 ;
   ly = 0.57 ;
   lw = 0.18 ;
   lh = 0.30 ;

   TLegend* tl_rms = new TLegend( lx, ly, lx+lw, ly+lh ) ;

   tl_rms -> AddEntry( tg_rms_dnn, "DNN" ) ;
   tl_rms -> AddEntry( tg_rms_e, "electron" ) ;
   tl_rms -> AddEntry( tg_rms_da, "DA" ) ;
   tl_rms -> AddEntry( tg_rms_h, "hadron" ) ;
   tl_rms -> AddEntry( tg_rms_is, "ISigma" ) ;
   tl_rms -> Draw() ;


   sprintf( save_fname, "reso-plots/%s-rms-%s-vs-y%s.pdf", experiment, var, cut_label ) ;
   can_rms -> SaveAs( save_fname ) ;

   sprintf( save_fname, "reso-plots/%s-rms-%s-vs-y%s.png", experiment, var, cut_label ) ;
   can_rms -> SaveAs( save_fname ) ;

   can -> Update() ;
   can -> Draw() ;
   gSystem -> ProcessEvents() ;





   //TCanvas* can_mean = (TCanvas*) gDirectory -> FindObject("can_mean") ;
   //if ( can_mean == 0x0 ) can_mean = new TCanvas("can_mean","", 950, 50, 900, can_size_y ) ;

   TCanvas* can_mean = new TCanvas("can_mean","", 950, 50, 900, can_size_y ) ;

   can_mean -> cd() ;

   hd_mean -> Draw() ;
   tg_mean_e->Draw("pl") ;
   tg_mean_da->Draw("pl") ;
   tg_mean_h->Draw("pl") ;
   tg_mean_is->Draw("pl") ;
   tg_mean_dnn->Draw("pl") ;

   gPad -> SetGridy(1) ;

   tt_title -> DrawTextNDC( 0.15, 0.93, plot_title ) ;

   lx = 0.79 ;
   ly = 0.65 ;
   lw = 0.18 ;
   lh = 0.30 ;


   TLegend* tl_mean = new TLegend( lx, ly, lx+lw, ly+lh ) ;

   tl_mean -> AddEntry( tg_mean_dnn, "DNN" ) ;
   tl_mean -> AddEntry( tg_mean_e, "electron" ) ;
   tl_mean -> AddEntry( tg_mean_da, "DA" ) ;
   tl_mean -> AddEntry( tg_mean_h, "hadron" ) ;
   tl_mean -> AddEntry( tg_mean_is, "ISigma" ) ;
   tl_mean -> Draw() ;


   sprintf( save_fname, "reso-plots/%s-mean-%s-vs-y%s.pdf", experiment, var, cut_label ) ;
   can_mean -> SaveAs( save_fname ) ;

   sprintf( save_fname, "reso-plots/%s-mean-%s-vs-y%s.png", experiment, var, cut_label ) ;
   can_mean -> SaveAs( save_fname ) ;

   can -> Update() ;
   can -> Draw() ;
   gSystem -> ProcessEvents() ;



   gDirectory->ls() ;

}












