


   TH1F* landau1( float mu=0, float sigma=0.05, float xmax = 4 ) {

      gDirectory -> Delete( "h*" ) ;


      gStyle -> SetOptStat( "mr" ) ;
      gStyle -> SetOptTitle(0) ;
      gStyle -> SetPadTopMargin(0.03) ;
      gStyle -> SetPadRightMargin(0.03) ;
      gStyle -> SetPadBottomMargin(0.15) ;
      gStyle -> SetPadLeftMargin(0.15) ;
      gStyle -> SetStatFormat( "5.3f" ) ;

      TRandom tran ;
      char hname[100] ;
      char htitle[100] ;
      sprintf(hname, "h_landau" ) ;
      sprintf(htitle, "Landau, mu = %5.2f, sigma = %5.2f", mu, sigma ) ;
      TH1F* hp = new TH1F( hname, htitle, 80, 0., xmax ) ;


      int nevts = 100000 ;
      for ( int i=0; i<100000; i++ ) {
         float rn = tran.Landau( mu, sigma ) ;
         hp->Fill( rn ) ;
      } // i

      hp -> Scale( 1./ (1.*nevts) ) ;

      hp -> SetXTitle("Energy (GeV)") ;
      hp -> SetYTitle("Arbitrary units") ;

      hp -> SetTitleSize( 0.05, "X") ;
      hp -> SetTitleSize( 0.05, "Y") ;
      hp -> SetLabelSize( 0.045, "X" ) ;
      hp -> SetLabelSize( 0.045, "Y" ) ;


      TCanvas* can = new TCanvas( "can", "", 50, 50, 900, 600 ) ;

      hp -> SetFillColor(30) ;
      hp -> Draw( "hist" ) ;
      hp -> SetStats(0) ;


   //---------


      gStyle -> SetPadTopMargin(0.00) ;
      gStyle -> SetPadBottomMargin(0.18) ;
      gStyle -> SetPadRightMargin(0.03) ;
      gStyle -> SetStatW(0.30) ;
      gStyle -> SetStatH(0.30) ;
      gStyle -> SetStatX(0.93) ;
      gStyle -> SetStatY(0.93) ;

      TH1F* hpc = (TH1F*) hp -> Clone( "hcopy" ) ;
      hpc -> SetStats(1) ;

      hpc -> SetYTitle("") ;
      hpc -> SetXTitle("") ;
      hpc -> SetLabelSize( 0.080, "X" ) ;
      hpc -> SetLabelSize( 0.080, "Y" ) ;
      hpc -> SetXTitle("Energy (GeV)") ;
      hpc -> SetYTitle("Arbitrary units") ;
      hpc -> SetTitleSize( 0.08, "X") ;
      hpc -> SetTitleSize( 0.08, "Y") ;

      //TPad* tp = new TPad( "tp", "", 0.30, 0.35, 0.88, 0.85 ) ;
      TPad* tp = new TPad( "tp", "", 0.30, 0.35, 0.93, 0.89 ) ;
      tp->Draw() ;
      tp->cd() ;
      hpc -> Draw("hist") ;
      gPad -> SetLogy(1) ;

      can -> SaveAs("noise.pdf") ;



      return hp ;

   }
