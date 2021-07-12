#define minitree_maker_cxx
#include "minitree_maker.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

float calc_dr( double phi1, double phi2, double eta1, double eta2 ) ;

void ptep_to_xyze( float pt, double eta, double phi, double m,
                   float& px, float& py, float& pz, float& e ) ;

void minitree_maker::Loop( bool verbose, int maxEvt )
{

   if (fChain == 0) return;

   TFile* tf_out = new TFile("mini-tree.root", "recreate") ;
   TTree* tt_out = new TTree("minitree", "Minimal flat TTree for hadronic Q2,x,y analysis" ) ;

   float beam_electron_energy ;
   float beam_proton_energy ;

   float gene_px ;
   float gene_py ;
   float gene_pz ;
   float gen_e_pt ;
   float gen_e_eta ;
   float gen_e_phi ;
   float gen_e_theta ;
   float gen_e_e ;

   float gen_Q2 ;
   float gen_y ;
   float gen_s ;
   float gen_x ;

   float HFS_px ;
   float HFS_py ;
   float HFS_pz ;
   float HFS_E  ;

   float rec_e_pt ;
   float rec_e_eta ;
   float rec_e_phi ;

   float e_px ;
   float e_py ;
   float e_pz ;

   float rec_e_theta ;
   float rec_e_e ;


   float Sigma ;
   float pth ;
   float Empz ; // this is Delta (eqn 7 in the paper)

   float y_sigma ;
   float Q2_sigma ;
   float Q2_esigma ;
   float x_sigma ;

   float y_esigma ;

   float y_e ;
   float Q2_e ;
   float x_e ;

   float y_h ;
   float Q2_h ;
   float x_h ;

   float y_da ;
   float Q2_da ;
   float x_da ;



   tt_out -> Branch( "Q2_e", &Q2_e, "Q2_e/F" ) ; // ***
   tt_out -> Branch( "Q2_sigma", &Q2_sigma, "Q2_sigma/F" ) ; // ***
   tt_out -> Branch( "Q2_esigma", &Q2_esigma, "Q2_esigma/F" ) ; // ***
   tt_out -> Branch( "Q2_da", &Q2_da, "Q2_da/F" ) ; // ***
   tt_out -> Branch( "gen_Q2", &gen_Q2, "gen_Q2/F" ) ; // ***
   tt_out -> Branch( "y_e", &y_e, "y_e/F" ) ; // ***
   tt_out -> Branch( "y_sigma", &y_sigma, "y_sigma/F" ) ; // ***
   tt_out -> Branch( "y_esigma", &y_esigma, "y_esigma/F" ) ; // ***
   tt_out -> Branch( "y_da", &y_da, "y_da/F" ) ; // ***
   tt_out -> Branch( "HFS_px", &HFS_px, "HFS_px/F" ) ; // ***
   tt_out -> Branch( "HFS_py", &HFS_py, "HFS_py/F" ) ; // ***
   tt_out -> Branch( "HFS_pz", &HFS_pz, "HFS_pz/F" ) ; // ***
   tt_out -> Branch( "gen_y", &gen_y, "gen_y/F" ) ; // ***
   tt_out -> Branch( "e_px", &e_px, "e_px/F" ) ; // ***
   tt_out -> Branch( "e_py", &e_py, "e_py/F" ) ; // ***
   tt_out -> Branch( "e_pz", &e_pz, "e_pz/F" ) ; // ***
   tt_out -> Branch( "gene_px", &gene_px, "gene_px/F" ) ; // ***
   tt_out -> Branch( "gene_py", &gene_py, "gene_py/F" ) ; // ***
   tt_out -> Branch( "gene_pz", &gene_pz, "gene_pz/F" ) ; // ***

   tt_out -> Branch( "Empz", &Empz, "Empz/F" ) ; // ***
   tt_out -> Branch( "pth", &pth, "pth/F" ) ; // ***


   tt_out -> Branch( "x_sigma", &x_sigma, "x_sigma/F" ) ; // +++
   tt_out -> Branch( "x_e", &x_e, "x_e/F" ) ; // +++
   tt_out -> Branch( "x_da", &x_da, "x_da/F" ) ; // +++

   tt_out -> Branch( "y_h", &y_h, "y_h/F" ) ; // +++
   tt_out -> Branch( "Q2_h", &Q2_h, "Q2_h/F" ) ; // +++
   tt_out -> Branch( "x_h", &x_h, "x_h/F" ) ; // +++

   tt_out -> Branch( "gen_s", &gen_s, "gen_s/F" ) ; // +++
   tt_out -> Branch( "gen_x", &gen_x, "gen_x/F" ) ; // ***

   tt_out -> Branch( "HFS_E", &HFS_E, "HFS_E/F" ) ; // +++



   gDirectory -> Delete( "h*" ) ;

   TH1F* h_gen_q2 = new TH1F("h_gen_g2","Generator Q2", 80, 0., 2000. ) ;
   TH1F* h_gen_x  = new TH1F("h_gen_x","Generator x", 80, -0.1, 1.1 ) ;
   TH1F* h_gen_y  = new TH1F("h_gen_y","Generator y", 80, -0.1, 1.1 ) ;

   TH2F* h_q2_sigma_vs_gen = new TH2F("h_q2_sigma_vs_gen","Q2: Sigma vs gen",80, 0., 2000., 80, 0., 2000.) ;
   TH2F* h_x_sigma_vs_gen = new TH2F("h_x_sigma_vs_gen","x: Sigma vs gen", 80, -0.1, 1.1 , 80, -0.1, 1.1 ) ;
   TH2F* h_y_sigma_vs_gen = new TH2F("h_y_sigma_vs_gen","y: Sigma vs gen", 80, -0.1, 1.1 , 80, -0.1, 1.1 ) ;

   TH1F* h_x_rec_over_true_y_01_05_sigma = new TH1F("h_x_rec_over_true_y_01_05_sigma","x/xtrue, ytrue [0.01,0.05], Sigma", 80, 0., 2.) ;
   TH1F* h_x_rec_over_true_y_05_10_sigma = new TH1F("h_x_rec_over_true_y_05_10_sigma","x/xtrue, ytrue [0.05,0.10], Sigma", 80, 0., 2.) ;
   TH1F* h_x_rec_over_true_y_10_20_sigma = new TH1F("h_x_rec_over_true_y_10_20_sigma","x/xtrue, ytrue [0.10,0.20], Sigma", 80, 0., 2.) ;
   TH1F* h_x_rec_over_true_y_20_50_sigma = new TH1F("h_x_rec_over_true_y_20_50_sigma","x/xtrue, ytrue [0.20,0.50], Sigma", 80, 0., 2.) ;
   TH1F* h_x_rec_over_true_y_50_80_sigma = new TH1F("h_x_rec_over_true_y_50_80_sigma","x/xtrue, ytrue [0.50,0.80], Sigma", 80, 0., 2.) ;

   TH1F* h_x_rec_over_true_y_01_05_e = new TH1F("h_x_rec_over_true_y_01_05_e","x/xtrue, ytrue [0.01,0.05], ele", 80, 0., 2.) ;
   TH1F* h_x_rec_over_true_y_05_10_e = new TH1F("h_x_rec_over_true_y_05_10_e","x/xtrue, ytrue [0.05,0.10], ele", 80, 0., 2.) ;
   TH1F* h_x_rec_over_true_y_10_20_e = new TH1F("h_x_rec_over_true_y_10_20_e","x/xtrue, ytrue [0.10,0.20], ele", 80, 0., 2.) ;
   TH1F* h_x_rec_over_true_y_20_50_e = new TH1F("h_x_rec_over_true_y_20_50_e","x/xtrue, ytrue [0.20,0.50], ele", 80, 0., 2.) ;
   TH1F* h_x_rec_over_true_y_50_80_e = new TH1F("h_x_rec_over_true_y_50_80_e","x/xtrue, ytrue [0.50,0.80], ele", 80, 0., 2.) ;

   TH1F* h_x_rec_over_true_y_01_05_h = new TH1F("h_x_rec_over_true_y_01_05_h","x/xtrue, ytrue [0.01,0.05], had", 80, 0., 2.) ;
   TH1F* h_x_rec_over_true_y_05_10_h = new TH1F("h_x_rec_over_true_y_05_10_h","x/xtrue, ytrue [0.05,0.10], had", 80, 0., 2.) ;
   TH1F* h_x_rec_over_true_y_10_20_h = new TH1F("h_x_rec_over_true_y_10_20_h","x/xtrue, ytrue [0.10,0.20], had", 80, 0., 2.) ;
   TH1F* h_x_rec_over_true_y_20_50_h = new TH1F("h_x_rec_over_true_y_20_50_h","x/xtrue, ytrue [0.20,0.50], had", 80, 0., 2.) ;
   TH1F* h_x_rec_over_true_y_50_80_h = new TH1F("h_x_rec_over_true_y_50_80_h","x/xtrue, ytrue [0.50,0.80], had", 80, 0., 2.) ;

   TH1F* h_x_rec_over_true_y_01_05_da = new TH1F("h_x_rec_over_true_y_01_05_da","x/xtrue, ytrue [0.01,0.05], DA", 80, 0., 2.) ;
   TH1F* h_x_rec_over_true_y_05_10_da = new TH1F("h_x_rec_over_true_y_05_10_da","x/xtrue, ytrue [0.05,0.10], DA", 80, 0., 2.) ;
   TH1F* h_x_rec_over_true_y_10_20_da = new TH1F("h_x_rec_over_true_y_10_20_da","x/xtrue, ytrue [0.10,0.20], DA", 80, 0., 2.) ;
   TH1F* h_x_rec_over_true_y_20_50_da = new TH1F("h_x_rec_over_true_y_20_50_da","x/xtrue, ytrue [0.20,0.50], DA", 80, 0., 2.) ;
   TH1F* h_x_rec_over_true_y_50_80_da = new TH1F("h_x_rec_over_true_y_50_80_da","x/xtrue, ytrue [0.50,0.80], DA", 80, 0., 2.) ;

   Long64_t nentries = fChain->GetEntries() ;

   printf("\n\n Found %lld entries in input file.\n\n", nentries ) ;

   if ( maxEvt > 0 ) { nentries = maxEvt ; }

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      int ei = jentry ;

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if ( !verbose && ei%100 == 0 ) {
         printf(" --- Event: %7d / %lld    %6.3f\r", ei, nentries, (1.*ei)/(1.*nentries) ) ;
         fflush(stdout) ;
      }

      if ( verbose ) printf("\n\n ----- Event : %d\n", ei ) ;


      gene_px = 0. ;
      gene_py = 0. ;
      gene_pz = 0. ;
      gen_e_pt = 0. ;
      gen_e_eta = 0. ;
      gen_e_phi = 0. ;
      gen_e_theta = 0. ;
      gen_e_e = -1. ;

      for ( int pi = 0; pi < Particle_ ; pi ++ ) {

         if ( Particle_PID[pi] == 11 && Particle_Status[pi] == 4 ) {
            if ( verbose ) printf("  %d : found beam electron.  pi = %d\n", ei, pi ) ;
            beam_electron_energy = Particle_E[pi] ;
         }
         if ( Particle_PID[pi] == 2212 && Particle_Status[pi] == 4 ) {
            if ( verbose ) printf("  %d : found beam proton.  pi = %d\n", ei, pi ) ;
            beam_proton_energy = Particle_E[pi] ;
         }
         if ( Particle_PID[pi] == 11 && Particle_Status[pi] == 1 && gen_e_e < 0. ) {
            if ( verbose ) printf("  %d : found scattered electron.  pi = %d, E = %7.1f\n", ei, pi, Particle_E[pi] ) ;
            gene_px = Particle_Px[pi] ;
            gene_py = Particle_Py[pi] ;
            gene_pz = Particle_Pz[pi] ;
            gen_e_pt = Particle_PT[pi] ;
            gen_e_eta = Particle_Eta[pi] ;
            gen_e_phi = Particle_Phi[pi] ;
            gen_e_e = Particle_E[pi] ;
            gen_e_theta = atan2( gen_e_pt, gene_pz ) ;
         }

      } // pi

      if ( beam_electron_energy < 0 ) {
         printf("\n\n *** %3d : can't find beam electron.\n\n", ei ) ;
         gSystem -> Exit(-1) ;
      }
      if ( beam_proton_energy < 0 ) {
         printf("\n\n *** %3d : can't find beam proton.\n\n", ei ) ;
         gSystem -> Exit(-1) ;
      }
      if ( gen_e_e < 0 ) {
         printf("\n\n *** %3d : can't find scattered electron.\n\n", ei ) ;
         gSystem -> Exit(-1) ;
      }

      gen_Q2 = 2. * beam_electron_energy * gen_e_e * ( 1. + cos( gen_e_theta )) ;

      gen_y = 1. - (gen_e_e * (1. - cos(gen_e_theta))) / (2. * beam_electron_energy) ;

      gen_s = 4. * beam_electron_energy * beam_proton_energy ;

      gen_x = gen_Q2 / ( gen_y * gen_s ) ;

      if ( verbose ) printf( "   %3d : generator Q2 = %7.1f , y = %6.3f , x = %6.3f   s = %7.1f\n",
        ei, gen_Q2, gen_y, gen_x, gen_s ) ;



      h_gen_q2 -> Fill( gen_Q2 ) ;
      h_gen_y -> Fill( gen_y ) ;
      h_gen_x -> Fill( gen_x ) ;



      HFS_px = 0. ;
      HFS_py = 0. ;
      HFS_pz = 0. ;
      HFS_E  = 0. ;

      for ( int i = 0; i < EFlowTrack_; i++ ) {

         if ( EFlowTrack_PID[i] == 11 ) {
            float dr = calc_dr( gen_e_phi, EFlowTrack_Phi[i], gen_e_eta, EFlowTrack_Eta[i] ) ;
            if ( verbose ) printf("  %3d : EFlowTrack %3d, pid=11, dR = %7.4f\n", ei, i, dr ) ;
            if ( dr < 0.05 ) continue ;
         }

         float px, py, pz, e ;

         ptep_to_xyze( EFlowTrack_PT[i], EFlowTrack_Eta[i], EFlowTrack_Phi[i], EFlowTrack_Mass[i],
                       px, py, pz, e ) ;

         if ( e < 0 ) continue ;

         if ( verbose ) printf( " %3d, trk %3d : pid = %5d  pt = %7.2f  px,py,pz = %7.2f, %7.2f, %7.2f   E = %7.2f\n",
           ei, i, EFlowTrack_PID[i], EFlowTrack_PT[i], px, py, pz, e ) ;

         HFS_px += px ;
         HFS_py += py ;
         HFS_pz += pz ;
         HFS_E  += e  ;

      } // i

      for ( int i = 0; i < EFlowPhoton_; i++ ) {

         float px, py, pz, e ;

         ptep_to_xyze( EFlowPhoton_ET[i], EFlowPhoton_Eta[i], EFlowPhoton_Phi[i], 0.,
                       px, py, pz, e ) ;

         if ( e < 0 ) continue ;

         if ( verbose ) printf( " %3d, trk %3d :  pt = %7.2f  px,py,pz = %7.2f, %7.2f, %7.2f   E = %7.2f\n",
           ei, i, EFlowPhoton_ET[i], px, py, pz, e ) ;

         HFS_px += px ;
         HFS_py += py ;
         HFS_pz += pz ;
         HFS_E  += e  ;

      } // i

      for ( int i = 0; i < EFlowNeutralHadron_; i++ ) {

         float px, py, pz, e ;

         ptep_to_xyze( EFlowNeutralHadron_ET[i], EFlowNeutralHadron_Eta[i], EFlowNeutralHadron_Phi[i], 0.,
                       px, py, pz, e ) ;

         if ( e < 0 ) continue ;

         if ( verbose ) printf( " %3d, trk %3d :  pt = %7.2f  px,py,pz = %7.2f, %7.2f, %7.2f   E = %7.2f\n",
           ei, i, EFlowNeutralHadron_ET[i], px, py, pz, e ) ;

         HFS_px += px ;
         HFS_py += py ;
         HFS_pz += pz ;
         HFS_E  += e  ;

      } // i

      if ( verbose ) printf( " %3d : sums :  px,py,pz = %7.2f, %7.2f, %7.2f   E = %7.2f\n", ei, HFS_px, HFS_py, HFS_pz, HFS_E ) ;

      float ef_sum_pt = sqrt( HFS_px * HFS_px + HFS_py * HFS_py ) ;



      //-- Get the reconstructed scattered electron

      if ( Electron_ <= 0 ) {
         if ( verbose ) printf(" \n\n *** no reconstructed electrons in event\n\n") ;
         continue ;
      }

      float drmin = 9999999. ;
      int best_ei = -1 ;
      for ( int i = 0; i < Electron_; i ++ ) {
         float dr = calc_dr( gen_e_phi, Electron_Phi[i], gen_e_eta, Electron_Eta[i] ) ;
         if ( verbose ) printf("  %3d : electron %3d : dr = %7.4f\n", ei, i, dr ) ;
         if ( dr < drmin ) {
            drmin = dr ;
            best_ei = i ;
         }
      } // i

      if ( best_ei < 0 ) {
         if ( verbose ) printf("\n\n *** can't find reco scattered electron.\n\n") ;
         continue ;
      }

      rec_e_pt = Electron_PT[best_ei] ;
      rec_e_eta = Electron_Eta[best_ei] ;
      rec_e_phi = Electron_Phi[best_ei] ;

      e_px = rec_e_pt * cos( rec_e_phi ) ;
      e_py = rec_e_pt * sin( rec_e_phi ) ;
      if ( fabs(sin( rec_e_theta )) > 1e-9 ) {
         e_pz = rec_e_pt * ( cos( rec_e_theta ) / sin( rec_e_theta ) ) ;
      } else {
         e_pz = 999. ;
      }

      rec_e_theta = 2. * atan( exp( -1. * rec_e_eta ) ) ;
      rec_e_e = rec_e_pt / sin( rec_e_theta ) ;

      if ( verbose ) printf(" %3d : reco electron,  pt = %7.2f (%7.2f)  e = %7.2f (%7.2f)  theta = %8.5f (%8.5f)\n",
        ei, rec_e_pt, gen_e_pt,  rec_e_e, gen_e_e,  rec_e_theta, gen_e_theta ) ;


      //-- Sigma method

      Sigma = HFS_E - HFS_pz ;

      Empz = Sigma + rec_e_e * (1. - cos(rec_e_theta)) ;

      pth = sqrt( HFS_px * HFS_px + HFS_py * HFS_py ) ;

      y_sigma = Sigma / ( Sigma + rec_e_e * (1. - cos(rec_e_theta)) ) ;

      Q2_sigma = pow( rec_e_e * sin(rec_e_theta) , 2 ) / ( 1. - y_sigma ) ;

      x_sigma = Q2_sigma / ( gen_s * y_sigma ) ;



      //-- Electron method

      y_e = 1. - (rec_e_e/beam_electron_energy) * pow( sin(rec_e_theta/2.), 2 ) ;

      Q2_e = 4. * rec_e_e * beam_electron_energy * pow( cos(rec_e_theta/2.), 2 ) ;

      x_e = Q2_e / ( gen_s * y_e ) ;



      //-- eSigma method

      y_esigma = 2. * Sigma * beam_electron_energy / pow( ( Sigma + rec_e_e * (1. - cos(rec_e_theta)) ), 2 ) ;

      Q2_esigma = Q2_e ;



      //-- hadrons method

      y_h = Sigma / ( 2. * beam_electron_energy ) ;

      Q2_h = pth * pth / ( 1. - y_h ) ;

      x_h = Q2_h / ( gen_s * y_h ) ;


      //-- double angle method

      float tan_gamma_over_2 = Sigma / pth ;

      y_da =  tan_gamma_over_2 / ( tan_gamma_over_2 + tan( rec_e_theta / 2. ) ) ;

      Q2_da = 4. * beam_electron_energy * beam_electron_energy * (1. / tan( rec_e_theta / 2. ))  /  ( tan_gamma_over_2 + tan( rec_e_theta / 2. ) ) ;

      x_da = Q2_da / ( gen_s * y_da ) ;


      if ( verbose ) {
         printf("\n") ;
         printf("     %3d :  Q2 :  true = %7.1f  sigma = %7.1f  e = %7.1f  h = %7.1f  da = %7.1f\n",
            ei, gen_Q2, Q2_sigma, Q2_e, Q2_h, Q2_da ) ;
         printf("     %3d :  x  :  true = %7.3f  sigma = %7.3f  e = %7.3f  h = %7.3f  da = %7.3f\n",
            ei, gen_x, x_sigma, x_e, x_h, x_da  ) ;
         printf("     %3d :  y  :  true = %7.3f  sigma = %7.3f  e = %7.3f  h = %7.3f  da = %7.3f\n",
            ei, gen_y, y_sigma, y_e, y_h, y_da  ) ;
      }


       h_q2_sigma_vs_gen -> Fill( gen_Q2, Q2_sigma ) ;
       h_x_sigma_vs_gen -> Fill( gen_x, x_sigma ) ;
       h_y_sigma_vs_gen -> Fill( gen_y, y_sigma ) ;

       if ( gen_y > 0.01 && gen_y <= 0.05 ) {
           h_x_rec_over_true_y_01_05_sigma -> Fill( x_sigma / gen_x ) ;
           h_x_rec_over_true_y_01_05_e -> Fill( x_e / gen_x ) ;
           h_x_rec_over_true_y_01_05_h -> Fill( x_h / gen_x ) ;
           h_x_rec_over_true_y_01_05_da -> Fill( x_da / gen_x ) ;
       }

       if ( gen_y > 0.05 && gen_y <= 0.10 ) {
           h_x_rec_over_true_y_05_10_sigma -> Fill( x_sigma / gen_x ) ;
           h_x_rec_over_true_y_05_10_e -> Fill( x_e / gen_x ) ;
           h_x_rec_over_true_y_05_10_h -> Fill( x_h / gen_x ) ;
           h_x_rec_over_true_y_05_10_da -> Fill( x_da / gen_x ) ;
       }

       if ( gen_y > 0.10 && gen_y <= 0.20 ) {
           h_x_rec_over_true_y_10_20_sigma -> Fill( x_sigma / gen_x ) ;
           h_x_rec_over_true_y_10_20_e -> Fill( x_e / gen_x ) ;
           h_x_rec_over_true_y_10_20_h -> Fill( x_h / gen_x ) ;
           h_x_rec_over_true_y_10_20_da -> Fill( x_da / gen_x ) ;
       }

       if ( gen_y > 0.20 && gen_y <= 0.50 ) {
           h_x_rec_over_true_y_20_50_sigma -> Fill( x_sigma / gen_x ) ;
           h_x_rec_over_true_y_20_50_e -> Fill( x_e / gen_x ) ;
           h_x_rec_over_true_y_20_50_h -> Fill( x_h / gen_x ) ;
           h_x_rec_over_true_y_20_50_da -> Fill( x_da / gen_x ) ;
       }

       if ( gen_y > 0.50 && gen_y <= 0.80 ) {
           h_x_rec_over_true_y_50_80_sigma -> Fill( x_sigma / gen_x ) ;
           h_x_rec_over_true_y_50_80_e -> Fill( x_e / gen_x ) ;
           h_x_rec_over_true_y_50_80_h -> Fill( x_h / gen_x ) ;
           h_x_rec_over_true_y_50_80_da -> Fill( x_da / gen_x ) ;
       }


       tt_out -> Fill() ;

   } // jentry

   printf("\n\n Done.\n\n") ;

   tt_out -> Write() ;

   tf_out -> Close() ;

} // Loop


//==========

float calc_dr( double phi1, double phi2, double eta1, double eta2 ) {
   float deta = fabs( eta1 - eta2 ) ;
   float dphi = phi1 - phi2 ;
   if ( dphi < -3.14159265 ) dphi += 2*3.14159265 ;
   if ( dphi >  3.14159265 ) dphi -= 2*3.14159265 ;
   return sqrt( deta*deta + dphi*dphi ) ;
}

//==========

void ptep_to_xyze( float pt, double eta, double phi, double m,
                   float& px, float& py, float& pz, float& e ) {

   px = 0. ;
   py = 0. ;
   pz = 0. ;
   e = -1. ;


   float theta = 2. * atan( exp( -1. * eta ) ) ;

   float st = sin( theta ) ;

   if ( fabs( st ) < 1.e-9 ) return ;

   float p = pt / st ;

   e = sqrt( p*p + m*m ) ;

   px = pt * cos( phi ) ;
   py = pt * sin( phi ) ;
   pz = p  * cos( theta ) ;

}






