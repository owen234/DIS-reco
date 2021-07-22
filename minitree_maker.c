#define minitree_maker_cxx
#include "minitree_maker.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

float calc_dr( double phi1, double phi2, double eta1, double eta2 ) ;

void ptep_to_xyze( float pt, double eta, double phi, double m,
                   float& px, float& py, float& pz, float& e ) ;

float eta_acceptance( float eta ) ;

void minitree_maker::Loop( bool verbose, int maxEvt )
{

   TRandom2 tran(12345) ;

   bool  towers_only_for_hfs = true ;

   float gen_HFS_max_eta = 4.0 ;
   //float gen_HFS_max_eta = 9.0 ;

   //bool useEtaTurnoff = true ;
   bool useEtaTurnoff = false ;

   //bool useJets = true ;
   bool useJets = false ;

  //--- nominal minimums
     float minTrkPt = 0.30 ; //--- YR says 150-400 MeV
     float minPhoE = 0.10 ; //--- YR says 50 MeV, use 100 MeV
     float minNHE = 0.50 ; //--- YR says 500 MeV



   //float maxEta = 2.4 ;
   //float maxEta = 2.7 ;
   //float maxEta = 3.3 ;
   //float maxEta = 4.0 ;
     float maxEta = 3.3 ;
   //float maxEta = 3.26 ;

   //float minTrkPt = 0.50 ;
   //float minPhoE = 0.10 ;
   //float minNHE = 0.50 ;

   //--- hardmin1
   //float minTrkPt = 0.50 ;
   //float minPhoE = 0.50 ;
   //float minNHE = 1.50 ;

   //--- hardmin2
   //float minTrkPt = 1.50 ;
   //float minPhoE = 1.50 ;
   //float minNHE = 4.50 ;



   bool includeExtraTTreeVars = true ;

   if (fChain == 0) return;

   TFile* tf_out = new TFile("mini-tree.root", "recreate") ;
   TTree* tt_out = new TTree("minitree", "Minimal flat TTree for hadronic Q2,x,y analysis" ) ;

   float beam_electron_energy ;
   float beam_proton_energy ;

   float gene_px ;
   float gene_py ;
   float gene_pz ;

   float gen_e_e ;
   float gen_e_pt ;
   float gen_e_eta ;
   float gen_e_phi ;
   float gen_e_theta ;

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

   float HFS_pt ;
   float HFS_eta ;
   float HFS_phi ;
   float HFS_theta ;
   float HFS_gamma ;

   float gen_HFS_px ;
   float gen_HFS_py ;
   float gen_HFS_pz ;
   float gen_HFS_e ;

   float gen_HFS_pt ;
   float gen_HFS_eta ;
   float gen_HFS_phi ;
   float gen_HFS_theta ;



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

   if ( includeExtraTTreeVars ) {
      tt_out -> Branch( "gen_e_e", &gen_e_e, "gen_e_e/F" ) ;
      tt_out -> Branch( "gen_e_pt", &gen_e_pt, "gen_e_pt/F" ) ;
      tt_out -> Branch( "gen_e_eta", &gen_e_eta, "gen_e_eta/F" ) ;
      tt_out -> Branch( "rec_e_eta", &rec_e_eta, "rec_e_eta/F" ) ;
      tt_out -> Branch( "gen_e_phi", &gen_e_phi, "gen_e_phi/F" ) ;
      tt_out -> Branch( "gen_e_theta", &gen_e_theta, "gen_e_theta/F" ) ;
      tt_out -> Branch( "HFS_pt", &HFS_pt, "HFS_pt/F" ) ;
      tt_out -> Branch( "HFS_eta", &HFS_eta, "HFS_eta/F" ) ;
      tt_out -> Branch( "HFS_phi", &HFS_phi, "HFS_phi/F" ) ;
      tt_out -> Branch( "HFS_theta", &HFS_theta, "HFS_theta/F" ) ;
      tt_out -> Branch( "HFS_gamma", &HFS_gamma, "HFS_gamma/F" ) ;
      tt_out -> Branch( "gen_HFS_px", &gen_HFS_px, "gen_HFS_px/F" ) ;
      tt_out -> Branch( "gen_HFS_py", &gen_HFS_py, "gen_HFS_py/F" ) ;
      tt_out -> Branch( "gen_HFS_pz", &gen_HFS_pz, "gen_HFS_pz/F" ) ;
      tt_out -> Branch( "gen_HFS_e", &gen_HFS_e, "gen_HFS_e/F" ) ;
      tt_out -> Branch( "gen_HFS_pt", &gen_HFS_pt, "gen_HFS_pt/F" ) ;
      tt_out -> Branch( "gen_HFS_eta", &gen_HFS_eta, "gen_HFS_eta/F" ) ;
      tt_out -> Branch( "gen_HFS_phi", &gen_HFS_phi, "gen_HFS_phi/F" ) ;
      tt_out -> Branch( "gen_HFS_theta", &gen_HFS_theta, "gen_HFS_theta/F" ) ;
      tt_out -> Branch( "beam_electron_energy", &beam_electron_energy, "beam_electron_energy/F" ) ;
      tt_out -> Branch( "beam_proton_energy", &beam_proton_energy, "beam_proton_energy/F" ) ;
   }


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

      gen_HFS_px = 0. ;
      gen_HFS_py = 0. ;
      gen_HFS_pz = 0. ;
      gen_HFS_e = 0. ;

      for ( int pi = 0; pi < Particle_ ; pi ++ ) {

         if ( Particle_PID[pi] == 11 && Particle_Status[pi] == 4 ) {
            if ( verbose ) printf("  %d : found beam electron.  pi = %d\n", ei, pi ) ;
            beam_electron_energy = Particle_E[pi] ;
         }
         if ( Particle_PID[pi] == 2212 && Particle_Status[pi] == 4 ) {
            if ( verbose ) printf("  %d : found beam proton.  pi = %d\n", ei, pi ) ;
            beam_proton_energy = Particle_E[pi] ;
         }
         if ( Particle_Status[pi] == 1 ) {

            if ( Particle_PID[pi] == 11 && gen_e_e < 0. ) {

               if ( verbose ) printf("  %d : found scattered electron.  pi = %d, E = %7.1f, eta,phi = %8.3f, %8.3f\n", 
                    ei, pi, Particle_E[pi], Particle_Eta[pi], Particle_Phi[pi] ) ;
               gene_px = Particle_Px[pi] ;
               gene_py = Particle_Py[pi] ;
               gene_pz = Particle_Pz[pi] ;
               gen_e_pt = Particle_PT[pi] ;
               gen_e_eta = Particle_Eta[pi] ;
               gen_e_phi = Particle_Phi[pi] ;
               gen_e_e = Particle_E[pi] ;
               gen_e_theta = atan2( gen_e_pt, gene_pz ) ;

            } else {

               if ( fabs(Particle_Eta[pi]) < gen_HFS_max_eta ) {
                  gen_HFS_px += Particle_Px[pi] ;
                  gen_HFS_py += Particle_Py[pi] ;
                  gen_HFS_pz += Particle_Pz[pi] ;
                  gen_HFS_e  += Particle_E[pi] ;
                  if ( verbose ) {
                     printf("    MC sum:  %3d :  %7d   :  px,py,pz = %7.2f, %7.2f, %7.2f    pt = %7.2f, eta = %8.3f, phi = %8.3f\n",
                       pi, Particle_PID[pi], Particle_Px[pi], Particle_Py[pi], Particle_Pz[pi], Particle_PT[pi], Particle_Eta[pi], Particle_Phi[pi] ) ;
                  }
               }

            }
         } // status 1 ?



      } // pi

      if ( verbose ) printf( " %3d : gen  sums :  px,py,pz = %7.2f, %7.2f, %7.2f   E = %7.2f\n", ei, gen_HFS_px, gen_HFS_py, gen_HFS_pz, gen_HFS_e ) ;

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


      float mcmatch_sum_px(0.) ;
      float mcmatch_sum_py(0.) ;
      float mcmatch_sum_pz(0.) ;



      if ( !useJets ) {


         if ( towers_only_for_hfs ) {

   //==== Towers +++++++++++++++++++++++++++++++++++++++++

            int ti_ele = -1 ;
            float min_dr_with_ele = 99999. ;

            for ( int i = 0; i < Tower_; i++ ) {

               float dr_with_ele = calc_dr( gen_e_phi, Tower_Phi[i], gen_e_eta, Tower_Eta[i] ) ;
               if ( dr_with_ele < min_dr_with_ele ) {
                  min_dr_with_ele = dr_with_ele ;
                  ti_ele = i ;
               }

            } // i

            if ( ti_ele < 0 ) { printf("\n\n *** can't find electron tower.\n\n") ; gSystem->Exit(-1) ; }

            if ( verbose ) printf("  %3d : ele tower %3d :  dr with electron = %8.4f,  E = %7.1f,  eta,phi = %8.3f, %8.3f\n",
              ei, ti_ele, min_dr_with_ele, Tower_E[ti_ele], Tower_Eta[ti_ele], Tower_Phi[ti_ele] ) ;

            for ( int i = 0; i < Tower_; i++ ) {

               if ( i == ti_ele ) continue ;

               if ( fabs( Tower_Eta[i] ) > maxEta ) continue ;

               float px, py, pz, e ;
               ptep_to_xyze( Tower_ET[i], Tower_Eta[i], Tower_Phi[i], 0.,
                             px, py, pz, e ) ;

               if ( e < 0 ) continue ;

               if ( verbose ) printf( " %3d, tower %3d :  pt = %7.2f  px,py,pz = %7.2f, %7.2f, %7.2f   E = %7.2f    eta = %8.3f, phi = %8.3f\n",
                 ei, i, Tower_ET[i], px, py, pz, e, Tower_Eta[i], Tower_Phi[i] ) ;

               if ( useEtaTurnoff ) {
                  float prob = eta_acceptance( Tower_Eta[i] ) ;
                  float rn = tran.Uniform() ;
                  if ( rn > prob ) continue ;
               }

               HFS_px += px ;
               HFS_py += py ;
               HFS_pz += pz ;
               HFS_E  += e  ;

            } // i

         } else {

   //==== EF Candidates +++++++++++++++++++++++++++++++++++++++++

            for ( int i = 0; i < EFlowTrack_; i++ ) {

               if ( EFlowTrack_PID[i] == 11 ) {
                  float dr = calc_dr( gen_e_phi, EFlowTrack_Phi[i], gen_e_eta, EFlowTrack_Eta[i] ) ;
                  if ( verbose ) printf("  %3d : EFlowTrack %3d, pid=11, dR = %7.4f\n", ei, i, dr ) ;
                  if ( dr < 0.05 ) continue ;
               }

               float px, py, pz, e ;

               if ( fabs(EFlowTrack_Eta[i]) > 4.0 ) continue ;

               if ( EFlowTrack_PT[i] < minTrkPt ) continue ;

               ptep_to_xyze( EFlowTrack_PT[i], EFlowTrack_Eta[i], EFlowTrack_Phi[i], EFlowTrack_Mass[i],
                             px, py, pz, e ) ;

               if ( e < 0 ) continue ;


               if ( verbose ) printf( " %3d, trk %3d : pid = %5d  pt = %7.2f  px,py,pz = %7.2f, %7.2f, %7.2f   E = %7.2f    eta = %8.3f, phi = %8.3f\n",
                 ei, i, EFlowTrack_PID[i], EFlowTrack_PT[i], px, py, pz, e, EFlowTrack_Eta[i], EFlowTrack_Phi[i] ) ;

               if ( verbose ) {
                  double min_mc_dr = 99999. ;
                  int mc_pi = -1 ;
                  for ( int pi = 0; pi < Particle_ ; pi ++ ) {
                     if  ( Particle_Status[pi] != 1 ) continue ;
                     if ( fabs(Particle_Charge[pi]) != 1 ) continue ;
                     float dr = calc_dr( EFlowTrack_Phi[i], Particle_Phi[pi], EFlowTrack_Eta[i], Particle_Eta[pi] ) ;
                     if ( dr < min_mc_dr ) {
                        min_mc_dr = dr ;
                        mc_pi = pi ;
                     }
                  } // pi
                  if ( mc_pi >= 0 && min_mc_dr < 0.05 ) {
                     printf( " %3d, trk %3d :                  MC match  px,py,pz = %7.2f, %7.2f, %7.2f,  dr = %8.4f, eta = %8.3f, phi = %8.3f\n",
                        ei, i, Particle_Px[mc_pi], Particle_Py[mc_pi], Particle_Pz[mc_pi], min_mc_dr, Particle_Eta[mc_pi], Particle_Phi[mc_pi] ) ;
                     mcmatch_sum_px += Particle_Px[mc_pi] ;
                     mcmatch_sum_py += Particle_Py[mc_pi] ;
                     mcmatch_sum_pz += Particle_Pz[mc_pi] ;
                  }
               }

               //if ( EFlowTrack_Eta[i] > maxEta ) continue ;
               if ( fabs(EFlowTrack_Eta[i]) > maxEta ) continue ;

               if ( useEtaTurnoff ) {
                  float prob = eta_acceptance( EFlowTrack_Eta[i] ) ;
                  float rn = tran.Uniform() ;
                  if ( rn > prob ) continue ;
               }

               HFS_px += px ;
               HFS_py += py ;
               HFS_pz += pz ;
               HFS_E  += e  ;

            } // i

            if ( verbose ) printf( "  --- after tracks:  HFS :  px,py,pz = %7.2f , %7.2f, %7.2f   E = %7.2f\n", HFS_px, HFS_py, HFS_pz, HFS_E ) ;





            for ( int i = 0; i < EFlowPhoton_; i++ ) {

               float px, py, pz, e ;

               if ( EFlowPhoton_Eta[i] > 4.0 ) continue ;

               if ( EFlowPhoton_E[i] < minPhoE ) continue ;

               ptep_to_xyze( EFlowPhoton_ET[i], EFlowPhoton_Eta[i], EFlowPhoton_Phi[i], 0.,
                             px, py, pz, e ) ;

               if ( e < 0 ) continue ;

               if ( verbose ) printf( " %3d, pho %3d :  pt = %7.2f  px,py,pz = %7.2f, %7.2f, %7.2f   E = %7.2f    eta = %8.3f, phi = %8.3f\n",
                 ei, i, EFlowPhoton_ET[i], px, py, pz, e, EFlowPhoton_Eta[i], EFlowPhoton_Phi[i] ) ;

               if ( verbose ) {
                  double min_mc_dr = 99999. ;
                  int mc_pi = -1 ;
                  for ( int pi = 0; pi < Particle_ ; pi ++ ) {
                     if  ( Particle_Status[pi] != 1 ) continue ;
                     if ( fabs(Particle_PID[pi]) != 22 ) continue ;
                     if ( fabs(Particle_Charge[pi]) != 0 ) continue ;
                     float dr = calc_dr( EFlowPhoton_Phi[i], Particle_Phi[pi], EFlowPhoton_Eta[i], Particle_Eta[pi] ) ;
                     if ( dr < min_mc_dr ) {
                        min_mc_dr = dr ;
                        mc_pi = pi ;
                     }
                  } // pi
                  if ( mc_pi >= 0 && min_mc_dr < 0.05 ) {
                     printf( " %3d, pho %3d :      MC match  px,py,pz = %7.2f, %7.2f, %7.2f,  dr = %8.4f, eta = %8.3f, phi = %8.3f\n",
                        ei, i, Particle_Px[mc_pi], Particle_Py[mc_pi], Particle_Pz[mc_pi], min_mc_dr, Particle_Eta[mc_pi], Particle_Phi[mc_pi] ) ;
                     mcmatch_sum_px += Particle_Px[mc_pi] ;
                     mcmatch_sum_py += Particle_Py[mc_pi] ;
                     mcmatch_sum_pz += Particle_Pz[mc_pi] ;
                  }
               }

               //if ( EFlowPhoton_Eta[i] > maxEta ) continue ;
               if ( fabs(EFlowPhoton_Eta[i]) > maxEta ) continue ;

               if ( useEtaTurnoff ) {
                  float prob = eta_acceptance( EFlowPhoton_Eta[i] ) ;
                  float rn = tran.Uniform() ;
                  if ( rn > prob ) continue ;
               }

               HFS_px += px ;
               HFS_py += py ;
               HFS_pz += pz ;
               HFS_E  += e  ;

            } // i


            if ( verbose ) printf( "  --- after photons:  HFS :  px,py,pz = %7.2f , %7.2f, %7.2f   E = %7.2f\n", HFS_px, HFS_py, HFS_pz, HFS_E ) ;





            for ( int i = 0; i < EFlowNeutralHadron_; i++ ) {

               float px, py, pz, e ;

               if ( fabs(EFlowNeutralHadron_Eta[i]) > 4.0 ) continue ;

               if ( EFlowNeutralHadron_E[i] < minNHE ) continue ;

               ptep_to_xyze( EFlowNeutralHadron_ET[i], EFlowNeutralHadron_Eta[i], EFlowNeutralHadron_Phi[i], 0.,
                             px, py, pz, e ) ;

               if ( e < 0 ) continue ;

               if ( verbose ) printf( " %3d, nh  %3d :  pt = %7.2f  px,py,pz = %7.2f, %7.2f, %7.2f   E = %7.2f    eta = %8.3f, phi = %8.3f\n",
                 ei, i, EFlowNeutralHadron_ET[i], px, py, pz, e, EFlowNeutralHadron_Eta[i], EFlowNeutralHadron_Phi[i] ) ;

               if ( verbose ) {
                  double min_mc_dr = 99999. ;
                  int mc_pi = -1 ;
                  for ( int pi = 0; pi < Particle_ ; pi ++ ) {
                     if  ( Particle_Status[pi] != 1 ) continue ;
                     if ( fabs(Particle_PID[pi]) == 22 ) continue ;
                     if ( fabs(Particle_Charge[pi]) != 0 ) continue ;
                     float dr = calc_dr( EFlowNeutralHadron_Phi[i], Particle_Phi[pi], EFlowNeutralHadron_Eta[i], Particle_Eta[pi] ) ;
                     if ( dr < min_mc_dr ) {
                        min_mc_dr = dr ;
                        mc_pi = pi ;
                     }
                  } // pi
                  if ( mc_pi >= 0 && min_mc_dr < 0.05 ) {
                     printf( " %3d, nh  %3d :      MC match  px,py,pz = %7.2f, %7.2f, %7.2f,  dr = %8.4f, eta = %8.3f, phi = %8.3f  PID = %d\n",
                        ei, i, Particle_Px[mc_pi], Particle_Py[mc_pi], Particle_Pz[mc_pi], min_mc_dr, Particle_Eta[mc_pi], Particle_Phi[mc_pi], Particle_PID[mc_pi] ) ;
                     mcmatch_sum_px += Particle_Px[mc_pi] ;
                     mcmatch_sum_py += Particle_Py[mc_pi] ;
                     mcmatch_sum_pz += Particle_Pz[mc_pi] ;
                  }
               }


               //if ( EFlowNeutralHadron_Eta[i] > maxEta ) continue ;
               if ( fabs(EFlowNeutralHadron_Eta[i]) > maxEta ) continue ;

               if ( useEtaTurnoff ) {
                  float prob = eta_acceptance( EFlowNeutralHadron_Eta[i] ) ;
                  float rn = tran.Uniform() ;
                  if ( rn > prob ) continue ;
               }

               HFS_px += px ;
               HFS_py += py ;
               HFS_pz += pz ;
               HFS_E  += e  ;

            } // i

            if ( verbose ) printf( "  --- after NH    :  HFS :  px,py,pz = %7.2f , %7.2f, %7.2f   E = %7.2f\n", HFS_px, HFS_py, HFS_pz, HFS_E ) ;


         }


      } else {

   //==== Jets +++++++++++++++++++++++++++++++++++++++++

         for ( int i = 0; i < Jet_; i++ ) {

            float px, py, pz, e ;

            if ( Jet_Eta[i] > 4.0 ) continue ;

            ptep_to_xyze( Jet_PT[i], Jet_Eta[i], Jet_Phi[i], Jet_Mass[i],
                          px, py, pz, e ) ;

            if ( e < 0 ) continue ;

            if ( verbose ) printf( " %3d, jet %3d :  pt = %7.2f, eta = %7.3f  px,py,pz = %7.2f, %7.2f, %7.2f   E = %7.2f\n",
              ei, i, Jet_PT[i], Jet_Eta[i],  px, py, pz, e ) ;

            if ( fabs(Jet_Eta[i]) > maxEta ) continue ;

            if ( useEtaTurnoff ) {
               float prob = eta_acceptance( Jet_Eta[i] ) ;
               float rn = tran.Uniform() ;
               if ( rn > prob ) continue ;
            }

            HFS_px += px ;
            HFS_py += py ;
            HFS_pz += pz ;
            HFS_E  += e  ;

         } // i

      } // useJets?


      if ( verbose ) printf( " %3d : reco  sums :  px,py,pz = %7.2f, %7.2f, %7.2f   E = %7.2f\n", ei, HFS_px, HFS_py, HFS_pz, HFS_E ) ;
      if ( verbose ) printf( " %3d : gen   sums :  px,py,pz = %7.2f, %7.2f, %7.2f   E = %7.2f\n", ei, gen_HFS_px, gen_HFS_py, gen_HFS_pz, gen_HFS_e ) ;
      if ( verbose ) printf( " %3d : match sums :  px,py,pz = %7.2f, %7.2f, %7.2f   \n", ei, mcmatch_sum_px, mcmatch_sum_py, mcmatch_sum_pz ) ;

      if ( HFS_E <= 0. ) {
         if ( verbose ) printf(" %3d :  empty HFS.  Skipping event.\n", ei ) ;
         continue ;
      }

      float ef_sum_pt = sqrt( HFS_px * HFS_px + HFS_py * HFS_py ) ;

      if ( includeExtraTTreeVars ) {
         HFS_pt = ef_sum_pt ;
         HFS_theta = atan2( HFS_pt, HFS_pz ) ;
         HFS_phi = atan2( HFS_py, HFS_px ) ;
         HFS_eta = -1. * log( tan( HFS_theta/2. ) ) ;


         gen_HFS_pt = sqrt( gen_HFS_px * gen_HFS_px + gen_HFS_py * gen_HFS_py ) ;
         gen_HFS_theta = atan2( gen_HFS_pt, gen_HFS_pz ) ;
         gen_HFS_phi = atan2( gen_HFS_py, gen_HFS_px ) ;
         gen_HFS_eta = -1. * log( tan( gen_HFS_theta/2. ) ) ;
      }




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

      ///e_px = rec_e_pt * cos( rec_e_phi ) ;
      ///e_py = rec_e_pt * sin( rec_e_phi ) ;
      ///if ( fabs(sin( rec_e_theta )) > 1e-9 ) {
      ///   e_pz = rec_e_pt * ( cos( rec_e_theta ) / sin( rec_e_theta ) ) ;
      ///} else {
      ///   e_pz = 999. ;
      ///}

      ptep_to_xyze( rec_e_pt, rec_e_eta, rec_e_phi, 0.000511,
                          e_px, e_py, e_pz, rec_e_e ) ;


      rec_e_theta = 2. * atan( exp( -1. * rec_e_eta ) ) ;
      float rec_e_theta_check = atan2( rec_e_pt, e_pz ) ;
      if (verbose) printf(" *** Check: 2. * atan( exp( -1. * rec_e_eta ) ) = %8.4f,  atan2( rec_e_pt, rec_e_pz ) = %8.4f\n", rec_e_theta, rec_e_theta_check ) ;
      ///rec_e_e = rec_e_pt / sin( rec_e_theta ) ;

      if ( verbose ) printf(" %3d : reco electron,  pt = %7.2f (%7.2f)  e = %7.2f (%7.2f)  theta = %8.5f (%8.5f)\n",
        ei, rec_e_pt, gen_e_pt,  rec_e_e, gen_e_e,  rec_e_theta, gen_e_theta ) ;
      if ( verbose ) printf(" %3d : reco electron,  px = %7.2f (%7.2f)  py = %7.2f (%7.2f)  pz = %7.2f (%7.2f)\n",
        ei, e_px, gene_px,  e_py, gene_py,  e_pz, gene_pz ) ;


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

      HFS_gamma = 2 * atan( tan_gamma_over_2 ) ;

      y_da =  tan_gamma_over_2 / ( tan_gamma_over_2 + tan( rec_e_theta / 2. ) ) ;

      Q2_da = 4. * beam_electron_energy * beam_electron_energy * (1. / tan( rec_e_theta / 2. ))  /  ( tan_gamma_over_2 + tan( rec_e_theta / 2. ) ) ;

      x_da = Q2_da / ( gen_s * y_da ) ;



      if ( verbose ) {
         printf("\n") ;

         float gamma = 2. * atan( tan_gamma_over_2 ) ;


         printf("     %3d :  HFS   px,py,pz = %7.2f, %7.2f, %7.2f    E = %7.2f   pt,eta,phi,theta = %7.2f, %8.4f, %8.4f, %8.4f\n",
               ei, HFS_px, HFS_py, HFS_pz, HFS_E, HFS_pt, HFS_eta, HFS_phi, HFS_theta ) ;
         printf("     %3d :  tan(gamma/2) = %8.4f,  gamma = %8.4f,  HFS_theta = %8.4f\n", ei, tan_gamma_over_2, gamma, HFS_theta ) ;

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


   ///// float theta = 2. * atan( exp( -1. * eta ) ) ;

   ///////float st = sin( theta ) ;

   ///////if ( fabs( st ) < 1.e-9 ) return ;

   ///////float p = pt / st ;

   /////e = sqrt( p*p + m*m ) ;

   ///////pz = p  * cos( theta ) ;

   px = pt * cos( phi ) ;
   py = pt * sin( phi ) ;
   pz = pt  * sinh( eta ) ;

   float p2 = px*px + py*py + pz*pz ;
   e = sqrt( p2 + m*m ) ;

   //////printf("  *** DEBUG:  px,py,pz = %7.2f, %7.2f, %7.2f   E = %7.2f\n", px, py, pz, e ) ;


}


//==========

float eta_acceptance( float eta ) {

   //--etf1
   //return 1./(1.+exp(-9.0*(3.3-eta))) ;

   //--etf2
   //return 1./(1.+exp(-5.0*(2.9-eta))) ;

   //--etf3
   return 1./(1.+exp(-5.0*(2.4-eta))) ;

} // eta_acceptance







