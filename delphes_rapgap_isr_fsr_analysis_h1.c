#define delphes_rapgap_isr_fsr_analysis_h1_cxx
#include "delphes_rapgap_isr_fsr_analysis_h1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

float calc_dr( double phi1, double phi2, double eta1, double eta2 ) ;

float calc_angle( double pt1, double eta1, double phi1,
                  double pt2, double eta2, double phi2 ) ;

void ptep_to_xyze( float pt, double eta, double phi, double m,
                   float& px, float& py, float& pz, float& e ) ;

void smear_theta( float pt0, float eta0, float phi0,
                  float& pts, float& etas, float& phis, bool verbose ) ;

TRandom3* tran ;

void delphes_rapgap_isr_fsr_analysis_h1::Loop( bool verbose, int maxEvts )
{

   tran = new TRandom3() ;

   const int nmeth = 9 ;

   float ele_min_pt = 5.0 ;

   float maxEta = 4.0 ;
   float minTrkPt = 0.15 ; //--- YR says 150-400 MeV
   float minPhoE = 0.05 ; //--- YR says 50 MeV
   float minNHE = 0.50 ; //--- YR says 500 MeV

   float add_noise = false ;
   //float add_noise = true ;

   float noise_landau_mu = 0.0 ;
   float noise_landau_sigma = 0.05 ;
   float noise_prob = 1.0 ;

   float use_electron_cand_not_track = true ;

   bool obs_recombine_fsr_gamma = false ;

   //bool smear_electron_theta = true ;
   bool smear_electron_theta = false ;

   bool save_only_dnn_inputs = true ;
   //bool save_only_dnn_inputs = false ;

   //bool recombine_fsr_gamma = true ;
   bool recombine_fsr_gamma = false ; // for v4

   //bool include_rad_gamma_in_hfs = true ;
   bool include_rad_gamma_in_hfs = false ; // for v4

   if (fChain == 0) return;

   //float escale = 1000. ;
   float escale = 1. ;

   char out_fname[1000] ;
   //if ( save_only_dnn_inputs ) {
   //   sprintf( out_fname, "dnn-inputs.root" ) ;
   //} else {
   //   sprintf( out_fname, "rad-tree.root" ) ;
   //}
   sprintf( out_fname, "%s", outfile_name ) ;

   TFile* tf_out = new TFile( out_fname, "recreate") ;
   TTree* tt_out = new TTree("minitree", "Minimal flat TTree for hadronic Q2,x,y analysis" ) ;

   float gen_y[nmeth] ;
   float gen_x[nmeth] ;
   float gen_Q2[nmeth] ;

   float obs_y[nmeth] ;
   float obs_x[nmeth] ;
   float obs_Q2[nmeth] ;

   char  meth_name[nmeth][100] ;


   float be0_e ;
   float be0_pt ;
   float be0_pz ;
   float be0_phi ;

   float bef_e ;
   float bef_pt ;
   float bef_pz ;
   float bef_phi ;

   float se0_e ;
   float se0_pt ;
   float se0_pz ;
   float se0_phi ;

   float sef_e ;
   float sef_pt ;
   float sef_pz ;
   float sef_phi ;
   float sef_eta ;

   float beam_e_e ;
   float beam_p_e ;

   float gen_e0_Q2 ;
   float gen_e0_y ;
   float gen_e0_x ;

   float gen_ef_Q2 ;
   float gen_ef_y ;
   float gen_ef_x ;

   bool  has_isr ;
   bool  has_fsr ;

   bool  gen_fsr_photon_recombined ;
   bool  obs_fsr_photon_recombined ;

   float gen_hfs_Sigma_no_rad  ;
   float gen_hfs_T_no_rad  ;
   float gen_hfs_E_no_rad  ;
   float gen_hfs_pz_no_rad  ;
   float gen_hfs_phi_no_rad  ;

   float gen_hfs_Sigma_with_rad  ;
   float gen_hfs_T_with_rad  ;
   float gen_hfs_E_with_rad  ;
   float gen_hfs_pz_with_rad  ;
   float gen_hfs_phi_with_rad  ;

   float gen_hfs_Sigma_no_rad_eta_lt_4  ;
   float gen_hfs_T_no_rad_eta_lt_4  ;
   float gen_hfs_E_no_rad_eta_lt_4  ;
   float gen_hfs_pz_no_rad_eta_lt_4  ;
   float gen_hfs_phi_no_rad_eta_lt_4  ;

   float gen_hfs_Sigma_with_rad_eta_lt_4  ;
   float gen_hfs_T_with_rad_eta_lt_4  ;
   float gen_hfs_E_with_rad_eta_lt_4  ;
   float gen_hfs_pz_with_rad_eta_lt_4  ;
   float gen_hfs_phi_with_rad_eta_lt_4  ;

   float from_tlv_gen_Q2 ;
   float from_tlv_gen_x ;
   float from_tlv_gen_y ;

   float fsr_gamma_dr ;
   float fsr_gamma_e ;
   float fsr_gamma_pt ;
   float fsr_gamma_pz ;
   float fsr_gamma_phi ;
   float fsr_gamma_angle ;
   float fsr_gamma_angle_wrt_beam ;
   float fsr_gamma_eta ;

   float isr_gamma_dr ;
   float isr_gamma_e ;
   float isr_gamma_pt ;
   float isr_gamma_pz ;
   float isr_gamma_phi ;
   float isr_gamma_angle ;
   float isr_gamma_angle_wrt_se ;
   float isr_gamma_eta ;


  //--- DNN inputs
   float gen_obs_e_e ;
   float gen_obs_e_pz ;
   float gen_obs_e_pt ;
   float gen_obs_e_phi ;
   float gen_obs_e_eta ;

   float gen_obs_hfs_e ;
   float gen_obs_hfs_pz ;
   float gen_obs_hfs_pt ;
   float gen_obs_hfs_phi ;

   float obs_e_e ;
   float obs_e_pz ;
   float obs_e_pt ;
   float obs_e_phi ;
   float obs_e_eta ;

   float obs_hfs_e ;
   float obs_hfs_pz ;
   float obs_hfs_pt ;
   float obs_hfs_phi ;

   float gen_obs_dphi_hfs_eta_lt_4 ;

   float obs_dphi ;


   int   n_pho ;
   vector<bool> pho_isr_match ;
   vector<bool> pho_fsr_match ;
   vector<float> pho_angle_wrt_ebeam ;
   vector<float> pho_angle_wrt_efinal ;
   vector<float> pho_e ;
   vector<float> pho_pt ;
   vector<float> pho_pz ;
   vector<float> pho_eta ;
   vector<float> pho_dphi_wrt_efinal ;


   float tower_sum_07 ;
   float tower_sum_20 ;
   float tower_sum_30 ;
   float tower_sum_40 ;

   int   n_towers_07 ;
   int   n_towers_20 ;
   int   n_towers_30 ;
   int   n_towers_40 ;


   float phi_pho_closest_to_ebeam ;
   float eta_pho_closest_to_ebeam ;
   float e_pho_closest_to_ebeam ;

   bool  e_is_correct_track ;



   char bname[100] ;



   if ( ! save_only_dnn_inputs ) {

      tt_out -> Branch( "be0_e", &be0_e, "be0_e/F" ) ;
      tt_out -> Branch( "be0_pt", &be0_pt, "be0_pt/F" ) ;
      tt_out -> Branch( "be0_pz", &be0_pz, "be0_pz/F" ) ;
      tt_out -> Branch( "be0_phi", &be0_phi, "be0_phi/F" ) ;

      tt_out -> Branch( "bef_e", &bef_e, "bef_e/F" ) ;
      tt_out -> Branch( "bef_pt", &bef_pt, "bef_pt/F" ) ;
      tt_out -> Branch( "bef_pz", &bef_pz, "bef_pz/F" ) ;
      tt_out -> Branch( "bef_phi", &bef_phi, "bef_phi/F" ) ;

      tt_out -> Branch( "se0_e", &se0_e, "se0_e/F" ) ;
      tt_out -> Branch( "se0_pt", &se0_pt, "se0_pt/F" ) ;
      tt_out -> Branch( "se0_pz", &se0_pz, "se0_pz/F" ) ;
      tt_out -> Branch( "se0_phi", &se0_phi, "se0_phi/F" ) ;

      tt_out -> Branch( "gen_e0_Q2", &gen_e0_Q2, "gen_e0_Q2/F" ) ;
      tt_out -> Branch( "gen_e0_y", &gen_e0_y, "gen_e0_y/F" ) ;
      tt_out -> Branch( "gen_e0_x", &gen_e0_x, "gen_e0_x/F" ) ;

      tt_out -> Branch( "gen_ef_Q2", &gen_ef_Q2, "gen_ef_Q2/F" ) ;
      tt_out -> Branch( "gen_ef_y", &gen_ef_y, "gen_ef_y/F" ) ;
      tt_out -> Branch( "gen_ef_x", &gen_ef_x, "gen_ef_x/F" ) ;

      tt_out -> Branch( "gen_hfs_Sigma_no_rad", &gen_hfs_Sigma_no_rad, "gen_hfs_Sigma_no_rad/F" ) ;
      tt_out -> Branch( "gen_hfs_T_no_rad", &gen_hfs_T_no_rad, "gen_hfs_T_no_rad/F" ) ;
      tt_out -> Branch( "gen_hfs_E_no_rad", &gen_hfs_E_no_rad, "gen_hfs_E_no_rad/F" ) ;
      tt_out -> Branch( "gen_hfs_pz_no_rad", &gen_hfs_pz_no_rad, "gen_hfs_pz_no_rad/F" ) ;
      tt_out -> Branch( "gen_hfs_phi_no_rad", &gen_hfs_phi_no_rad, "gen_hfs_phi_no_rad/F" ) ;

      tt_out -> Branch( "gen_hfs_Sigma_with_rad", &gen_hfs_Sigma_with_rad, "gen_hfs_Sigma_with_rad/F" ) ;
      tt_out -> Branch( "gen_hfs_T_with_rad", &gen_hfs_T_with_rad, "gen_hfs_T_with_rad/F" ) ;
      tt_out -> Branch( "gen_hfs_E_with_rad", &gen_hfs_E_with_rad, "gen_hfs_E_with_rad/F" ) ;
      tt_out -> Branch( "gen_hfs_pz_with_rad", &gen_hfs_pz_with_rad, "gen_hfs_pz_with_rad/F" ) ;
      tt_out -> Branch( "gen_hfs_phi_with_rad", &gen_hfs_phi_with_rad, "gen_hfs_phi_with_rad/F" ) ;

      tt_out -> Branch( "gen_hfs_Sigma_no_rad_eta_lt_4", &gen_hfs_Sigma_no_rad_eta_lt_4, "gen_hfs_Sigma_no_rad_eta_lt_4/F" ) ;
      tt_out -> Branch( "gen_hfs_T_no_rad_eta_lt_4", &gen_hfs_T_no_rad_eta_lt_4, "gen_hfs_T_no_rad_eta_lt_4/F" ) ;
      tt_out -> Branch( "gen_hfs_E_no_rad_eta_lt_4", &gen_hfs_E_no_rad_eta_lt_4, "gen_hfs_E_no_rad_eta_lt_4/F" ) ;
      tt_out -> Branch( "gen_hfs_pz_no_rad_eta_lt_4", &gen_hfs_pz_no_rad_eta_lt_4, "gen_hfs_pz_no_rad_eta_lt_4/F" ) ;
      tt_out -> Branch( "gen_hfs_phi_no_rad_eta_lt_4", &gen_hfs_phi_no_rad_eta_lt_4, "gen_hfs_phi_no_rad_eta_lt_4/F" ) ;

      tt_out -> Branch( "gen_hfs_Sigma_with_rad_eta_lt_4", &gen_hfs_Sigma_with_rad_eta_lt_4, "gen_hfs_Sigma_with_rad_eta_lt_4/F" ) ;
      tt_out -> Branch( "gen_hfs_T_with_rad_eta_lt_4", &gen_hfs_T_with_rad_eta_lt_4, "gen_hfs_T_with_rad_eta_lt_4/F" ) ;
      tt_out -> Branch( "gen_hfs_E_with_rad_eta_lt_4", &gen_hfs_E_with_rad_eta_lt_4, "gen_hfs_E_with_rad_eta_lt_4/F" ) ;
      tt_out -> Branch( "gen_hfs_pz_with_rad_eta_lt_4", &gen_hfs_pz_with_rad_eta_lt_4, "gen_hfs_pz_with_rad_eta_lt_4/F" ) ;
      tt_out -> Branch( "gen_hfs_phi_with_rad_eta_lt_4", &gen_hfs_phi_with_rad_eta_lt_4, "gen_hfs_phi_with_rad_eta_lt_4/F" ) ;

      tt_out -> Branch( "fsr_gamma_dr", &fsr_gamma_dr, "fsr_gamma_dr/F" ) ;
      tt_out -> Branch( "fsr_gamma_e", &fsr_gamma_e, "fsr_gamma_e/F" ) ;
      tt_out -> Branch( "fsr_gamma_pt", &fsr_gamma_pt, "fsr_gamma_pt/F" ) ;
      tt_out -> Branch( "fsr_gamma_pz", &fsr_gamma_pz, "fsr_gamma_pz/F" ) ;
      tt_out -> Branch( "fsr_gamma_phi", &fsr_gamma_phi, "fsr_gamma_phi/F" ) ;
      tt_out -> Branch( "fsr_gamma_angle", &fsr_gamma_angle, "fsr_gamma_angle/F" ) ;
      tt_out -> Branch( "fsr_gamma_angle_wrt_beam", &fsr_gamma_angle_wrt_beam, "fsr_gamma_angle_wrt_beam/F" ) ;
      tt_out -> Branch( "fsr_gamma_eta", &fsr_gamma_eta, "fsr_gamma_eta/F" ) ;

      tt_out -> Branch( "isr_gamma_dr", &isr_gamma_dr, "isr_gamma_dr/F" ) ;
      tt_out -> Branch( "isr_gamma_e", &isr_gamma_e, "isr_gamma_e/F" ) ;
      tt_out -> Branch( "isr_gamma_pt", &isr_gamma_pt, "isr_gamma_pt/F" ) ;
      tt_out -> Branch( "isr_gamma_pz", &isr_gamma_pz, "isr_gamma_pz/F" ) ;
      tt_out -> Branch( "isr_gamma_phi", &isr_gamma_phi, "isr_gamma_phi/F" ) ;
      tt_out -> Branch( "isr_gamma_angle", &isr_gamma_angle, "isr_gamma_angle/F" ) ;
      tt_out -> Branch( "isr_gamma_angle_wrt_se", &isr_gamma_angle_wrt_se, "isr_gamma_angle_wrt_se/F" ) ;
      tt_out -> Branch( "isr_gamma_eta", &isr_gamma_eta, "isr_gamma_eta/F" ) ;

      tt_out -> Branch( "n_pho", &n_pho, "n_pho/I" ) ;
      tt_out -> Branch( "pho_isr_match", &pho_isr_match ) ;
      tt_out -> Branch( "pho_fsr_match", &pho_fsr_match ) ;
      tt_out -> Branch( "pho_angle_wrt_ebeam", &pho_angle_wrt_ebeam ) ;
      tt_out -> Branch( "pho_angle_wrt_efinal", &pho_angle_wrt_efinal ) ;
      tt_out -> Branch( "pho_e", &pho_e ) ;
      tt_out -> Branch( "pho_pt", &pho_pt ) ;
      tt_out -> Branch( "pho_pz", &pho_pz ) ;
      tt_out -> Branch( "pho_eta", &pho_eta ) ;
      tt_out -> Branch( "pho_dphi_wrt_efinal", &pho_dphi_wrt_efinal ) ;

   } // not save only dnn inputs?

   tt_out -> Branch( "sef_e", &sef_e, "sef_e/F" ) ;
   tt_out -> Branch( "sef_pt", &sef_pt, "sef_pt/F" ) ;
   tt_out -> Branch( "sef_pz", &sef_pz, "sef_pz/F" ) ;
   tt_out -> Branch( "sef_phi", &sef_phi, "sef_phi/F" ) ;
   tt_out -> Branch( "sef_eta", &sef_eta, "sef_eta/F" ) ;

   tt_out -> Branch( "gen_hfs_Sigma_no_rad", &gen_hfs_Sigma_no_rad, "gen_hfs_Sigma_no_rad/F" ) ;
   tt_out -> Branch( "gen_hfs_T_no_rad", &gen_hfs_T_no_rad, "gen_hfs_T_no_rad/F" ) ;
   tt_out -> Branch( "gen_hfs_E_no_rad", &gen_hfs_E_no_rad, "gen_hfs_E_no_rad/F" ) ;
   tt_out -> Branch( "gen_hfs_pz_no_rad", &gen_hfs_pz_no_rad, "gen_hfs_pz_no_rad/F" ) ;
   tt_out -> Branch( "gen_hfs_phi_no_rad", &gen_hfs_phi_no_rad, "gen_hfs_phi_no_rad/F" ) ;

   tt_out -> Branch( "tower_sum_07", &tower_sum_07, "tower_sum_07/F" ) ;
   tt_out -> Branch( "tower_sum_20", &tower_sum_20, "tower_sum_20/F" ) ;
   tt_out -> Branch( "tower_sum_30", &tower_sum_30, "tower_sum_30/F" ) ;
   tt_out -> Branch( "tower_sum_40", &tower_sum_40, "tower_sum_40/F" ) ;
   tt_out -> Branch( "n_towers_07", &n_towers_07, "n_towers_07/I" ) ;
   tt_out -> Branch( "n_towers_20", &n_towers_20, "n_towers_20/I" ) ;
   tt_out -> Branch( "n_towers_30", &n_towers_30, "n_towers_30/I" ) ;
   tt_out -> Branch( "n_towers_40", &n_towers_40, "n_towers_40/I" ) ;

   tt_out -> Branch( "phi_pho_closest_to_ebeam", &phi_pho_closest_to_ebeam, "phi_pho_closest_to_ebeam/F" ) ;
   tt_out -> Branch( "eta_pho_closest_to_ebeam", &eta_pho_closest_to_ebeam, "eta_pho_closest_to_ebeam/F" ) ;
   tt_out -> Branch( "e_pho_closest_to_ebeam", &e_pho_closest_to_ebeam, "e_pho_closest_to_ebeam/F" ) ;


   tt_out -> Branch( "beam_e_e", &beam_e_e, "beam_e_e/F" ) ;
   tt_out -> Branch( "beam_p_e", &beam_p_e, "beam_p_e/F" ) ;


   sprintf( bname, "gen_y[%d]/F", nmeth ) ;
   tt_out -> Branch( "gen_y", &gen_y, bname ) ;

   sprintf( bname, "gen_x[%d]/F", nmeth ) ;
   tt_out -> Branch( "gen_x", &gen_x, bname ) ;

   sprintf( bname, "gen_Q2[%d]/F", nmeth ) ;
   tt_out -> Branch( "gen_Q2", &gen_Q2, bname ) ;



   sprintf( bname, "obs_y[%d]/F", nmeth ) ;
   tt_out -> Branch( "obs_y", &obs_y, bname ) ;

   sprintf( bname, "obs_x[%d]/F", nmeth ) ;
   tt_out -> Branch( "obs_x", &obs_x, bname ) ;

   sprintf( bname, "obs_Q2[%d]/F", nmeth ) ;
   tt_out -> Branch( "obs_Q2", &obs_Q2, bname ) ;


   tt_out -> Branch( "has_isr", &has_isr, "has_isr/B" ) ;
   tt_out -> Branch( "has_fsr", &has_fsr, "has_fsr/B" ) ;
   tt_out -> Branch( "gen_fsr_photon_recombined", &gen_fsr_photon_recombined, "gen_fsr_photon_recombined/B" ) ;
   tt_out -> Branch( "obs_fsr_photon_recombined", &obs_fsr_photon_recombined, "obs_fsr_photon_recombined/B" ) ;

   tt_out -> Branch( "from_tlv_gen_Q2", &from_tlv_gen_Q2, "from_tlv_gen_Q2/F" ) ;
   tt_out -> Branch( "from_tlv_gen_x", &from_tlv_gen_x, "from_tlv_gen_x/F" ) ;
   tt_out -> Branch( "from_tlv_gen_y", &from_tlv_gen_y, "from_tlv_gen_y/F" ) ;

   tt_out -> Branch( "gen_obs_e_e", &gen_obs_e_e, "gen_obs_e_e/F" ) ;
   tt_out -> Branch( "gen_obs_e_pz", &gen_obs_e_pz, "gen_obs_e_pz/F" ) ;
   tt_out -> Branch( "gen_obs_e_pt", &gen_obs_e_pt, "gen_obs_e_pt/F" ) ;
   tt_out -> Branch( "gen_obs_e_phi", &gen_obs_e_phi, "gen_obs_e_phi/F" ) ;
   tt_out -> Branch( "gen_obs_e_eta", &gen_obs_e_eta, "gen_obs_e_eta/F" ) ;

   tt_out -> Branch( "gen_obs_hfs_e", &gen_obs_hfs_e, "gen_obs_hfs_e/F" ) ;
   tt_out -> Branch( "gen_obs_hfs_pz", &gen_obs_hfs_pz, "gen_obs_hfs_pz/F" ) ;
   tt_out -> Branch( "gen_obs_hfs_pt", &gen_obs_hfs_pt, "gen_obs_hfs_pt/F" ) ;
   tt_out -> Branch( "gen_obs_hfs_phi", &gen_obs_hfs_phi, "gen_obs_hfs_phi/F" ) ;

   tt_out -> Branch( "gen_obs_dphi_hfs_eta_lt_4", &gen_obs_dphi_hfs_eta_lt_4, "gen_obs_dphi_hfs_eta_lt_4/F" ) ;
   tt_out -> Branch( "obs_dphi", &obs_dphi, "obs_dphi/F" ) ;

   tt_out -> Branch( "obs_e_e", &obs_e_e, "obs_e_e/F" ) ;
   tt_out -> Branch( "obs_e_pz", &obs_e_pz, "obs_e_pz/F" ) ;
   tt_out -> Branch( "obs_e_pt", &obs_e_pt, "obs_e_pt/F" ) ;
   tt_out -> Branch( "obs_e_phi", &obs_e_phi, "obs_e_phi/F" ) ;
   tt_out -> Branch( "obs_e_eta", &obs_e_eta, "obs_e_eta/F" ) ;

   tt_out -> Branch( "obs_hfs_e", &obs_hfs_e, "obs_hfs_e/F" ) ;
   tt_out -> Branch( "obs_hfs_pz", &obs_hfs_pz, "obs_hfs_pz/F" ) ;
   tt_out -> Branch( "obs_hfs_pt", &obs_hfs_pt, "obs_hfs_pt/F" ) ;
   tt_out -> Branch( "obs_hfs_phi", &obs_hfs_phi, "obs_hfs_phi/F" ) ;

   tt_out -> Branch( "e_is_correct_track", &e_is_correct_track, "e_is_correct_track/B" ) ;




   TDatabasePDG* pdg = new TDatabasePDG() ;
   pdg->ReadPDGTable() ;

   Long64_t nentries = fChain->GetEntries();

   int last_evt = nentries ;
   if ( maxEvts >= 0 ) last_evt = maxEvts ;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<last_evt;jentry++) {

      int ei = jentry ;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;


      if ( !verbose && ei%100 == 0 ) {
         printf(" --- Event: %7d / %lld    %6.3f\r", ei, nentries, (1.*ei)/(1.*nentries) ) ;
         fflush(stdout) ;
      }



    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //-- Gen-level-based analysis here.

      if ( verbose ) {
         printf("\n\n\n   %4d :  n particles = %d\n", ei, Particle_ ) ;
         for ( int pi = 0; pi < Particle_; pi ++ ) {
            TParticlePDG* pdg_particle = pdg -> GetParticle( Particle_PID[pi] ) ;
            char pdg_name[100];
            float pdg_charge(0.) ;
            if ( pdg_particle != 0x0 ) {
               sprintf( pdg_name, "%s", pdg_particle->GetName() ) ;
               pdg_charge = pdg_particle->Charge() / (3.) ;
            }
            printf("   %4d : %3d :  PID %10d  %15s  status = %3d  M1 %3d M2 %3d  D1 %3d D2 %3d   Pt %15.6f  Pz %15.6f  Eta %9.5f  Phi %9.5f ", ei, pi,
            Particle_PID[pi], pdg_name, Particle_Status[pi],
            Particle_M1[pi],
            Particle_M2[pi],
            Particle_D1[pi],
            Particle_D2[pi],
            escale * Particle_PT[pi],
            escale * Particle_Pz[pi],
            Particle_Eta[pi],
            Particle_Phi[pi]
            ) ;
            if ( Particle_PID[pi] == 22 ) {
               int mom_pid = 0 ;
               int mom_ind = Particle_M1[pi] ;
               if ( mom_ind >= 0 ) {
                  mom_pid = Particle_PID[mom_ind] ;
               }
               if ( mom_pid != 111 ) {
                  TParticlePDG* mom_pdg_particle = pdg -> GetParticle( mom_pid ) ;
                  char mom_pdg_name[100];
                  if ( mom_pdg_particle != 0x0 ) {
                     sprintf( mom_pdg_name, "%s", mom_pdg_particle->GetName() ) ;
                  }
                  printf(" **** %s", mom_pdg_name ) ;
               }
            }
            if ( abs(Particle_PID[pi]) == 11 ) {
               bool conversion = false ;
               bool decay = false ;
               int mom_pid = 0 ;
               int mom_ind = Particle_M1[pi] ;
               if ( mom_ind >= 0 ) {
                  mom_pid = Particle_PID[mom_ind] ;
               }
               if ( mom_pid == 22 ) conversion = true ;
               if ( ( mom_pid != 0 && abs(mom_pid) != 11) && ! conversion ) decay = true ;
               if ( ! conversion && ! decay ) {
                  if ( Particle_Status[pi] == 1 ) {
                     printf(" -- final scattered electron") ;
                  } else if ( Particle_Status[pi] == 4 ) {
                     printf(" -- beam electron") ;
                  } else {
                     printf(" -- electron") ;
                  }
               }
            }
            printf("\n") ;
         } // pi
         //continue ;
      } // end verbose.



      int be0_pi = -1 ;
      int bef_pi = -1 ;

      int se0_pi = -1 ;
      int sef_pi = -1 ;
      int rad_gam_pi = -1 ;

      int post_isr_beam_e_pi = -1 ;
      int pre_fsr_scattered_e_pi = -1 ;
      int beam_p_pi = -1 ;

      int fsr_gamma_pi = -1 ;
      int isr_gamma_pi = -1 ;

      has_isr = false ;
      has_fsr = false ;
      gen_fsr_photon_recombined = false ;

      for ( int pi = 0; pi < Particle_; pi ++ ) {

         if ( Particle_PID[pi] == 2212 && Particle_Status[pi] == 4 ) {
            beam_p_e = escale * Particle_E[pi] ;
            beam_p_pi = pi ;
            continue ;
         }

         if ( Particle_PID[pi] == 23 ) {
            if ( Particle_M1[pi] >= 0 ) {
               if ( abs(Particle_PID[ Particle_M1[pi] ]) == 11 ) {
                  post_isr_beam_e_pi = Particle_M1[pi] ;
               }
            }
         }

         if ( Particle_PID[pi] == 22 ) {
            int mom_pi = Particle_M1[pi] ;
            int mom_pid = 0 ;
            int mom_status = 0 ;
            if ( mom_pi >= 0 ) {
               mom_pid = Particle_PID[mom_pi] ;
               mom_status = Particle_Status[mom_pi] ;
            }
            if ( abs(mom_pid) == 11 ) {
               rad_gam_pi = pi ;
               if ( mom_status == 4 ) {
                  has_isr = true ;
                  isr_gamma_pi = pi ;
               } else {
                  has_fsr = true ;
                  pre_fsr_scattered_e_pi = mom_pi ;
                  fsr_gamma_pi = pi ;
               }
            }
         }

         if ( abs(Particle_PID[pi]) != 11 )  continue ;

         if ( Particle_Status[pi] == 4 ) {
            be0_pi = pi ;
            beam_e_e = escale * Particle_E[pi] ;
         }

         if ( Particle_Status[pi] == 11 ) se0_pi = pi ;
         if ( Particle_Status[pi] ==  1 ) {
            if ( Particle_M1[pi] > 0 ) {
               int mom_pi = Particle_M1[pi] ;
               if ( Particle_Status[mom_pi] < 4 ) continue ;
            }
            sef_pi = pi ;
         }

      } // pi

      if ( !has_fsr && pre_fsr_scattered_e_pi < 0 ) {
         pre_fsr_scattered_e_pi = sef_pi ;
      }





      gen_hfs_Sigma_no_rad = 0. ;
      gen_hfs_E_no_rad = 0. ;
      gen_hfs_pz_no_rad = 0. ;
      float gen_hfs_px_no_rad = 0. ;
      float gen_hfs_py_no_rad = 0. ;

      gen_hfs_Sigma_with_rad = 0. ;
      gen_hfs_E_with_rad = 0. ;
      gen_hfs_pz_with_rad = 0. ;
      float gen_hfs_px_with_rad = 0. ;
      float gen_hfs_py_with_rad = 0. ;

      gen_hfs_Sigma_no_rad_eta_lt_4 = 0. ;
      gen_hfs_E_no_rad_eta_lt_4 = 0. ;
      gen_hfs_pz_no_rad_eta_lt_4 = 0. ;
      float gen_hfs_px_no_rad_eta_lt_4 = 0. ;
      float gen_hfs_py_no_rad_eta_lt_4 = 0. ;

      gen_hfs_Sigma_with_rad_eta_lt_4 = 0. ;
      gen_hfs_E_with_rad_eta_lt_4 = 0. ;
      gen_hfs_pz_with_rad_eta_lt_4 = 0. ;
      float gen_hfs_px_with_rad_eta_lt_4 = 0. ;
      float gen_hfs_py_with_rad_eta_lt_4 = 0. ;

      for ( int pi = 0; pi < Particle_; pi ++ ) {

         if ( pi == sef_pi ) continue ;

         if ( Particle_Status[pi] != 1 ) continue ;

         gen_hfs_Sigma_with_rad += escale * ( Particle_E[pi] - Particle_Pz[pi] ) ;
         gen_hfs_px_with_rad += escale * Particle_Px[pi] ;
         gen_hfs_py_with_rad += escale * Particle_Py[pi] ;
         gen_hfs_pz_with_rad += escale * Particle_Pz[pi] ;
         gen_hfs_E_with_rad  += escale * Particle_E[pi] ;

         if ( fabs( Particle_Eta[pi] ) < 4.0 ) {
            gen_hfs_Sigma_with_rad_eta_lt_4 += escale * ( Particle_E[pi] - Particle_Pz[pi] ) ;
            gen_hfs_px_with_rad_eta_lt_4 += escale * Particle_Px[pi] ;
            gen_hfs_py_with_rad_eta_lt_4 += escale * Particle_Py[pi] ;
            gen_hfs_pz_with_rad_eta_lt_4 += escale * Particle_Pz[pi] ;
            gen_hfs_E_with_rad_eta_lt_4  += escale * Particle_E[pi] ;
         }

         if ( Particle_PID[pi] == 22 ) {
            int mom_pi = Particle_M1[pi] ;
            int mom_pid = 0 ;
            int mom_status = 0 ;
            if ( mom_pi >= 0 ) {
               mom_pid = Particle_PID[mom_pi] ;
               mom_status = Particle_Status[mom_pi] ;
            }
            if ( abs(mom_pid) == 11 ) continue ;
         }

         gen_hfs_Sigma_no_rad += escale * ( Particle_E[pi] - Particle_Pz[pi] ) ;
         gen_hfs_px_no_rad += escale * Particle_Px[pi] ;
         gen_hfs_py_no_rad += escale * Particle_Py[pi] ;
         gen_hfs_pz_no_rad += escale * Particle_Pz[pi] ;
         gen_hfs_E_no_rad  += escale * Particle_E[pi] ;

         if ( fabs( Particle_Eta[pi] ) < 4.0 ) {
            gen_hfs_Sigma_no_rad_eta_lt_4 += escale * ( Particle_E[pi] - Particle_Pz[pi] ) ;
            gen_hfs_px_no_rad_eta_lt_4 += escale * Particle_Px[pi] ;
            gen_hfs_py_no_rad_eta_lt_4 += escale * Particle_Py[pi] ;
            gen_hfs_pz_no_rad_eta_lt_4 += escale * Particle_Pz[pi] ;
            gen_hfs_E_no_rad_eta_lt_4  += escale * Particle_E[pi] ;
         }

         if ( verbose ) {
            TParticlePDG* pdg_particle = pdg -> GetParticle( Particle_PID[pi] ) ;
            char pdg_name[100];
            if ( pdg_particle != 0x0 ) {
               sprintf( pdg_name, "%s", pdg_particle->GetName() ) ;
            }
            char eta_tag[10] ;
            sprintf( eta_tag, " " ) ;
            if ( fabs( Particle_Eta[pi] ) > 4. ) sprintf( eta_tag, "*" ) ;
            printf("  HFS : %3d :  %15s   E = %9.4f , pz = %9.3f , px = %9.3f , py = %9.3f , eta = %9.5f %s  E-pz %9.5f\n", pi, pdg_name,
              escale * Particle_E[pi], escale * Particle_Pz[pi], escale * Particle_Px[pi], escale * Particle_Py[pi], Particle_Eta[pi], eta_tag, 
              escale * Particle_E[pi] - escale * Particle_Pz[pi] ) ;
         }

      } // pi

      gen_hfs_T_no_rad = sqrt( gen_hfs_px_no_rad * gen_hfs_px_no_rad + gen_hfs_py_no_rad * gen_hfs_py_no_rad ) ;
      gen_hfs_T_with_rad = sqrt( gen_hfs_px_with_rad * gen_hfs_px_with_rad + gen_hfs_py_with_rad * gen_hfs_py_with_rad ) ;
      gen_hfs_T_with_rad_eta_lt_4 = sqrt( gen_hfs_px_with_rad_eta_lt_4 * gen_hfs_px_with_rad_eta_lt_4 + gen_hfs_py_with_rad_eta_lt_4 * gen_hfs_py_with_rad_eta_lt_4 ) ;
      gen_hfs_T_no_rad_eta_lt_4 = sqrt( gen_hfs_px_no_rad_eta_lt_4 * gen_hfs_px_no_rad_eta_lt_4 + gen_hfs_py_no_rad_eta_lt_4 * gen_hfs_py_no_rad_eta_lt_4 ) ;

      gen_hfs_phi_no_rad = atan2( gen_hfs_py_no_rad, gen_hfs_pz_no_rad ) ;
      gen_hfs_phi_with_rad = atan2( gen_hfs_py_with_rad, gen_hfs_pz_with_rad ) ;
      gen_hfs_phi_with_rad_eta_lt_4 = atan2( gen_hfs_py_with_rad_eta_lt_4, gen_hfs_pz_with_rad_eta_lt_4 ) ;
      gen_hfs_phi_no_rad_eta_lt_4 = atan2( gen_hfs_py_no_rad_eta_lt_4, gen_hfs_pz_no_rad_eta_lt_4 ) ;




      if ( se0_pi < 0 ) {
         //printf("\n\n *** can't find initial scattered electron.\n\n") ;
         //gSystem -> Exit(-1) ;
         se0_pi = sef_pi ;
      }
      if ( sef_pi < 0 ) {
         printf("\n\n *** can't find final scattered electron.\n\n") ;
         ///gSystem -> Exit(-1) ;
         continue ;
      }

      if ( has_isr ) {
         bef_pi = se0_pi ;
         se0_pi = sef_pi ;
      }

      if ( bef_pi < 0 ) {
         bef_pi = be0_pi ;
      }




      if ( has_fsr && fsr_gamma_pi >= 0 ) {
         float dphi = Particle_Phi[fsr_gamma_pi] - Particle_Phi[sef_pi] ;
         if ( dphi < - 3.14159265 ) dphi += 2 * 3.14159265 ;
         if ( dphi >   3.14159265 ) dphi -= 2 * 3.14159265 ;
         float deta = Particle_Eta[fsr_gamma_pi] - Particle_Eta[sef_pi] ;
         fsr_gamma_dr = sqrt( dphi*dphi + deta*deta ) ;
         fsr_gamma_e  = escale * Particle_E[fsr_gamma_pi] ;
         fsr_gamma_pt = escale * Particle_PT[fsr_gamma_pi] ;
         fsr_gamma_pz = escale * Particle_Pz[fsr_gamma_pi] ;
         fsr_gamma_phi = Particle_Phi[fsr_gamma_pi] ;
         fsr_gamma_eta = Particle_Eta[fsr_gamma_pi] ;
         float pgamma = sqrt( Particle_Px[fsr_gamma_pi] * Particle_Px[fsr_gamma_pi] + Particle_Py[fsr_gamma_pi] * Particle_Py[fsr_gamma_pi] + Particle_Pz[fsr_gamma_pi] * Particle_Pz[fsr_gamma_pi] ) ;
         float pe     = sqrt( Particle_Px[sef_pi] * Particle_Px[sef_pi] + Particle_Py[sef_pi] * Particle_Py[sef_pi] + Particle_Pz[sef_pi] * Particle_Pz[sef_pi] ) ;
         float cos_angle = ( Particle_Px[fsr_gamma_pi] * Particle_Px[sef_pi]
                           + Particle_Py[fsr_gamma_pi] * Particle_Py[sef_pi]
                           + Particle_Pz[fsr_gamma_pi] * Particle_Pz[sef_pi] ) /
                           ( pgamma * pe ) ;
         fsr_gamma_angle = acos( cos_angle ) ;
               cos_angle = ( Particle_Px[fsr_gamma_pi] * Particle_Px[be0_pi]
                           + Particle_Py[fsr_gamma_pi] * Particle_Py[be0_pi]
                           + Particle_Pz[fsr_gamma_pi] * Particle_Pz[be0_pi] ) /
                           ( pgamma * fabs(Particle_Pz[be0_pi]) ) ;
         fsr_gamma_angle_wrt_beam = acos( cos_angle ) ;
         //if ( verbose ) { printf(" *** debug FSR : cos_angle = %9.5f , angle = %9.5f , P1 %9.5f, P2 %9.5f\n", cos_angle, fsr_gamma_angle, pgamma, pe ) ; }
      } else {
         fsr_gamma_dr = -1. ;
         fsr_gamma_e = 0. ;
         fsr_gamma_pt = 0. ;
         fsr_gamma_pz = 0. ;
         fsr_gamma_phi = -9. ;
         fsr_gamma_angle = -1. ;
      }


      if ( has_isr && isr_gamma_pi >= 0 ) {
         float dphi = Particle_Phi[isr_gamma_pi] - Particle_Phi[bef_pi] ;
         if ( dphi < - 3.14159265 ) dphi += 2 * 3.14159265 ;
         if ( dphi >   3.14159265 ) dphi -= 2 * 3.14159265 ;
         float deta = Particle_Eta[isr_gamma_pi] - Particle_Eta[bef_pi] ;
         //if (verbose) printf("  *** debug ISR :  phi gamma = %8.4f, phi bef = %8.4f,  dphi = %8.4f   eta gamma = %8.4f, eta bef = %8.4f,  deta = %8.4f\n",
          // Particle_Phi[isr_gamma_pi], Particle_Phi[bef_pi], dphi,   Particle_Eta[isr_gamma_pi], Particle_Eta[bef_pi], deta ) ;
         isr_gamma_dr = sqrt( dphi*dphi + deta*deta ) ;
         isr_gamma_e  = escale * Particle_E[isr_gamma_pi] ;
         isr_gamma_pt = escale * Particle_PT[isr_gamma_pi] ;
         isr_gamma_pz = escale * Particle_Pz[isr_gamma_pi] ;
         isr_gamma_phi = Particle_Phi[isr_gamma_pi] ;
         isr_gamma_eta = Particle_Eta[isr_gamma_pi] ;
         float pgamma = sqrt( Particle_Px[isr_gamma_pi] * Particle_Px[isr_gamma_pi] + Particle_Py[isr_gamma_pi] * Particle_Py[isr_gamma_pi] + Particle_Pz[isr_gamma_pi] * Particle_Pz[isr_gamma_pi] ) ;
         float pe     = sqrt( Particle_Px[bef_pi] * Particle_Px[bef_pi] + Particle_Py[bef_pi] * Particle_Py[bef_pi] + Particle_Pz[bef_pi] * Particle_Pz[bef_pi] ) ;
         float cos_angle = ( Particle_Px[isr_gamma_pi] * Particle_Px[bef_pi]
                           + Particle_Py[isr_gamma_pi] * Particle_Py[bef_pi]
                           + Particle_Pz[isr_gamma_pi] * Particle_Pz[bef_pi] ) /
                           ( pgamma * pe ) ;
         isr_gamma_angle = acos( cos_angle ) ;
         float pse     = sqrt( Particle_Px[sef_pi] * Particle_Px[sef_pi] + Particle_Py[sef_pi] * Particle_Py[sef_pi] + Particle_Pz[sef_pi] * Particle_Pz[sef_pi] ) ;
               cos_angle = ( Particle_Px[isr_gamma_pi] * Particle_Px[sef_pi]
                           + Particle_Py[isr_gamma_pi] * Particle_Py[sef_pi]
                           + Particle_Pz[isr_gamma_pi] * Particle_Pz[sef_pi] ) /
                           ( pgamma * pse ) ;
         isr_gamma_angle_wrt_se = acos( cos_angle ) ;
         //if ( verbose ) { printf(" *** debug ISR : cos_angle = %9.5f , angle = %9.5f , P1 %9.5f, P2 %9.5f\n", cos_angle, isr_gamma_angle, Particle_P[isr_gamma_pi], Particle_P[bef_pi] ) ; }
      } else {
         isr_gamma_dr = -1. ;
         isr_gamma_e = 0. ;
         isr_gamma_pt = 0. ;
         isr_gamma_pz = 0. ;
         isr_gamma_phi = -9. ;
         isr_gamma_eta = -9. ;
         isr_gamma_angle = -1. ;
      }









      be0_e   = escale * Particle_E[be0_pi] ;
      be0_pt  = escale * Particle_PT[be0_pi] ;
      be0_pz  = escale * Particle_Pz[be0_pi] ;
      be0_phi = Particle_Phi[be0_pi] ;

      bef_e   = escale * Particle_E[bef_pi] ;
      bef_pt  = escale * Particle_PT[bef_pi] ;
      bef_pz  = escale * Particle_Pz[bef_pi] ;
      bef_phi = Particle_Phi[bef_pi] ;

      se0_e   = escale * Particle_E[se0_pi] ;
      se0_pt  = escale * Particle_PT[se0_pi] ;
      se0_pz  = escale * Particle_Pz[se0_pi] ;
      se0_phi = Particle_Phi[se0_pi] ;

      sef_e   = escale * Particle_E[sef_pi] ;
      sef_pt  = escale * Particle_PT[sef_pi] ;
      sef_pz  = escale * Particle_Pz[sef_pi] ;
      sef_phi = Particle_Phi[sef_pi] ;
      sef_eta = Particle_Eta[sef_pi] ;

      gen_e0_Q2 = 2. * beam_e_e * ( se0_e + se0_pz ) ;
      gen_e0_y  = 1. - ( se0_e - se0_pz ) / ( 2. * beam_e_e ) ;
      float gen_e0_s = 4. * beam_e_e * beam_p_e ;
      gen_e0_x = gen_e0_Q2 / ( gen_e0_y * gen_e0_s ) ;

      gen_ef_Q2 = 2. * beam_e_e * ( sef_e + sef_pz ) ;
      gen_ef_y  = 1. - ( sef_e - sef_pz ) / ( 2. * beam_e_e ) ;
      float gen_ef_s = 4. * beam_e_e * beam_p_e ;
      gen_ef_x = gen_ef_Q2 / ( gen_ef_y * gen_ef_s ) ;

      if ( verbose ) {
         if ( rad_gam_pi >= 0 ) {
            char ptype[100] ;
            if ( has_isr ) sprintf( ptype, "isr" ) ;
            if ( has_fsr ) sprintf( ptype, "fsr" ) ;
            printf(" %4d : %3d  %s radiated photon  pt = %9.3f\n", ei, rad_gam_pi, ptype, escale * Particle_PT[rad_gam_pi] ) ;
            if ( has_fsr ) {
               printf(" %4d : FSR photon DeltaR = %8.5f,  E = %9.5f,  pt = %9.5f,  pz = %9.5f,  phi = %8.4f\n",
                  ei, fsr_gamma_dr, fsr_gamma_e, fsr_gamma_pt, fsr_gamma_pz, fsr_gamma_phi ) ;
            }
            if ( has_isr ) {
               printf(" %4d : ISR photon DeltaR = %8.5f,  E = %9.5f,  pt = %9.5f,  pz = %9.5f,  phi = %8.4f\n",
                  ei, isr_gamma_dr, isr_gamma_e, isr_gamma_pt, isr_gamma_pz, isr_gamma_phi ) ;
            }
         }
         printf(" %4d : %3d be0_pz = %9.3f , %3d bef_pz = %9.3f , diff = %9.3f\n", ei, be0_pi, be0_pz, bef_pi, bef_pz, be0_pz-bef_pz ) ;
         printf(" %4d : %3d se0_pt = %9.3f , %3d sef_pt = %9.3f , diff = %9.3f\n", ei, se0_pi, se0_pt, sef_pi, sef_pt, se0_pt-sef_pt ) ;
         printf(" %4d :     HFS pt = %9.3f , HFS Sigma = %9.3f\n", ei, gen_hfs_T_no_rad, gen_hfs_Sigma_no_rad ) ;
         printf(" %4d :  beam_p_pi = %3d   post_isr_beam_e_pi = %3d   pre_fsr_scattered_e_pi = %3d\n", ei, beam_p_pi, post_isr_beam_e_pi, pre_fsr_scattered_e_pi ) ;
      }



     //--- Calculate everything using post ISR / pre FSR 4-vectors

      TLorentzVector tlv_p( Particle_Px[beam_p_pi], Particle_Py[beam_p_pi], Particle_Pz[beam_p_pi], Particle_E[beam_p_pi] ) ;
      TLorentzVector tlv_l( Particle_Px[post_isr_beam_e_pi], Particle_Py[post_isr_beam_e_pi], Particle_Pz[post_isr_beam_e_pi], Particle_E[post_isr_beam_e_pi] ) ;
      TLorentzVector tlv_lprime( Particle_Px[pre_fsr_scattered_e_pi], Particle_Py[pre_fsr_scattered_e_pi], Particle_Pz[pre_fsr_scattered_e_pi], Particle_E[pre_fsr_scattered_e_pi] ) ;

      tlv_p = escale * tlv_p ;
      tlv_l = escale * tlv_l ;
      tlv_lprime = escale * tlv_lprime ;

      TLorentzVector tlv_q = tlv_l - tlv_lprime ;

      from_tlv_gen_Q2 = -1 * tlv_q * tlv_q ;
      from_tlv_gen_y = ( tlv_p * tlv_q ) / ( tlv_p * tlv_l ) ;
      from_tlv_gen_x = -1 * ( tlv_q * tlv_q ) / ( 2. * tlv_p * tlv_q ) ;



      {

         int im ;

         float e ;
         float e0 ;
         float theta ;
         float Sigma ;
         float T ;
         float tan_gamma_over_2 ;
         float ep ;

         e = sef_e ;
         e0 = be0_e ;
         theta = atan2( sef_pt, sef_pz ) ;
         Sigma = gen_hfs_Sigma_no_rad ;
         T = gen_hfs_T_no_rad ;
         tan_gamma_over_2 = Sigma / T ;
         ep = beam_p_e ;

         if ( has_fsr && recombine_fsr_gamma ) {
            if ( fsr_gamma_dr * 180. / 3.14159265 < 4.0 ) {
               float gam_px = fsr_gamma_pt * cos( fsr_gamma_phi ) ;
               float gam_py = fsr_gamma_pt * sin( fsr_gamma_phi ) ;
               float sef_px = sef_pt * cos( sef_phi ) ;
               float sef_py = sef_pt * sin( sef_phi ) ;
               float recomb_px = sef_px + gam_px ;
               float recomb_py = sef_py + gam_py ;
               float recomb_pz = sef_pz + fsr_gamma_pz ;
               float recomb_pt = sqrt( recomb_px * recomb_px + recomb_py * recomb_py ) ;
               if (verbose) printf("   FSR, pre recomb  E = %10.5f , theta = %10.6f   ", sef_e, theta ) ;
               e = sef_e + fsr_gamma_e ;
               theta = atan2( recomb_pt, recomb_pz ) ;
               if (verbose) printf("   post recomb  E = %10.5f , theta = %10.6f   ", sef_e, theta ) ;
               if (verbose) printf("   se0  E = %10.5f , theta = %10.6f\n", se0_e, atan2( se0_pt, se0_pz ) ) ;
               gen_fsr_photon_recombined = true ;
            }
         }

         if ( include_rad_gamma_in_hfs ) {
            if ( (has_isr && fabs(isr_gamma_eta)<4 ) || ( has_fsr && fsr_gamma_dr * 180. / 3.14159265 >= 4.0 && fabs(fsr_gamma_eta)<4 ) ) {
               Sigma = gen_hfs_Sigma_with_rad ;
               T = gen_hfs_T_with_rad ;
               tan_gamma_over_2 = Sigma / T ;
            }
         }


        //-- Method  0 : electron
         im = 0 ;
         sprintf( meth_name[im], "0 : electron (E0,E,theta)" ) ;

         gen_y[im] = 1. - (e/e0) * pow( sin(theta/2.), 2 ) ;

         gen_Q2[im] = 4. * e0 * e * pow( cos(theta/2.), 2 ) ;

         gen_x[im] = ( e0 * e * pow( cos(theta/2.), 2 ) ) / ( ep * ( e0 - e * pow( sin(theta/2), 2 ) ) ) ;




        //-- Method  1 : unpub1
         im = 1 ;
         sprintf( meth_name[im], "1 : unpub1 (E0,E,Sigma)" ) ;

         gen_y[im] = Sigma / ( 2. * e0 ) ; // yh

         gen_Q2[im] = 4. * e0 * e - 4. * e0 * e0 + 2 * e0 * Sigma ;

         gen_x[im] = gen_Q2[im] / ( 2. * ep * Sigma ) ;





        //-- Method  2 : unpub2
         im = 2 ;
         sprintf( meth_name[im], "2 : unpub2 (E0,theta,Sigma)" ) ;

         gen_y[im] = Sigma / ( 2. * e0 ) ; // yh

         gen_Q2[im] = 2. * e0 * ( 2. * e0 - Sigma ) / pow( tan(theta/2.), 2. ) ;

         gen_x[im] = gen_Q2[im] / ( 2. * ep * Sigma ) ;




        //-- Method  3 : DA
         im = 3 ;
         sprintf( meth_name[im], "3 : DA (E0,theta,gamma)" ) ;

         gen_y[im] = tan_gamma_over_2 / ( tan_gamma_over_2 + tan(theta/2.) ) ;

         gen_Q2[im] = 4. * e0 * e0 * (1./tan(theta/2.)) / ( tan_gamma_over_2 + tan(theta/2.) ) ;

         gen_x[im] = gen_Q2[im] / ( 4. * ep * e0 * gen_y[im] ) ;




        //-- Method  4 : h
         im = 4 ;
         sprintf( meth_name[im], "4 : hadron (E0,Sigma,T)" ) ;

         gen_y[im] = Sigma / ( 2. * e0 ) ;

         gen_Q2[im] = T * T / ( 1. - gen_y[im] ) ;

         gen_x[im] = gen_Q2[im] / ( 2. * ep * Sigma ) ;




        //-- Method  5 : ISigma
         im = 5 ;
         sprintf( meth_name[im], "5 : ISigma (E,theta,Sigma)" ) ;


         gen_y[im] = Sigma / ( Sigma + e * ( 1. - cos(theta) ) ) ;

         gen_Q2[im] = e * e * pow( sin(theta), 2 ) / ( 1. - gen_y[im] ) ;

         gen_x[im] = e * pow( cos(theta/2.), 2 ) / ( ep * gen_y[im] ) ;




        //-- Method  6 : IDA
         im = 6 ;
         sprintf( meth_name[im], "6 : IDA (E,theta,gamma)" ) ;

         gen_y[im] = tan_gamma_over_2 / ( tan_gamma_over_2 + tan(theta/2.) ) ; // yDA

         //gen_Q2[im] = e * e * tan(theta/2.) * ( tan_gamma_over_2 + tan(theta/2.) ) / ( 1./tan(theta/2.) + tan(theta/2.) ) ; // *** this is wrong somehow
         gen_Q2[im] = 4. * e * e * (1./tan(theta/2.)) * ( tan_gamma_over_2 + tan(theta/2.) ) / pow( (1./tan(theta/2.) + tan(theta/2.)), 2 ) ;

         gen_x[im] = e * pow( cos(theta/2.), 2. ) / ( ep * gen_y[im] ) ;





        //-- Method  7 : unpub3
         im = 7 ;
         sprintf( meth_name[im], "7 : unpub3 (theta,Sigma,gamma)" ) ;

         gen_y[im] = tan_gamma_over_2 / ( tan_gamma_over_2 + tan(theta/2.) ) ; // yDA

         gen_Q2[im] = pow( Sigma / tan_gamma_over_2, 2 ) / ( 1. - gen_y[im] ) ;

         gen_x[im] = e * pow( cos(theta/2.), 2. ) / ( ep * gen_y[im] ) ;





        //-- Method  8 : eSigma
         im = 8 ;
         sprintf( meth_name[im], "8 : eSigma (E0,E,Sigma,theta)" ) ;

         gen_y[im] = 2. * e0 * Sigma / pow( ( Sigma + e*(1.-cos(theta)) ) , 2 ) ;

         gen_Q2[im] = 4. * e0 * e * pow( cos(theta/2.), 2 ) ; // Q2e

         gen_x[im] = ( e * pow( cos(theta/2), 2 ) + (e * e/(2.*Sigma)) * pow( sin(theta), 2 ) ) / ep ; // *** fixed typo in note



      }

      if ( verbose ) {

         printf("                                   " ) ;
         for ( int i=0; i<nmeth; i++ ) { printf("   %6d   ", i) ; }
         printf("\n") ;


         printf("  %4d :  gen_y   : ", ei ) ;
         printf("  %9.5f  | ", from_tlv_gen_y ) ;
         for ( int i=0; i<nmeth; i++ ) {
            printf("  %9.5f ", gen_y[i] ) ;
         } // i
         printf("\n") ;

         printf("  %4d :  gen_x   : ", ei ) ;
         printf("  %9.5f  | ", from_tlv_gen_x ) ;
         for ( int i=0; i<nmeth; i++ ) {
            printf("  %9.5f ", gen_x[i] ) ;
         } // i
         printf("\n") ;

         printf("  %4d :  gen_Q2  : ", ei ) ;
         printf("  %9.1f  | ", from_tlv_gen_Q2 ) ;
         for ( int i=0; i<nmeth; i++ ) {
            printf("  %9.1f ", gen_Q2[i] ) ;
         } // i
         printf("\n") ;


      }





     //-- default values for no ISR or FSR radiation.

      gen_obs_e_e = sef_e ;
      gen_obs_e_pz = sef_pz ;
      gen_obs_e_pt = sef_pt ;
      gen_obs_e_phi = sef_phi ;
      gen_obs_e_eta = sef_eta ;

      gen_obs_hfs_e = gen_hfs_E_no_rad ;
      gen_obs_hfs_pz = gen_hfs_pz_no_rad ;
      gen_obs_hfs_pt = gen_hfs_T_no_rad ;
      gen_obs_hfs_phi = gen_hfs_phi_no_rad_eta_lt_4 ;

      if ( has_fsr && recombine_fsr_gamma ) {
         if ( fsr_gamma_dr * 180. / 3.14159265 < 4.0 ) {
           //-- recombine FSR photon with scattered electron.
            float gam_px = fsr_gamma_pt * cos( fsr_gamma_phi ) ;
            float gam_py = fsr_gamma_pt * sin( fsr_gamma_phi ) ;
            float sef_px = sef_pt * cos( sef_phi ) ;
            float sef_py = sef_pt * sin( sef_phi ) ;
            float recomb_px = sef_px + gam_px ;
            float recomb_py = sef_py + gam_py ;
            float recomb_pz = sef_pz + fsr_gamma_pz ;
            float recomb_pt = sqrt( recomb_px * recomb_px + recomb_py * recomb_py ) ;
            gen_obs_e_e = sef_e + fsr_gamma_e ;
            gen_obs_e_pz = recomb_pz ;
            gen_obs_e_pt = recomb_pt ;
            gen_obs_e_phi = atan2( recomb_py, recomb_px ) ;
            float theta = atan2( recomb_pt, recomb_pz ) ;
            gen_obs_e_eta = -1 * log( tan( theta/2. ) ) ;
         }
      }

      gen_obs_dphi_hfs_eta_lt_4 = gen_obs_e_phi - gen_obs_hfs_phi ;

      if ( include_rad_gamma_in_hfs ) {
         if ( (has_isr && fabs(isr_gamma_eta)<4 ) || ( has_fsr && fsr_gamma_dr * 180. / 3.14159265 >= 4.0 && fabs(fsr_gamma_eta)<4 ) ) {
            //-- include radiated photon in HFS.
            gen_obs_hfs_e = gen_hfs_E_with_rad ;
            gen_obs_hfs_pz = gen_hfs_pz_with_rad ;
            gen_obs_hfs_pt = gen_hfs_T_with_rad ;
            gen_obs_hfs_phi = gen_hfs_phi_with_rad_eta_lt_4 ;
            gen_obs_dphi_hfs_eta_lt_4 = gen_obs_e_phi - gen_obs_hfs_phi ;
         }
      }

     //-- define dphi distribution that is centered at pi.
      if ( gen_obs_dphi_hfs_eta_lt_4 < -3.14159265 ) gen_obs_dphi_hfs_eta_lt_4 += 2 * 3.14159265 ;
      if ( gen_obs_dphi_hfs_eta_lt_4 >  3.14159265 ) gen_obs_dphi_hfs_eta_lt_4 -= 2 * 3.14159265 ;
      if ( gen_obs_dphi_hfs_eta_lt_4 < 0 ) gen_obs_dphi_hfs_eta_lt_4 += 2 * 3.14159265 ;






    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //-- Reconstruction-based analysis here.

      if ( verbose ) printf("\n\n Begin reconstruction-level section for event %d\n", ei ) ;


      float HFS_px = 0. ;
      float HFS_py = 0. ;
      float HFS_pz = 0. ;
      float HFS_E  = 0. ;


      float mcmatch_sum_px(0.) ;
      float mcmatch_sum_py(0.) ;
      float mcmatch_sum_pz(0.) ;
      float mcmatch_sum_E(0.) ;

      obs_fsr_photon_recombined = false ;

   //==== EF Candidates

    //-- Find the scattered electron in the energy flow tracks.

      int cheat_ele_efti(-1) ;
      int rad_gam_efgi(-1) ;
      float cheat_ele_iso(0.) ;

      for ( int i = 0; i < EFlowTrack_; i++ ) {

         if ( abs(EFlowTrack_PID[i]) != 11 ) continue ;

         float dr = calc_dr( sef_phi, EFlowTrack_Phi[i], sef_eta, EFlowTrack_Eta[i] ) ;

         if ( dr > 0.05 ) continue ;

         cheat_ele_efti = i ;

         if ( verbose ) printf("\n  %3d : electron EFlowTrack %3d, pid=11, dR = %7.4f   Pt = %9.4f,  Eta = %9.4f,  phi = %9.4f, charge=%d\n", ei, i, dr,
               EFlowTrack_PT[i], EFlowTrack_Eta[i], EFlowTrack_Phi[i], EFlowTrack_Charge[i] ) ;

         if ( verbose ) printf("\n Begin checking things within DeltaR of 0.5 around electron\n") ;

         for ( int j = 0; j < EFlowTrack_; j++ ) {
            if ( j == i ) continue ;
            float trk_dr = calc_dr( EFlowTrack_Phi[i], EFlowTrack_Phi[j], EFlowTrack_Eta[i], EFlowTrack_Eta[j] ) ;
            if ( trk_dr < 0.5 ) {
               cheat_ele_iso += EFlowTrack_PT[j] / EFlowTrack_PT[i] ;
               if ( verbose ) {
                  printf("             trk %3d in iso 0.5 cone:  dR = %9.4f,  Pt = %9.4f,  Eta = %9.4f,  phi = %9.4f\n",
                      j, trk_dr, EFlowTrack_PT[j], EFlowTrack_Eta[j], EFlowTrack_Phi[j] ) ;
               }
            }
         } // j
         for ( int j = 0; j < EFlowNeutralHadron_; j++ ) {
            float nh_dr = calc_dr( EFlowTrack_Phi[i], EFlowNeutralHadron_Phi[j], EFlowTrack_Eta[i], EFlowNeutralHadron_Eta[j] ) ;
            if ( nh_dr < 0.5 ) {
               cheat_ele_iso += EFlowNeutralHadron_ET[j] / EFlowTrack_PT[i] ;
               if ( verbose ) {
                  printf("             NH  %3d in iso 0.5 cone:  dR = %9.4f,  Pt = %9.4f,  Eta = %9.4f,  phi = %9.4f\n",
                      j, nh_dr, EFlowNeutralHadron_ET[j], EFlowNeutralHadron_Eta[j], EFlowNeutralHadron_Phi[j] ) ;
               }
            }
         } // j
         for ( int j = 0; j < EFlowPhoton_; j++ ) {
            float g_dr = calc_dr( EFlowTrack_Phi[i], EFlowPhoton_Phi[j], EFlowTrack_Eta[i], EFlowPhoton_Eta[j] ) ;
            if ( g_dr < 0.5 ) {
               if ( verbose ) {
                  printf("         photon  %3d in iso 0.5 cone:  dR = %9.4f,  Pt = %9.4f,  Eta = %9.4f,  phi = %9.4f\n",
                      j, g_dr, EFlowPhoton_ET[j], EFlowPhoton_Eta[j], EFlowPhoton_Phi[j] ) ;
               }
            }
            if ( g_dr * 180 / 3.14159265 < 4. ) {  // photon within 4 degrees?
               rad_gam_efgi = j ;
               if ( verbose ) printf("              found FSR photon:  dR = %9.5f ,   E = %9.5f,  Pt = %9.5f,  Eta = %9.5f,  phi = %9.5f\n",
                    g_dr, EFlowPhoton_E[j], EFlowPhoton_ET[j], EFlowPhoton_Eta[j], EFlowPhoton_Phi[j] ) ;
            }

         } // j
         if ( verbose ) printf(" End checking things within DeltaR of 0.5 around electron\n") ;

      } // i

      int ele_efti = cheat_ele_efti ;







//   //-- Try and do an honest scattered electron selection.


//    float ele_cand_highest_pt(0.) ;
//    float ele_cand_iso(0.) ;
//    int   ele_cand_efti(-1) ;
//    float ele_cand_emfr(0.) ;

//    for ( int i = 0; i < EFlowTrack_; i++ ) {

//       if ( EFlowTrack_Charge[i] != -1 ) continue ;

//       if ( EFlowTrack_PT[i] < ele_min_pt ) {
//          if ( verbose && i == cheat_ele_efti ) printf(" *** rejecting cheat electron because pt is too low: %9.3f\n", EFlowTrack_PT[i] ) ;
//          continue ;
//       }

//       float iso_sum(0.) ;
//       for ( int j = 0; j < EFlowTrack_; j++ ) {
//          if ( j == i ) continue ;
//          float trk_dr = calc_dr( EFlowTrack_Phi[i], EFlowTrack_Phi[j], EFlowTrack_Eta[i], EFlowTrack_Eta[j] ) ;
//          if ( trk_dr < 0.5 ) {
//             iso_sum += EFlowTrack_PT[j] / EFlowTrack_PT[i] ;
//          }
//       } // j
//       for ( int j = 0; j < EFlowNeutralHadron_; j++ ) {
//          float nh_dr = calc_dr( EFlowTrack_Phi[i], EFlowNeutralHadron_Phi[j], EFlowTrack_Eta[i], EFlowNeutralHadron_Eta[j] ) ;
//          if ( nh_dr < 0.5 ) {
//             iso_sum += EFlowNeutralHadron_ET[j] / EFlowTrack_PT[i] ;
//          }
//       } // j

//       if ( iso_sum > 0.20 ) {
//          if ( verbose && i == cheat_ele_efti ) printf(" *** rejecting cheat electron with iso : %9.3f\n", iso_sum ) ;
//          continue ;
//       }

//       int tower_index(-1) ;
//       float tower_mindr(999.) ;
//       for ( int j = 0; j < Tower_; j++ ) {
//          float tow_dr = calc_dr( EFlowTrack_Phi[i], Tower_Phi[j], EFlowTrack_Eta[i], Tower_Eta[j] ) ;
//          if ( tow_dr < tower_mindr ) {
//             tower_index = j ;
//             tower_mindr = tow_dr ;
//          }
//       }
//       if ( tower_index < 0 ) {
//          if (  verbose && i == cheat_ele_efti ) printf(" *** rejecting cheat electron because no matched tower.\n" ) ;
//          continue ;
//       }
//       float tower_em_frac(0.) ;
//       if ( Tower_Eem[tower_index] + Tower_Ehad[tower_index] > 0 ) tower_em_frac = Tower_Eem[tower_index] / ( Tower_Eem[tower_index] + Tower_Ehad[tower_index] ) ;

//       if ( tower_em_frac < 0.80 ) {
//          if ( verbose && i == cheat_ele_efti ) printf(" *** rejecting cheat electron because emfr = %9.3f.\n", tower_em_frac ) ;
//          continue ;
//       }

//       if ( EFlowTrack_PT[i] > ele_cand_highest_pt ) {
//          ele_cand_highest_pt = EFlowTrack_PT[i] ;
//          ele_cand_efti = i ;
//          ele_cand_iso = iso_sum ;
//          ele_cand_emfr = tower_em_frac ;
//       }

//    } // i

//    if ( verbose ) {
//       printf("   Honest electron selection:  track %d , Pt = %9.2f, iso = %9.3f  emfr = %9.3f ", ele_cand_efti, ele_cand_highest_pt, ele_cand_iso, ele_cand_emfr  ) ;
//       if ( ele_cand_efti != cheat_ele_efti ) printf(" *** does not agree with cheat!  %d  iso %9.3f", cheat_ele_efti, cheat_ele_iso ) ;
//       printf("\n") ;
//    }


//    ele_efti = ele_cand_efti ;  //-- use honest selction, not cheat.

//    if ( ele_cand_efti == cheat_ele_efti ) {
//       e_is_correct_track = true ;
//    } else {
//       e_is_correct_track = false ;
//    }




      if ( verbose ) printf("\n") ;

    //-- Now, go through tracks for HFS.

      for ( int i = 0; i < EFlowTrack_; i++ ) {

         if ( i == ele_efti ) continue ;

         float px, py, pz, e ;

         if ( fabs(EFlowTrack_Eta[i]) > 4.0 ) continue ;

         if ( EFlowTrack_PT[i] < minTrkPt ) continue ;

         ptep_to_xyze( EFlowTrack_PT[i], EFlowTrack_Eta[i], EFlowTrack_Phi[i], EFlowTrack_Mass[i],
                       px, py, pz, e ) ;

         if ( e < 0 ) continue ;

         if ( fabs(EFlowTrack_Eta[i]) > maxEta ) continue ;

         HFS_px += px ;
         HFS_py += py ;
         HFS_pz += pz ;
         HFS_E  += e  ;

         if ( verbose ) printf( " %3d, trk %3d : pid = %5d  pt = %7.2f  px,py,pz = %7.2f, %7.2f, %7.2f   E = %7.2f    eta = %8.3f, phi = %8.3f,  E-pz  %9.5f   sum E-pz  %9.5f\n",
           ei, i, EFlowTrack_PID[i], EFlowTrack_PT[i], px, py, pz, e, EFlowTrack_Eta[i], EFlowTrack_Phi[i], e-pz, HFS_E-HFS_pz ) ;

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
               mcmatch_sum_px += Particle_Px[mc_pi] ;
               mcmatch_sum_py += Particle_Py[mc_pi] ;
               mcmatch_sum_pz += Particle_Pz[mc_pi] ;
               mcmatch_sum_E += Particle_E[mc_pi] ;
               printf( " %3d, trk %3d :                  MC match  px,py,pz = %7.2f, %7.2f, %7.2f,  dr = %8.4f, eta = %8.3f, phi = %8.3f,  E-pz  %9.5f   sum E-pz  %9.5f mc\n",
                  ei, i, Particle_Px[mc_pi], Particle_Py[mc_pi], Particle_Pz[mc_pi], min_mc_dr, Particle_Eta[mc_pi], Particle_Phi[mc_pi], Particle_E[mc_pi]-Particle_Pz[mc_pi],
                  mcmatch_sum_E-mcmatch_sum_pz) ;
            }
         }



      } // i

      if ( verbose ) printf( "  --- after tracks:  HFS :  px,py,pz = %7.2f , %7.2f, %7.2f   E = %7.2f  E-pz = %9.5f\n", HFS_px, HFS_py, HFS_pz, HFS_E, HFS_E-HFS_pz ) ;





      for ( int i = 0; i < EFlowPhoton_; i++ ) {

         float px, py, pz, e ;

         if ( EFlowPhoton_Eta[i] > 4.0 ) continue ;

         if ( EFlowPhoton_E[i] < minPhoE ) continue ;

         if ( i == rad_gam_efgi ) {
            if ( verbose ) printf("   skipping photon %d (FSR photon)\n", i ) ;
            continue ;
         }

         ptep_to_xyze( EFlowPhoton_ET[i], EFlowPhoton_Eta[i], EFlowPhoton_Phi[i], 0.,
                       px, py, pz, e ) ;

         if ( e < 0 ) continue ;

         if ( fabs(EFlowPhoton_Eta[i]) > maxEta ) continue ;

         float dr_wrt_se = calc_dr( EFlowPhoton_Phi[i], sef_phi , EFlowPhoton_Eta[i], sef_eta  ) ;

         if ( dr_wrt_se < 0.5 ) continue ; //--- reject any photons close to scattered electron for HFS.

         HFS_px += px ;
         HFS_py += py ;
         HFS_pz += pz ;
         HFS_E  += e  ;


         if ( verbose ) printf( " %3d, pho %3d :  pt = %7.2f  px,py,pz = %7.2f, %7.2f, %7.2f   E = %7.2f    eta = %8.3f, phi = %8.3f,  E-pz  %9.5f   sum E-pz  %9.5f  dR wrt se %9.5f\n",
           ei, i, EFlowPhoton_ET[i], px, py, pz, e, EFlowPhoton_Eta[i], EFlowPhoton_Phi[i], e-pz, HFS_E-HFS_pz,  dr_wrt_se ) ;

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
               mcmatch_sum_px += Particle_Px[mc_pi] ;
               mcmatch_sum_py += Particle_Py[mc_pi] ;
               mcmatch_sum_pz += Particle_Pz[mc_pi] ;
               mcmatch_sum_E += Particle_E[mc_pi] ;
               printf( " %3d, pho %3d :      MC match  px,py,pz = %7.2f, %7.2f, %7.2f,  dr = %8.4f, eta = %8.3f, phi = %8.3f,  E-pz  %9.5f   sum E-pz  %9.5f mc\n",
                  ei, i, Particle_Px[mc_pi], Particle_Py[mc_pi], Particle_Pz[mc_pi], min_mc_dr, Particle_Eta[mc_pi], Particle_Phi[mc_pi], Particle_E[mc_pi]-Particle_Pz[mc_pi],
                  mcmatch_sum_E-mcmatch_sum_pz) ;
            }
         }

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


         if ( fabs(EFlowNeutralHadron_Eta[i]) > maxEta ) continue ;


         HFS_px += px ;
         HFS_py += py ;
         HFS_pz += pz ;
         HFS_E  += e  ;

      } // i

      if ( add_noise ) {
         float urn = tran->Uniform() ;
         if ( urn < noise_prob ) {

            float grn_px = tran->Landau( noise_landau_mu, noise_landau_sigma  ) ;
            float grn_py = tran->Landau( noise_landau_mu, noise_landau_sigma  ) ;
            float grn_pz = tran->Landau( noise_landau_mu, noise_landau_sigma  ) ;

            if ( tran->Uniform() < 0.5 ) { grn_px = -1*grn_px ; }
            if ( tran->Uniform() < 0.5 ) { grn_py = -1*grn_py ; }
            if ( tran->Uniform() < 0.5 ) { grn_pz = -1*grn_pz ; }

            float grn_e = sqrt( grn_px*grn_px + grn_py*grn_py + grn_pz*grn_pz ) ;

            HFS_px += grn_px ;
            HFS_py += grn_py ;
            HFS_pz += grn_pz ;
            HFS_E += grn_e ;
         }
      }


      obs_hfs_pz = HFS_pz ;
      obs_hfs_e  = HFS_E ;
      obs_hfs_phi = atan2( HFS_py, HFS_px ) ;
      obs_hfs_pt = sqrt( HFS_px * HFS_px  +  HFS_py * HFS_py ) ;

      if ( verbose ) printf( "  --- after NH    :  HFS :  px,py,pz = %7.2f , %7.2f, %7.2f   E = %7.2f\n", HFS_px, HFS_py, HFS_pz, HFS_E ) ;

      if ( verbose ) printf( "\n  reco HFS :  Sigma  %9.4f (%9.4f)   T  %9.4f (%9.4f)\n\n",
            (obs_hfs_e-obs_hfs_pz),  (gen_obs_hfs_e-gen_obs_hfs_pz),
            obs_hfs_pt, gen_obs_hfs_pt ) ;


      if ( use_electron_cand_not_track ) {

      //-- Get the reconstructed scattered electron

         if ( Electron_ <= 0 ) {
            if ( verbose ) printf(" \n\n *** no reconstructed electrons in event\n\n") ;
            continue ;
         }

         float drmin = 9999999. ;
         int best_ei = -1 ;
         for ( int i = 0; i < Electron_; i ++ ) {
            float dr = calc_dr( sef_phi, Electron_Phi[i], sef_eta, Electron_Eta[i] ) ;
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


         obs_e_pt = Electron_PT[best_ei] ;
         obs_e_eta = Electron_Eta[best_ei] ;
         obs_e_phi = Electron_Phi[best_ei] ;



        float obs_e_px ;
        float obs_e_py ;
        ptep_to_xyze( obs_e_pt, obs_e_eta, obs_e_phi, 0.000511,
                            obs_e_px, obs_e_py, obs_e_pz, obs_e_e ) ;

        if ( verbose ) {
           printf("\n Electron collection\n") ;
           for ( int i = 0; i < Electron_; i ++ ) {
              printf("   Electron collection     electron:  pt = %9.4f (%9.4f)   phi = %9.4f (%9.4f)   eta = %9.5f (%9.5f)\n",
                 Electron_PT[i], sef_pt,  Electron_Phi[i], sef_phi,  Electron_Eta[i], sef_eta ) ;
           } // i
        }


     } else {


        if ( ele_efti < 0 ) {
           if ( verbose ) printf("\n\n *** can't find reco scattered electron track.\n\n") ;
           continue ;
        }

        obs_e_pt = EFlowTrack_PT[ele_efti] ;
        obs_e_eta = EFlowTrack_Eta[ele_efti] ;
        obs_e_phi = EFlowTrack_Phi[ele_efti] ;
        if ( smear_electron_theta ) {
           float pt0 = obs_e_pt ;
           float eta0 = obs_e_eta ;
           float phi0 = obs_e_phi ;
           smear_theta( pt0, eta0, phi0, obs_e_pt, obs_e_eta, obs_e_phi, verbose ) ;
        }


        float obs_e_px ;
        float obs_e_py ;
        ptep_to_xyze( obs_e_pt, obs_e_eta, obs_e_phi, 0.000511,
                            obs_e_px, obs_e_py, obs_e_pz, obs_e_e ) ;

        if ( rad_gam_efgi >= 0 && obs_recombine_fsr_gamma ) {

           float gam_px, gam_py, gam_pz, gam_e ;
           ptep_to_xyze( EFlowPhoton_ET[rad_gam_efgi], EFlowPhoton_Eta[rad_gam_efgi], EFlowPhoton_Phi[rad_gam_efgi], 0.0,
                            gam_px, gam_py, gam_pz, gam_e ) ;

           obs_e_px += gam_px ;
           obs_e_py += gam_py ;
           obs_e_pz += gam_pz ;
           obs_e_e += gam_e ;

           obs_e_pt = sqrt( obs_e_px * obs_e_px  +  obs_e_py * obs_e_py ) ;
           obs_e_phi = atan2( obs_e_py, obs_e_px ) ;
           obs_e_eta = -1. * log( tan(   atan2( obs_e_pt, obs_e_pz ) / 2. ) ) ;

           obs_fsr_photon_recombined = true ;

        }

     }

     if ( verbose ) printf("\n   reconstructed scattered electron:  pt = %9.4f (%9.4f)   phi = %9.4f (%9.4f)   eta = %9.5f (%9.5f)\n\n",
         obs_e_pt, sef_pt,  obs_e_phi, sef_phi, obs_e_eta, sef_eta ) ;


      obs_dphi = obs_e_phi - obs_hfs_phi ;

     //-- define dphi distribution that is centered at pi.
      if ( obs_dphi < -3.14159265 ) obs_dphi += 2 * 3.14159265 ;
      if ( obs_dphi >  3.14159265 ) obs_dphi -= 2 * 3.14159265 ;
      if ( obs_dphi < 0 ) obs_dphi += 2 * 3.14159265 ;




     //-- save photons


      n_pho = 0 ;
      pho_isr_match.clear() ;
      pho_fsr_match.clear() ;
      pho_angle_wrt_ebeam.clear() ;
      pho_angle_wrt_efinal.clear() ;
      pho_e.clear() ;
      pho_pt.clear() ;
      pho_pz.clear() ;
      pho_eta.clear() ;
      pho_dphi_wrt_efinal.clear() ;

      phi_pho_closest_to_ebeam = 0. ;
      eta_pho_closest_to_ebeam = 9. ;
      e_pho_closest_to_ebeam = -1. ;

      for ( int pi=0; pi<Photon_; pi++ ) {

         pho_e.emplace_back( Photon_E[pi] ) ;
         pho_pt.emplace_back( Photon_PT[pi] ) ;
         pho_eta.emplace_back( Photon_Eta[pi] ) ;

         float theta = 2. * atan( exp( -1.* Photon_Eta[pi] ) ) ;
         float pz = Photon_E[pi] * cos( theta ) ;
         float px = Photon_PT[pi] * cos( Photon_Phi[pi] ) ;
         float py = Photon_PT[pi] * sin( Photon_Phi[pi] ) ;

         pho_pz.emplace_back( pz ) ;

         float dr_wrt_sef = calc_dr( sef_phi, Photon_Phi[pi], sef_eta, Photon_Eta[pi] ) ;

         if ( dr_wrt_sef > 0.5 && Photon_Eta[pi] < eta_pho_closest_to_ebeam ) {
            phi_pho_closest_to_ebeam = Photon_Phi[pi] ;
            eta_pho_closest_to_ebeam = Photon_Eta[pi] ;
            e_pho_closest_to_ebeam = Photon_E[pi] ;
            if ( Photon_Phi[pi] > 3.2 ) {
               printf("\n\n %5d : *** BAD Photon_Phi??? : %9.5f\n\n", ei, Photon_Phi[pi] ) ;
            }
         }

         float cos_angle_wrt_beam = -1. * beam_e_e * pz / (beam_e_e * Photon_E[pi]) ;

         pho_angle_wrt_ebeam.emplace_back( acos( cos_angle_wrt_beam ) ) ;

         float angle_wrt_efinal = calc_angle( sef_pt, sef_eta, sef_phi,
                                              Photon_PT[pi], Photon_Eta[pi], Photon_Phi[pi] ) ;
         pho_angle_wrt_efinal.emplace_back( angle_wrt_efinal ) ;

         float isr_angle(9.) ;
         if ( has_isr ) isr_angle = calc_angle( isr_gamma_pt, isr_gamma_eta, isr_gamma_phi,
                                                Photon_PT[pi], Photon_Eta[pi], Photon_Phi[pi] ) ;

         bool isr_match = false ;
         if ( isr_angle < 0.05 ) isr_match = true ;

         pho_isr_match.push_back( isr_match ) ;



         float fsr_angle(9.) ;
         if ( has_fsr ) fsr_angle = calc_angle( fsr_gamma_pt, fsr_gamma_eta, fsr_gamma_phi,
                                                Photon_PT[pi], Photon_Eta[pi], Photon_Phi[pi] ) ;

         bool fsr_match = false ;
         if ( fsr_angle < 0.05 ) fsr_match = true ;

         pho_fsr_match.push_back( fsr_match ) ;



         float dphi = sef_phi - Photon_Phi[pi] ;
         if ( dphi >  3.14159265 ) dphi -= 2 * 3.14159265 ;
         if ( dphi < -3.14159265 ) dphi += 2 * 3.14159265 ;
         pho_dphi_wrt_efinal.emplace_back( dphi ) ;

         if ( verbose ) {
            printf("  %5d : pho %3d :  E %9.5f  pt %9.5f  pz  %9.5f  eta %9.5f  phi %9.5f  angle(ebeam)  %9.5f  angle(efinal)  %9.5f",
               ei, pi, Photon_E[pi], Photon_PT[pi], pz, Photon_Eta[pi], Photon_Phi[pi], acos( cos_angle_wrt_beam ), angle_wrt_efinal ) ;
            if ( isr_match ) printf(" *** ISR match : %9.5f", isr_angle ) ;
            if ( fsr_match ) printf(" *** FSR match : %9.5f", fsr_angle ) ;
            printf("\n") ;
         }

         n_pho ++ ;

      } // pi



     //-- Calculate tower sum around final scattered electron direction

      tower_sum_07 = 0. ;
      tower_sum_20 = 0. ;
      tower_sum_30 = 0. ;
      tower_sum_40 = 0. ;

      n_towers_07 = 0 ;
      n_towers_20 = 0 ;
      n_towers_30 = 0 ;
      n_towers_40 = 0 ;

      for ( int ti=0; ti<Tower_; ti++ ) {

         if ( Tower_Eem[ti] < 0.05 ) continue ;

         float dr = calc_dr( Tower_Phi[ti], sef_phi,  Tower_Eta[ti], sef_eta ) ;

         if ( dr < 0.07 ) { tower_sum_07 += Tower_Eem[ti] ; n_towers_07 ++ ; }
         if ( dr < 0.20 ) { tower_sum_20 += Tower_Eem[ti] ; n_towers_20 ++ ; }
         if ( dr < 0.30 ) { tower_sum_30 += Tower_Eem[ti] ; n_towers_30 ++ ; }
         if ( dr < 0.40 ) { tower_sum_40 += Tower_Eem[ti] ; n_towers_40 ++ ; }

         if ( verbose && dr < 0.40 ) {
            printf("  %5d :  tower  %3d :  dr = %9.5f  E %9.5f  Eem %9.5f\n", ei, ti, dr, Tower_E[ti], Tower_Eem[ti] ) ;
         }

      } // ti

      if ( verbose ) {
         printf("\n") ;
         printf("  %5d :  tower sum for dr<0.07 :  %9.5f, sef E = %9.5f,  ratio  %9.5f\n", ei, tower_sum_07, sef_e, tower_sum_07/sef_e ) ;
         printf("  %5d :  tower sum for dr<0.20 :  %9.5f, sef E = %9.5f,  ratio  %9.5f\n", ei, tower_sum_20, sef_e, tower_sum_20/sef_e ) ;
         printf("  %5d :  tower sum for dr<0.30 :  %9.5f, sef E = %9.5f,  ratio  %9.5f\n", ei, tower_sum_30, sef_e, tower_sum_30/sef_e ) ;
         printf("  %5d :  tower sum for dr<0.40 :  %9.5f, sef E = %9.5f,  ratio  %9.5f\n", ei, tower_sum_40, sef_e, tower_sum_40/sef_e ) ;
         printf("\n") ;
      }



      {

         int im ;

         float e ;
         float e0 ;
         float theta ;
         float Sigma ;
         float T ;
         float tan_gamma_over_2 ;
         float ep ;

         e = obs_e_e ;
         e0 = be0_e ;
         theta = atan2( obs_e_pt, obs_e_pz ) ;
         Sigma = obs_hfs_e - obs_hfs_pz  ;
         T = obs_hfs_pt ;
         tan_gamma_over_2 = Sigma / T ;
         ep = beam_p_e ;



        //-- Method  0 : electron
         im = 0 ;

         obs_y[im] = 1. - (e/e0) * pow( sin(theta/2.), 2 ) ;

         obs_Q2[im] = 4. * e0 * e * pow( cos(theta/2.), 2 ) ;

         obs_x[im] = ( e0 * e * pow( cos(theta/2.), 2 ) ) / ( ep * ( e0 - e * pow( sin(theta/2), 2 ) ) ) ;




        //-- Method  1 : unpub1
         im = 1 ;

         obs_y[im] = Sigma / ( 2. * e0 ) ; // yh

         obs_Q2[im] = 4. * e0 * e - 4. * e0 * e0 + 2 * e0 * Sigma ;

         obs_x[im] = obs_Q2[im] / ( 2. * ep * Sigma ) ;





        //-- Method  2 : unpub2
         im = 2 ;

         obs_y[im] = Sigma / ( 2. * e0 ) ; // yh

         obs_Q2[im] = 2. * e0 * ( 2. * e0 - Sigma ) / pow( tan(theta/2.), 2. ) ;

         obs_x[im] = obs_Q2[im] / ( 2. * ep * Sigma ) ;




        //-- Method  3 : DA
         im = 3 ;

         obs_y[im] = tan_gamma_over_2 / ( tan_gamma_over_2 + tan(theta/2.) ) ;

         obs_Q2[im] = 4. * e0 * e0 * (1./tan(theta/2.)) / ( tan_gamma_over_2 + tan(theta/2.) ) ;

         obs_x[im] = obs_Q2[im] / ( 4. * ep * e0 * obs_y[im] ) ;




        //-- Method  4 : h
         im = 4 ;

         obs_y[im] = Sigma / ( 2. * e0 ) ;

         obs_Q2[im] = T * T / ( 1. - obs_y[im] ) ;

         obs_x[im] = obs_Q2[im] / ( 2. * ep * Sigma ) ;




        //-- Method  5 : ISigma
         im = 5 ;


         obs_y[im] = Sigma / ( Sigma + e * ( 1. - cos(theta) ) ) ;

         obs_Q2[im] = e * e * pow( sin(theta), 2 ) / ( 1. - obs_y[im] ) ;

         obs_x[im] = e * pow( cos(theta/2.), 2 ) / ( ep * obs_y[im] ) ;




        //-- Method  6 : IDA
         im = 6 ;

         obs_y[im] = tan_gamma_over_2 / ( tan_gamma_over_2 + tan(theta/2.) ) ; // yDA

         //obs_Q2[im] = e * e * tan(theta/2.) * ( tan_gamma_over_2 + tan(theta/2.) ) / ( 1./tan(theta/2.) + tan(theta/2.) ) ; // *** this is wrong somehow
         obs_Q2[im] = 4. * e * e * (1./tan(theta/2.)) * ( tan_gamma_over_2 + tan(theta/2.) ) / pow( (1./tan(theta/2.) + tan(theta/2.)), 2 ) ;

         obs_x[im] = e * pow( cos(theta/2.), 2. ) / ( ep * obs_y[im] ) ;





        //-- Method  7 : unpub3
         im = 7 ;

         obs_y[im] = tan_gamma_over_2 / ( tan_gamma_over_2 + tan(theta/2.) ) ; // yDA

         obs_Q2[im] = pow( Sigma / tan_gamma_over_2, 2 ) / ( 1. - obs_y[im] ) ;

         obs_x[im] = e * pow( cos(theta/2.), 2. ) / ( ep * obs_y[im] ) ;





        //-- Method  8 : eSigma
         im = 8 ;

         obs_y[im] = 2. * e0 * Sigma / pow( ( Sigma + e*(1.-cos(theta)) ) , 2 ) ;

         obs_Q2[im] = 4. * e0 * e * pow( cos(theta/2.), 2 ) ; // Q2e

         obs_x[im] = ( e * pow( cos(theta/2), 2 ) + (e * e/(2.*Sigma)) * pow( sin(theta), 2 ) ) / ep ; // *** fixed typo in note



      }



      if ( verbose ) {

         printf("                                   " ) ;
         for ( int i=0; i<nmeth; i++ ) { printf("   %6d   ", i) ; }
         printf("\n") ;


         printf("  %4d :  obs_y   : ", ei ) ;
         printf("  %9.5f  | ", from_tlv_gen_y ) ;
         for ( int i=0; i<nmeth; i++ ) {
            printf("  %9.5f ", obs_y[i] ) ;
         } // i
         printf("\n") ;

         printf("  %4d :  obs_x   : ", ei ) ;
         printf("  %9.5f  | ", from_tlv_gen_x ) ;
         for ( int i=0; i<nmeth; i++ ) {
            printf("  %9.5f ", obs_x[i] ) ;
         } // i
         printf("\n") ;

         printf("  %4d :  obs_Q2  : ", ei ) ;
         printf("  %9.1f  | ", from_tlv_gen_Q2 ) ;
         for ( int i=0; i<nmeth; i++ ) {
            printf("  %9.1f ", obs_Q2[i] ) ;
         } // i
         printf("\n") ;


         ///float reso_x_da = obs_x[3] / from_tlv_gen_x ;
         ///if ( !has_isr && !has_fsr && from_tlv_gen_y > 0.01 && from_tlv_gen_y < 0.05 && reso_x_da < 0.5 ) {
         ///   if (verbose) printf("\n\n  *** %5d :  check this event:  y_true = %9.5f    x/x_true = %9.5f\n\n", ei, from_tlv_gen_y, reso_x_da ) ;
         ///}

         ///if ( !has_isr && !has_fsr && from_tlv_gen_y > 0.5 && obs_x[0] / from_tlv_gen_x > 1.1 ) {
         ///   if (verbose) printf("\n\n  *** %5d :  check this event:  obs_x[0] / from_tlv_gen_x = %9.5f / %9.5f = %9.5f\n\n", ei, obs_x[0], from_tlv_gen_x, obs_x[0] / from_tlv_gen_x ) ;
         ///}

      }




      tt_out -> Fill() ;

   } // jentry

   printf("\n\n Done.\n\n") ;

   tt_out -> Write() ;

   if ( !save_only_dnn_inputs ) {
      TH1F* h_meth_names = new TH1F( "h_meth_names", "Method names", nmeth, 0.5, 0.5+nmeth ) ;
      for ( int im=0; im<nmeth; im++ ) {  h_meth_names -> GetXaxis() -> SetBinLabel( im+1, meth_name[im] ) ; }
      h_meth_names -> LabelsOption("v") ;
      h_meth_names -> Write() ;
   }


   tf_out -> Close() ;

   printf("  Saved results in %s\n\n", out_fname ) ;

}


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

   px = pt * cos( phi ) ;
   py = pt * sin( phi ) ;
   pz = pt  * sinh( eta ) ;

   float p2 = px*px + py*py + pz*pz ;
   e = sqrt( p2 + m*m ) ;

}


//==========

float calc_angle( double pt1, double eta1, double phi1,
                  double pt2, double eta2, double phi2 ) {

   float px1, py1, pz1, p1 ;
   ptep_to_xyze( pt1, eta1, phi1, 0., px1, py1, pz1, p1 ) ;

   float px2, py2, pz2, p2 ;
   ptep_to_xyze( pt2, eta2, phi2, 0., px2, py2, pz2, p2 ) ;

   float dot = px1 * px2  +  py1 * py2  +  pz1 * pz2 ;

   float cos_angle = dot / (p1 * p2) ;

   return acos( cos_angle ) ;

}

//==========

void smear_theta( float pt0, float eta0, float phi0,
                  float& pts, float& etas, float& phis, bool verbose ) {

   float theta0 = 2 * atan( exp( -1*eta0 ) ) ;
   float p0 = pt0 / sin( theta0 ) ;

  //-- parameterization of theta resolution from slides here: https://indico.bnl.gov/event/12293/
  //   See draw_trk_angle_reso.c
  //
   float sigma_theta = ((0.048839+0.001190*eta0+0.011905*eta0*eta0)/(0.048839))*(0.439387/p0+0.409963/(p0*p0)-0.227344/(p0*p0*p0)) / 1000. ;

   float dtheta = tran->Gaus(0., sigma_theta) ;

   float thetas = theta0 + dtheta ;

   pts = p0 * sin(thetas) ;

   etas = -1 * log( tan( thetas/2. ) ) ;

   phis = phi0 ;

   if ( verbose ) {
      printf("\n") ;
      printf("   smear_theta:  p0 = %9.3f  eta0 = %9.5f  sigma_theta = %10.6f (mrad)  dtheta = %10.6f (mrad)\n", p0, eta0, 1000*sigma_theta, 1000*dtheta ) ;
      printf("   smear_theta:  pt  %9.5f (%9.5f)   eta %9.5f (%9.5f)\n", pt0, pts, eta0, etas ) ;
      printf("\n") ;
   }

}

//==========


