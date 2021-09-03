#define delphes_rapgap_isr_fsr_analysis_cxx
#include "delphes_rapgap_isr_fsr_analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void delphes_rapgap_isr_fsr_analysis::Loop( bool verbose, int maxEvts )
{

   const int nmeth = 9 ;

   if (fChain == 0) return;

   float escale = 1000. ;

   TFile* tf_out = new TFile("rad-tree.root", "recreate") ;
   TTree* tt_out = new TTree("minitree", "Minimal flat TTree for hadronic Q2,x,y analysis" ) ;

   float gen_y[nmeth] ;
   float gen_x[nmeth] ;
   float gen_Q2[nmeth] ;

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

   float gen_hfs_Sigma_no_rad  ;
   float gen_hfs_T_no_rad  ;

   float gen_hfs_Sigma_with_rad  ;
   float gen_hfs_T_with_rad  ;

   float gen_hfs_Sigma_with_rad_eta_lt_4  ;
   float gen_hfs_T_with_rad_eta_lt_4  ;

   char bname[100] ;


   sprintf( bname, "gen_y[%d]/F", nmeth ) ;
   tt_out -> Branch( "gen_y", &gen_y, bname ) ;

   sprintf( bname, "gen_x[%d]/F", nmeth ) ;
   tt_out -> Branch( "gen_x", &gen_x, bname ) ;

   sprintf( bname, "gen_Q2[%d]/F", nmeth ) ;
   tt_out -> Branch( "gen_Q2", &gen_Q2, bname ) ;

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

   tt_out -> Branch( "sef_e", &sef_e, "sef_e/F" ) ;
   tt_out -> Branch( "sef_pt", &sef_pt, "sef_pt/F" ) ;
   tt_out -> Branch( "sef_pz", &sef_pz, "sef_pz/F" ) ;
   tt_out -> Branch( "sef_phi", &sef_phi, "sef_phi/F" ) ;

   tt_out -> Branch( "beam_e_e", &beam_e_e, "beam_e_e/F" ) ;
   tt_out -> Branch( "beam_p_e", &beam_p_e, "beam_p_e/F" ) ;

   tt_out -> Branch( "gen_e0_Q2", &gen_e0_Q2, "gen_e0_Q2/F" ) ;
   tt_out -> Branch( "gen_e0_y", &gen_e0_y, "gen_e0_y/F" ) ;
   tt_out -> Branch( "gen_e0_x", &gen_e0_x, "gen_e0_x/F" ) ;

   tt_out -> Branch( "gen_ef_Q2", &gen_ef_Q2, "gen_ef_Q2/F" ) ;
   tt_out -> Branch( "gen_ef_y", &gen_ef_y, "gen_ef_y/F" ) ;
   tt_out -> Branch( "gen_ef_x", &gen_ef_x, "gen_ef_x/F" ) ;

   tt_out -> Branch( "has_isr", &has_isr, "has_isr/B" ) ;
   tt_out -> Branch( "has_fsr", &has_fsr, "has_fsr/B" ) ;

   tt_out -> Branch( "gen_hfs_Sigma_no_rad", &gen_hfs_Sigma_no_rad, "gen_hfs_Sigma_no_rad/F" ) ;
   tt_out -> Branch( "gen_hfs_T_no_rad", &gen_hfs_T_no_rad, "gen_hfs_T_no_rad/F" ) ;

   tt_out -> Branch( "gen_hfs_Sigma_with_rad", &gen_hfs_Sigma_with_rad, "gen_hfs_Sigma_with_rad/F" ) ;
   tt_out -> Branch( "gen_hfs_T_with_rad", &gen_hfs_T_with_rad, "gen_hfs_T_with_rad/F" ) ;

   tt_out -> Branch( "gen_hfs_Sigma_with_rad_eta_lt_4", &gen_hfs_Sigma_with_rad_eta_lt_4, "gen_hfs_Sigma_with_rad_eta_lt_4/F" ) ;
   tt_out -> Branch( "gen_hfs_T_with_rad_eta_lt_4", &gen_hfs_T_with_rad_eta_lt_4, "gen_hfs_T_with_rad_eta_lt_4/F" ) ;

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
            printf("   %4d : %3d :  PID %10d  %15s  status = %3d  M1 %3d M2 %3d  D1 %3d D2 %3d   Pt %15.6f  Pz %15.6f  Phi %9.5f ", ei, pi,
            Particle_PID[pi], pdg_name, Particle_Status[pi],
            Particle_M1[pi],
            Particle_M2[pi],
            Particle_D1[pi],
            Particle_D2[pi],
            escale * Particle_PT[pi],
            escale * Particle_Pz[pi],
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

      has_isr = false ;
      has_fsr = false ;

      for ( int pi = 0; pi < Particle_; pi ++ ) {

         if ( Particle_PID[pi] == 2212 && Particle_Status[pi] == 4 ) {
            beam_p_e = escale * Particle_E[pi] ;
            continue ;
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
               } else {
                  has_fsr = true ;
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





      gen_hfs_Sigma_no_rad = 0. ;
      float gen_hfs_px_no_rad = 0. ;
      float gen_hfs_py_no_rad = 0. ;

      gen_hfs_Sigma_with_rad = 0. ;
      float gen_hfs_px_with_rad = 0. ;
      float gen_hfs_py_with_rad = 0. ;

      gen_hfs_Sigma_with_rad_eta_lt_4 = 0. ;
      float gen_hfs_px_with_rad_eta_lt_4 = 0. ;
      float gen_hfs_py_with_rad_eta_lt_4 = 0. ;

      for ( int pi = 0; pi < Particle_; pi ++ ) {

         if ( pi == sef_pi ) continue ;

         if ( Particle_Status[pi] != 1 ) continue ;

         gen_hfs_Sigma_with_rad += escale * ( Particle_E[pi] - Particle_Pz[pi] ) ;

         gen_hfs_px_with_rad += escale * Particle_Px[pi] ;
         gen_hfs_py_with_rad += escale * Particle_Py[pi] ;

         if ( fabs( Particle_Eta[pi] ) < 4.0 ) {
            gen_hfs_Sigma_with_rad_eta_lt_4 += escale * ( Particle_E[pi] - Particle_Pz[pi] ) ;
            gen_hfs_px_with_rad_eta_lt_4 += escale * Particle_Px[pi] ;
            gen_hfs_py_with_rad_eta_lt_4 += escale * Particle_Py[pi] ;
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
         if ( verbose ) { printf("  HFS : %3d :  E = %9.4f , pz = %9.3f , px = %9.3f , py = %9.3f\n", pi,
           escale * Particle_E[pi], escale * Particle_Pz[pi], escale * Particle_Px[pi], escale * Particle_Py[pi] ) ; }

      } // pi

      gen_hfs_T_no_rad = sqrt( gen_hfs_px_no_rad * gen_hfs_px_no_rad + gen_hfs_py_no_rad * gen_hfs_py_no_rad ) ;
      gen_hfs_T_with_rad = sqrt( gen_hfs_px_with_rad * gen_hfs_px_with_rad + gen_hfs_py_with_rad * gen_hfs_py_with_rad ) ;
      gen_hfs_T_with_rad_eta_lt_4 = sqrt( gen_hfs_px_with_rad_eta_lt_4 * gen_hfs_px_with_rad_eta_lt_4 + gen_hfs_py_with_rad_eta_lt_4 * gen_hfs_py_with_rad_eta_lt_4 ) ;




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
         }
         printf(" %4d : %3d be0_pz = %9.3f , %3d bef_pz = %9.3f , diff = %9.3f\n", ei, be0_pi, be0_pz, bef_pi, bef_pz, be0_pz-bef_pz ) ;
         printf(" %4d : %3d se0_pt = %9.3f , %3d sef_pt = %9.3f , diff = %9.3f\n", ei, se0_pi, se0_pt, sef_pi, sef_pt, se0_pt-sef_pt ) ;
         printf(" %4d :     HFS pt = %9.3f , HFS Sigma = %9.3f\n", ei, gen_hfs_T_no_rad, gen_hfs_Sigma_no_rad ) ;
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

         e = sef_e ;
         e0 = be0_e ;
         theta = atan2( sef_pt, sef_pz ) ;
         Sigma = gen_hfs_Sigma_no_rad ;
         T = gen_hfs_T_no_rad ;
         tan_gamma_over_2 = Sigma / T ;
         ep = beam_p_e ;


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

         printf("                    " ) ;
         for ( int i=0; i<nmeth; i++ ) { printf("   %6d   ", i) ; }
         printf("\n") ;


         printf("  %4d :  gen_y   : ", ei ) ;
         for ( int i=0; i<nmeth; i++ ) {
            printf("  %9.5f ", gen_y[i] ) ;
         } // i
         printf("\n") ;

         printf("  %4d :  gen_x   : ", ei ) ;
         for ( int i=0; i<nmeth; i++ ) {
            printf("  %9.5f ", gen_x[i] ) ;
         } // i
         printf("\n") ;

         printf("  %4d :  gen_Q2  : ", ei ) ;
         for ( int i=0; i<nmeth; i++ ) {
            printf("  %9.1f ", gen_Q2[i] ) ;
         } // i
         printf("\n") ;


      }


      tt_out -> Fill() ;

   } // jentry

   printf("\n\n Done.\n\n") ;

   TH1F* h_meth_names = new TH1F( "h_meth_names", "Method names", nmeth, 0.5, 0.5+nmeth ) ;
   for ( int im=0; im<nmeth; im++ ) {  h_meth_names -> GetXaxis() -> SetBinLabel( im+1, meth_name[im] ) ; }
   h_meth_names -> LabelsOption("v") ;

   tt_out -> Write() ;
   h_meth_names -> Write() ;

   tf_out -> Close() ;

}







