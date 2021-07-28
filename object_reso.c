#define object_reso_cxx
#include "object_reso.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

float calc_dr( double phi1, double phi2, double eta1, double eta2 ) ;

void ptep_to_xyze( float pt, double eta, double phi, double m,
                   float& px, float& py, float& pz, float& e ) ;

void object_reso::Loop( bool verbose, int maxEvt )
{

   TRandom2 tran(12345) ;


   bool includeExtraTTreeVars = true ;

   if (fChain == 0) return;

   TFile* tf_out = new TFile("mini-tree.root", "recreate") ;
   TTree* tt_out = new TTree("minitree", "Minimal flat TTree for hadronic Q2,x,y analysis" ) ;



   float gen_pt ;
   float gen_eta ;
   float rec_pt ;
   float rec_eta ;
   float match_dr ;

   tt_out -> Branch( "gen_pt", &gen_pt, "gen_pt/F" ) ;
   tt_out -> Branch( "gen_eta", &gen_eta, "gen_eta/F" ) ;
   tt_out -> Branch( "rec_pt", &rec_pt, "rec_pt/F" ) ;
   tt_out -> Branch( "rec_eta", &rec_eta, "rec_eta/F" ) ;
   tt_out -> Branch( "match_dr", &match_dr, "match_dr/F" ) ;



   gDirectory -> Delete( "h*" ) ;


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



    //-----------------------------

///   for ( int i = 0; i < EFlowTrack_; i++ ) {


///      if ( verbose ) printf( " %3d, trk %3d : pid = %5d  pt = %9.4f  eta = %9.5f, phi = %9.5f\n",
///        ei, i, EFlowTrack_PID[i], EFlowTrack_PT[i], EFlowTrack_Eta[i], EFlowTrack_Phi[i] ) ;

///      double min_mc_dr = 99999. ;
///      int mc_pi = -1 ;
///      for ( int pi = 0; pi < Particle_ ; pi ++ ) {
///         if  ( Particle_Status[pi] != 1 ) continue ;
///         if ( fabs(Particle_Charge[pi]) != 1 ) continue ;
///         float dr = calc_dr( EFlowTrack_Phi[i], Particle_Phi[pi], EFlowTrack_Eta[i], Particle_Eta[pi] ) ;
///         if ( dr < min_mc_dr ) {
///            min_mc_dr = dr ;
///            mc_pi = pi ;
///         }
///      } // pi
///      if ( mc_pi < 0 ) continue ;
///      if ( min_mc_dr > 0.05 ) continue ;

///      if ( verbose ) {
///         float ptratio = Track_PT[i] / Particle_PT[mc_pi] ;
///         printf( " %3d, trk %3d :  MC match    pt = %9.4f  eta = %9.5f, phi = %9.5f   dr = %11.7f   rec pt / gen pt = %8.4f\n",
///            ei, i, Particle_PT[mc_pi], Particle_Eta[mc_pi], Particle_Phi[mc_pi], min_mc_dr, ptratio ) ;
///         printf("\n") ;
///      }

///      rec_pt = EFlowTrack_PT[i] ;
///      rec_eta = EFlowTrack_Eta[i] ;
///      gen_pt = Particle_PT[mc_pi] ;
///      gen_eta = Particle_Eta[mc_pi] ;
///      match_dr = min_mc_dr ;

///      tt_out -> Fill() ;

///   } // i

    //-----------------------------

      for ( int i = 0; i < Track_; i++ ) {


         if ( verbose ) printf( " %3d, trk %3d : pid = %5d  pt = %9.4f  eta = %9.5f, phi = %9.5f\n",
           ei, i, Track_PID[i], Track_PT[i], Track_Eta[i], Track_Phi[i] ) ;

         double min_mc_dr = 99999. ;
         int mc_pi = -1 ;
         for ( int pi = 0; pi < Particle_ ; pi ++ ) {
            if  ( Particle_Status[pi] != 1 ) continue ;
            if ( fabs(Particle_Charge[pi]) != 1 ) continue ;
            float dr = calc_dr( Track_Phi[i], Particle_Phi[pi], Track_Eta[i], Particle_Eta[pi] ) ;
            if ( dr < min_mc_dr ) {
               min_mc_dr = dr ;
               mc_pi = pi ;
            }
         } // pi
         if ( mc_pi < 0 ) continue ;
         if ( min_mc_dr > 0.05 ) continue ;

         if ( verbose ) {
            float ptratio = Track_PT[i] / Particle_PT[mc_pi] ;
            printf( " %3d, trk %3d :  MC match    pt = %9.4f  eta = %9.5f, phi = %9.5f   dr = %11.7f   rec pt / gen pt = %8.4f\n",
               ei, i, Particle_PT[mc_pi], Particle_Eta[mc_pi], Particle_Phi[mc_pi], min_mc_dr, ptratio ) ;
            printf("\n") ;
         }

         rec_pt = Track_PT[i] ;
         rec_eta = Track_Eta[i] ;
         gen_pt = Particle_PT[mc_pi] ;
         gen_eta = Particle_Eta[mc_pi] ;
         match_dr = min_mc_dr ;

         tt_out -> Fill() ;

      } // i


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






