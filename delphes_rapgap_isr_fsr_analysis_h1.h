//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep  2 07:53:17 2021 by ROOT version 6.24/01
// from TTree Delphes/Analysis tree
// found on file: /Volumes/Ext_2020_08/dis-reco-work/work-2021-09-02/athena-rapgap/athena-rapgap001.root
//////////////////////////////////////////////////////////

#ifndef delphes_rapgap_isr_fsr_analysis_h1_h
#define delphes_rapgap_isr_fsr_analysis_h1_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"

class delphes_rapgap_isr_fsr_analysis_h1 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
// static constexpr Int_t kMaxEvent = 1;
// static constexpr Int_t kMaxWeight = 1;
// static constexpr Int_t kMaxParticle = 155;
// static constexpr Int_t kMaxBeamSpot = 1;
// static constexpr Int_t kMaxTrack = 1;
// static constexpr Int_t kMaxTower = 2;
// static constexpr Int_t kMaxEFlowTrack = 1;
// static constexpr Int_t kMaxEFlowPhoton = 2;
// static constexpr Int_t kMaxEFlowNeutralHadron = 1;
// static constexpr Int_t kMaxmRICHTrack = 1;
// static constexpr Int_t kMaxbarrelDIRCTrack = 1;
// static constexpr Int_t kMaxdualRICHagTrack = 1;
// static constexpr Int_t kMaxdualRICHcfTrack = 1;
// static constexpr Int_t kMaxGenJet = 1;
// static constexpr Int_t kMaxGenMissingET = 1;
// static constexpr Int_t kMaxJet = 1;
// static constexpr Int_t kMaxElectron = 1;
// static constexpr Int_t kMaxPhoton = 2;
// static constexpr Int_t kMaxMissingET = 1;
// static constexpr Int_t kMaxScalarHT = 1;

   static constexpr Int_t kMaxEvent = 10;
   static constexpr Int_t kMaxWeight = 10;
   static constexpr Int_t kMaxParticle = 555;
   static constexpr Int_t kMaxBeamSpot = 10;
   static constexpr Int_t kMaxTrack = 10;
   static constexpr Int_t kMaxTower = 12;
   static constexpr Int_t kMaxEFlowTrack = 10;
   static constexpr Int_t kMaxEFlowPhoton = 12;
   static constexpr Int_t kMaxEFlowNeutralHadron = 10;
   static constexpr Int_t kMaxmRICHTrack = 10;
   static constexpr Int_t kMaxbarrelDIRCTrack = 10;
   static constexpr Int_t kMaxdualRICHagTrack = 10;
   static constexpr Int_t kMaxdualRICHcfTrack = 10;
   static constexpr Int_t kMaxGenJet = 10;
   static constexpr Int_t kMaxGenMissingET = 10;
   static constexpr Int_t kMaxJet = 10;
   static constexpr Int_t kMaxElectron = 10;
   static constexpr Int_t kMaxPhoton = 52;
   static constexpr Int_t kMaxMissingET = 10;
   static constexpr Int_t kMaxScalarHT = 10;

   // Declaration of leaf types
   Int_t           Event_;
   UInt_t          Event_fUniqueID[kMaxEvent];   //[Event_]
   UInt_t          Event_fBits[kMaxEvent];   //[Event_]
   Long64_t        Event_Number[kMaxEvent];   //[Event_]
   Float_t         Event_ReadTime[kMaxEvent];   //[Event_]
   Float_t         Event_ProcTime[kMaxEvent];   //[Event_]
   Int_t           Event_ProcessID[kMaxEvent];   //[Event_]
   Int_t           Event_MPI[kMaxEvent];   //[Event_]
   Float_t         Event_Weight[kMaxEvent];   //[Event_]
   Float_t         Event_CrossSection[kMaxEvent];   //[Event_]
   Float_t         Event_CrossSectionError[kMaxEvent];   //[Event_]
   Float_t         Event_Scale[kMaxEvent];   //[Event_]
   Float_t         Event_AlphaQED[kMaxEvent];   //[Event_]
   Float_t         Event_AlphaQCD[kMaxEvent];   //[Event_]
   Int_t           Event_ID1[kMaxEvent];   //[Event_]
   Int_t           Event_ID2[kMaxEvent];   //[Event_]
   Float_t         Event_X1[kMaxEvent];   //[Event_]
   Float_t         Event_X2[kMaxEvent];   //[Event_]
   Float_t         Event_ScalePDF[kMaxEvent];   //[Event_]
   Float_t         Event_PDF1[kMaxEvent];   //[Event_]
   Float_t         Event_PDF2[kMaxEvent];   //[Event_]
   Int_t           Event_size;
   Int_t           Weight_;
   UInt_t          Weight_fUniqueID[kMaxWeight];   //[Weight_]
   UInt_t          Weight_fBits[kMaxWeight];   //[Weight_]
   Float_t         Weight_Weight[kMaxWeight];   //[Weight_]
   Int_t           Weight_size;
   Int_t           Particle_;
   UInt_t          Particle_fUniqueID[kMaxParticle];   //[Particle_]
   UInt_t          Particle_fBits[kMaxParticle];   //[Particle_]
   Int_t           Particle_PID[kMaxParticle];   //[Particle_]
   Int_t           Particle_Status[kMaxParticle];   //[Particle_]
   Int_t           Particle_IsPU[kMaxParticle];   //[Particle_]
   Int_t           Particle_M1[kMaxParticle];   //[Particle_]
   Int_t           Particle_M2[kMaxParticle];   //[Particle_]
   Int_t           Particle_D1[kMaxParticle];   //[Particle_]
   Int_t           Particle_D2[kMaxParticle];   //[Particle_]
   Int_t           Particle_Charge[kMaxParticle];   //[Particle_]
   Float_t         Particle_Mass[kMaxParticle];   //[Particle_]
   Float_t         Particle_E[kMaxParticle];   //[Particle_]
   Float_t         Particle_Px[kMaxParticle];   //[Particle_]
   Float_t         Particle_Py[kMaxParticle];   //[Particle_]
   Float_t         Particle_Pz[kMaxParticle];   //[Particle_]
   Float_t         Particle_P[kMaxParticle];   //[Particle_]
   Float_t         Particle_PT[kMaxParticle];   //[Particle_]
   Float_t         Particle_Eta[kMaxParticle];   //[Particle_]
   Float_t         Particle_Phi[kMaxParticle];   //[Particle_]
   Float_t         Particle_Rapidity[kMaxParticle];   //[Particle_]
   Float_t         Particle_T[kMaxParticle];   //[Particle_]
   Float_t         Particle_X[kMaxParticle];   //[Particle_]
   Float_t         Particle_Y[kMaxParticle];   //[Particle_]
   Float_t         Particle_Z[kMaxParticle];   //[Particle_]
   Int_t           Particle_size;
   Int_t           BeamSpot_;
   UInt_t          BeamSpot_fUniqueID[kMaxBeamSpot];   //[BeamSpot_]
   UInt_t          BeamSpot_fBits[kMaxBeamSpot];   //[BeamSpot_]
   Int_t           BeamSpot_PID[kMaxBeamSpot];   //[BeamSpot_]
   Int_t           BeamSpot_Status[kMaxBeamSpot];   //[BeamSpot_]
   Int_t           BeamSpot_IsPU[kMaxBeamSpot];   //[BeamSpot_]
   Int_t           BeamSpot_M1[kMaxBeamSpot];   //[BeamSpot_]
   Int_t           BeamSpot_M2[kMaxBeamSpot];   //[BeamSpot_]
   Int_t           BeamSpot_D1[kMaxBeamSpot];   //[BeamSpot_]
   Int_t           BeamSpot_D2[kMaxBeamSpot];   //[BeamSpot_]
   Int_t           BeamSpot_Charge[kMaxBeamSpot];   //[BeamSpot_]
   Float_t         BeamSpot_Mass[kMaxBeamSpot];   //[BeamSpot_]
   Float_t         BeamSpot_E[kMaxBeamSpot];   //[BeamSpot_]
   Float_t         BeamSpot_Px[kMaxBeamSpot];   //[BeamSpot_]
   Float_t         BeamSpot_Py[kMaxBeamSpot];   //[BeamSpot_]
   Float_t         BeamSpot_Pz[kMaxBeamSpot];   //[BeamSpot_]
   Float_t         BeamSpot_P[kMaxBeamSpot];   //[BeamSpot_]
   Float_t         BeamSpot_PT[kMaxBeamSpot];   //[BeamSpot_]
   Float_t         BeamSpot_Eta[kMaxBeamSpot];   //[BeamSpot_]
   Float_t         BeamSpot_Phi[kMaxBeamSpot];   //[BeamSpot_]
   Float_t         BeamSpot_Rapidity[kMaxBeamSpot];   //[BeamSpot_]
   Float_t         BeamSpot_T[kMaxBeamSpot];   //[BeamSpot_]
   Float_t         BeamSpot_X[kMaxBeamSpot];   //[BeamSpot_]
   Float_t         BeamSpot_Y[kMaxBeamSpot];   //[BeamSpot_]
   Float_t         BeamSpot_Z[kMaxBeamSpot];   //[BeamSpot_]
   Int_t           BeamSpot_size;
   Int_t           Track_;
   UInt_t          Track_fUniqueID[kMaxTrack];   //[Track_]
   UInt_t          Track_fBits[kMaxTrack];   //[Track_]
   Int_t           Track_PID[kMaxTrack];   //[Track_]
   Int_t           Track_Charge[kMaxTrack];   //[Track_]
   Float_t         Track_P[kMaxTrack];   //[Track_]
   Float_t         Track_PT[kMaxTrack];   //[Track_]
   Float_t         Track_Eta[kMaxTrack];   //[Track_]
   Float_t         Track_Phi[kMaxTrack];   //[Track_]
   Float_t         Track_CtgTheta[kMaxTrack];   //[Track_]
   Float_t         Track_C[kMaxTrack];   //[Track_]
   Float_t         Track_Mass[kMaxTrack];   //[Track_]
   Float_t         Track_EtaOuter[kMaxTrack];   //[Track_]
   Float_t         Track_PhiOuter[kMaxTrack];   //[Track_]
   Float_t         Track_T[kMaxTrack];   //[Track_]
   Float_t         Track_X[kMaxTrack];   //[Track_]
   Float_t         Track_Y[kMaxTrack];   //[Track_]
   Float_t         Track_Z[kMaxTrack];   //[Track_]
   Float_t         Track_TOuter[kMaxTrack];   //[Track_]
   Float_t         Track_XOuter[kMaxTrack];   //[Track_]
   Float_t         Track_YOuter[kMaxTrack];   //[Track_]
   Float_t         Track_ZOuter[kMaxTrack];   //[Track_]
   Float_t         Track_Xd[kMaxTrack];   //[Track_]
   Float_t         Track_Yd[kMaxTrack];   //[Track_]
   Float_t         Track_Zd[kMaxTrack];   //[Track_]
   Float_t         Track_L[kMaxTrack];   //[Track_]
   Float_t         Track_D0[kMaxTrack];   //[Track_]
   Float_t         Track_DZ[kMaxTrack];   //[Track_]
   Float_t         Track_Nclusters[kMaxTrack];   //[Track_]
   Float_t         Track_dNdx[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorP[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorPT[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorPhi[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorCtgTheta[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorT[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorD0[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorDZ[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorC[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorD0Phi[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorD0C[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorD0DZ[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorD0CtgTheta[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorPhiC[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorPhiDZ[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorPhiCtgTheta[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorCDZ[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorCCtgTheta[kMaxTrack];   //[Track_]
   Float_t         Track_ErrorDZCtgTheta[kMaxTrack];   //[Track_]
   TRef            Track_Particle[kMaxTrack];
   Int_t           Track_VertexIndex[kMaxTrack];   //[Track_]
   Int_t           Track_size;
   Int_t           Tower_;
   UInt_t          Tower_fUniqueID[kMaxTower];   //[Tower_]
   UInt_t          Tower_fBits[kMaxTower];   //[Tower_]
   Float_t         Tower_ET[kMaxTower];   //[Tower_]
   Float_t         Tower_Eta[kMaxTower];   //[Tower_]
   Float_t         Tower_Phi[kMaxTower];   //[Tower_]
   Float_t         Tower_E[kMaxTower];   //[Tower_]
   Float_t         Tower_T[kMaxTower];   //[Tower_]
   Int_t           Tower_NTimeHits[kMaxTower];   //[Tower_]
   Float_t         Tower_Eem[kMaxTower];   //[Tower_]
   Float_t         Tower_Ehad[kMaxTower];   //[Tower_]
   Float_t         Tower_Edges[kMaxTower][4];   //[Tower_]
   TRefArray       Tower_Particles[kMaxTower];
   Int_t           Tower_size;
   Int_t           EFlowTrack_;
   UInt_t          EFlowTrack_fUniqueID[kMaxEFlowTrack];   //[EFlowTrack_]
   UInt_t          EFlowTrack_fBits[kMaxEFlowTrack];   //[EFlowTrack_]
   Int_t           EFlowTrack_PID[kMaxEFlowTrack];   //[EFlowTrack_]
   Int_t           EFlowTrack_Charge[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_P[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_PT[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Eta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Phi[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_CtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_C[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Mass[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_EtaOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_PhiOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_T[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_X[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Y[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Z[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_TOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_XOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_YOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ZOuter[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Xd[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Yd[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Zd[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_L[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_D0[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_DZ[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_Nclusters[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_dNdx[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorP[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorPT[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorPhi[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorCtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorT[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorD0[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorDZ[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorC[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorD0Phi[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorD0C[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorD0DZ[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorD0CtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorPhiC[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorPhiDZ[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorPhiCtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorCDZ[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorCCtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   Float_t         EFlowTrack_ErrorDZCtgTheta[kMaxEFlowTrack];   //[EFlowTrack_]
   TRef            EFlowTrack_Particle[kMaxEFlowTrack];
   Int_t           EFlowTrack_VertexIndex[kMaxEFlowTrack];   //[EFlowTrack_]
   Int_t           EFlowTrack_size;
   Int_t           EFlowPhoton_;
   UInt_t          EFlowPhoton_fUniqueID[kMaxEFlowPhoton];   //[EFlowPhoton_]
   UInt_t          EFlowPhoton_fBits[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Float_t         EFlowPhoton_ET[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Float_t         EFlowPhoton_Eta[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Float_t         EFlowPhoton_Phi[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Float_t         EFlowPhoton_E[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Float_t         EFlowPhoton_T[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Int_t           EFlowPhoton_NTimeHits[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Float_t         EFlowPhoton_Eem[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Float_t         EFlowPhoton_Ehad[kMaxEFlowPhoton];   //[EFlowPhoton_]
   Float_t         EFlowPhoton_Edges[kMaxEFlowPhoton][4];   //[EFlowPhoton_]
   TRefArray       EFlowPhoton_Particles[kMaxEFlowPhoton];
   Int_t           EFlowPhoton_size;
   Int_t           EFlowNeutralHadron_;
   UInt_t          EFlowNeutralHadron_fUniqueID[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   UInt_t          EFlowNeutralHadron_fBits[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Float_t         EFlowNeutralHadron_ET[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Float_t         EFlowNeutralHadron_Eta[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Float_t         EFlowNeutralHadron_Phi[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Float_t         EFlowNeutralHadron_E[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Float_t         EFlowNeutralHadron_T[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Int_t           EFlowNeutralHadron_NTimeHits[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Float_t         EFlowNeutralHadron_Eem[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Float_t         EFlowNeutralHadron_Ehad[kMaxEFlowNeutralHadron];   //[EFlowNeutralHadron_]
   Float_t         EFlowNeutralHadron_Edges[kMaxEFlowNeutralHadron][4];   //[EFlowNeutralHadron_]
   TRefArray       EFlowNeutralHadron_Particles[kMaxEFlowNeutralHadron];
   Int_t           EFlowNeutralHadron_size;
   Int_t           mRICHTrack_;
   UInt_t          mRICHTrack_fUniqueID[kMaxmRICHTrack];   //[mRICHTrack_]
   UInt_t          mRICHTrack_fBits[kMaxmRICHTrack];   //[mRICHTrack_]
   Int_t           mRICHTrack_PID[kMaxmRICHTrack];   //[mRICHTrack_]
   Int_t           mRICHTrack_Charge[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_P[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_PT[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_Eta[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_Phi[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_CtgTheta[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_C[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_Mass[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_EtaOuter[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_PhiOuter[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_T[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_X[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_Y[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_Z[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_TOuter[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_XOuter[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_YOuter[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ZOuter[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_Xd[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_Yd[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_Zd[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_L[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_D0[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_DZ[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_Nclusters[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_dNdx[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorP[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorPT[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorPhi[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorCtgTheta[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorT[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorD0[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorDZ[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorC[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorD0Phi[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorD0C[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorD0DZ[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorD0CtgTheta[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorPhiC[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorPhiDZ[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorPhiCtgTheta[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorCDZ[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorCCtgTheta[kMaxmRICHTrack];   //[mRICHTrack_]
   Float_t         mRICHTrack_ErrorDZCtgTheta[kMaxmRICHTrack];   //[mRICHTrack_]
   TRef            mRICHTrack_Particle[kMaxmRICHTrack];
   Int_t           mRICHTrack_VertexIndex[kMaxmRICHTrack];   //[mRICHTrack_]
   Int_t           mRICHTrack_size;
   Int_t           barrelDIRCTrack_;
   UInt_t          barrelDIRCTrack_fUniqueID[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   UInt_t          barrelDIRCTrack_fBits[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Int_t           barrelDIRCTrack_PID[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Int_t           barrelDIRCTrack_Charge[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_P[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_PT[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_Eta[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_Phi[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_CtgTheta[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_C[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_Mass[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_EtaOuter[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_PhiOuter[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_T[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_X[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_Y[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_Z[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_TOuter[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_XOuter[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_YOuter[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ZOuter[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_Xd[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_Yd[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_Zd[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_L[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_D0[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_DZ[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_Nclusters[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_dNdx[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorP[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorPT[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorPhi[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorCtgTheta[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorT[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorD0[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorDZ[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorC[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorD0Phi[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorD0C[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorD0DZ[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorD0CtgTheta[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorPhiC[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorPhiDZ[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorPhiCtgTheta[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorCDZ[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorCCtgTheta[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Float_t         barrelDIRCTrack_ErrorDZCtgTheta[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   TRef            barrelDIRCTrack_Particle[kMaxbarrelDIRCTrack];
   Int_t           barrelDIRCTrack_VertexIndex[kMaxbarrelDIRCTrack];   //[barrelDIRCTrack_]
   Int_t           barrelDIRCTrack_size;
   Int_t           dualRICHagTrack_;
   UInt_t          dualRICHagTrack_fUniqueID[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   UInt_t          dualRICHagTrack_fBits[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Int_t           dualRICHagTrack_PID[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Int_t           dualRICHagTrack_Charge[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_P[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_PT[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_Eta[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_Phi[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_CtgTheta[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_C[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_Mass[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_EtaOuter[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_PhiOuter[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_T[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_X[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_Y[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_Z[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_TOuter[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_XOuter[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_YOuter[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ZOuter[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_Xd[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_Yd[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_Zd[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_L[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_D0[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_DZ[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_Nclusters[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_dNdx[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorP[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorPT[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorPhi[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorCtgTheta[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorT[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorD0[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorDZ[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorC[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorD0Phi[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorD0C[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorD0DZ[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorD0CtgTheta[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorPhiC[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorPhiDZ[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorPhiCtgTheta[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorCDZ[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorCCtgTheta[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Float_t         dualRICHagTrack_ErrorDZCtgTheta[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   TRef            dualRICHagTrack_Particle[kMaxdualRICHagTrack];
   Int_t           dualRICHagTrack_VertexIndex[kMaxdualRICHagTrack];   //[dualRICHagTrack_]
   Int_t           dualRICHagTrack_size;
   Int_t           dualRICHcfTrack_;
   UInt_t          dualRICHcfTrack_fUniqueID[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   UInt_t          dualRICHcfTrack_fBits[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Int_t           dualRICHcfTrack_PID[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Int_t           dualRICHcfTrack_Charge[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_P[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_PT[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_Eta[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_Phi[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_CtgTheta[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_C[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_Mass[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_EtaOuter[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_PhiOuter[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_T[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_X[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_Y[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_Z[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_TOuter[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_XOuter[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_YOuter[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ZOuter[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_Xd[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_Yd[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_Zd[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_L[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_D0[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_DZ[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_Nclusters[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_dNdx[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorP[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorPT[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorPhi[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorCtgTheta[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorT[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorD0[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorDZ[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorC[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorD0Phi[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorD0C[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorD0DZ[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorD0CtgTheta[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorPhiC[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorPhiDZ[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorPhiCtgTheta[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorCDZ[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorCCtgTheta[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Float_t         dualRICHcfTrack_ErrorDZCtgTheta[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   TRef            dualRICHcfTrack_Particle[kMaxdualRICHcfTrack];
   Int_t           dualRICHcfTrack_VertexIndex[kMaxdualRICHcfTrack];   //[dualRICHcfTrack_]
   Int_t           dualRICHcfTrack_size;
   Int_t           GenJet_;
   UInt_t          GenJet_fUniqueID[kMaxGenJet];   //[GenJet_]
   UInt_t          GenJet_fBits[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_PT[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_Eta[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_Phi[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_T[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_Mass[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_DeltaEta[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_DeltaPhi[kMaxGenJet];   //[GenJet_]
   UInt_t          GenJet_Flavor[kMaxGenJet];   //[GenJet_]
   UInt_t          GenJet_FlavorAlgo[kMaxGenJet];   //[GenJet_]
   UInt_t          GenJet_FlavorPhys[kMaxGenJet];   //[GenJet_]
   UInt_t          GenJet_BTag[kMaxGenJet];   //[GenJet_]
   UInt_t          GenJet_BTagAlgo[kMaxGenJet];   //[GenJet_]
   UInt_t          GenJet_BTagPhys[kMaxGenJet];   //[GenJet_]
   UInt_t          GenJet_TauTag[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_TauWeight[kMaxGenJet];   //[GenJet_]
   Int_t           GenJet_Charge[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_EhadOverEem[kMaxGenJet];   //[GenJet_]
   Int_t           GenJet_NCharged[kMaxGenJet];   //[GenJet_]
   Int_t           GenJet_NNeutrals[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_NeutralEnergyFraction[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_ChargedEnergyFraction[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_Beta[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_BetaStar[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_MeanSqDeltaR[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_PTD[kMaxGenJet];   //[GenJet_]
   Float_t         GenJet_FracPt[kMaxGenJet][5];   //[GenJet_]
   Float_t         GenJet_Tau[kMaxGenJet][5];   //[GenJet_]
   TLorentzVector  GenJet_SoftDroppedJet[kMaxGenJet];
   TLorentzVector  GenJet_SoftDroppedSubJet1[kMaxGenJet];
   TLorentzVector  GenJet_SoftDroppedSubJet2[kMaxGenJet];
   TLorentzVector  GenJet_TrimmedP4[5][kMaxGenJet];
   TLorentzVector  GenJet_PrunedP4[5][kMaxGenJet];
   TLorentzVector  GenJet_SoftDroppedP4[5][kMaxGenJet];
   Int_t           GenJet_NSubJetsTrimmed[kMaxGenJet];   //[GenJet_]
   Int_t           GenJet_NSubJetsPruned[kMaxGenJet];   //[GenJet_]
   Int_t           GenJet_NSubJetsSoftDropped[kMaxGenJet];   //[GenJet_]
   Double_t        GenJet_ExclYmerge23[kMaxGenJet];   //[GenJet_]
   Double_t        GenJet_ExclYmerge34[kMaxGenJet];   //[GenJet_]
   Double_t        GenJet_ExclYmerge45[kMaxGenJet];   //[GenJet_]
   Double_t        GenJet_ExclYmerge56[kMaxGenJet];   //[GenJet_]
   TRefArray       GenJet_Constituents[kMaxGenJet];
   TRefArray       GenJet_Particles[kMaxGenJet];
   TLorentzVector  GenJet_Area[kMaxGenJet];
   Int_t           GenJet_size;
   Int_t           GenMissingET_;
   UInt_t          GenMissingET_fUniqueID[kMaxGenMissingET];   //[GenMissingET_]
   UInt_t          GenMissingET_fBits[kMaxGenMissingET];   //[GenMissingET_]
   Float_t         GenMissingET_MET[kMaxGenMissingET];   //[GenMissingET_]
   Float_t         GenMissingET_Eta[kMaxGenMissingET];   //[GenMissingET_]
   Float_t         GenMissingET_Phi[kMaxGenMissingET];   //[GenMissingET_]
   Int_t           GenMissingET_size;
   Int_t           Jet_;
   UInt_t          Jet_fUniqueID[kMaxJet];   //[Jet_]
   UInt_t          Jet_fBits[kMaxJet];   //[Jet_]
   Float_t         Jet_PT[kMaxJet];   //[Jet_]
   Float_t         Jet_Eta[kMaxJet];   //[Jet_]
   Float_t         Jet_Phi[kMaxJet];   //[Jet_]
   Float_t         Jet_T[kMaxJet];   //[Jet_]
   Float_t         Jet_Mass[kMaxJet];   //[Jet_]
   Float_t         Jet_DeltaEta[kMaxJet];   //[Jet_]
   Float_t         Jet_DeltaPhi[kMaxJet];   //[Jet_]
   UInt_t          Jet_Flavor[kMaxJet];   //[Jet_]
   UInt_t          Jet_FlavorAlgo[kMaxJet];   //[Jet_]
   UInt_t          Jet_FlavorPhys[kMaxJet];   //[Jet_]
   UInt_t          Jet_BTag[kMaxJet];   //[Jet_]
   UInt_t          Jet_BTagAlgo[kMaxJet];   //[Jet_]
   UInt_t          Jet_BTagPhys[kMaxJet];   //[Jet_]
   UInt_t          Jet_TauTag[kMaxJet];   //[Jet_]
   Float_t         Jet_TauWeight[kMaxJet];   //[Jet_]
   Int_t           Jet_Charge[kMaxJet];   //[Jet_]
   Float_t         Jet_EhadOverEem[kMaxJet];   //[Jet_]
   Int_t           Jet_NCharged[kMaxJet];   //[Jet_]
   Int_t           Jet_NNeutrals[kMaxJet];   //[Jet_]
   Float_t         Jet_NeutralEnergyFraction[kMaxJet];   //[Jet_]
   Float_t         Jet_ChargedEnergyFraction[kMaxJet];   //[Jet_]
   Float_t         Jet_Beta[kMaxJet];   //[Jet_]
   Float_t         Jet_BetaStar[kMaxJet];   //[Jet_]
   Float_t         Jet_MeanSqDeltaR[kMaxJet];   //[Jet_]
   Float_t         Jet_PTD[kMaxJet];   //[Jet_]
   Float_t         Jet_FracPt[kMaxJet][5];   //[Jet_]
   Float_t         Jet_Tau[kMaxJet][5];   //[Jet_]
   TLorentzVector  Jet_SoftDroppedJet[kMaxJet];
   TLorentzVector  Jet_SoftDroppedSubJet1[kMaxJet];
   TLorentzVector  Jet_SoftDroppedSubJet2[kMaxJet];
   TLorentzVector  Jet_TrimmedP4[5][kMaxJet];
   TLorentzVector  Jet_PrunedP4[5][kMaxJet];
   TLorentzVector  Jet_SoftDroppedP4[5][kMaxJet];
   Int_t           Jet_NSubJetsTrimmed[kMaxJet];   //[Jet_]
   Int_t           Jet_NSubJetsPruned[kMaxJet];   //[Jet_]
   Int_t           Jet_NSubJetsSoftDropped[kMaxJet];   //[Jet_]
   Double_t        Jet_ExclYmerge23[kMaxJet];   //[Jet_]
   Double_t        Jet_ExclYmerge34[kMaxJet];   //[Jet_]
   Double_t        Jet_ExclYmerge45[kMaxJet];   //[Jet_]
   Double_t        Jet_ExclYmerge56[kMaxJet];   //[Jet_]
   TRefArray       Jet_Constituents[kMaxJet];
   TRefArray       Jet_Particles[kMaxJet];
   TLorentzVector  Jet_Area[kMaxJet];
   Int_t           Jet_size;
   Int_t           Electron_;
   UInt_t          Electron_fUniqueID[kMaxElectron];   //[Electron_]
   UInt_t          Electron_fBits[kMaxElectron];   //[Electron_]
   Float_t         Electron_PT[kMaxElectron];   //[Electron_]
   Float_t         Electron_Eta[kMaxElectron];   //[Electron_]
   Float_t         Electron_Phi[kMaxElectron];   //[Electron_]
   Float_t         Electron_T[kMaxElectron];   //[Electron_]
   Int_t           Electron_Charge[kMaxElectron];   //[Electron_]
   Float_t         Electron_EhadOverEem[kMaxElectron];   //[Electron_]
   TRef            Electron_Particle[kMaxElectron];
   Float_t         Electron_IsolationVar[kMaxElectron];   //[Electron_]
   Float_t         Electron_IsolationVarRhoCorr[kMaxElectron];   //[Electron_]
   Float_t         Electron_SumPtCharged[kMaxElectron];   //[Electron_]
   Float_t         Electron_SumPtNeutral[kMaxElectron];   //[Electron_]
   Float_t         Electron_SumPtChargedPU[kMaxElectron];   //[Electron_]
   Float_t         Electron_SumPt[kMaxElectron];   //[Electron_]
   Float_t         Electron_D0[kMaxElectron];   //[Electron_]
   Float_t         Electron_DZ[kMaxElectron];   //[Electron_]
   Float_t         Electron_ErrorD0[kMaxElectron];   //[Electron_]
   Float_t         Electron_ErrorDZ[kMaxElectron];   //[Electron_]
   Int_t           Electron_size;
   Int_t           Photon_;
   UInt_t          Photon_fUniqueID[kMaxPhoton];   //[Photon_]
   UInt_t          Photon_fBits[kMaxPhoton];   //[Photon_]
   Float_t         Photon_PT[kMaxPhoton];   //[Photon_]
   Float_t         Photon_Eta[kMaxPhoton];   //[Photon_]
   Float_t         Photon_Phi[kMaxPhoton];   //[Photon_]
   Float_t         Photon_E[kMaxPhoton];   //[Photon_]
   Float_t         Photon_T[kMaxPhoton];   //[Photon_]
   Float_t         Photon_EhadOverEem[kMaxPhoton];   //[Photon_]
   TRefArray       Photon_Particles[kMaxPhoton];
   Float_t         Photon_IsolationVar[kMaxPhoton];   //[Photon_]
   Float_t         Photon_IsolationVarRhoCorr[kMaxPhoton];   //[Photon_]
   Float_t         Photon_SumPtCharged[kMaxPhoton];   //[Photon_]
   Float_t         Photon_SumPtNeutral[kMaxPhoton];   //[Photon_]
   Float_t         Photon_SumPtChargedPU[kMaxPhoton];   //[Photon_]
   Float_t         Photon_SumPt[kMaxPhoton];   //[Photon_]
   Int_t           Photon_Status[kMaxPhoton];   //[Photon_]
   Int_t           Photon_size;
   Int_t           MissingET_;
   UInt_t          MissingET_fUniqueID[kMaxMissingET];   //[MissingET_]
   UInt_t          MissingET_fBits[kMaxMissingET];   //[MissingET_]
   Float_t         MissingET_MET[kMaxMissingET];   //[MissingET_]
   Float_t         MissingET_Eta[kMaxMissingET];   //[MissingET_]
   Float_t         MissingET_Phi[kMaxMissingET];   //[MissingET_]
   Int_t           MissingET_size;
   Int_t           ScalarHT_;
   UInt_t          ScalarHT_fUniqueID[kMaxScalarHT];   //[ScalarHT_]
   UInt_t          ScalarHT_fBits[kMaxScalarHT];   //[ScalarHT_]
   Float_t         ScalarHT_HT[kMaxScalarHT];   //[ScalarHT_]
   Int_t           ScalarHT_size;

   // List of branches
   TBranch        *b_Event_;   //!
   TBranch        *b_Event_fUniqueID;   //!
   TBranch        *b_Event_fBits;   //!
   TBranch        *b_Event_Number;   //!
   TBranch        *b_Event_ReadTime;   //!
   TBranch        *b_Event_ProcTime;   //!
   TBranch        *b_Event_ProcessID;   //!
   TBranch        *b_Event_MPI;   //!
   TBranch        *b_Event_Weight;   //!
   TBranch        *b_Event_CrossSection;   //!
   TBranch        *b_Event_CrossSectionError;   //!
   TBranch        *b_Event_Scale;   //!
   TBranch        *b_Event_AlphaQED;   //!
   TBranch        *b_Event_AlphaQCD;   //!
   TBranch        *b_Event_ID1;   //!
   TBranch        *b_Event_ID2;   //!
   TBranch        *b_Event_X1;   //!
   TBranch        *b_Event_X2;   //!
   TBranch        *b_Event_ScalePDF;   //!
   TBranch        *b_Event_PDF1;   //!
   TBranch        *b_Event_PDF2;   //!
   TBranch        *b_Event_size;   //!
   TBranch        *b_Weight_;   //!
   TBranch        *b_Weight_fUniqueID;   //!
   TBranch        *b_Weight_fBits;   //!
   TBranch        *b_Weight_Weight;   //!
   TBranch        *b_Weight_size;   //!
   TBranch        *b_Particle_;   //!
   TBranch        *b_Particle_fUniqueID;   //!
   TBranch        *b_Particle_fBits;   //!
   TBranch        *b_Particle_PID;   //!
   TBranch        *b_Particle_Status;   //!
   TBranch        *b_Particle_IsPU;   //!
   TBranch        *b_Particle_M1;   //!
   TBranch        *b_Particle_M2;   //!
   TBranch        *b_Particle_D1;   //!
   TBranch        *b_Particle_D2;   //!
   TBranch        *b_Particle_Charge;   //!
   TBranch        *b_Particle_Mass;   //!
   TBranch        *b_Particle_E;   //!
   TBranch        *b_Particle_Px;   //!
   TBranch        *b_Particle_Py;   //!
   TBranch        *b_Particle_Pz;   //!
   TBranch        *b_Particle_P;   //!
   TBranch        *b_Particle_PT;   //!
   TBranch        *b_Particle_Eta;   //!
   TBranch        *b_Particle_Phi;   //!
   TBranch        *b_Particle_Rapidity;   //!
   TBranch        *b_Particle_T;   //!
   TBranch        *b_Particle_X;   //!
   TBranch        *b_Particle_Y;   //!
   TBranch        *b_Particle_Z;   //!
   TBranch        *b_Particle_size;   //!
   TBranch        *b_BeamSpot_;   //!
   TBranch        *b_BeamSpot_fUniqueID;   //!
   TBranch        *b_BeamSpot_fBits;   //!
   TBranch        *b_BeamSpot_PID;   //!
   TBranch        *b_BeamSpot_Status;   //!
   TBranch        *b_BeamSpot_IsPU;   //!
   TBranch        *b_BeamSpot_M1;   //!
   TBranch        *b_BeamSpot_M2;   //!
   TBranch        *b_BeamSpot_D1;   //!
   TBranch        *b_BeamSpot_D2;   //!
   TBranch        *b_BeamSpot_Charge;   //!
   TBranch        *b_BeamSpot_Mass;   //!
   TBranch        *b_BeamSpot_E;   //!
   TBranch        *b_BeamSpot_Px;   //!
   TBranch        *b_BeamSpot_Py;   //!
   TBranch        *b_BeamSpot_Pz;   //!
   TBranch        *b_BeamSpot_P;   //!
   TBranch        *b_BeamSpot_PT;   //!
   TBranch        *b_BeamSpot_Eta;   //!
   TBranch        *b_BeamSpot_Phi;   //!
   TBranch        *b_BeamSpot_Rapidity;   //!
   TBranch        *b_BeamSpot_T;   //!
   TBranch        *b_BeamSpot_X;   //!
   TBranch        *b_BeamSpot_Y;   //!
   TBranch        *b_BeamSpot_Z;   //!
   TBranch        *b_BeamSpot_size;   //!
   TBranch        *b_Track_;   //!
   TBranch        *b_Track_fUniqueID;   //!
   TBranch        *b_Track_fBits;   //!
   TBranch        *b_Track_PID;   //!
   TBranch        *b_Track_Charge;   //!
   TBranch        *b_Track_P;   //!
   TBranch        *b_Track_PT;   //!
   TBranch        *b_Track_Eta;   //!
   TBranch        *b_Track_Phi;   //!
   TBranch        *b_Track_CtgTheta;   //!
   TBranch        *b_Track_C;   //!
   TBranch        *b_Track_Mass;   //!
   TBranch        *b_Track_EtaOuter;   //!
   TBranch        *b_Track_PhiOuter;   //!
   TBranch        *b_Track_T;   //!
   TBranch        *b_Track_X;   //!
   TBranch        *b_Track_Y;   //!
   TBranch        *b_Track_Z;   //!
   TBranch        *b_Track_TOuter;   //!
   TBranch        *b_Track_XOuter;   //!
   TBranch        *b_Track_YOuter;   //!
   TBranch        *b_Track_ZOuter;   //!
   TBranch        *b_Track_Xd;   //!
   TBranch        *b_Track_Yd;   //!
   TBranch        *b_Track_Zd;   //!
   TBranch        *b_Track_L;   //!
   TBranch        *b_Track_D0;   //!
   TBranch        *b_Track_DZ;   //!
   TBranch        *b_Track_Nclusters;   //!
   TBranch        *b_Track_dNdx;   //!
   TBranch        *b_Track_ErrorP;   //!
   TBranch        *b_Track_ErrorPT;   //!
   TBranch        *b_Track_ErrorPhi;   //!
   TBranch        *b_Track_ErrorCtgTheta;   //!
   TBranch        *b_Track_ErrorT;   //!
   TBranch        *b_Track_ErrorD0;   //!
   TBranch        *b_Track_ErrorDZ;   //!
   TBranch        *b_Track_ErrorC;   //!
   TBranch        *b_Track_ErrorD0Phi;   //!
   TBranch        *b_Track_ErrorD0C;   //!
   TBranch        *b_Track_ErrorD0DZ;   //!
   TBranch        *b_Track_ErrorD0CtgTheta;   //!
   TBranch        *b_Track_ErrorPhiC;   //!
   TBranch        *b_Track_ErrorPhiDZ;   //!
   TBranch        *b_Track_ErrorPhiCtgTheta;   //!
   TBranch        *b_Track_ErrorCDZ;   //!
   TBranch        *b_Track_ErrorCCtgTheta;   //!
   TBranch        *b_Track_ErrorDZCtgTheta;   //!
   TBranch        *b_Track_Particle;   //!
   TBranch        *b_Track_VertexIndex;   //!
   TBranch        *b_Track_size;   //!
   TBranch        *b_Tower_;   //!
   TBranch        *b_Tower_fUniqueID;   //!
   TBranch        *b_Tower_fBits;   //!
   TBranch        *b_Tower_ET;   //!
   TBranch        *b_Tower_Eta;   //!
   TBranch        *b_Tower_Phi;   //!
   TBranch        *b_Tower_E;   //!
   TBranch        *b_Tower_T;   //!
   TBranch        *b_Tower_NTimeHits;   //!
   TBranch        *b_Tower_Eem;   //!
   TBranch        *b_Tower_Ehad;   //!
   TBranch        *b_Tower_Edges;   //!
   TBranch        *b_Tower_Particles;   //!
   TBranch        *b_Tower_size;   //!
   TBranch        *b_EFlowTrack_;   //!
   TBranch        *b_EFlowTrack_fUniqueID;   //!
   TBranch        *b_EFlowTrack_fBits;   //!
   TBranch        *b_EFlowTrack_PID;   //!
   TBranch        *b_EFlowTrack_Charge;   //!
   TBranch        *b_EFlowTrack_P;   //!
   TBranch        *b_EFlowTrack_PT;   //!
   TBranch        *b_EFlowTrack_Eta;   //!
   TBranch        *b_EFlowTrack_Phi;   //!
   TBranch        *b_EFlowTrack_CtgTheta;   //!
   TBranch        *b_EFlowTrack_C;   //!
   TBranch        *b_EFlowTrack_Mass;   //!
   TBranch        *b_EFlowTrack_EtaOuter;   //!
   TBranch        *b_EFlowTrack_PhiOuter;   //!
   TBranch        *b_EFlowTrack_T;   //!
   TBranch        *b_EFlowTrack_X;   //!
   TBranch        *b_EFlowTrack_Y;   //!
   TBranch        *b_EFlowTrack_Z;   //!
   TBranch        *b_EFlowTrack_TOuter;   //!
   TBranch        *b_EFlowTrack_XOuter;   //!
   TBranch        *b_EFlowTrack_YOuter;   //!
   TBranch        *b_EFlowTrack_ZOuter;   //!
   TBranch        *b_EFlowTrack_Xd;   //!
   TBranch        *b_EFlowTrack_Yd;   //!
   TBranch        *b_EFlowTrack_Zd;   //!
   TBranch        *b_EFlowTrack_L;   //!
   TBranch        *b_EFlowTrack_D0;   //!
   TBranch        *b_EFlowTrack_DZ;   //!
   TBranch        *b_EFlowTrack_Nclusters;   //!
   TBranch        *b_EFlowTrack_dNdx;   //!
   TBranch        *b_EFlowTrack_ErrorP;   //!
   TBranch        *b_EFlowTrack_ErrorPT;   //!
   TBranch        *b_EFlowTrack_ErrorPhi;   //!
   TBranch        *b_EFlowTrack_ErrorCtgTheta;   //!
   TBranch        *b_EFlowTrack_ErrorT;   //!
   TBranch        *b_EFlowTrack_ErrorD0;   //!
   TBranch        *b_EFlowTrack_ErrorDZ;   //!
   TBranch        *b_EFlowTrack_ErrorC;   //!
   TBranch        *b_EFlowTrack_ErrorD0Phi;   //!
   TBranch        *b_EFlowTrack_ErrorD0C;   //!
   TBranch        *b_EFlowTrack_ErrorD0DZ;   //!
   TBranch        *b_EFlowTrack_ErrorD0CtgTheta;   //!
   TBranch        *b_EFlowTrack_ErrorPhiC;   //!
   TBranch        *b_EFlowTrack_ErrorPhiDZ;   //!
   TBranch        *b_EFlowTrack_ErrorPhiCtgTheta;   //!
   TBranch        *b_EFlowTrack_ErrorCDZ;   //!
   TBranch        *b_EFlowTrack_ErrorCCtgTheta;   //!
   TBranch        *b_EFlowTrack_ErrorDZCtgTheta;   //!
   TBranch        *b_EFlowTrack_Particle;   //!
   TBranch        *b_EFlowTrack_VertexIndex;   //!
   TBranch        *b_EFlowTrack_size;   //!
   TBranch        *b_EFlowPhoton_;   //!
   TBranch        *b_EFlowPhoton_fUniqueID;   //!
   TBranch        *b_EFlowPhoton_fBits;   //!
   TBranch        *b_EFlowPhoton_ET;   //!
   TBranch        *b_EFlowPhoton_Eta;   //!
   TBranch        *b_EFlowPhoton_Phi;   //!
   TBranch        *b_EFlowPhoton_E;   //!
   TBranch        *b_EFlowPhoton_T;   //!
   TBranch        *b_EFlowPhoton_NTimeHits;   //!
   TBranch        *b_EFlowPhoton_Eem;   //!
   TBranch        *b_EFlowPhoton_Ehad;   //!
   TBranch        *b_EFlowPhoton_Edges;   //!
   TBranch        *b_EFlowPhoton_Particles;   //!
   TBranch        *b_EFlowPhoton_size;   //!
   TBranch        *b_EFlowNeutralHadron_;   //!
   TBranch        *b_EFlowNeutralHadron_fUniqueID;   //!
   TBranch        *b_EFlowNeutralHadron_fBits;   //!
   TBranch        *b_EFlowNeutralHadron_ET;   //!
   TBranch        *b_EFlowNeutralHadron_Eta;   //!
   TBranch        *b_EFlowNeutralHadron_Phi;   //!
   TBranch        *b_EFlowNeutralHadron_E;   //!
   TBranch        *b_EFlowNeutralHadron_T;   //!
   TBranch        *b_EFlowNeutralHadron_NTimeHits;   //!
   TBranch        *b_EFlowNeutralHadron_Eem;   //!
   TBranch        *b_EFlowNeutralHadron_Ehad;   //!
   TBranch        *b_EFlowNeutralHadron_Edges;   //!
   TBranch        *b_EFlowNeutralHadron_Particles;   //!
   TBranch        *b_EFlowNeutralHadron_size;   //!
   TBranch        *b_mRICHTrack_;   //!
   TBranch        *b_mRICHTrack_fUniqueID;   //!
   TBranch        *b_mRICHTrack_fBits;   //!
   TBranch        *b_mRICHTrack_PID;   //!
   TBranch        *b_mRICHTrack_Charge;   //!
   TBranch        *b_mRICHTrack_P;   //!
   TBranch        *b_mRICHTrack_PT;   //!
   TBranch        *b_mRICHTrack_Eta;   //!
   TBranch        *b_mRICHTrack_Phi;   //!
   TBranch        *b_mRICHTrack_CtgTheta;   //!
   TBranch        *b_mRICHTrack_C;   //!
   TBranch        *b_mRICHTrack_Mass;   //!
   TBranch        *b_mRICHTrack_EtaOuter;   //!
   TBranch        *b_mRICHTrack_PhiOuter;   //!
   TBranch        *b_mRICHTrack_T;   //!
   TBranch        *b_mRICHTrack_X;   //!
   TBranch        *b_mRICHTrack_Y;   //!
   TBranch        *b_mRICHTrack_Z;   //!
   TBranch        *b_mRICHTrack_TOuter;   //!
   TBranch        *b_mRICHTrack_XOuter;   //!
   TBranch        *b_mRICHTrack_YOuter;   //!
   TBranch        *b_mRICHTrack_ZOuter;   //!
   TBranch        *b_mRICHTrack_Xd;   //!
   TBranch        *b_mRICHTrack_Yd;   //!
   TBranch        *b_mRICHTrack_Zd;   //!
   TBranch        *b_mRICHTrack_L;   //!
   TBranch        *b_mRICHTrack_D0;   //!
   TBranch        *b_mRICHTrack_DZ;   //!
   TBranch        *b_mRICHTrack_Nclusters;   //!
   TBranch        *b_mRICHTrack_dNdx;   //!
   TBranch        *b_mRICHTrack_ErrorP;   //!
   TBranch        *b_mRICHTrack_ErrorPT;   //!
   TBranch        *b_mRICHTrack_ErrorPhi;   //!
   TBranch        *b_mRICHTrack_ErrorCtgTheta;   //!
   TBranch        *b_mRICHTrack_ErrorT;   //!
   TBranch        *b_mRICHTrack_ErrorD0;   //!
   TBranch        *b_mRICHTrack_ErrorDZ;   //!
   TBranch        *b_mRICHTrack_ErrorC;   //!
   TBranch        *b_mRICHTrack_ErrorD0Phi;   //!
   TBranch        *b_mRICHTrack_ErrorD0C;   //!
   TBranch        *b_mRICHTrack_ErrorD0DZ;   //!
   TBranch        *b_mRICHTrack_ErrorD0CtgTheta;   //!
   TBranch        *b_mRICHTrack_ErrorPhiC;   //!
   TBranch        *b_mRICHTrack_ErrorPhiDZ;   //!
   TBranch        *b_mRICHTrack_ErrorPhiCtgTheta;   //!
   TBranch        *b_mRICHTrack_ErrorCDZ;   //!
   TBranch        *b_mRICHTrack_ErrorCCtgTheta;   //!
   TBranch        *b_mRICHTrack_ErrorDZCtgTheta;   //!
   TBranch        *b_mRICHTrack_Particle;   //!
   TBranch        *b_mRICHTrack_VertexIndex;   //!
   TBranch        *b_mRICHTrack_size;   //!
   TBranch        *b_barrelDIRCTrack_;   //!
   TBranch        *b_barrelDIRCTrack_fUniqueID;   //!
   TBranch        *b_barrelDIRCTrack_fBits;   //!
   TBranch        *b_barrelDIRCTrack_PID;   //!
   TBranch        *b_barrelDIRCTrack_Charge;   //!
   TBranch        *b_barrelDIRCTrack_P;   //!
   TBranch        *b_barrelDIRCTrack_PT;   //!
   TBranch        *b_barrelDIRCTrack_Eta;   //!
   TBranch        *b_barrelDIRCTrack_Phi;   //!
   TBranch        *b_barrelDIRCTrack_CtgTheta;   //!
   TBranch        *b_barrelDIRCTrack_C;   //!
   TBranch        *b_barrelDIRCTrack_Mass;   //!
   TBranch        *b_barrelDIRCTrack_EtaOuter;   //!
   TBranch        *b_barrelDIRCTrack_PhiOuter;   //!
   TBranch        *b_barrelDIRCTrack_T;   //!
   TBranch        *b_barrelDIRCTrack_X;   //!
   TBranch        *b_barrelDIRCTrack_Y;   //!
   TBranch        *b_barrelDIRCTrack_Z;   //!
   TBranch        *b_barrelDIRCTrack_TOuter;   //!
   TBranch        *b_barrelDIRCTrack_XOuter;   //!
   TBranch        *b_barrelDIRCTrack_YOuter;   //!
   TBranch        *b_barrelDIRCTrack_ZOuter;   //!
   TBranch        *b_barrelDIRCTrack_Xd;   //!
   TBranch        *b_barrelDIRCTrack_Yd;   //!
   TBranch        *b_barrelDIRCTrack_Zd;   //!
   TBranch        *b_barrelDIRCTrack_L;   //!
   TBranch        *b_barrelDIRCTrack_D0;   //!
   TBranch        *b_barrelDIRCTrack_DZ;   //!
   TBranch        *b_barrelDIRCTrack_Nclusters;   //!
   TBranch        *b_barrelDIRCTrack_dNdx;   //!
   TBranch        *b_barrelDIRCTrack_ErrorP;   //!
   TBranch        *b_barrelDIRCTrack_ErrorPT;   //!
   TBranch        *b_barrelDIRCTrack_ErrorPhi;   //!
   TBranch        *b_barrelDIRCTrack_ErrorCtgTheta;   //!
   TBranch        *b_barrelDIRCTrack_ErrorT;   //!
   TBranch        *b_barrelDIRCTrack_ErrorD0;   //!
   TBranch        *b_barrelDIRCTrack_ErrorDZ;   //!
   TBranch        *b_barrelDIRCTrack_ErrorC;   //!
   TBranch        *b_barrelDIRCTrack_ErrorD0Phi;   //!
   TBranch        *b_barrelDIRCTrack_ErrorD0C;   //!
   TBranch        *b_barrelDIRCTrack_ErrorD0DZ;   //!
   TBranch        *b_barrelDIRCTrack_ErrorD0CtgTheta;   //!
   TBranch        *b_barrelDIRCTrack_ErrorPhiC;   //!
   TBranch        *b_barrelDIRCTrack_ErrorPhiDZ;   //!
   TBranch        *b_barrelDIRCTrack_ErrorPhiCtgTheta;   //!
   TBranch        *b_barrelDIRCTrack_ErrorCDZ;   //!
   TBranch        *b_barrelDIRCTrack_ErrorCCtgTheta;   //!
   TBranch        *b_barrelDIRCTrack_ErrorDZCtgTheta;   //!
   TBranch        *b_barrelDIRCTrack_Particle;   //!
   TBranch        *b_barrelDIRCTrack_VertexIndex;   //!
   TBranch        *b_barrelDIRCTrack_size;   //!
   TBranch        *b_dualRICHagTrack_;   //!
   TBranch        *b_dualRICHagTrack_fUniqueID;   //!
   TBranch        *b_dualRICHagTrack_fBits;   //!
   TBranch        *b_dualRICHagTrack_PID;   //!
   TBranch        *b_dualRICHagTrack_Charge;   //!
   TBranch        *b_dualRICHagTrack_P;   //!
   TBranch        *b_dualRICHagTrack_PT;   //!
   TBranch        *b_dualRICHagTrack_Eta;   //!
   TBranch        *b_dualRICHagTrack_Phi;   //!
   TBranch        *b_dualRICHagTrack_CtgTheta;   //!
   TBranch        *b_dualRICHagTrack_C;   //!
   TBranch        *b_dualRICHagTrack_Mass;   //!
   TBranch        *b_dualRICHagTrack_EtaOuter;   //!
   TBranch        *b_dualRICHagTrack_PhiOuter;   //!
   TBranch        *b_dualRICHagTrack_T;   //!
   TBranch        *b_dualRICHagTrack_X;   //!
   TBranch        *b_dualRICHagTrack_Y;   //!
   TBranch        *b_dualRICHagTrack_Z;   //!
   TBranch        *b_dualRICHagTrack_TOuter;   //!
   TBranch        *b_dualRICHagTrack_XOuter;   //!
   TBranch        *b_dualRICHagTrack_YOuter;   //!
   TBranch        *b_dualRICHagTrack_ZOuter;   //!
   TBranch        *b_dualRICHagTrack_Xd;   //!
   TBranch        *b_dualRICHagTrack_Yd;   //!
   TBranch        *b_dualRICHagTrack_Zd;   //!
   TBranch        *b_dualRICHagTrack_L;   //!
   TBranch        *b_dualRICHagTrack_D0;   //!
   TBranch        *b_dualRICHagTrack_DZ;   //!
   TBranch        *b_dualRICHagTrack_Nclusters;   //!
   TBranch        *b_dualRICHagTrack_dNdx;   //!
   TBranch        *b_dualRICHagTrack_ErrorP;   //!
   TBranch        *b_dualRICHagTrack_ErrorPT;   //!
   TBranch        *b_dualRICHagTrack_ErrorPhi;   //!
   TBranch        *b_dualRICHagTrack_ErrorCtgTheta;   //!
   TBranch        *b_dualRICHagTrack_ErrorT;   //!
   TBranch        *b_dualRICHagTrack_ErrorD0;   //!
   TBranch        *b_dualRICHagTrack_ErrorDZ;   //!
   TBranch        *b_dualRICHagTrack_ErrorC;   //!
   TBranch        *b_dualRICHagTrack_ErrorD0Phi;   //!
   TBranch        *b_dualRICHagTrack_ErrorD0C;   //!
   TBranch        *b_dualRICHagTrack_ErrorD0DZ;   //!
   TBranch        *b_dualRICHagTrack_ErrorD0CtgTheta;   //!
   TBranch        *b_dualRICHagTrack_ErrorPhiC;   //!
   TBranch        *b_dualRICHagTrack_ErrorPhiDZ;   //!
   TBranch        *b_dualRICHagTrack_ErrorPhiCtgTheta;   //!
   TBranch        *b_dualRICHagTrack_ErrorCDZ;   //!
   TBranch        *b_dualRICHagTrack_ErrorCCtgTheta;   //!
   TBranch        *b_dualRICHagTrack_ErrorDZCtgTheta;   //!
   TBranch        *b_dualRICHagTrack_Particle;   //!
   TBranch        *b_dualRICHagTrack_VertexIndex;   //!
   TBranch        *b_dualRICHagTrack_size;   //!
   TBranch        *b_dualRICHcfTrack_;   //!
   TBranch        *b_dualRICHcfTrack_fUniqueID;   //!
   TBranch        *b_dualRICHcfTrack_fBits;   //!
   TBranch        *b_dualRICHcfTrack_PID;   //!
   TBranch        *b_dualRICHcfTrack_Charge;   //!
   TBranch        *b_dualRICHcfTrack_P;   //!
   TBranch        *b_dualRICHcfTrack_PT;   //!
   TBranch        *b_dualRICHcfTrack_Eta;   //!
   TBranch        *b_dualRICHcfTrack_Phi;   //!
   TBranch        *b_dualRICHcfTrack_CtgTheta;   //!
   TBranch        *b_dualRICHcfTrack_C;   //!
   TBranch        *b_dualRICHcfTrack_Mass;   //!
   TBranch        *b_dualRICHcfTrack_EtaOuter;   //!
   TBranch        *b_dualRICHcfTrack_PhiOuter;   //!
   TBranch        *b_dualRICHcfTrack_T;   //!
   TBranch        *b_dualRICHcfTrack_X;   //!
   TBranch        *b_dualRICHcfTrack_Y;   //!
   TBranch        *b_dualRICHcfTrack_Z;   //!
   TBranch        *b_dualRICHcfTrack_TOuter;   //!
   TBranch        *b_dualRICHcfTrack_XOuter;   //!
   TBranch        *b_dualRICHcfTrack_YOuter;   //!
   TBranch        *b_dualRICHcfTrack_ZOuter;   //!
   TBranch        *b_dualRICHcfTrack_Xd;   //!
   TBranch        *b_dualRICHcfTrack_Yd;   //!
   TBranch        *b_dualRICHcfTrack_Zd;   //!
   TBranch        *b_dualRICHcfTrack_L;   //!
   TBranch        *b_dualRICHcfTrack_D0;   //!
   TBranch        *b_dualRICHcfTrack_DZ;   //!
   TBranch        *b_dualRICHcfTrack_Nclusters;   //!
   TBranch        *b_dualRICHcfTrack_dNdx;   //!
   TBranch        *b_dualRICHcfTrack_ErrorP;   //!
   TBranch        *b_dualRICHcfTrack_ErrorPT;   //!
   TBranch        *b_dualRICHcfTrack_ErrorPhi;   //!
   TBranch        *b_dualRICHcfTrack_ErrorCtgTheta;   //!
   TBranch        *b_dualRICHcfTrack_ErrorT;   //!
   TBranch        *b_dualRICHcfTrack_ErrorD0;   //!
   TBranch        *b_dualRICHcfTrack_ErrorDZ;   //!
   TBranch        *b_dualRICHcfTrack_ErrorC;   //!
   TBranch        *b_dualRICHcfTrack_ErrorD0Phi;   //!
   TBranch        *b_dualRICHcfTrack_ErrorD0C;   //!
   TBranch        *b_dualRICHcfTrack_ErrorD0DZ;   //!
   TBranch        *b_dualRICHcfTrack_ErrorD0CtgTheta;   //!
   TBranch        *b_dualRICHcfTrack_ErrorPhiC;   //!
   TBranch        *b_dualRICHcfTrack_ErrorPhiDZ;   //!
   TBranch        *b_dualRICHcfTrack_ErrorPhiCtgTheta;   //!
   TBranch        *b_dualRICHcfTrack_ErrorCDZ;   //!
   TBranch        *b_dualRICHcfTrack_ErrorCCtgTheta;   //!
   TBranch        *b_dualRICHcfTrack_ErrorDZCtgTheta;   //!
   TBranch        *b_dualRICHcfTrack_Particle;   //!
   TBranch        *b_dualRICHcfTrack_VertexIndex;   //!
   TBranch        *b_dualRICHcfTrack_size;   //!
   TBranch        *b_GenJet_;   //!
   TBranch        *b_GenJet_fUniqueID;   //!
   TBranch        *b_GenJet_fBits;   //!
   TBranch        *b_GenJet_PT;   //!
   TBranch        *b_GenJet_Eta;   //!
   TBranch        *b_GenJet_Phi;   //!
   TBranch        *b_GenJet_T;   //!
   TBranch        *b_GenJet_Mass;   //!
   TBranch        *b_GenJet_DeltaEta;   //!
   TBranch        *b_GenJet_DeltaPhi;   //!
   TBranch        *b_GenJet_Flavor;   //!
   TBranch        *b_GenJet_FlavorAlgo;   //!
   TBranch        *b_GenJet_FlavorPhys;   //!
   TBranch        *b_GenJet_BTag;   //!
   TBranch        *b_GenJet_BTagAlgo;   //!
   TBranch        *b_GenJet_BTagPhys;   //!
   TBranch        *b_GenJet_TauTag;   //!
   TBranch        *b_GenJet_TauWeight;   //!
   TBranch        *b_GenJet_Charge;   //!
   TBranch        *b_GenJet_EhadOverEem;   //!
   TBranch        *b_GenJet_NCharged;   //!
   TBranch        *b_GenJet_NNeutrals;   //!
   TBranch        *b_GenJet_NeutralEnergyFraction;   //!
   TBranch        *b_GenJet_ChargedEnergyFraction;   //!
   TBranch        *b_GenJet_Beta;   //!
   TBranch        *b_GenJet_BetaStar;   //!
   TBranch        *b_GenJet_MeanSqDeltaR;   //!
   TBranch        *b_GenJet_PTD;   //!
   TBranch        *b_GenJet_FracPt;   //!
   TBranch        *b_GenJet_Tau;   //!
   TBranch        *b_GenJet_SoftDroppedJet;   //!
   TBranch        *b_GenJet_SoftDroppedSubJet1;   //!
   TBranch        *b_GenJet_SoftDroppedSubJet2;   //!
   TBranch        *b_GenJet_TrimmedP4;   //!
   TBranch        *b_GenJet_PrunedP4;   //!
   TBranch        *b_GenJet_SoftDroppedP4;   //!
   TBranch        *b_GenJet_NSubJetsTrimmed;   //!
   TBranch        *b_GenJet_NSubJetsPruned;   //!
   TBranch        *b_GenJet_NSubJetsSoftDropped;   //!
   TBranch        *b_GenJet_ExclYmerge23;   //!
   TBranch        *b_GenJet_ExclYmerge34;   //!
   TBranch        *b_GenJet_ExclYmerge45;   //!
   TBranch        *b_GenJet_ExclYmerge56;   //!
   TBranch        *b_GenJet_Constituents;   //!
   TBranch        *b_GenJet_Particles;   //!
   TBranch        *b_GenJet_Area;   //!
   TBranch        *b_GenJet_size;   //!
   TBranch        *b_GenMissingET_;   //!
   TBranch        *b_GenMissingET_fUniqueID;   //!
   TBranch        *b_GenMissingET_fBits;   //!
   TBranch        *b_GenMissingET_MET;   //!
   TBranch        *b_GenMissingET_Eta;   //!
   TBranch        *b_GenMissingET_Phi;   //!
   TBranch        *b_GenMissingET_size;   //!
   TBranch        *b_Jet_;   //!
   TBranch        *b_Jet_fUniqueID;   //!
   TBranch        *b_Jet_fBits;   //!
   TBranch        *b_Jet_PT;   //!
   TBranch        *b_Jet_Eta;   //!
   TBranch        *b_Jet_Phi;   //!
   TBranch        *b_Jet_T;   //!
   TBranch        *b_Jet_Mass;   //!
   TBranch        *b_Jet_DeltaEta;   //!
   TBranch        *b_Jet_DeltaPhi;   //!
   TBranch        *b_Jet_Flavor;   //!
   TBranch        *b_Jet_FlavorAlgo;   //!
   TBranch        *b_Jet_FlavorPhys;   //!
   TBranch        *b_Jet_BTag;   //!
   TBranch        *b_Jet_BTagAlgo;   //!
   TBranch        *b_Jet_BTagPhys;   //!
   TBranch        *b_Jet_TauTag;   //!
   TBranch        *b_Jet_TauWeight;   //!
   TBranch        *b_Jet_Charge;   //!
   TBranch        *b_Jet_EhadOverEem;   //!
   TBranch        *b_Jet_NCharged;   //!
   TBranch        *b_Jet_NNeutrals;   //!
   TBranch        *b_Jet_NeutralEnergyFraction;   //!
   TBranch        *b_Jet_ChargedEnergyFraction;   //!
   TBranch        *b_Jet_Beta;   //!
   TBranch        *b_Jet_BetaStar;   //!
   TBranch        *b_Jet_MeanSqDeltaR;   //!
   TBranch        *b_Jet_PTD;   //!
   TBranch        *b_Jet_FracPt;   //!
   TBranch        *b_Jet_Tau;   //!
   TBranch        *b_Jet_SoftDroppedJet;   //!
   TBranch        *b_Jet_SoftDroppedSubJet1;   //!
   TBranch        *b_Jet_SoftDroppedSubJet2;   //!
   TBranch        *b_Jet_TrimmedP4;   //!
   TBranch        *b_Jet_PrunedP4;   //!
   TBranch        *b_Jet_SoftDroppedP4;   //!
   TBranch        *b_Jet_NSubJetsTrimmed;   //!
   TBranch        *b_Jet_NSubJetsPruned;   //!
   TBranch        *b_Jet_NSubJetsSoftDropped;   //!
   TBranch        *b_Jet_ExclYmerge23;   //!
   TBranch        *b_Jet_ExclYmerge34;   //!
   TBranch        *b_Jet_ExclYmerge45;   //!
   TBranch        *b_Jet_ExclYmerge56;   //!
   TBranch        *b_Jet_Constituents;   //!
   TBranch        *b_Jet_Particles;   //!
   TBranch        *b_Jet_Area;   //!
   TBranch        *b_Jet_size;   //!
   TBranch        *b_Electron_;   //!
   TBranch        *b_Electron_fUniqueID;   //!
   TBranch        *b_Electron_fBits;   //!
   TBranch        *b_Electron_PT;   //!
   TBranch        *b_Electron_Eta;   //!
   TBranch        *b_Electron_Phi;   //!
   TBranch        *b_Electron_T;   //!
   TBranch        *b_Electron_Charge;   //!
   TBranch        *b_Electron_EhadOverEem;   //!
   TBranch        *b_Electron_Particle;   //!
   TBranch        *b_Electron_IsolationVar;   //!
   TBranch        *b_Electron_IsolationVarRhoCorr;   //!
   TBranch        *b_Electron_SumPtCharged;   //!
   TBranch        *b_Electron_SumPtNeutral;   //!
   TBranch        *b_Electron_SumPtChargedPU;   //!
   TBranch        *b_Electron_SumPt;   //!
   TBranch        *b_Electron_D0;   //!
   TBranch        *b_Electron_DZ;   //!
   TBranch        *b_Electron_ErrorD0;   //!
   TBranch        *b_Electron_ErrorDZ;   //!
   TBranch        *b_Electron_size;   //!
   TBranch        *b_Photon_;   //!
   TBranch        *b_Photon_fUniqueID;   //!
   TBranch        *b_Photon_fBits;   //!
   TBranch        *b_Photon_PT;   //!
   TBranch        *b_Photon_Eta;   //!
   TBranch        *b_Photon_Phi;   //!
   TBranch        *b_Photon_E;   //!
   TBranch        *b_Photon_T;   //!
   TBranch        *b_Photon_EhadOverEem;   //!
   TBranch        *b_Photon_Particles;   //!
   TBranch        *b_Photon_IsolationVar;   //!
   TBranch        *b_Photon_IsolationVarRhoCorr;   //!
   TBranch        *b_Photon_SumPtCharged;   //!
   TBranch        *b_Photon_SumPtNeutral;   //!
   TBranch        *b_Photon_SumPtChargedPU;   //!
   TBranch        *b_Photon_SumPt;   //!
   TBranch        *b_Photon_Status;   //!
   TBranch        *b_Photon_size;   //!
   TBranch        *b_MissingET_;   //!
   TBranch        *b_MissingET_fUniqueID;   //!
   TBranch        *b_MissingET_fBits;   //!
   TBranch        *b_MissingET_MET;   //!
   TBranch        *b_MissingET_Eta;   //!
   TBranch        *b_MissingET_Phi;   //!
   TBranch        *b_MissingET_size;   //!
   TBranch        *b_ScalarHT_;   //!
   TBranch        *b_ScalarHT_fUniqueID;   //!
   TBranch        *b_ScalarHT_fBits;   //!
   TBranch        *b_ScalarHT_HT;   //!
   TBranch        *b_ScalarHT_size;   //!

   char outfile_name[1000] ;

   delphes_rapgap_isr_fsr_analysis_h1( const char* infilename = "/Volumes/Ext_2020_08/dis-reco-work/work-2021-10-09/h1-rapgap1/h1-rapgap-q2min-200-80200.root",
                                    const char* outfile = "h1-rad-tree.root" );
   virtual ~delphes_rapgap_isr_fsr_analysis_h1();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(bool verbose=true, int maxEvts=10);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef delphes_rapgap_isr_fsr_analysis_h1_cxx
//delphes_rapgap_isr_fsr_analysis_h1::delphes_rapgap_isr_fsr_analysis_h1(TTree *tree) : fChain(0) 
//{
//// if parameter tree is not specified (or zero), connect the file
//// used to generate this class and read the Tree.
//   if (tree == 0) {
//      //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-02/athena-rapgap/athena-rapgap001.root");
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-02/athena-rapgap/athena-rapgap-q2min-200-000.root");
//      if (!f || !f->IsOpen()) {
//         //f = new TFile("/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-02/athena-rapgap/athena-rapgap001.root");
//         f = new TFile("/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-02/athena-rapgap/athena-rapgap-q2min-200-000.root");
//      }
//      f->GetObject("Delphes",tree);
//
//   }
//   Init(tree);
//}
delphes_rapgap_isr_fsr_analysis_h1::delphes_rapgap_isr_fsr_analysis_h1( const char* infilename, const char* outfile ) : fChain(0) 
{
   sprintf( outfile_name, "%s", outfile ) ;
   TChain* ch = new TChain("Delphes") ;

   if ( strcmp( infilename, "chunk1" ) == 0 ) {

      sprintf( outfile_name, "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/dnn-inputs-chunk1.root" ) ;

      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20200.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20201.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20202.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20203.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20204.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20205.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20206.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20207.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20208.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20209.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20210.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20211.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20212.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20213.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20220.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20221.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20222.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20223.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20224.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20225.root" ) ;


   } else if ( strcmp( infilename, "chunk2" ) == 0 ) {

      sprintf( outfile_name, "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/dnn-inputs-chunk2.root" ) ;

      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20226.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20227.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20228.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20229.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20230.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20231.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20232.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20233.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20234.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20240.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20241.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20242.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20243.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20244.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20245.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20246.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20247.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20248.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20249.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20250.root" ) ;

   } else if ( strcmp( infilename, "chunk3" ) == 0 ) {

      sprintf( outfile_name, "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/dnn-inputs-chunk3.root" ) ;

      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20251.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20252.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20253.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20254.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20260.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20261.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20262.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20263.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20264.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20265.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20266.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20267.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20268.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20269.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20270.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20271.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20272.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20280.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20281.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20282.root" ) ;

   } else if ( strcmp( infilename, "chunk4" ) == 0 ) {

      sprintf( outfile_name, "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/dnn-inputs-chunk4.root" ) ;

      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20283.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20284.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20285.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20286.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20287.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20288.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20289.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20290.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20291.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20292.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20293.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20295.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20300.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20301.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20302.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20303.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20304.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20305.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20306.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20307.root" ) ;

   } else if ( strcmp( infilename, "chunk5" ) == 0 ) {

      sprintf( outfile_name, "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/dnn-inputs-chunk5.root" ) ;

      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20308.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20309.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20310.root" ) ;
      ch -> Add( "/Volumes/Ext_2020_08/dis-reco-work/work-2021-09-14/athena-rapgap1/athena-rapgap-q2min-200-20311.root" ) ;

   } else {

      ch -> Add( infilename ) ;

   }
   int entries_check = ch -> GetEntries() ;
   if ( entries_check <= 0 ) {
      printf("\n\n *** file %s seems to be missing the Delphes TTree.\n\n", infilename ) ;
      gSystem -> Exit(-1) ;
   }
   printf("\n\n File %s has %d events.\n\n", infilename, entries_check ) ;
   printf(" Output file will be : %s\n\n", outfile_name ) ;

   TTree* tree = ch ;

   Init(tree);
}

delphes_rapgap_isr_fsr_analysis_h1::~delphes_rapgap_isr_fsr_analysis_h1()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t delphes_rapgap_isr_fsr_analysis_h1::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t delphes_rapgap_isr_fsr_analysis_h1::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void delphes_rapgap_isr_fsr_analysis_h1::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event_, &b_Event_);
   fChain->SetBranchAddress("Event.fUniqueID", Event_fUniqueID, &b_Event_fUniqueID);
   fChain->SetBranchAddress("Event.fBits", Event_fBits, &b_Event_fBits);
   fChain->SetBranchAddress("Event.Number", Event_Number, &b_Event_Number);
   fChain->SetBranchAddress("Event.ReadTime", Event_ReadTime, &b_Event_ReadTime);
   fChain->SetBranchAddress("Event.ProcTime", Event_ProcTime, &b_Event_ProcTime);
   fChain->SetBranchAddress("Event.ProcessID", Event_ProcessID, &b_Event_ProcessID);
   fChain->SetBranchAddress("Event.MPI", Event_MPI, &b_Event_MPI);
   fChain->SetBranchAddress("Event.Weight", Event_Weight, &b_Event_Weight);
   fChain->SetBranchAddress("Event.CrossSection", Event_CrossSection, &b_Event_CrossSection);
   fChain->SetBranchAddress("Event.CrossSectionError", Event_CrossSectionError, &b_Event_CrossSectionError);
   fChain->SetBranchAddress("Event.Scale", Event_Scale, &b_Event_Scale);
   fChain->SetBranchAddress("Event.AlphaQED", Event_AlphaQED, &b_Event_AlphaQED);
   fChain->SetBranchAddress("Event.AlphaQCD", Event_AlphaQCD, &b_Event_AlphaQCD);
   fChain->SetBranchAddress("Event.ID1", Event_ID1, &b_Event_ID1);
   fChain->SetBranchAddress("Event.ID2", Event_ID2, &b_Event_ID2);
   fChain->SetBranchAddress("Event.X1", Event_X1, &b_Event_X1);
   fChain->SetBranchAddress("Event.X2", Event_X2, &b_Event_X2);
   fChain->SetBranchAddress("Event.ScalePDF", Event_ScalePDF, &b_Event_ScalePDF);
   fChain->SetBranchAddress("Event.PDF1", Event_PDF1, &b_Event_PDF1);
   fChain->SetBranchAddress("Event.PDF2", Event_PDF2, &b_Event_PDF2);
   fChain->SetBranchAddress("Event_size", &Event_size, &b_Event_size);
   fChain->SetBranchAddress("Weight", &Weight_, &b_Weight_);
   fChain->SetBranchAddress("Weight.fUniqueID", Weight_fUniqueID, &b_Weight_fUniqueID);
   fChain->SetBranchAddress("Weight.fBits", Weight_fBits, &b_Weight_fBits);
   fChain->SetBranchAddress("Weight.Weight", Weight_Weight, &b_Weight_Weight);
   fChain->SetBranchAddress("Weight_size", &Weight_size, &b_Weight_size);
   fChain->SetBranchAddress("Particle", &Particle_, &b_Particle_);
   fChain->SetBranchAddress("Particle.fUniqueID", Particle_fUniqueID, &b_Particle_fUniqueID);
   fChain->SetBranchAddress("Particle.fBits", Particle_fBits, &b_Particle_fBits);
   fChain->SetBranchAddress("Particle.PID", Particle_PID, &b_Particle_PID);
   fChain->SetBranchAddress("Particle.Status", Particle_Status, &b_Particle_Status);
   fChain->SetBranchAddress("Particle.IsPU", Particle_IsPU, &b_Particle_IsPU);
   fChain->SetBranchAddress("Particle.M1", Particle_M1, &b_Particle_M1);
   fChain->SetBranchAddress("Particle.M2", Particle_M2, &b_Particle_M2);
   fChain->SetBranchAddress("Particle.D1", Particle_D1, &b_Particle_D1);
   fChain->SetBranchAddress("Particle.D2", Particle_D2, &b_Particle_D2);
   fChain->SetBranchAddress("Particle.Charge", Particle_Charge, &b_Particle_Charge);
   fChain->SetBranchAddress("Particle.Mass", Particle_Mass, &b_Particle_Mass);
   fChain->SetBranchAddress("Particle.E", Particle_E, &b_Particle_E);
   fChain->SetBranchAddress("Particle.Px", Particle_Px, &b_Particle_Px);
   fChain->SetBranchAddress("Particle.Py", Particle_Py, &b_Particle_Py);
   fChain->SetBranchAddress("Particle.Pz", Particle_Pz, &b_Particle_Pz);
   fChain->SetBranchAddress("Particle.P", Particle_P, &b_Particle_P);
   fChain->SetBranchAddress("Particle.PT", Particle_PT, &b_Particle_PT);
   fChain->SetBranchAddress("Particle.Eta", Particle_Eta, &b_Particle_Eta);
   fChain->SetBranchAddress("Particle.Phi", Particle_Phi, &b_Particle_Phi);
   fChain->SetBranchAddress("Particle.Rapidity", Particle_Rapidity, &b_Particle_Rapidity);
   fChain->SetBranchAddress("Particle.T", Particle_T, &b_Particle_T);
   fChain->SetBranchAddress("Particle.X", Particle_X, &b_Particle_X);
   fChain->SetBranchAddress("Particle.Y", Particle_Y, &b_Particle_Y);
   fChain->SetBranchAddress("Particle.Z", Particle_Z, &b_Particle_Z);
   fChain->SetBranchAddress("Particle_size", &Particle_size, &b_Particle_size);
   fChain->SetBranchAddress("BeamSpot", &BeamSpot_, &b_BeamSpot_);
   fChain->SetBranchAddress("BeamSpot.fUniqueID", BeamSpot_fUniqueID, &b_BeamSpot_fUniqueID);
   fChain->SetBranchAddress("BeamSpot.fBits", BeamSpot_fBits, &b_BeamSpot_fBits);
   fChain->SetBranchAddress("BeamSpot.PID", BeamSpot_PID, &b_BeamSpot_PID);
   fChain->SetBranchAddress("BeamSpot.Status", BeamSpot_Status, &b_BeamSpot_Status);
   fChain->SetBranchAddress("BeamSpot.IsPU", BeamSpot_IsPU, &b_BeamSpot_IsPU);
   fChain->SetBranchAddress("BeamSpot.M1", BeamSpot_M1, &b_BeamSpot_M1);
   fChain->SetBranchAddress("BeamSpot.M2", BeamSpot_M2, &b_BeamSpot_M2);
   fChain->SetBranchAddress("BeamSpot.D1", BeamSpot_D1, &b_BeamSpot_D1);
   fChain->SetBranchAddress("BeamSpot.D2", BeamSpot_D2, &b_BeamSpot_D2);
   fChain->SetBranchAddress("BeamSpot.Charge", BeamSpot_Charge, &b_BeamSpot_Charge);
   fChain->SetBranchAddress("BeamSpot.Mass", BeamSpot_Mass, &b_BeamSpot_Mass);
   fChain->SetBranchAddress("BeamSpot.E", BeamSpot_E, &b_BeamSpot_E);
   fChain->SetBranchAddress("BeamSpot.Px", BeamSpot_Px, &b_BeamSpot_Px);
   fChain->SetBranchAddress("BeamSpot.Py", BeamSpot_Py, &b_BeamSpot_Py);
   fChain->SetBranchAddress("BeamSpot.Pz", BeamSpot_Pz, &b_BeamSpot_Pz);
   fChain->SetBranchAddress("BeamSpot.P", BeamSpot_P, &b_BeamSpot_P);
   fChain->SetBranchAddress("BeamSpot.PT", BeamSpot_PT, &b_BeamSpot_PT);
   fChain->SetBranchAddress("BeamSpot.Eta", BeamSpot_Eta, &b_BeamSpot_Eta);
   fChain->SetBranchAddress("BeamSpot.Phi", BeamSpot_Phi, &b_BeamSpot_Phi);
   fChain->SetBranchAddress("BeamSpot.Rapidity", BeamSpot_Rapidity, &b_BeamSpot_Rapidity);
   fChain->SetBranchAddress("BeamSpot.T", BeamSpot_T, &b_BeamSpot_T);
   fChain->SetBranchAddress("BeamSpot.X", BeamSpot_X, &b_BeamSpot_X);
   fChain->SetBranchAddress("BeamSpot.Y", BeamSpot_Y, &b_BeamSpot_Y);
   fChain->SetBranchAddress("BeamSpot.Z", BeamSpot_Z, &b_BeamSpot_Z);
   fChain->SetBranchAddress("BeamSpot_size", &BeamSpot_size, &b_BeamSpot_size);
   fChain->SetBranchAddress("Track", &Track_, &b_Track_);
   fChain->SetBranchAddress("Track.fUniqueID", &Track_fUniqueID, &b_Track_fUniqueID);
   fChain->SetBranchAddress("Track.fBits", &Track_fBits, &b_Track_fBits);
   fChain->SetBranchAddress("Track.PID", &Track_PID, &b_Track_PID);
   fChain->SetBranchAddress("Track.Charge", &Track_Charge, &b_Track_Charge);
   fChain->SetBranchAddress("Track.P", &Track_P, &b_Track_P);
   fChain->SetBranchAddress("Track.PT", &Track_PT, &b_Track_PT);
   fChain->SetBranchAddress("Track.Eta", &Track_Eta, &b_Track_Eta);
   fChain->SetBranchAddress("Track.Phi", &Track_Phi, &b_Track_Phi);
   fChain->SetBranchAddress("Track.CtgTheta", &Track_CtgTheta, &b_Track_CtgTheta);
   fChain->SetBranchAddress("Track.C", &Track_C, &b_Track_C);
   fChain->SetBranchAddress("Track.Mass", &Track_Mass, &b_Track_Mass);
   fChain->SetBranchAddress("Track.EtaOuter", &Track_EtaOuter, &b_Track_EtaOuter);
   fChain->SetBranchAddress("Track.PhiOuter", &Track_PhiOuter, &b_Track_PhiOuter);
   fChain->SetBranchAddress("Track.T", &Track_T, &b_Track_T);
   fChain->SetBranchAddress("Track.X", &Track_X, &b_Track_X);
   fChain->SetBranchAddress("Track.Y", &Track_Y, &b_Track_Y);
   fChain->SetBranchAddress("Track.Z", &Track_Z, &b_Track_Z);
   fChain->SetBranchAddress("Track.TOuter", &Track_TOuter, &b_Track_TOuter);
   fChain->SetBranchAddress("Track.XOuter", &Track_XOuter, &b_Track_XOuter);
   fChain->SetBranchAddress("Track.YOuter", &Track_YOuter, &b_Track_YOuter);
   fChain->SetBranchAddress("Track.ZOuter", &Track_ZOuter, &b_Track_ZOuter);
   fChain->SetBranchAddress("Track.Xd", &Track_Xd, &b_Track_Xd);
   fChain->SetBranchAddress("Track.Yd", &Track_Yd, &b_Track_Yd);
   fChain->SetBranchAddress("Track.Zd", &Track_Zd, &b_Track_Zd);
   fChain->SetBranchAddress("Track.L", &Track_L, &b_Track_L);
   fChain->SetBranchAddress("Track.D0", &Track_D0, &b_Track_D0);
   fChain->SetBranchAddress("Track.DZ", &Track_DZ, &b_Track_DZ);
   fChain->SetBranchAddress("Track.Nclusters", &Track_Nclusters, &b_Track_Nclusters);
   fChain->SetBranchAddress("Track.dNdx", &Track_dNdx, &b_Track_dNdx);
   fChain->SetBranchAddress("Track.ErrorP", &Track_ErrorP, &b_Track_ErrorP);
   fChain->SetBranchAddress("Track.ErrorPT", &Track_ErrorPT, &b_Track_ErrorPT);
   fChain->SetBranchAddress("Track.ErrorPhi", &Track_ErrorPhi, &b_Track_ErrorPhi);
   fChain->SetBranchAddress("Track.ErrorCtgTheta", &Track_ErrorCtgTheta, &b_Track_ErrorCtgTheta);
   fChain->SetBranchAddress("Track.ErrorT", &Track_ErrorT, &b_Track_ErrorT);
   fChain->SetBranchAddress("Track.ErrorD0", &Track_ErrorD0, &b_Track_ErrorD0);
   fChain->SetBranchAddress("Track.ErrorDZ", &Track_ErrorDZ, &b_Track_ErrorDZ);
   fChain->SetBranchAddress("Track.ErrorC", &Track_ErrorC, &b_Track_ErrorC);
   fChain->SetBranchAddress("Track.ErrorD0Phi", &Track_ErrorD0Phi, &b_Track_ErrorD0Phi);
   fChain->SetBranchAddress("Track.ErrorD0C", &Track_ErrorD0C, &b_Track_ErrorD0C);
   fChain->SetBranchAddress("Track.ErrorD0DZ", &Track_ErrorD0DZ, &b_Track_ErrorD0DZ);
   fChain->SetBranchAddress("Track.ErrorD0CtgTheta", &Track_ErrorD0CtgTheta, &b_Track_ErrorD0CtgTheta);
   fChain->SetBranchAddress("Track.ErrorPhiC", &Track_ErrorPhiC, &b_Track_ErrorPhiC);
   fChain->SetBranchAddress("Track.ErrorPhiDZ", &Track_ErrorPhiDZ, &b_Track_ErrorPhiDZ);
   fChain->SetBranchAddress("Track.ErrorPhiCtgTheta", &Track_ErrorPhiCtgTheta, &b_Track_ErrorPhiCtgTheta);
   fChain->SetBranchAddress("Track.ErrorCDZ", &Track_ErrorCDZ, &b_Track_ErrorCDZ);
   fChain->SetBranchAddress("Track.ErrorCCtgTheta", &Track_ErrorCCtgTheta, &b_Track_ErrorCCtgTheta);
   fChain->SetBranchAddress("Track.ErrorDZCtgTheta", &Track_ErrorDZCtgTheta, &b_Track_ErrorDZCtgTheta);
   fChain->SetBranchAddress("Track.Particle", &Track_Particle, &b_Track_Particle);
   fChain->SetBranchAddress("Track.VertexIndex", &Track_VertexIndex, &b_Track_VertexIndex);
   fChain->SetBranchAddress("Track_size", &Track_size, &b_Track_size);
   fChain->SetBranchAddress("Tower", &Tower_, &b_Tower_);
   fChain->SetBranchAddress("Tower.fUniqueID", Tower_fUniqueID, &b_Tower_fUniqueID);
   fChain->SetBranchAddress("Tower.fBits", Tower_fBits, &b_Tower_fBits);
   fChain->SetBranchAddress("Tower.ET", Tower_ET, &b_Tower_ET);
   fChain->SetBranchAddress("Tower.Eta", Tower_Eta, &b_Tower_Eta);
   fChain->SetBranchAddress("Tower.Phi", Tower_Phi, &b_Tower_Phi);
   fChain->SetBranchAddress("Tower.E", Tower_E, &b_Tower_E);
   fChain->SetBranchAddress("Tower.T", Tower_T, &b_Tower_T);
   fChain->SetBranchAddress("Tower.NTimeHits", Tower_NTimeHits, &b_Tower_NTimeHits);
   fChain->SetBranchAddress("Tower.Eem", Tower_Eem, &b_Tower_Eem);
   fChain->SetBranchAddress("Tower.Ehad", Tower_Ehad, &b_Tower_Ehad);
   fChain->SetBranchAddress("Tower.Edges[4]", Tower_Edges, &b_Tower_Edges);
   fChain->SetBranchAddress("Tower.Particles", Tower_Particles, &b_Tower_Particles);
   fChain->SetBranchAddress("Tower_size", &Tower_size, &b_Tower_size);
   fChain->SetBranchAddress("EFlowTrack", &EFlowTrack_, &b_EFlowTrack_);
   fChain->SetBranchAddress("EFlowTrack.fUniqueID", &EFlowTrack_fUniqueID, &b_EFlowTrack_fUniqueID);
   fChain->SetBranchAddress("EFlowTrack.fBits", &EFlowTrack_fBits, &b_EFlowTrack_fBits);
   fChain->SetBranchAddress("EFlowTrack.PID", &EFlowTrack_PID, &b_EFlowTrack_PID);
   fChain->SetBranchAddress("EFlowTrack.Charge", &EFlowTrack_Charge, &b_EFlowTrack_Charge);
   fChain->SetBranchAddress("EFlowTrack.P", &EFlowTrack_P, &b_EFlowTrack_P);
   fChain->SetBranchAddress("EFlowTrack.PT", &EFlowTrack_PT, &b_EFlowTrack_PT);
   fChain->SetBranchAddress("EFlowTrack.Eta", &EFlowTrack_Eta, &b_EFlowTrack_Eta);
   fChain->SetBranchAddress("EFlowTrack.Phi", &EFlowTrack_Phi, &b_EFlowTrack_Phi);
   fChain->SetBranchAddress("EFlowTrack.CtgTheta", &EFlowTrack_CtgTheta, &b_EFlowTrack_CtgTheta);
   fChain->SetBranchAddress("EFlowTrack.C", &EFlowTrack_C, &b_EFlowTrack_C);
   fChain->SetBranchAddress("EFlowTrack.Mass", &EFlowTrack_Mass, &b_EFlowTrack_Mass);
   fChain->SetBranchAddress("EFlowTrack.EtaOuter", &EFlowTrack_EtaOuter, &b_EFlowTrack_EtaOuter);
   fChain->SetBranchAddress("EFlowTrack.PhiOuter", &EFlowTrack_PhiOuter, &b_EFlowTrack_PhiOuter);
   fChain->SetBranchAddress("EFlowTrack.T", &EFlowTrack_T, &b_EFlowTrack_T);
   fChain->SetBranchAddress("EFlowTrack.X", &EFlowTrack_X, &b_EFlowTrack_X);
   fChain->SetBranchAddress("EFlowTrack.Y", &EFlowTrack_Y, &b_EFlowTrack_Y);
   fChain->SetBranchAddress("EFlowTrack.Z", &EFlowTrack_Z, &b_EFlowTrack_Z);
   fChain->SetBranchAddress("EFlowTrack.TOuter", &EFlowTrack_TOuter, &b_EFlowTrack_TOuter);
   fChain->SetBranchAddress("EFlowTrack.XOuter", &EFlowTrack_XOuter, &b_EFlowTrack_XOuter);
   fChain->SetBranchAddress("EFlowTrack.YOuter", &EFlowTrack_YOuter, &b_EFlowTrack_YOuter);
   fChain->SetBranchAddress("EFlowTrack.ZOuter", &EFlowTrack_ZOuter, &b_EFlowTrack_ZOuter);
   fChain->SetBranchAddress("EFlowTrack.Xd", &EFlowTrack_Xd, &b_EFlowTrack_Xd);
   fChain->SetBranchAddress("EFlowTrack.Yd", &EFlowTrack_Yd, &b_EFlowTrack_Yd);
   fChain->SetBranchAddress("EFlowTrack.Zd", &EFlowTrack_Zd, &b_EFlowTrack_Zd);
   fChain->SetBranchAddress("EFlowTrack.L", &EFlowTrack_L, &b_EFlowTrack_L);
   fChain->SetBranchAddress("EFlowTrack.D0", &EFlowTrack_D0, &b_EFlowTrack_D0);
   fChain->SetBranchAddress("EFlowTrack.DZ", &EFlowTrack_DZ, &b_EFlowTrack_DZ);
   fChain->SetBranchAddress("EFlowTrack.Nclusters", &EFlowTrack_Nclusters, &b_EFlowTrack_Nclusters);
   fChain->SetBranchAddress("EFlowTrack.dNdx", &EFlowTrack_dNdx, &b_EFlowTrack_dNdx);
   fChain->SetBranchAddress("EFlowTrack.ErrorP", &EFlowTrack_ErrorP, &b_EFlowTrack_ErrorP);
   fChain->SetBranchAddress("EFlowTrack.ErrorPT", &EFlowTrack_ErrorPT, &b_EFlowTrack_ErrorPT);
   fChain->SetBranchAddress("EFlowTrack.ErrorPhi", &EFlowTrack_ErrorPhi, &b_EFlowTrack_ErrorPhi);
   fChain->SetBranchAddress("EFlowTrack.ErrorCtgTheta", &EFlowTrack_ErrorCtgTheta, &b_EFlowTrack_ErrorCtgTheta);
   fChain->SetBranchAddress("EFlowTrack.ErrorT", &EFlowTrack_ErrorT, &b_EFlowTrack_ErrorT);
   fChain->SetBranchAddress("EFlowTrack.ErrorD0", &EFlowTrack_ErrorD0, &b_EFlowTrack_ErrorD0);
   fChain->SetBranchAddress("EFlowTrack.ErrorDZ", &EFlowTrack_ErrorDZ, &b_EFlowTrack_ErrorDZ);
   fChain->SetBranchAddress("EFlowTrack.ErrorC", &EFlowTrack_ErrorC, &b_EFlowTrack_ErrorC);
   fChain->SetBranchAddress("EFlowTrack.ErrorD0Phi", &EFlowTrack_ErrorD0Phi, &b_EFlowTrack_ErrorD0Phi);
   fChain->SetBranchAddress("EFlowTrack.ErrorD0C", &EFlowTrack_ErrorD0C, &b_EFlowTrack_ErrorD0C);
   fChain->SetBranchAddress("EFlowTrack.ErrorD0DZ", &EFlowTrack_ErrorD0DZ, &b_EFlowTrack_ErrorD0DZ);
   fChain->SetBranchAddress("EFlowTrack.ErrorD0CtgTheta", &EFlowTrack_ErrorD0CtgTheta, &b_EFlowTrack_ErrorD0CtgTheta);
   fChain->SetBranchAddress("EFlowTrack.ErrorPhiC", &EFlowTrack_ErrorPhiC, &b_EFlowTrack_ErrorPhiC);
   fChain->SetBranchAddress("EFlowTrack.ErrorPhiDZ", &EFlowTrack_ErrorPhiDZ, &b_EFlowTrack_ErrorPhiDZ);
   fChain->SetBranchAddress("EFlowTrack.ErrorPhiCtgTheta", &EFlowTrack_ErrorPhiCtgTheta, &b_EFlowTrack_ErrorPhiCtgTheta);
   fChain->SetBranchAddress("EFlowTrack.ErrorCDZ", &EFlowTrack_ErrorCDZ, &b_EFlowTrack_ErrorCDZ);
   fChain->SetBranchAddress("EFlowTrack.ErrorCCtgTheta", &EFlowTrack_ErrorCCtgTheta, &b_EFlowTrack_ErrorCCtgTheta);
   fChain->SetBranchAddress("EFlowTrack.ErrorDZCtgTheta", &EFlowTrack_ErrorDZCtgTheta, &b_EFlowTrack_ErrorDZCtgTheta);
   fChain->SetBranchAddress("EFlowTrack.Particle", &EFlowTrack_Particle, &b_EFlowTrack_Particle);
   fChain->SetBranchAddress("EFlowTrack.VertexIndex", &EFlowTrack_VertexIndex, &b_EFlowTrack_VertexIndex);
   fChain->SetBranchAddress("EFlowTrack_size", &EFlowTrack_size, &b_EFlowTrack_size);
   fChain->SetBranchAddress("EFlowPhoton", &EFlowPhoton_, &b_EFlowPhoton_);
   fChain->SetBranchAddress("EFlowPhoton.fUniqueID", EFlowPhoton_fUniqueID, &b_EFlowPhoton_fUniqueID);
   fChain->SetBranchAddress("EFlowPhoton.fBits", EFlowPhoton_fBits, &b_EFlowPhoton_fBits);
   fChain->SetBranchAddress("EFlowPhoton.ET", EFlowPhoton_ET, &b_EFlowPhoton_ET);
   fChain->SetBranchAddress("EFlowPhoton.Eta", EFlowPhoton_Eta, &b_EFlowPhoton_Eta);
   fChain->SetBranchAddress("EFlowPhoton.Phi", EFlowPhoton_Phi, &b_EFlowPhoton_Phi);
   fChain->SetBranchAddress("EFlowPhoton.E", EFlowPhoton_E, &b_EFlowPhoton_E);
   fChain->SetBranchAddress("EFlowPhoton.T", EFlowPhoton_T, &b_EFlowPhoton_T);
   fChain->SetBranchAddress("EFlowPhoton.NTimeHits", EFlowPhoton_NTimeHits, &b_EFlowPhoton_NTimeHits);
   fChain->SetBranchAddress("EFlowPhoton.Eem", EFlowPhoton_Eem, &b_EFlowPhoton_Eem);
   fChain->SetBranchAddress("EFlowPhoton.Ehad", EFlowPhoton_Ehad, &b_EFlowPhoton_Ehad);
   fChain->SetBranchAddress("EFlowPhoton.Edges[4]", EFlowPhoton_Edges, &b_EFlowPhoton_Edges);
   fChain->SetBranchAddress("EFlowPhoton.Particles", EFlowPhoton_Particles, &b_EFlowPhoton_Particles);
   fChain->SetBranchAddress("EFlowPhoton_size", &EFlowPhoton_size, &b_EFlowPhoton_size);
   fChain->SetBranchAddress("EFlowNeutralHadron", &EFlowNeutralHadron_, &b_EFlowNeutralHadron_);
   fChain->SetBranchAddress("EFlowNeutralHadron.fUniqueID", EFlowNeutralHadron_fUniqueID, &b_EFlowNeutralHadron_fUniqueID);
   fChain->SetBranchAddress("EFlowNeutralHadron.fBits", EFlowNeutralHadron_fBits, &b_EFlowNeutralHadron_fBits);
   fChain->SetBranchAddress("EFlowNeutralHadron.ET", EFlowNeutralHadron_ET, &b_EFlowNeutralHadron_ET);
   fChain->SetBranchAddress("EFlowNeutralHadron.Eta", EFlowNeutralHadron_Eta, &b_EFlowNeutralHadron_Eta);
   fChain->SetBranchAddress("EFlowNeutralHadron.Phi", EFlowNeutralHadron_Phi, &b_EFlowNeutralHadron_Phi);
   fChain->SetBranchAddress("EFlowNeutralHadron.E", EFlowNeutralHadron_E, &b_EFlowNeutralHadron_E);
   fChain->SetBranchAddress("EFlowNeutralHadron.T", EFlowNeutralHadron_T, &b_EFlowNeutralHadron_T);
   fChain->SetBranchAddress("EFlowNeutralHadron.NTimeHits", EFlowNeutralHadron_NTimeHits, &b_EFlowNeutralHadron_NTimeHits);
   fChain->SetBranchAddress("EFlowNeutralHadron.Eem", EFlowNeutralHadron_Eem, &b_EFlowNeutralHadron_Eem);
   fChain->SetBranchAddress("EFlowNeutralHadron.Ehad", EFlowNeutralHadron_Ehad, &b_EFlowNeutralHadron_Ehad);
   fChain->SetBranchAddress("EFlowNeutralHadron.Edges[4]", EFlowNeutralHadron_Edges, &b_EFlowNeutralHadron_Edges);
   fChain->SetBranchAddress("EFlowNeutralHadron.Particles", EFlowNeutralHadron_Particles, &b_EFlowNeutralHadron_Particles);
   fChain->SetBranchAddress("EFlowNeutralHadron_size", &EFlowNeutralHadron_size, &b_EFlowNeutralHadron_size);
   fChain->SetBranchAddress("mRICHTrack", &mRICHTrack_, &b_mRICHTrack_);
   fChain->SetBranchAddress("mRICHTrack.fUniqueID", &mRICHTrack_fUniqueID, &b_mRICHTrack_fUniqueID);
   fChain->SetBranchAddress("mRICHTrack.fBits", &mRICHTrack_fBits, &b_mRICHTrack_fBits);
   fChain->SetBranchAddress("mRICHTrack.PID", &mRICHTrack_PID, &b_mRICHTrack_PID);
   fChain->SetBranchAddress("mRICHTrack.Charge", &mRICHTrack_Charge, &b_mRICHTrack_Charge);
   fChain->SetBranchAddress("mRICHTrack.P", &mRICHTrack_P, &b_mRICHTrack_P);
   fChain->SetBranchAddress("mRICHTrack.PT", &mRICHTrack_PT, &b_mRICHTrack_PT);
   fChain->SetBranchAddress("mRICHTrack.Eta", &mRICHTrack_Eta, &b_mRICHTrack_Eta);
   fChain->SetBranchAddress("mRICHTrack.Phi", &mRICHTrack_Phi, &b_mRICHTrack_Phi);
   fChain->SetBranchAddress("mRICHTrack.CtgTheta", &mRICHTrack_CtgTheta, &b_mRICHTrack_CtgTheta);
   fChain->SetBranchAddress("mRICHTrack.C", &mRICHTrack_C, &b_mRICHTrack_C);
   fChain->SetBranchAddress("mRICHTrack.Mass", &mRICHTrack_Mass, &b_mRICHTrack_Mass);
   fChain->SetBranchAddress("mRICHTrack.EtaOuter", &mRICHTrack_EtaOuter, &b_mRICHTrack_EtaOuter);
   fChain->SetBranchAddress("mRICHTrack.PhiOuter", &mRICHTrack_PhiOuter, &b_mRICHTrack_PhiOuter);
   fChain->SetBranchAddress("mRICHTrack.T", &mRICHTrack_T, &b_mRICHTrack_T);
   fChain->SetBranchAddress("mRICHTrack.X", &mRICHTrack_X, &b_mRICHTrack_X);
   fChain->SetBranchAddress("mRICHTrack.Y", &mRICHTrack_Y, &b_mRICHTrack_Y);
   fChain->SetBranchAddress("mRICHTrack.Z", &mRICHTrack_Z, &b_mRICHTrack_Z);
   fChain->SetBranchAddress("mRICHTrack.TOuter", &mRICHTrack_TOuter, &b_mRICHTrack_TOuter);
   fChain->SetBranchAddress("mRICHTrack.XOuter", &mRICHTrack_XOuter, &b_mRICHTrack_XOuter);
   fChain->SetBranchAddress("mRICHTrack.YOuter", &mRICHTrack_YOuter, &b_mRICHTrack_YOuter);
   fChain->SetBranchAddress("mRICHTrack.ZOuter", &mRICHTrack_ZOuter, &b_mRICHTrack_ZOuter);
   fChain->SetBranchAddress("mRICHTrack.Xd", &mRICHTrack_Xd, &b_mRICHTrack_Xd);
   fChain->SetBranchAddress("mRICHTrack.Yd", &mRICHTrack_Yd, &b_mRICHTrack_Yd);
   fChain->SetBranchAddress("mRICHTrack.Zd", &mRICHTrack_Zd, &b_mRICHTrack_Zd);
   fChain->SetBranchAddress("mRICHTrack.L", &mRICHTrack_L, &b_mRICHTrack_L);
   fChain->SetBranchAddress("mRICHTrack.D0", &mRICHTrack_D0, &b_mRICHTrack_D0);
   fChain->SetBranchAddress("mRICHTrack.DZ", &mRICHTrack_DZ, &b_mRICHTrack_DZ);
   fChain->SetBranchAddress("mRICHTrack.Nclusters", &mRICHTrack_Nclusters, &b_mRICHTrack_Nclusters);
   fChain->SetBranchAddress("mRICHTrack.dNdx", &mRICHTrack_dNdx, &b_mRICHTrack_dNdx);
   fChain->SetBranchAddress("mRICHTrack.ErrorP", &mRICHTrack_ErrorP, &b_mRICHTrack_ErrorP);
   fChain->SetBranchAddress("mRICHTrack.ErrorPT", &mRICHTrack_ErrorPT, &b_mRICHTrack_ErrorPT);
   fChain->SetBranchAddress("mRICHTrack.ErrorPhi", &mRICHTrack_ErrorPhi, &b_mRICHTrack_ErrorPhi);
   fChain->SetBranchAddress("mRICHTrack.ErrorCtgTheta", &mRICHTrack_ErrorCtgTheta, &b_mRICHTrack_ErrorCtgTheta);
   fChain->SetBranchAddress("mRICHTrack.ErrorT", &mRICHTrack_ErrorT, &b_mRICHTrack_ErrorT);
   fChain->SetBranchAddress("mRICHTrack.ErrorD0", &mRICHTrack_ErrorD0, &b_mRICHTrack_ErrorD0);
   fChain->SetBranchAddress("mRICHTrack.ErrorDZ", &mRICHTrack_ErrorDZ, &b_mRICHTrack_ErrorDZ);
   fChain->SetBranchAddress("mRICHTrack.ErrorC", &mRICHTrack_ErrorC, &b_mRICHTrack_ErrorC);
   fChain->SetBranchAddress("mRICHTrack.ErrorD0Phi", &mRICHTrack_ErrorD0Phi, &b_mRICHTrack_ErrorD0Phi);
   fChain->SetBranchAddress("mRICHTrack.ErrorD0C", &mRICHTrack_ErrorD0C, &b_mRICHTrack_ErrorD0C);
   fChain->SetBranchAddress("mRICHTrack.ErrorD0DZ", &mRICHTrack_ErrorD0DZ, &b_mRICHTrack_ErrorD0DZ);
   fChain->SetBranchAddress("mRICHTrack.ErrorD0CtgTheta", &mRICHTrack_ErrorD0CtgTheta, &b_mRICHTrack_ErrorD0CtgTheta);
   fChain->SetBranchAddress("mRICHTrack.ErrorPhiC", &mRICHTrack_ErrorPhiC, &b_mRICHTrack_ErrorPhiC);
   fChain->SetBranchAddress("mRICHTrack.ErrorPhiDZ", &mRICHTrack_ErrorPhiDZ, &b_mRICHTrack_ErrorPhiDZ);
   fChain->SetBranchAddress("mRICHTrack.ErrorPhiCtgTheta", &mRICHTrack_ErrorPhiCtgTheta, &b_mRICHTrack_ErrorPhiCtgTheta);
   fChain->SetBranchAddress("mRICHTrack.ErrorCDZ", &mRICHTrack_ErrorCDZ, &b_mRICHTrack_ErrorCDZ);
   fChain->SetBranchAddress("mRICHTrack.ErrorCCtgTheta", &mRICHTrack_ErrorCCtgTheta, &b_mRICHTrack_ErrorCCtgTheta);
   fChain->SetBranchAddress("mRICHTrack.ErrorDZCtgTheta", &mRICHTrack_ErrorDZCtgTheta, &b_mRICHTrack_ErrorDZCtgTheta);
   fChain->SetBranchAddress("mRICHTrack.Particle", &mRICHTrack_Particle, &b_mRICHTrack_Particle);
   fChain->SetBranchAddress("mRICHTrack.VertexIndex", &mRICHTrack_VertexIndex, &b_mRICHTrack_VertexIndex);
   fChain->SetBranchAddress("mRICHTrack_size", &mRICHTrack_size, &b_mRICHTrack_size);
   fChain->SetBranchAddress("barrelDIRCTrack", &barrelDIRCTrack_, &b_barrelDIRCTrack_);
   fChain->SetBranchAddress("barrelDIRCTrack.fUniqueID", &barrelDIRCTrack_fUniqueID, &b_barrelDIRCTrack_fUniqueID);
   fChain->SetBranchAddress("barrelDIRCTrack.fBits", &barrelDIRCTrack_fBits, &b_barrelDIRCTrack_fBits);
   fChain->SetBranchAddress("barrelDIRCTrack.PID", &barrelDIRCTrack_PID, &b_barrelDIRCTrack_PID);
   fChain->SetBranchAddress("barrelDIRCTrack.Charge", &barrelDIRCTrack_Charge, &b_barrelDIRCTrack_Charge);
   fChain->SetBranchAddress("barrelDIRCTrack.P", &barrelDIRCTrack_P, &b_barrelDIRCTrack_P);
   fChain->SetBranchAddress("barrelDIRCTrack.PT", &barrelDIRCTrack_PT, &b_barrelDIRCTrack_PT);
   fChain->SetBranchAddress("barrelDIRCTrack.Eta", &barrelDIRCTrack_Eta, &b_barrelDIRCTrack_Eta);
   fChain->SetBranchAddress("barrelDIRCTrack.Phi", &barrelDIRCTrack_Phi, &b_barrelDIRCTrack_Phi);
   fChain->SetBranchAddress("barrelDIRCTrack.CtgTheta", &barrelDIRCTrack_CtgTheta, &b_barrelDIRCTrack_CtgTheta);
   fChain->SetBranchAddress("barrelDIRCTrack.C", &barrelDIRCTrack_C, &b_barrelDIRCTrack_C);
   fChain->SetBranchAddress("barrelDIRCTrack.Mass", &barrelDIRCTrack_Mass, &b_barrelDIRCTrack_Mass);
   fChain->SetBranchAddress("barrelDIRCTrack.EtaOuter", &barrelDIRCTrack_EtaOuter, &b_barrelDIRCTrack_EtaOuter);
   fChain->SetBranchAddress("barrelDIRCTrack.PhiOuter", &barrelDIRCTrack_PhiOuter, &b_barrelDIRCTrack_PhiOuter);
   fChain->SetBranchAddress("barrelDIRCTrack.T", &barrelDIRCTrack_T, &b_barrelDIRCTrack_T);
   fChain->SetBranchAddress("barrelDIRCTrack.X", &barrelDIRCTrack_X, &b_barrelDIRCTrack_X);
   fChain->SetBranchAddress("barrelDIRCTrack.Y", &barrelDIRCTrack_Y, &b_barrelDIRCTrack_Y);
   fChain->SetBranchAddress("barrelDIRCTrack.Z", &barrelDIRCTrack_Z, &b_barrelDIRCTrack_Z);
   fChain->SetBranchAddress("barrelDIRCTrack.TOuter", &barrelDIRCTrack_TOuter, &b_barrelDIRCTrack_TOuter);
   fChain->SetBranchAddress("barrelDIRCTrack.XOuter", &barrelDIRCTrack_XOuter, &b_barrelDIRCTrack_XOuter);
   fChain->SetBranchAddress("barrelDIRCTrack.YOuter", &barrelDIRCTrack_YOuter, &b_barrelDIRCTrack_YOuter);
   fChain->SetBranchAddress("barrelDIRCTrack.ZOuter", &barrelDIRCTrack_ZOuter, &b_barrelDIRCTrack_ZOuter);
   fChain->SetBranchAddress("barrelDIRCTrack.Xd", &barrelDIRCTrack_Xd, &b_barrelDIRCTrack_Xd);
   fChain->SetBranchAddress("barrelDIRCTrack.Yd", &barrelDIRCTrack_Yd, &b_barrelDIRCTrack_Yd);
   fChain->SetBranchAddress("barrelDIRCTrack.Zd", &barrelDIRCTrack_Zd, &b_barrelDIRCTrack_Zd);
   fChain->SetBranchAddress("barrelDIRCTrack.L", &barrelDIRCTrack_L, &b_barrelDIRCTrack_L);
   fChain->SetBranchAddress("barrelDIRCTrack.D0", &barrelDIRCTrack_D0, &b_barrelDIRCTrack_D0);
   fChain->SetBranchAddress("barrelDIRCTrack.DZ", &barrelDIRCTrack_DZ, &b_barrelDIRCTrack_DZ);
   fChain->SetBranchAddress("barrelDIRCTrack.Nclusters", &barrelDIRCTrack_Nclusters, &b_barrelDIRCTrack_Nclusters);
   fChain->SetBranchAddress("barrelDIRCTrack.dNdx", &barrelDIRCTrack_dNdx, &b_barrelDIRCTrack_dNdx);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorP", &barrelDIRCTrack_ErrorP, &b_barrelDIRCTrack_ErrorP);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorPT", &barrelDIRCTrack_ErrorPT, &b_barrelDIRCTrack_ErrorPT);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorPhi", &barrelDIRCTrack_ErrorPhi, &b_barrelDIRCTrack_ErrorPhi);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorCtgTheta", &barrelDIRCTrack_ErrorCtgTheta, &b_barrelDIRCTrack_ErrorCtgTheta);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorT", &barrelDIRCTrack_ErrorT, &b_barrelDIRCTrack_ErrorT);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorD0", &barrelDIRCTrack_ErrorD0, &b_barrelDIRCTrack_ErrorD0);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorDZ", &barrelDIRCTrack_ErrorDZ, &b_barrelDIRCTrack_ErrorDZ);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorC", &barrelDIRCTrack_ErrorC, &b_barrelDIRCTrack_ErrorC);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorD0Phi", &barrelDIRCTrack_ErrorD0Phi, &b_barrelDIRCTrack_ErrorD0Phi);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorD0C", &barrelDIRCTrack_ErrorD0C, &b_barrelDIRCTrack_ErrorD0C);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorD0DZ", &barrelDIRCTrack_ErrorD0DZ, &b_barrelDIRCTrack_ErrorD0DZ);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorD0CtgTheta", &barrelDIRCTrack_ErrorD0CtgTheta, &b_barrelDIRCTrack_ErrorD0CtgTheta);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorPhiC", &barrelDIRCTrack_ErrorPhiC, &b_barrelDIRCTrack_ErrorPhiC);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorPhiDZ", &barrelDIRCTrack_ErrorPhiDZ, &b_barrelDIRCTrack_ErrorPhiDZ);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorPhiCtgTheta", &barrelDIRCTrack_ErrorPhiCtgTheta, &b_barrelDIRCTrack_ErrorPhiCtgTheta);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorCDZ", &barrelDIRCTrack_ErrorCDZ, &b_barrelDIRCTrack_ErrorCDZ);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorCCtgTheta", &barrelDIRCTrack_ErrorCCtgTheta, &b_barrelDIRCTrack_ErrorCCtgTheta);
   fChain->SetBranchAddress("barrelDIRCTrack.ErrorDZCtgTheta", &barrelDIRCTrack_ErrorDZCtgTheta, &b_barrelDIRCTrack_ErrorDZCtgTheta);
   fChain->SetBranchAddress("barrelDIRCTrack.Particle", &barrelDIRCTrack_Particle, &b_barrelDIRCTrack_Particle);
   fChain->SetBranchAddress("barrelDIRCTrack.VertexIndex", &barrelDIRCTrack_VertexIndex, &b_barrelDIRCTrack_VertexIndex);
   fChain->SetBranchAddress("barrelDIRCTrack_size", &barrelDIRCTrack_size, &b_barrelDIRCTrack_size);
   fChain->SetBranchAddress("dualRICHagTrack", &dualRICHagTrack_, &b_dualRICHagTrack_);
   fChain->SetBranchAddress("dualRICHagTrack.fUniqueID", &dualRICHagTrack_fUniqueID, &b_dualRICHagTrack_fUniqueID);
   fChain->SetBranchAddress("dualRICHagTrack.fBits", &dualRICHagTrack_fBits, &b_dualRICHagTrack_fBits);
   fChain->SetBranchAddress("dualRICHagTrack.PID", &dualRICHagTrack_PID, &b_dualRICHagTrack_PID);
   fChain->SetBranchAddress("dualRICHagTrack.Charge", &dualRICHagTrack_Charge, &b_dualRICHagTrack_Charge);
   fChain->SetBranchAddress("dualRICHagTrack.P", &dualRICHagTrack_P, &b_dualRICHagTrack_P);
   fChain->SetBranchAddress("dualRICHagTrack.PT", &dualRICHagTrack_PT, &b_dualRICHagTrack_PT);
   fChain->SetBranchAddress("dualRICHagTrack.Eta", &dualRICHagTrack_Eta, &b_dualRICHagTrack_Eta);
   fChain->SetBranchAddress("dualRICHagTrack.Phi", &dualRICHagTrack_Phi, &b_dualRICHagTrack_Phi);
   fChain->SetBranchAddress("dualRICHagTrack.CtgTheta", &dualRICHagTrack_CtgTheta, &b_dualRICHagTrack_CtgTheta);
   fChain->SetBranchAddress("dualRICHagTrack.C", &dualRICHagTrack_C, &b_dualRICHagTrack_C);
   fChain->SetBranchAddress("dualRICHagTrack.Mass", &dualRICHagTrack_Mass, &b_dualRICHagTrack_Mass);
   fChain->SetBranchAddress("dualRICHagTrack.EtaOuter", &dualRICHagTrack_EtaOuter, &b_dualRICHagTrack_EtaOuter);
   fChain->SetBranchAddress("dualRICHagTrack.PhiOuter", &dualRICHagTrack_PhiOuter, &b_dualRICHagTrack_PhiOuter);
   fChain->SetBranchAddress("dualRICHagTrack.T", &dualRICHagTrack_T, &b_dualRICHagTrack_T);
   fChain->SetBranchAddress("dualRICHagTrack.X", &dualRICHagTrack_X, &b_dualRICHagTrack_X);
   fChain->SetBranchAddress("dualRICHagTrack.Y", &dualRICHagTrack_Y, &b_dualRICHagTrack_Y);
   fChain->SetBranchAddress("dualRICHagTrack.Z", &dualRICHagTrack_Z, &b_dualRICHagTrack_Z);
   fChain->SetBranchAddress("dualRICHagTrack.TOuter", &dualRICHagTrack_TOuter, &b_dualRICHagTrack_TOuter);
   fChain->SetBranchAddress("dualRICHagTrack.XOuter", &dualRICHagTrack_XOuter, &b_dualRICHagTrack_XOuter);
   fChain->SetBranchAddress("dualRICHagTrack.YOuter", &dualRICHagTrack_YOuter, &b_dualRICHagTrack_YOuter);
   fChain->SetBranchAddress("dualRICHagTrack.ZOuter", &dualRICHagTrack_ZOuter, &b_dualRICHagTrack_ZOuter);
   fChain->SetBranchAddress("dualRICHagTrack.Xd", &dualRICHagTrack_Xd, &b_dualRICHagTrack_Xd);
   fChain->SetBranchAddress("dualRICHagTrack.Yd", &dualRICHagTrack_Yd, &b_dualRICHagTrack_Yd);
   fChain->SetBranchAddress("dualRICHagTrack.Zd", &dualRICHagTrack_Zd, &b_dualRICHagTrack_Zd);
   fChain->SetBranchAddress("dualRICHagTrack.L", &dualRICHagTrack_L, &b_dualRICHagTrack_L);
   fChain->SetBranchAddress("dualRICHagTrack.D0", &dualRICHagTrack_D0, &b_dualRICHagTrack_D0);
   fChain->SetBranchAddress("dualRICHagTrack.DZ", &dualRICHagTrack_DZ, &b_dualRICHagTrack_DZ);
   fChain->SetBranchAddress("dualRICHagTrack.Nclusters", &dualRICHagTrack_Nclusters, &b_dualRICHagTrack_Nclusters);
   fChain->SetBranchAddress("dualRICHagTrack.dNdx", &dualRICHagTrack_dNdx, &b_dualRICHagTrack_dNdx);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorP", &dualRICHagTrack_ErrorP, &b_dualRICHagTrack_ErrorP);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorPT", &dualRICHagTrack_ErrorPT, &b_dualRICHagTrack_ErrorPT);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorPhi", &dualRICHagTrack_ErrorPhi, &b_dualRICHagTrack_ErrorPhi);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorCtgTheta", &dualRICHagTrack_ErrorCtgTheta, &b_dualRICHagTrack_ErrorCtgTheta);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorT", &dualRICHagTrack_ErrorT, &b_dualRICHagTrack_ErrorT);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorD0", &dualRICHagTrack_ErrorD0, &b_dualRICHagTrack_ErrorD0);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorDZ", &dualRICHagTrack_ErrorDZ, &b_dualRICHagTrack_ErrorDZ);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorC", &dualRICHagTrack_ErrorC, &b_dualRICHagTrack_ErrorC);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorD0Phi", &dualRICHagTrack_ErrorD0Phi, &b_dualRICHagTrack_ErrorD0Phi);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorD0C", &dualRICHagTrack_ErrorD0C, &b_dualRICHagTrack_ErrorD0C);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorD0DZ", &dualRICHagTrack_ErrorD0DZ, &b_dualRICHagTrack_ErrorD0DZ);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorD0CtgTheta", &dualRICHagTrack_ErrorD0CtgTheta, &b_dualRICHagTrack_ErrorD0CtgTheta);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorPhiC", &dualRICHagTrack_ErrorPhiC, &b_dualRICHagTrack_ErrorPhiC);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorPhiDZ", &dualRICHagTrack_ErrorPhiDZ, &b_dualRICHagTrack_ErrorPhiDZ);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorPhiCtgTheta", &dualRICHagTrack_ErrorPhiCtgTheta, &b_dualRICHagTrack_ErrorPhiCtgTheta);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorCDZ", &dualRICHagTrack_ErrorCDZ, &b_dualRICHagTrack_ErrorCDZ);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorCCtgTheta", &dualRICHagTrack_ErrorCCtgTheta, &b_dualRICHagTrack_ErrorCCtgTheta);
   fChain->SetBranchAddress("dualRICHagTrack.ErrorDZCtgTheta", &dualRICHagTrack_ErrorDZCtgTheta, &b_dualRICHagTrack_ErrorDZCtgTheta);
   fChain->SetBranchAddress("dualRICHagTrack.Particle", &dualRICHagTrack_Particle, &b_dualRICHagTrack_Particle);
   fChain->SetBranchAddress("dualRICHagTrack.VertexIndex", &dualRICHagTrack_VertexIndex, &b_dualRICHagTrack_VertexIndex);
   fChain->SetBranchAddress("dualRICHagTrack_size", &dualRICHagTrack_size, &b_dualRICHagTrack_size);
   fChain->SetBranchAddress("dualRICHcfTrack", &dualRICHcfTrack_, &b_dualRICHcfTrack_);
   fChain->SetBranchAddress("dualRICHcfTrack.fUniqueID", &dualRICHcfTrack_fUniqueID, &b_dualRICHcfTrack_fUniqueID);
   fChain->SetBranchAddress("dualRICHcfTrack.fBits", &dualRICHcfTrack_fBits, &b_dualRICHcfTrack_fBits);
   fChain->SetBranchAddress("dualRICHcfTrack.PID", &dualRICHcfTrack_PID, &b_dualRICHcfTrack_PID);
   fChain->SetBranchAddress("dualRICHcfTrack.Charge", &dualRICHcfTrack_Charge, &b_dualRICHcfTrack_Charge);
   fChain->SetBranchAddress("dualRICHcfTrack.P", &dualRICHcfTrack_P, &b_dualRICHcfTrack_P);
   fChain->SetBranchAddress("dualRICHcfTrack.PT", &dualRICHcfTrack_PT, &b_dualRICHcfTrack_PT);
   fChain->SetBranchAddress("dualRICHcfTrack.Eta", &dualRICHcfTrack_Eta, &b_dualRICHcfTrack_Eta);
   fChain->SetBranchAddress("dualRICHcfTrack.Phi", &dualRICHcfTrack_Phi, &b_dualRICHcfTrack_Phi);
   fChain->SetBranchAddress("dualRICHcfTrack.CtgTheta", &dualRICHcfTrack_CtgTheta, &b_dualRICHcfTrack_CtgTheta);
   fChain->SetBranchAddress("dualRICHcfTrack.C", &dualRICHcfTrack_C, &b_dualRICHcfTrack_C);
   fChain->SetBranchAddress("dualRICHcfTrack.Mass", &dualRICHcfTrack_Mass, &b_dualRICHcfTrack_Mass);
   fChain->SetBranchAddress("dualRICHcfTrack.EtaOuter", &dualRICHcfTrack_EtaOuter, &b_dualRICHcfTrack_EtaOuter);
   fChain->SetBranchAddress("dualRICHcfTrack.PhiOuter", &dualRICHcfTrack_PhiOuter, &b_dualRICHcfTrack_PhiOuter);
   fChain->SetBranchAddress("dualRICHcfTrack.T", &dualRICHcfTrack_T, &b_dualRICHcfTrack_T);
   fChain->SetBranchAddress("dualRICHcfTrack.X", &dualRICHcfTrack_X, &b_dualRICHcfTrack_X);
   fChain->SetBranchAddress("dualRICHcfTrack.Y", &dualRICHcfTrack_Y, &b_dualRICHcfTrack_Y);
   fChain->SetBranchAddress("dualRICHcfTrack.Z", &dualRICHcfTrack_Z, &b_dualRICHcfTrack_Z);
   fChain->SetBranchAddress("dualRICHcfTrack.TOuter", &dualRICHcfTrack_TOuter, &b_dualRICHcfTrack_TOuter);
   fChain->SetBranchAddress("dualRICHcfTrack.XOuter", &dualRICHcfTrack_XOuter, &b_dualRICHcfTrack_XOuter);
   fChain->SetBranchAddress("dualRICHcfTrack.YOuter", &dualRICHcfTrack_YOuter, &b_dualRICHcfTrack_YOuter);
   fChain->SetBranchAddress("dualRICHcfTrack.ZOuter", &dualRICHcfTrack_ZOuter, &b_dualRICHcfTrack_ZOuter);
   fChain->SetBranchAddress("dualRICHcfTrack.Xd", &dualRICHcfTrack_Xd, &b_dualRICHcfTrack_Xd);
   fChain->SetBranchAddress("dualRICHcfTrack.Yd", &dualRICHcfTrack_Yd, &b_dualRICHcfTrack_Yd);
   fChain->SetBranchAddress("dualRICHcfTrack.Zd", &dualRICHcfTrack_Zd, &b_dualRICHcfTrack_Zd);
   fChain->SetBranchAddress("dualRICHcfTrack.L", &dualRICHcfTrack_L, &b_dualRICHcfTrack_L);
   fChain->SetBranchAddress("dualRICHcfTrack.D0", &dualRICHcfTrack_D0, &b_dualRICHcfTrack_D0);
   fChain->SetBranchAddress("dualRICHcfTrack.DZ", &dualRICHcfTrack_DZ, &b_dualRICHcfTrack_DZ);
   fChain->SetBranchAddress("dualRICHcfTrack.Nclusters", &dualRICHcfTrack_Nclusters, &b_dualRICHcfTrack_Nclusters);
   fChain->SetBranchAddress("dualRICHcfTrack.dNdx", &dualRICHcfTrack_dNdx, &b_dualRICHcfTrack_dNdx);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorP", &dualRICHcfTrack_ErrorP, &b_dualRICHcfTrack_ErrorP);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorPT", &dualRICHcfTrack_ErrorPT, &b_dualRICHcfTrack_ErrorPT);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorPhi", &dualRICHcfTrack_ErrorPhi, &b_dualRICHcfTrack_ErrorPhi);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorCtgTheta", &dualRICHcfTrack_ErrorCtgTheta, &b_dualRICHcfTrack_ErrorCtgTheta);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorT", &dualRICHcfTrack_ErrorT, &b_dualRICHcfTrack_ErrorT);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorD0", &dualRICHcfTrack_ErrorD0, &b_dualRICHcfTrack_ErrorD0);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorDZ", &dualRICHcfTrack_ErrorDZ, &b_dualRICHcfTrack_ErrorDZ);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorC", &dualRICHcfTrack_ErrorC, &b_dualRICHcfTrack_ErrorC);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorD0Phi", &dualRICHcfTrack_ErrorD0Phi, &b_dualRICHcfTrack_ErrorD0Phi);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorD0C", &dualRICHcfTrack_ErrorD0C, &b_dualRICHcfTrack_ErrorD0C);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorD0DZ", &dualRICHcfTrack_ErrorD0DZ, &b_dualRICHcfTrack_ErrorD0DZ);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorD0CtgTheta", &dualRICHcfTrack_ErrorD0CtgTheta, &b_dualRICHcfTrack_ErrorD0CtgTheta);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorPhiC", &dualRICHcfTrack_ErrorPhiC, &b_dualRICHcfTrack_ErrorPhiC);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorPhiDZ", &dualRICHcfTrack_ErrorPhiDZ, &b_dualRICHcfTrack_ErrorPhiDZ);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorPhiCtgTheta", &dualRICHcfTrack_ErrorPhiCtgTheta, &b_dualRICHcfTrack_ErrorPhiCtgTheta);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorCDZ", &dualRICHcfTrack_ErrorCDZ, &b_dualRICHcfTrack_ErrorCDZ);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorCCtgTheta", &dualRICHcfTrack_ErrorCCtgTheta, &b_dualRICHcfTrack_ErrorCCtgTheta);
   fChain->SetBranchAddress("dualRICHcfTrack.ErrorDZCtgTheta", &dualRICHcfTrack_ErrorDZCtgTheta, &b_dualRICHcfTrack_ErrorDZCtgTheta);
   fChain->SetBranchAddress("dualRICHcfTrack.Particle", &dualRICHcfTrack_Particle, &b_dualRICHcfTrack_Particle);
   fChain->SetBranchAddress("dualRICHcfTrack.VertexIndex", &dualRICHcfTrack_VertexIndex, &b_dualRICHcfTrack_VertexIndex);
   fChain->SetBranchAddress("dualRICHcfTrack_size", &dualRICHcfTrack_size, &b_dualRICHcfTrack_size);
   fChain->SetBranchAddress("GenJet", &GenJet_, &b_GenJet_);
   fChain->SetBranchAddress("GenJet.fUniqueID", &GenJet_fUniqueID, &b_GenJet_fUniqueID);
   fChain->SetBranchAddress("GenJet.fBits", &GenJet_fBits, &b_GenJet_fBits);
   fChain->SetBranchAddress("GenJet.PT", &GenJet_PT, &b_GenJet_PT);
   fChain->SetBranchAddress("GenJet.Eta", &GenJet_Eta, &b_GenJet_Eta);
   fChain->SetBranchAddress("GenJet.Phi", &GenJet_Phi, &b_GenJet_Phi);
   fChain->SetBranchAddress("GenJet.T", &GenJet_T, &b_GenJet_T);
   fChain->SetBranchAddress("GenJet.Mass", &GenJet_Mass, &b_GenJet_Mass);
   fChain->SetBranchAddress("GenJet.DeltaEta", &GenJet_DeltaEta, &b_GenJet_DeltaEta);
   fChain->SetBranchAddress("GenJet.DeltaPhi", &GenJet_DeltaPhi, &b_GenJet_DeltaPhi);
   fChain->SetBranchAddress("GenJet.Flavor", &GenJet_Flavor, &b_GenJet_Flavor);
   fChain->SetBranchAddress("GenJet.FlavorAlgo", &GenJet_FlavorAlgo, &b_GenJet_FlavorAlgo);
   fChain->SetBranchAddress("GenJet.FlavorPhys", &GenJet_FlavorPhys, &b_GenJet_FlavorPhys);
   fChain->SetBranchAddress("GenJet.BTag", &GenJet_BTag, &b_GenJet_BTag);
   fChain->SetBranchAddress("GenJet.BTagAlgo", &GenJet_BTagAlgo, &b_GenJet_BTagAlgo);
   fChain->SetBranchAddress("GenJet.BTagPhys", &GenJet_BTagPhys, &b_GenJet_BTagPhys);
   fChain->SetBranchAddress("GenJet.TauTag", &GenJet_TauTag, &b_GenJet_TauTag);
   fChain->SetBranchAddress("GenJet.TauWeight", &GenJet_TauWeight, &b_GenJet_TauWeight);
   fChain->SetBranchAddress("GenJet.Charge", &GenJet_Charge, &b_GenJet_Charge);
   fChain->SetBranchAddress("GenJet.EhadOverEem", &GenJet_EhadOverEem, &b_GenJet_EhadOverEem);
   fChain->SetBranchAddress("GenJet.NCharged", &GenJet_NCharged, &b_GenJet_NCharged);
   fChain->SetBranchAddress("GenJet.NNeutrals", &GenJet_NNeutrals, &b_GenJet_NNeutrals);
   fChain->SetBranchAddress("GenJet.NeutralEnergyFraction", &GenJet_NeutralEnergyFraction, &b_GenJet_NeutralEnergyFraction);
   fChain->SetBranchAddress("GenJet.ChargedEnergyFraction", &GenJet_ChargedEnergyFraction, &b_GenJet_ChargedEnergyFraction);
   fChain->SetBranchAddress("GenJet.Beta", &GenJet_Beta, &b_GenJet_Beta);
   fChain->SetBranchAddress("GenJet.BetaStar", &GenJet_BetaStar, &b_GenJet_BetaStar);
   fChain->SetBranchAddress("GenJet.MeanSqDeltaR", &GenJet_MeanSqDeltaR, &b_GenJet_MeanSqDeltaR);
   fChain->SetBranchAddress("GenJet.PTD", &GenJet_PTD, &b_GenJet_PTD);
   fChain->SetBranchAddress("GenJet.FracPt[5]", &GenJet_FracPt, &b_GenJet_FracPt);
   fChain->SetBranchAddress("GenJet.Tau[5]", &GenJet_Tau, &b_GenJet_Tau);
   fChain->SetBranchAddress("GenJet.SoftDroppedJet", &GenJet_SoftDroppedJet, &b_GenJet_SoftDroppedJet);
   fChain->SetBranchAddress("GenJet.SoftDroppedSubJet1", &GenJet_SoftDroppedSubJet1, &b_GenJet_SoftDroppedSubJet1);
   fChain->SetBranchAddress("GenJet.SoftDroppedSubJet2", &GenJet_SoftDroppedSubJet2, &b_GenJet_SoftDroppedSubJet2);
   fChain->SetBranchAddress("GenJet.TrimmedP4[5]", &GenJet_TrimmedP4, &b_GenJet_TrimmedP4);
   fChain->SetBranchAddress("GenJet.PrunedP4[5]", &GenJet_PrunedP4, &b_GenJet_PrunedP4);
   fChain->SetBranchAddress("GenJet.SoftDroppedP4[5]", &GenJet_SoftDroppedP4, &b_GenJet_SoftDroppedP4);
   fChain->SetBranchAddress("GenJet.NSubJetsTrimmed", &GenJet_NSubJetsTrimmed, &b_GenJet_NSubJetsTrimmed);
   fChain->SetBranchAddress("GenJet.NSubJetsPruned", &GenJet_NSubJetsPruned, &b_GenJet_NSubJetsPruned);
   fChain->SetBranchAddress("GenJet.NSubJetsSoftDropped", &GenJet_NSubJetsSoftDropped, &b_GenJet_NSubJetsSoftDropped);
   fChain->SetBranchAddress("GenJet.ExclYmerge23", &GenJet_ExclYmerge23, &b_GenJet_ExclYmerge23);
   fChain->SetBranchAddress("GenJet.ExclYmerge34", &GenJet_ExclYmerge34, &b_GenJet_ExclYmerge34);
   fChain->SetBranchAddress("GenJet.ExclYmerge45", &GenJet_ExclYmerge45, &b_GenJet_ExclYmerge45);
   fChain->SetBranchAddress("GenJet.ExclYmerge56", &GenJet_ExclYmerge56, &b_GenJet_ExclYmerge56);
   fChain->SetBranchAddress("GenJet.Constituents", &GenJet_Constituents, &b_GenJet_Constituents);
   fChain->SetBranchAddress("GenJet.Particles", &GenJet_Particles, &b_GenJet_Particles);
   fChain->SetBranchAddress("GenJet.Area", &GenJet_Area, &b_GenJet_Area);
   fChain->SetBranchAddress("GenJet_size", &GenJet_size, &b_GenJet_size);
   fChain->SetBranchAddress("GenMissingET", &GenMissingET_, &b_GenMissingET_);
   fChain->SetBranchAddress("GenMissingET.fUniqueID", GenMissingET_fUniqueID, &b_GenMissingET_fUniqueID);
   fChain->SetBranchAddress("GenMissingET.fBits", GenMissingET_fBits, &b_GenMissingET_fBits);
   fChain->SetBranchAddress("GenMissingET.MET", GenMissingET_MET, &b_GenMissingET_MET);
   fChain->SetBranchAddress("GenMissingET.Eta", GenMissingET_Eta, &b_GenMissingET_Eta);
   fChain->SetBranchAddress("GenMissingET.Phi", GenMissingET_Phi, &b_GenMissingET_Phi);
   fChain->SetBranchAddress("GenMissingET_size", &GenMissingET_size, &b_GenMissingET_size);
   fChain->SetBranchAddress("Jet", &Jet_, &b_Jet_);
   fChain->SetBranchAddress("Jet.fUniqueID", &Jet_fUniqueID, &b_Jet_fUniqueID);
   fChain->SetBranchAddress("Jet.fBits", &Jet_fBits, &b_Jet_fBits);
   fChain->SetBranchAddress("Jet.PT", &Jet_PT, &b_Jet_PT);
   fChain->SetBranchAddress("Jet.Eta", &Jet_Eta, &b_Jet_Eta);
   fChain->SetBranchAddress("Jet.Phi", &Jet_Phi, &b_Jet_Phi);
   fChain->SetBranchAddress("Jet.T", &Jet_T, &b_Jet_T);
   fChain->SetBranchAddress("Jet.Mass", &Jet_Mass, &b_Jet_Mass);
   fChain->SetBranchAddress("Jet.DeltaEta", &Jet_DeltaEta, &b_Jet_DeltaEta);
   fChain->SetBranchAddress("Jet.DeltaPhi", &Jet_DeltaPhi, &b_Jet_DeltaPhi);
   fChain->SetBranchAddress("Jet.Flavor", &Jet_Flavor, &b_Jet_Flavor);
   fChain->SetBranchAddress("Jet.FlavorAlgo", &Jet_FlavorAlgo, &b_Jet_FlavorAlgo);
   fChain->SetBranchAddress("Jet.FlavorPhys", &Jet_FlavorPhys, &b_Jet_FlavorPhys);
   fChain->SetBranchAddress("Jet.BTag", &Jet_BTag, &b_Jet_BTag);
   fChain->SetBranchAddress("Jet.BTagAlgo", &Jet_BTagAlgo, &b_Jet_BTagAlgo);
   fChain->SetBranchAddress("Jet.BTagPhys", &Jet_BTagPhys, &b_Jet_BTagPhys);
   fChain->SetBranchAddress("Jet.TauTag", &Jet_TauTag, &b_Jet_TauTag);
   fChain->SetBranchAddress("Jet.TauWeight", &Jet_TauWeight, &b_Jet_TauWeight);
   fChain->SetBranchAddress("Jet.Charge", &Jet_Charge, &b_Jet_Charge);
   fChain->SetBranchAddress("Jet.EhadOverEem", &Jet_EhadOverEem, &b_Jet_EhadOverEem);
   fChain->SetBranchAddress("Jet.NCharged", &Jet_NCharged, &b_Jet_NCharged);
   fChain->SetBranchAddress("Jet.NNeutrals", &Jet_NNeutrals, &b_Jet_NNeutrals);
   fChain->SetBranchAddress("Jet.NeutralEnergyFraction", &Jet_NeutralEnergyFraction, &b_Jet_NeutralEnergyFraction);
   fChain->SetBranchAddress("Jet.ChargedEnergyFraction", &Jet_ChargedEnergyFraction, &b_Jet_ChargedEnergyFraction);
   fChain->SetBranchAddress("Jet.Beta", &Jet_Beta, &b_Jet_Beta);
   fChain->SetBranchAddress("Jet.BetaStar", &Jet_BetaStar, &b_Jet_BetaStar);
   fChain->SetBranchAddress("Jet.MeanSqDeltaR", &Jet_MeanSqDeltaR, &b_Jet_MeanSqDeltaR);
   fChain->SetBranchAddress("Jet.PTD", &Jet_PTD, &b_Jet_PTD);
   fChain->SetBranchAddress("Jet.FracPt[5]", &Jet_FracPt, &b_Jet_FracPt);
   fChain->SetBranchAddress("Jet.Tau[5]", &Jet_Tau, &b_Jet_Tau);
   fChain->SetBranchAddress("Jet.SoftDroppedJet", &Jet_SoftDroppedJet, &b_Jet_SoftDroppedJet);
   fChain->SetBranchAddress("Jet.SoftDroppedSubJet1", &Jet_SoftDroppedSubJet1, &b_Jet_SoftDroppedSubJet1);
   fChain->SetBranchAddress("Jet.SoftDroppedSubJet2", &Jet_SoftDroppedSubJet2, &b_Jet_SoftDroppedSubJet2);
   fChain->SetBranchAddress("Jet.TrimmedP4[5]", &Jet_TrimmedP4, &b_Jet_TrimmedP4);
   fChain->SetBranchAddress("Jet.PrunedP4[5]", &Jet_PrunedP4, &b_Jet_PrunedP4);
   fChain->SetBranchAddress("Jet.SoftDroppedP4[5]", &Jet_SoftDroppedP4, &b_Jet_SoftDroppedP4);
   fChain->SetBranchAddress("Jet.NSubJetsTrimmed", &Jet_NSubJetsTrimmed, &b_Jet_NSubJetsTrimmed);
   fChain->SetBranchAddress("Jet.NSubJetsPruned", &Jet_NSubJetsPruned, &b_Jet_NSubJetsPruned);
   fChain->SetBranchAddress("Jet.NSubJetsSoftDropped", &Jet_NSubJetsSoftDropped, &b_Jet_NSubJetsSoftDropped);
   fChain->SetBranchAddress("Jet.ExclYmerge23", &Jet_ExclYmerge23, &b_Jet_ExclYmerge23);
   fChain->SetBranchAddress("Jet.ExclYmerge34", &Jet_ExclYmerge34, &b_Jet_ExclYmerge34);
   fChain->SetBranchAddress("Jet.ExclYmerge45", &Jet_ExclYmerge45, &b_Jet_ExclYmerge45);
   fChain->SetBranchAddress("Jet.ExclYmerge56", &Jet_ExclYmerge56, &b_Jet_ExclYmerge56);
   fChain->SetBranchAddress("Jet.Constituents", &Jet_Constituents, &b_Jet_Constituents);
   fChain->SetBranchAddress("Jet.Particles", &Jet_Particles, &b_Jet_Particles);
   fChain->SetBranchAddress("Jet.Area", &Jet_Area, &b_Jet_Area);
   fChain->SetBranchAddress("Jet_size", &Jet_size, &b_Jet_size);
   fChain->SetBranchAddress("Electron", &Electron_, &b_Electron_);
   fChain->SetBranchAddress("Electron.fUniqueID", &Electron_fUniqueID, &b_Electron_fUniqueID);
   fChain->SetBranchAddress("Electron.fBits", &Electron_fBits, &b_Electron_fBits);
   fChain->SetBranchAddress("Electron.PT", &Electron_PT, &b_Electron_PT);
   fChain->SetBranchAddress("Electron.Eta", &Electron_Eta, &b_Electron_Eta);
   fChain->SetBranchAddress("Electron.Phi", &Electron_Phi, &b_Electron_Phi);
   fChain->SetBranchAddress("Electron.T", &Electron_T, &b_Electron_T);
   fChain->SetBranchAddress("Electron.Charge", &Electron_Charge, &b_Electron_Charge);
   fChain->SetBranchAddress("Electron.EhadOverEem", &Electron_EhadOverEem, &b_Electron_EhadOverEem);
   fChain->SetBranchAddress("Electron.Particle", &Electron_Particle, &b_Electron_Particle);
   fChain->SetBranchAddress("Electron.IsolationVar", &Electron_IsolationVar, &b_Electron_IsolationVar);
   fChain->SetBranchAddress("Electron.IsolationVarRhoCorr", &Electron_IsolationVarRhoCorr, &b_Electron_IsolationVarRhoCorr);
   fChain->SetBranchAddress("Electron.SumPtCharged", &Electron_SumPtCharged, &b_Electron_SumPtCharged);
   fChain->SetBranchAddress("Electron.SumPtNeutral", &Electron_SumPtNeutral, &b_Electron_SumPtNeutral);
   fChain->SetBranchAddress("Electron.SumPtChargedPU", &Electron_SumPtChargedPU, &b_Electron_SumPtChargedPU);
   fChain->SetBranchAddress("Electron.SumPt", &Electron_SumPt, &b_Electron_SumPt);
   fChain->SetBranchAddress("Electron.D0", &Electron_D0, &b_Electron_D0);
   fChain->SetBranchAddress("Electron.DZ", &Electron_DZ, &b_Electron_DZ);
   fChain->SetBranchAddress("Electron.ErrorD0", &Electron_ErrorD0, &b_Electron_ErrorD0);
   fChain->SetBranchAddress("Electron.ErrorDZ", &Electron_ErrorDZ, &b_Electron_ErrorDZ);
   fChain->SetBranchAddress("Electron_size", &Electron_size, &b_Electron_size);
   fChain->SetBranchAddress("Photon", &Photon_, &b_Photon_);
   fChain->SetBranchAddress("Photon.fUniqueID", Photon_fUniqueID, &b_Photon_fUniqueID);
   fChain->SetBranchAddress("Photon.fBits", Photon_fBits, &b_Photon_fBits);
   fChain->SetBranchAddress("Photon.PT", Photon_PT, &b_Photon_PT);
   fChain->SetBranchAddress("Photon.Eta", Photon_Eta, &b_Photon_Eta);
   fChain->SetBranchAddress("Photon.Phi", Photon_Phi, &b_Photon_Phi);
   fChain->SetBranchAddress("Photon.E", Photon_E, &b_Photon_E);
   fChain->SetBranchAddress("Photon.T", Photon_T, &b_Photon_T);
   fChain->SetBranchAddress("Photon.EhadOverEem", Photon_EhadOverEem, &b_Photon_EhadOverEem);
   fChain->SetBranchAddress("Photon.Particles", Photon_Particles, &b_Photon_Particles);
   fChain->SetBranchAddress("Photon.IsolationVar", Photon_IsolationVar, &b_Photon_IsolationVar);
   fChain->SetBranchAddress("Photon.IsolationVarRhoCorr", Photon_IsolationVarRhoCorr, &b_Photon_IsolationVarRhoCorr);
   fChain->SetBranchAddress("Photon.SumPtCharged", Photon_SumPtCharged, &b_Photon_SumPtCharged);
   fChain->SetBranchAddress("Photon.SumPtNeutral", Photon_SumPtNeutral, &b_Photon_SumPtNeutral);
   fChain->SetBranchAddress("Photon.SumPtChargedPU", Photon_SumPtChargedPU, &b_Photon_SumPtChargedPU);
   fChain->SetBranchAddress("Photon.SumPt", Photon_SumPt, &b_Photon_SumPt);
   fChain->SetBranchAddress("Photon.Status", Photon_Status, &b_Photon_Status);
   fChain->SetBranchAddress("Photon_size", &Photon_size, &b_Photon_size);
   fChain->SetBranchAddress("MissingET", &MissingET_, &b_MissingET_);
   fChain->SetBranchAddress("MissingET.fUniqueID", MissingET_fUniqueID, &b_MissingET_fUniqueID);
   fChain->SetBranchAddress("MissingET.fBits", MissingET_fBits, &b_MissingET_fBits);
   fChain->SetBranchAddress("MissingET.MET", MissingET_MET, &b_MissingET_MET);
   fChain->SetBranchAddress("MissingET.Eta", MissingET_Eta, &b_MissingET_Eta);
   fChain->SetBranchAddress("MissingET.Phi", MissingET_Phi, &b_MissingET_Phi);
   fChain->SetBranchAddress("MissingET_size", &MissingET_size, &b_MissingET_size);
   fChain->SetBranchAddress("ScalarHT", &ScalarHT_, &b_ScalarHT_);
   fChain->SetBranchAddress("ScalarHT.fUniqueID", ScalarHT_fUniqueID, &b_ScalarHT_fUniqueID);
   fChain->SetBranchAddress("ScalarHT.fBits", ScalarHT_fBits, &b_ScalarHT_fBits);
   fChain->SetBranchAddress("ScalarHT.HT", ScalarHT_HT, &b_ScalarHT_HT);
   fChain->SetBranchAddress("ScalarHT_size", &ScalarHT_size, &b_ScalarHT_size);
   Notify();
}

Bool_t delphes_rapgap_isr_fsr_analysis_h1::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void delphes_rapgap_isr_fsr_analysis_h1::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t delphes_rapgap_isr_fsr_analysis_h1::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef delphes_rapgap_isr_fsr_analysis_h1_cxx
