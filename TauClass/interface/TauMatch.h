//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 19 14:41:29 2022 by ROOT version 6.24/06
// from TTree triplet/triplet
// found on file: final.root
//////////////////////////////////////////////////////////

#ifndef TauMatch_h
#define TauMatch_h

#define PI_MASS 0.1396 // GeV
#define RHO_MASS 0.77 // GeV

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

using namespace std;

class TauMatch {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Long64_t         evt;
   Long64_t         run;
   Float_t         pi1_pt;
   Float_t         pi2_pt;
   Float_t         pi3_pt;
   Float_t         pi1_eta;
   Float_t         pi2_eta;
   Float_t         pi3_eta;
   Float_t         pi1_phi;
   Float_t         pi2_phi;
   Float_t         pi3_phi;
   Float_t         tau_pt;
   Float_t         pi1r_pt;
   Float_t         pi2r_pt;
   Float_t         pi3r_pt;
   Float_t         pi1r_eta;
   Float_t         pi2r_eta;
   Float_t         pi3r_eta;
   Float_t         pi1r_phi;
   Float_t         pi2r_phi;
   Float_t         pi3r_phi;
   Int_t           pi1r_q;
   Int_t           pi2r_q;
   Int_t           pi3r_q;
   Float_t         pi1r_doca3d;
   Float_t         pi2r_doca3d;
   Float_t         pi3r_doca3d;
   Float_t         pi1r_doca3de;
   Float_t         pi2r_doca3de;
   Float_t         pi3r_doca3de;
   Float_t         pi1r_doca2d;
   Float_t         pi2r_doca2d;
   Float_t         pi3r_doca2d;
   Float_t         pi1r_doca2de;
   Float_t         pi2r_doca2de;
   Float_t         pi3r_doca2de;
   Float_t         pi1r_dz;
   Float_t         pi2r_dz;
   Float_t         pi3r_dz;
   Bool_t          flag;
   Bool_t          dup_flag;
   Int_t           dups_count;

   // List of branches
   TBranch        *b_pi1_pt;   //!
   TBranch        *b_pi2_pt;   //!
   TBranch        *b_pi3_pt;   //!
   TBranch        *b_pi1_eta;   //!
   TBranch        *b_pi2_eta;   //!
   TBranch        *b_pi3_eta;   //!
   TBranch        *b_pi1_phi;   //!
   TBranch        *b_pi2_phi;   //!
   TBranch        *b_pi3_phi;   //!
   TBranch        *b_tau_pt;   //!
   TBranch        *b_pi1r_pt;   //!
   TBranch        *b_pi2r_pt;   //!
   TBranch        *b_pi3r_pt;   //!
   TBranch        *b_pi1r_eta;   //!
   TBranch        *b_pi2r_eta;   //!
   TBranch        *b_pi3r_eta;   //!
   TBranch        *b_pi1r_phi;   //!
   TBranch        *b_pi2r_phi;   //!
   TBranch        *b_pi3r_phi;   //!
   TBranch        *b_pi1r_q;   //!
   TBranch        *b_pi2r_q;   //!
   TBranch        *b_pi3r_q;   //!
   TBranch        *b_pi1r_doca3d;   //!
   TBranch        *b_pi2r_doca3d;   //!
   TBranch        *b_pi3r_doca3d;   //!
   TBranch        *b_pi1r_doca3de;   //!
   TBranch        *b_pi2r_doca3de;   //!
   TBranch        *b_pi3r_doca3de;   //!
   TBranch        *b_pi1r_doca2d;   //!
   TBranch        *b_pi2r_doca2d;   //!
   TBranch        *b_pi3r_doca2d;   //!
   TBranch        *b_pi1r_doca2de;   //!
   TBranch        *b_pi2r_doca2de;   //!
   TBranch        *b_pi3r_doca2de;   //!
   TBranch        *b_pi1r_dz;   //!
   TBranch        *b_pi2r_dz;   //!
   TBranch        *b_pi3r_dz;   //!
   TBranch        *b_flag;   //!
   TBranch        *b_dup_flag;   //!
   TBranch        *b_dups_count;   //!

   TString        FileNameIn;
   TFile          *FileIn=NULL;
   TFile          *FileOut=NULL;
   Float_t        pir_pt[3];
   Float_t        pir_eta[3];
   Float_t        pir_phi[3];
   Int_t          pir_q[3];
   Float_t        pir_doca3d[3]; // same as doca, sry :)
   Float_t        pir_doca3de[3];
   Float_t        pir_doca2d[3]; // same as doca, sry :)
   Float_t        pir_doca2de[3];
   Float_t        pir_doca2ds[3];
   Float_t        pir_dz[3];
   Float_t        weighteddR;
   Float_t        mpipi;
   Float_t        mrho;
   Float_t        deta12;
   Float_t        dphi12;
   Float_t        dR12;
   Float_t        deta23;
   Float_t        dphi23;
   Float_t        dR23;
   Float_t        deta31;
   Float_t        dphi31;
   Float_t        dR31;
   Float_t        vis_mass;
   Float_t        w_1;
   TTree          *tree1 = NULL;

   TauMatch(TTree *tree=0);
   TauMatch(const char *file);
   virtual ~TauMatch();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     InitOut();
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual float    deltaPhi(float phi1, float phi2);
   virtual Float_t  getDeltaR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2);
};

#endif

#ifdef TauMatch_cxx
TauMatch::TauMatch(TTree *tree) : fChain(0) 
{
   if(FileOut) delete FileOut;
   FileOut = new TFile("vars.root", "RECREATE");
   tree1=NULL;
   InitOut();

// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      if(FileIn) delete FileIn;
      FileIn = (TFile*)gROOT->GetListOfFiles()->FindObject("final.root");
      if (!FileIn || !FileIn->IsOpen()) {
         FileIn = new TFile("final.root");
      }
      FileIn->GetObject("triplet",tree);

   }
   Init(tree);
}

TauMatch::TauMatch(const char *file)
{
   if(FileOut) delete FileOut;
   FileOut = new TFile("vars.root", "RECREATE");
   tree1=NULL;
   InitOut();

   if(FileIn) delete FileIn;
   FileNameIn = file;
   FileIn = new TFile(file);
   TTree *tree=NULL;
   FileIn->GetObject("triplet",tree);

   Init(tree);
}

TauMatch::~TauMatch()
{
   delete tree1;
   delete FileIn;
   delete FileOut;

   // if (!fChain) return;
   // delete fChain->GetCurrentFile(); // seg fault sometimes since file already deleted
}

Int_t TauMatch::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TauMatch::LoadTree(Long64_t entry)
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

void TauMatch::InitOut()
{
   if (tree1) delete tree1;

   tree1 = new TTree("data_vars","data_vars");
   // TTree *tree1 = new TTree("reco_triplet","reco_triplet");
   tree1->Branch("evt", &evt, "evt/L");
   tree1->Branch("run", &run, "run/L");
   tree1->Branch("pi1r_pt", &pir_pt[0], "pi1r_pt/F");
   tree1->Branch("pi2r_pt", &pir_pt[1], "pi2r_pt/F");
   tree1->Branch("pi3r_pt", &pir_pt[2], "pi3r_pt/F");
   tree1->Branch("pi1r_eta", &pir_eta[0], "pi1r_eta/F");
   tree1->Branch("pi2r_eta", &pir_eta[1], "pi2r_eta/F");
   tree1->Branch("pi3r_eta", &pir_eta[2], "pi3r_eta/F");
   tree1->Branch("pi1r_phi", &pir_phi[0], "pi1r_phi/F");
   tree1->Branch("pi2r_phi", &pir_phi[1], "pi2r_phi/F");
   tree1->Branch("pi3r_phi", &pir_phi[2], "pi3r_phi/F");
   tree1->Branch("pi1r_q", &pir_q[0], "pi1r_q/I");
   tree1->Branch("pi2r_q", &pir_q[1], "pi2r_q/I");
   tree1->Branch("pi3r_q", &pir_q[2], "pi3r_q/I");
   tree1->Branch("pi1r_doca3d", &pir_doca3d[0], "pi1r_doca3d/F"); // Jbypsi vetrex to pion vertex
   tree1->Branch("pi2r_doca3d", &pir_doca3d[1], "pi2r_doca3d/F");
   tree1->Branch("pi3r_doca3d", &pir_doca3d[2], "pi3r_doca3d/F");
   tree1->Branch("pi1r_doca3de", &pir_doca3de[0], "pi1r_doca3de/F"); // Jbypsi vetrex to pion vertex
   tree1->Branch("pi2r_doca3de", &pir_doca3de[1], "pi2r_doca3de/F");
   tree1->Branch("pi3r_doca3de", &pir_doca3de[2], "pi3r_doca3de/F");
   tree1->Branch("pi1r_doca2d", &pir_doca2d[0], "pi1r_doca2d/F"); // Jbypsi vetrex to pion vertex
   tree1->Branch("pi2r_doca2d", &pir_doca2d[1], "pi2r_doca2d/F");
   tree1->Branch("pi3r_doca2d", &pir_doca2d[2], "pi3r_doca2d/F");
   tree1->Branch("pi1r_doca2de", &pir_doca2de[0], "pi1r_doca2de/F"); // Jbypsi vetrex to pion vertex
   tree1->Branch("pi2r_doca2de", &pir_doca2de[1], "pi2r_doca2de/F");
   tree1->Branch("pi3r_doca2de", &pir_doca2de[2], "pi3r_doca2de/F");
   tree1->Branch("pi1r_doca2ds", &pir_doca2ds[0], "pi1r_doca2ds/F"); // Jbypsi vetrex to pion vertex
   tree1->Branch("pi2r_doca2ds", &pir_doca2ds[1], "pi2r_doca2ds/F");
   tree1->Branch("pi3r_doca2ds", &pir_doca2ds[2], "pi3r_doca2ds/F");
   tree1->Branch("pi1r_dz", &pir_dz[0], "pi1r_dz/F");
   tree1->Branch("pi2r_dz", &pir_dz[1], "pi2r_dz/F");
   tree1->Branch("pi3r_dz", &pir_dz[2], "pi3r_dz/F");
   tree1->Branch("weighteddR", &weighteddR, "weighteddR/F");
   tree1->Branch("mpipi", &mpipi, "mpipi/F");
   tree1->Branch("mrho", &mrho, "mrho/F");
   tree1->Branch("vis_mass", &vis_mass, "vis_mass/F");
   tree1->Branch("deta12", &deta12, "deta12/F");
   tree1->Branch("dphi12", &dphi12, "dphi12/F");
   tree1->Branch("dR12", &dR12, "dR12/F");
   tree1->Branch("deta23", &deta23, "deta23/F");
   tree1->Branch("dphi23", &dphi23, "dphi23/F");
   tree1->Branch("dR23", &dR23, "dR23/F");
   tree1->Branch("deta31", &deta31, "deta31/F");
   tree1->Branch("dphi31", &dphi31, "dphi31/F");
   tree1->Branch("dR31", &dR31, "dR31/F");
   tree1->Branch("flag", &flag, "flag/O");
   tree1->Branch("dup_flag", &dup_flag, "dup_flag/O");
   tree1->Branch("dups_count", &dups_count, "dups_count/I");
   tree1->Branch("w_1", &w_1, "w_1/F");
}

void TauMatch::Init(TTree *tree)
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

   fChain->SetBranchAddress("pi1_pt", &pi1_pt, &b_pi1_pt);
   fChain->SetBranchAddress("pi2_pt", &pi2_pt, &b_pi2_pt);
   fChain->SetBranchAddress("pi3_pt", &pi3_pt, &b_pi3_pt);
   fChain->SetBranchAddress("pi1_eta", &pi1_eta, &b_pi1_eta);
   fChain->SetBranchAddress("pi2_eta", &pi2_eta, &b_pi2_eta);
   fChain->SetBranchAddress("pi3_eta", &pi3_eta, &b_pi3_eta);
   fChain->SetBranchAddress("pi1_phi", &pi1_phi, &b_pi1_phi);
   fChain->SetBranchAddress("pi2_phi", &pi2_phi, &b_pi2_phi);
   fChain->SetBranchAddress("pi3_phi", &pi3_phi, &b_pi3_phi);
   fChain->SetBranchAddress("tau_pt", &tau_pt, &b_tau_pt);
   fChain->SetBranchAddress("pi1r_pt", &pi1r_pt, &b_pi1r_pt);
   fChain->SetBranchAddress("pi2r_pt", &pi2r_pt, &b_pi2r_pt);
   fChain->SetBranchAddress("pi3r_pt", &pi3r_pt, &b_pi3r_pt);
   fChain->SetBranchAddress("pi1r_eta", &pi1r_eta, &b_pi1r_eta);
   fChain->SetBranchAddress("pi2r_eta", &pi2r_eta, &b_pi2r_eta);
   fChain->SetBranchAddress("pi3r_eta", &pi3r_eta, &b_pi3r_eta);
   fChain->SetBranchAddress("pi1r_phi", &pi1r_phi, &b_pi1r_phi);
   fChain->SetBranchAddress("pi2r_phi", &pi2r_phi, &b_pi2r_phi);
   fChain->SetBranchAddress("pi3r_phi", &pi3r_phi, &b_pi3r_phi);
   fChain->SetBranchAddress("pi1r_q", &pi1r_q, &b_pi1r_q);
   fChain->SetBranchAddress("pi2r_q", &pi2r_q, &b_pi2r_q);
   fChain->SetBranchAddress("pi3r_q", &pi3r_q, &b_pi3r_q);
   fChain->SetBranchAddress("pi1r_doca3d", &pi1r_doca3d, &b_pi1r_doca3d);
   fChain->SetBranchAddress("pi2r_doca3d", &pi2r_doca3d, &b_pi2r_doca3d);
   fChain->SetBranchAddress("pi3r_doca3d", &pi3r_doca3d, &b_pi3r_doca3d);
   fChain->SetBranchAddress("pi1r_doca3de", &pi1r_doca3de, &b_pi1r_doca3de);
   fChain->SetBranchAddress("pi2r_doca3de", &pi2r_doca3de, &b_pi2r_doca3de);
   fChain->SetBranchAddress("pi3r_doca3de", &pi3r_doca3de, &b_pi3r_doca3de);
   fChain->SetBranchAddress("pi1r_doca2d", &pi1r_doca2d, &b_pi1r_doca2d);
   fChain->SetBranchAddress("pi2r_doca2d", &pi2r_doca2d, &b_pi2r_doca2d);
   fChain->SetBranchAddress("pi3r_doca2d", &pi3r_doca2d, &b_pi3r_doca2d);
   fChain->SetBranchAddress("pi1r_doca2de", &pi1r_doca2de, &b_pi1r_doca2de);
   fChain->SetBranchAddress("pi2r_doca2de", &pi2r_doca2de, &b_pi2r_doca2de);
   fChain->SetBranchAddress("pi3r_doca2de", &pi3r_doca2de, &b_pi3r_doca2de);
   fChain->SetBranchAddress("pi1r_dz", &pi1r_dz, &b_pi1r_dz);
   fChain->SetBranchAddress("pi2r_dz", &pi2r_dz, &b_pi2r_dz);
   fChain->SetBranchAddress("pi3r_dz", &pi3r_dz, &b_pi3r_dz);
   fChain->SetBranchAddress("flag", &flag, &b_flag);
   fChain->SetBranchAddress("dup_flag", &dup_flag, &b_dup_flag);
   fChain->SetBranchAddress("dups_count", &dups_count, &b_dups_count);
   Notify();
}

Bool_t TauMatch::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TauMatch::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TauMatch::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TauMatch_cxx
