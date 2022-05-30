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
   Float_t         pi1r_dz;
   Float_t         pi2r_dz;
   Float_t         pi3r_dz;
   Bool_t          flag;

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
   TBranch        *b_pi1r_dz;   //!
   TBranch        *b_pi2r_dz;   //!
   TBranch        *b_pi3r_dz;   //!
   TBranch        *b_flag;   //!

   TauMatch(TTree *tree=0);
   TauMatch(const char *file);
   virtual ~TauMatch();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
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
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("final.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("final.root");
      }
      f->GetObject("triplet",tree);

   }
   Init(tree);
}

TauMatch::TauMatch(const char *file)
{
   TFile *FileIn = new TFile(file);
   TTree *tree=NULL;
   FileIn->GetObject("triplet",tree);

   Init(tree);
}

TauMatch::~TauMatch()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
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
   fChain->SetBranchAddress("pi1r_dz", &pi1r_dz, &b_pi1r_dz);
   fChain->SetBranchAddress("pi2r_dz", &pi2r_dz, &b_pi2r_dz);
   fChain->SetBranchAddress("pi3r_dz", &pi3r_dz, &b_pi3r_dz);
   fChain->SetBranchAddress("flag", &flag, &b_flag);
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
