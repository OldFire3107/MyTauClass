#define TauMatch_cxx
#include "TauMatch.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TMath.h>
#include <iostream>

void TauMatch::Loop()
{
//   In a ROOT session, you can do:
//      root> .L TauMatch.C
//      root> TauMatch t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

// Ploting histograms
   TH1F *h1ratiosig = new TH1F("h1ratiosig", "pt_pi1_matched/pt_tau", 30, 0, 1);
   TH1F *h1ratiobg  = new TH1F("h1ratiobg", "pt_pi1_unmatched/pt_tau", 30, 0, 1);

   TH1F *h1deta = new TH1F("h1deta", "deta", 30, 0, 0.01);
   TH1F *h1dphi  = new TH1F("h1dphi", "dphi", 30, 0, 0.01);
   TH1F *h1dR = new TH1F("h1dR", "dR", 30, 0, 0.015);
   TH1F *h1dpT  = new TH1F("h1dpT", "dpT", 30, 0, 0.1);
   
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(jentry%50000==0) cout << "processed " << jentry << endl;
   
      Float_t pi_pt[3];
      pi_pt[0] = pi1r_pt;
      pi_pt[1] = pi2r_pt;
      pi_pt[2] = pi3r_pt;
      sort(pi_pt, pi_pt+3);

      Float_t dpT_max = 0;
      Float_t dR_max = 0;
      Float_t deta_max = 0;
      Float_t dphi_max = 0;
      
      Float_t dpT, dR, deta, dphi;
      
      dpT_max = abs(pi1r_pt-pi1_pt);
      deta_max = abs(pi1r_eta-pi1_eta);
      dphi_max = abs(pi1r_phi-pi1_phi);
      dR_max = TMath::Sqrt(deta_max*deta_max + dphi_max*dphi_max);

      dpT = abs(pi2r_pt-pi2_pt);
      deta = abs(pi2r_eta-pi2_eta);
      dphi = abs(pi2r_phi-pi2_phi);
      dR = TMath::Sqrt(deta*deta + dphi*dphi);

      if (dpT_max < dpT) dpT_max = dpT;
      if (deta_max < deta) deta_max = deta;
      if (dphi_max < dphi) dphi_max = dphi;
      if (dR_max < dR) dR_max = dR;

      dpT = abs(pi3r_pt-pi3_pt);
      deta = abs(pi3r_eta-pi3_eta);
      dphi = abs(pi3r_phi-pi3_phi);
      dR = TMath::Sqrt(deta*deta + dphi*dphi);

      if (dpT_max < dpT) dpT_max = dpT;
      if (deta_max < deta) deta_max = deta;
      if (dphi_max < dphi) dphi_max = dphi;
      if (dR_max < dR) dR_max = dR;


      if(flag)
      {
         h1ratiosig->Fill(pi_pt[0]/tau_pt);
         h1dpT->Fill(dpT_max);
         h1deta->Fill(deta_max);
         h1dphi->Fill(dphi_max);
         h1dR->Fill(dR_max);
      }
      else
      {
         h1ratiobg->Fill(pi_pt[1]/tau_pt);
      }
   }


/// Ploting section
   TCanvas *can1 = new TCanvas("can6", "Ratios_signbg", 300,20,1000,750);
   can1->Divide(1,2);
   can1->cd(1);
   h1ratiosig->Draw();
   can1->cd(2);
   h1ratiobg->Draw();
   can1->SaveAs("pTRatios.pdf");

   TCanvas *can2 = new TCanvas("can2","Difference",200,10,1400,1000);
   can2->Divide(2,2);
   can2->cd(1);
   h1dpT->Draw();
   can2->cd(2);
   h1deta->Draw();
   can2->cd(3);
   h1dphi->Draw();
   can2->cd(4);
   h1dR->Draw();
   can2->SaveAs("DeltaStuff.pdf");
}
