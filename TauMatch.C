#define TauMatch_cxx
#include "TauMatch.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TMath.h>
#include <TLorentzVector.h>
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
   TH1F *h1ratiosig1 = new TH1F("h1ratiosig1", "pt_pi1_matched/pt_tau", 30, 0, 1);
   TH1F *h1ratiobg1  = new TH1F("h1ratiobg1", "pt_pi1_unmatched/pt_tau", 30, 0, 1);

   TH1F *h1ratiosig2 = new TH1F("h1ratiosig2", "pt_pi2_matched/pt_tau", 30, 0, 1);
   TH1F *h1ratiobg2  = new TH1F("h1ratiobg2", "pt_pi2_unmatched/pt_tau", 30, 0, 1);

   TH1F *h1ratiosig3 = new TH1F("h1ratiosig3", "pt_pi3_matched/pt_tau", 30, 0, 1);
   TH1F *h1ratiobg3  = new TH1F("h1ratiobg3", "pt_pi3_unmatched/pt_tau", 30, 0, 1);

   TH1F *h1ratiosig4 = new TH1F("h1ratiosig4", "pt_sum/pt_tau", 30, 0, 2);
   TH1F *h1ratiobg4  = new TH1F("h1ratiobg4", "pt_sum/pt_tau", 30, 0, 2);

   TH1F *h1deta = new TH1F("h1deta", "deta abs", 30, 0, 0.01);
   TH1F *h1dphi  = new TH1F("h1dphi", "dphi abs", 30, 0, 0.01);
   TH1F *h1dR = new TH1F("h1dR", "dR", 30, 0, 0.015);
   TH1F *h1dpT  = new TH1F("h1dpT", "dpT", 30, 0, 0.1);

   TH1F *h1vismasssig = new TH1F("h1vismasssig", "vis_mass", 20, 0,4);
   TH1F *h1vismassbg = new TH1F("h1vismassbg", "vis_mass", 20, 0,4);
   
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
      Float_t pi_ptg[3];
      pi_ptg[0] = pi1_pt;
      pi_ptg[1] = pi2_pt;
      pi_ptg[2] = pi3_pt;
      sort(pi_pt, pi_pt+3);
      sort(pi_ptg, pi_ptg+3);


      TLorentzVector P[3];
      P[0].SetPtEtaPhiM(pi1r_pt, pi1r_eta, pi1r_phi, PI_MASS);
      P[1].SetPtEtaPhiM(pi2r_pt, pi2r_eta, pi2r_phi, PI_MASS);
      P[2].SetPtEtaPhiM(pi3r_pt, pi3r_eta, pi3r_phi, PI_MASS);
      TLorentzVector Invar = P[0] + P[1] + P[2];
      Double_t vis_mass = Invar.Mag();
      Double_t Pt_tot = Invar.Pt();

      Float_t dR_max = 0;
      Float_t deta_max = 0;
      Float_t dphi_max = 0;
      
      Float_t dR, deta, dphi;
      
      Float_t dphitemp = pi1r_phi-pi1_phi;
      if(dphitemp > 2*TMath::Pi()) dphitemp -= 2*TMath::Pi();
      if(dphitemp < -2*TMath::Pi()) dphitemp += 2*TMath::Pi();
      deta_max = abs(pi1r_eta-pi1_eta);
      dphi_max = abs(dphitemp);
      dR_max = TMath::Sqrt(deta_max*deta_max + dphi_max*dphi_max);

      dphitemp = pi2r_phi-pi2_phi;
      if(dphitemp > 2*TMath::Pi()) dphitemp -= 2*TMath::Pi();
      if(dphitemp < -2*TMath::Pi()) dphitemp += 2*TMath::Pi();
      deta = abs(pi2r_eta-pi2_eta);
      dphi = abs(dphitemp);
      dR = TMath::Sqrt(deta*deta + dphi*dphi);

      if (deta_max < deta) deta_max = deta;
      if (dphi_max < dphi) dphi_max = dphi;
      if (dR_max < dR) dR_max = dR;

      dphitemp = pi3r_phi-pi3_phi;
      if(dphitemp > 2*TMath::Pi()) dphitemp -= 2*TMath::Pi();
      if(dphitemp < -2*TMath::Pi()) dphitemp += 2*TMath::Pi();
      deta = abs(pi3r_eta-pi3_eta);
      dphi = abs(dphitemp);
      dR = TMath::Sqrt(deta*deta + dphi*dphi);

      if (deta_max < deta) deta_max = deta;
      if (dphi_max < dphi) dphi_max = dphi;
      if (dR_max < dR) dR_max = dR;


      if(flag)
      {
         h1ratiosig1->Fill(pi_ptg[2]/Pt_tot);
         h1ratiosig2->Fill(pi_ptg[1]/Pt_tot);
         h1ratiosig3->Fill(pi_ptg[0]/Pt_tot);
         h1ratiosig4->Fill(Pt_tot/tau_pt);
         h1vismasssig->Fill(vis_mass);
         h1dpT->Fill((pi1_pt-pi1r_pt)/pi1_pt);
         h1dpT->Fill((pi2_pt-pi2r_pt)/pi2_pt);
         h1dpT->Fill((pi3_pt-pi3r_pt)/pi3_pt);
         h1deta->Fill(deta_max);
         h1dphi->Fill(dphi_max);
         h1dR->Fill(dR_max);
      }
      else
      {
         h1ratiobg1->Fill(pi_pt[2]/Pt_tot);
         h1ratiobg2->Fill(pi_pt[1]/Pt_tot);
         h1ratiobg3->Fill(pi_pt[0]/Pt_tot);
         h1ratiobg4->Fill(Pt_tot/tau_pt);
         h1vismassbg->Fill(vis_mass);
      }
   }


/// Ploting section
   TCanvas *can1 = new TCanvas("can1", "Ratios1_signbg_vis", 300,20,1000,750);
   can1->Divide(1,2);
   can1->cd(1);
   h1ratiosig1->Draw();
   can1->cd(2);
   h1ratiobg1->Draw();
   can1->SaveAs("pT1Ratios_vis.pdf");

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

   TCanvas *can3 = new TCanvas("can3", "Ratios2_signbg_vis", 300,20,1000,750);
   can3->Divide(1,2);
   can3->cd(1);
   h1ratiosig2->Draw();
   can3->cd(2);
   h1ratiobg2->Draw();
   can3->SaveAs("pT2Ratios_vis.pdf");

   TCanvas *can4 = new TCanvas("can4", "Ratios3_signbg_vis", 300,20,1000,750);
   can4->Divide(1,2);
   can4->cd(1);
   h1ratiosig3->Draw();
   can4->cd(2);
   h1ratiobg3->Draw();
   can4->SaveAs("pT3Ratios_vis.pdf");

   TCanvas *can5 = new TCanvas("can5", "Ratios4_signbg_vis", 300,20,1000,750);
   can5->Divide(1,2);
   can5->cd(1);
   h1ratiosig4->Draw();
   can5->cd(2);
   h1ratiobg4->Draw();
   can5->SaveAs("pT4Ratios_vis.pdf");

   TCanvas *can6 = new TCanvas("can6", "vis_mass", 300,20,1000,750);
   can6->Divide(1,2);
   can6->cd(1);
   h1vismasssig->Draw();
   can6->cd(2);
   h1vismassbg->Draw();
   can6->SaveAs("vis_mass.pdf");
}
