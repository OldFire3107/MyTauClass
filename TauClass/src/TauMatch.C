#define TauMatch_cxx
#include "MyTauClass/TauClass/interface/TauMatch.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>

using namespace std;

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

   TH1F *h1ratiosig4 = new TH1F("h1ratiosig4", "pt_sum_matched/pt_tau", 30, 0, 2);
   TH1F *h1ratiobg4  = new TH1F("h1ratiobg4", "pt_sum_unmatched/pt_tau", 30, 0, 2);

   TH1F *h1ratiosig5 = new TH1F("h1ratiosig5", "pt_pi1/pt_pt2 matched", 30, 0, 3);
   TH1F *h1ratiobg5 = new TH1F("h1ratiobg5", "pt_pi1/pt_pt2 unmatched", 30, 0, 3);

   TH1F *h1ratiosig6 = new TH1F("h1ratiosig6", "pt_pi2/pt_pt3 matched", 30, 0, 3);
   TH1F *h1ratiobg6 = new TH1F("h1ratiobg6", "pt_pi2/pt_pt3 unmatched", 30, 0, 3);

   TH1F *h1deta = new TH1F("h1deta", "deta abs", 30, 0, 0.01);
   TH1F *h1dphi  = new TH1F("h1dphi", "dphi abs", 30, 0, 0.01);
   TH1F *h1dR = new TH1F("h1dR", "dR", 30, 0, 0.015);
   TH1F *h1dpT  = new TH1F("h1dpT", "dpT reso", 30, -0.1, 0.1);

   TH1F *h1vismasssig = new TH1F("h1vismasssig", "vis_mass_matched", 20, 0,4);
   TH1F *h1vismassbg = new TH1F("h1vismassbg", "vis_mass_background", 20, 0,4);

   TH1F *h1weightpTdRsig = new TH1F("h1weightpTdRsig", "dRij*(Pti+Ptj)/2Ptot sig", 30, 0, 3);
   TH1F *h1weightpTdRbg = new TH1F("h1weightpTdRbg", "dRij*(Pti+Ptj)/2Ptot bg", 30, 0, 3);

   TH1F *h1mrhosig = new TH1F("h1mrhosig", "mpipi closest to rho matched", 30, 0, 1.4);
   TH1F *h1mrhobg = new TH1F("h1mrhobg", "mpipi closest to rho unmatched", 30, 0, 1.4);

   TH1F *h1mpipisig = new TH1F("h1mpipisig", "mpipi matched", 30, 0, 1.6);
   TH1F *h1mpipisbg = new TH1F("h1mpipisbg", "mpipi matched", 30, 0, 1.6);
   TH1F *h1mpipidbg = new TH1F("h1mpipisdg", "mpipi matched", 30, 0, 1.6);

   TH1F *h1drsig = new TH1F("h1drsig", "dr 12 matched", 30, 0, 2);
   TH1F *h1drbg = new TH1F("h1drbg", "dr 12 unmatched", 30, 0, 2);

   TH1F *h1detasig = new TH1F("h1detesig", "deta 12 matched", 30, 0, 1.5);
   TH1F *h1detabg = new TH1F("h1detebg", "deta 12 unmatched", 30, 0, 1.5);

   TH1F *h1dphisig = new TH1F("h1dphisig", "dphi 12 matched", 30, 0, 1.5);
   TH1F *h1dphibg = new TH1F("h1dphibg", "dphi 12 unmatched", 30, 0, 1.5);

   TH1F *h1dzpi1sig = new TH1F("h1dzpi1sig", "dz pi1 matched", 30, -0.05, 0.05);
   TH1F *h1dzpi1bg = new TH1F("h1dzpi1bg", "dz pi1 unmatched", 30, -0.05, 0.05);
   TH1F *h1doca3dpi1sig = new TH1F("h1doca3dpi1sig", "doca3d pi1 matched", 30, -0.05, 0.05);
   TH1F *h1doca3dpi1bg = new TH1F("h1doca3dpi1bg", "doca3d pi1 unmatched", 30, -0.05, 0.05);

   TH2F *h2drvsptsig = new TH2F("h2drvsptsignal", "dr vs pt signal", 50, 1, 10, 25, 0, 2);
   TH2F *h2drvsptbg = new TH2F("h2drvsptbg", "dr vs pt bg", 50, 1, 10, 25, 0, 2);
   // TH2F *h2dr12vsptsig = new TH2F("h2dr12vsptsignal", "dr12 vs pt signal", 50, 1, 10, 25, 0, 2);
   // TH2F *h2dr12vsptbg = new TH2F("h2dr12vsptbg", "dr12 vs pt bg", 50, 1, 10, 25, 0, 2);
   // TH2F *h2dr23vsptsig = new TH2F("h2dr23vsptsignal", "dr23 vs pt signal", 50, 1, 10, 25, 0, 2);
   // TH2F *h2dr23vsptbg = new TH2F("h2dr23vsptbg", "dr23 vs pt bg", 50, 1, 10, 25, 0, 2);
   // TH2F *h2dr31vsptsig = new TH2F("h2dr31vsptsignal", "dr31 vs pt signal", 50, 1, 10, 25, 0, 2);
   // TH2F *h2dr31vsptbg = new TH2F("h2dr31vsptbg", "dr31 vs pt bg", 50, 1, 10, 25, 0, 2);
   TH2F *h2detavsptsig = new TH2F("h2detavsptsig", "deta vs pt sig", 50, 1, 10, 25, 0, 1.5);
   TH2F *h2detavsptbg = new TH2F("h2detavsptbg", "deta vs pt bg", 50, 1, 10, 25, 0, 1.5);

   TH2F *h2dphivsptsig = new TH2F("h2dphivsptsig", "dphi vs pt sig", 50, 1, 10, 25, 0, 1.5);
   TH2F *h2dphivsptbg = new TH2F("h2dphivsptbg", "dphi vs pt bg", 50, 1, 10, 25, 0, 1.5);

   TH2F *h2resoptratio = new TH2F("h2resoptratio", "dpt resolution vs pt", 30, 0, 10, 30, -0.1, 0.1);
   
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
      // sort(pi_pt, pi_pt+3);
      sort(pi_ptg, pi_ptg+3);


      TLorentzVector P[3];
      P[0].SetPtEtaPhiM(pi1r_pt, pi1r_eta, pi1r_phi, PI_MASS);
      P[1].SetPtEtaPhiM(pi2r_pt, pi2r_eta, pi2r_phi, PI_MASS);
      P[2].SetPtEtaPhiM(pi3r_pt, pi3r_eta, pi3r_phi, PI_MASS);
      TLorentzVector Invar = P[0] + P[1] + P[2];
      Double_t vis_mass = Invar.M();
      Double_t Pt_tot = Invar.Pt();

      Float_t dR_max = 0;
      Float_t deta_max = 0;
      Float_t dphi_max = 0;
      
      Float_t dR, deta, dphi;
      
      deta_max = abs(pi1r_eta-pi1_eta);
      dphi_max = deltaPhi(pi1r_phi, pi1_phi);
      dR_max = TMath::Sqrt(deta_max*deta_max + dphi_max*dphi_max);

      deta = abs(pi2r_eta-pi2_eta);
      dphi = deltaPhi(pi2r_phi, pi2_phi);
      dR = TMath::Sqrt(deta*deta + dphi*dphi);

      if (deta_max < deta) deta_max = deta;
      if (dphi_max < dphi) dphi_max = dphi;
      if (dR_max < dR) dR_max = dR;

      deta = abs(pi3r_eta-pi3_eta);
      dphi = deltaPhi(pi3r_phi, pi3_phi);
      dR = TMath::Sqrt(deta*deta + dphi*dphi);

      if (deta_max < deta) deta_max = deta;
      if (dphi_max < dphi) dphi_max = dphi;
      if (dR_max < dR) dR_max = dR;

      Float_t b_dr_max = 0;
      Float_t b_eta_max = 0;
      Float_t b_phi_max = 0;

      Float_t dr12, dr23, dr31;
      Float_t deta12, deta23, deta31;
      Float_t dphi12, dphi23, dphi31;

      dphi12 = deltaPhi(pi1r_phi, pi2r_phi);
      deta12 = abs(pi1r_eta - pi2r_eta);
      dr12 = getDeltaR(pi1r_eta, pi1r_phi, pi2r_eta, pi2r_phi);
      b_dr_max = dr12 > b_dr_max ? dr12 : b_dr_max;
      b_eta_max = deta12 > b_eta_max ? deta12 : b_eta_max;
      b_phi_max = dphi12 > b_phi_max ? dphi12 : b_phi_max;
      dphi23 = deltaPhi(pi3r_phi, pi2r_phi);
      deta23 = abs(pi3r_eta - pi2r_eta);
      dr23 = getDeltaR(pi3r_eta, pi3r_phi, pi2r_eta, pi2r_phi);
      b_dr_max = dr23 > b_dr_max ? dr23 : b_dr_max;
      b_eta_max = deta23 > b_eta_max ? deta23 : b_eta_max;
      b_phi_max = dphi23 > b_phi_max ? dphi23 : b_phi_max;
      dphi31 = deltaPhi(pi1r_phi, pi3r_phi);
      deta31 = abs(pi1r_eta - pi3r_eta);
      dr31 = getDeltaR(pi1r_eta, pi1r_phi, pi3r_eta, pi3r_phi);
      b_dr_max = dr31 > b_dr_max ? dr31 : b_dr_max;
      b_eta_max = deta31 > b_eta_max ? deta31 : b_eta_max;
      b_phi_max = dphi31 > b_phi_max ? dphi31 : b_phi_max;


      TLorentzVector P12 = P[0] + P[1];
      TLorentzVector P23 = P[1] + P[2];
      TLorentzVector P31 = P[2] + P[0];

      Float_t weighteddR = dr12*P12.Pt() + dr23*P23.Pt() + dr31*P31.Pt();
      weighteddR /= 2*Pt_tot;

      Float_t mpipi12 = P12.M();
      Float_t mpipi23 = P23.M();
      Float_t mpipi31 = P31.M();

      Float_t mrho, mindiff = 10;
      if(abs(mindiff) > abs(mpipi12 - RHO_MASS) && pi1r_q * pi2r_q < 0) 
         mindiff = mpipi12 - RHO_MASS;
      if(abs(mindiff) > abs(mpipi23 - RHO_MASS) && pi2r_q * pi3r_q < 0) 
         mindiff = mpipi23 - RHO_MASS;
      if(abs(mindiff) > abs(mpipi31 - RHO_MASS) && pi3r_q * pi1r_q < 0) 
         mindiff = mpipi31 - RHO_MASS;

      mrho = mindiff + RHO_MASS;
   
      if(flag)
      {
         h1ratiosig1->Fill(pi_pt[0]/Pt_tot);
         h1ratiosig2->Fill(pi_pt[1]/Pt_tot);
         h1ratiosig3->Fill(pi_pt[2]/Pt_tot);
         h1ratiosig4->Fill(Pt_tot/tau_pt);
         h1vismasssig->Fill(vis_mass);
         h1dpT->Fill((pi1_pt-pi1r_pt)/pi1_pt);
         h1dpT->Fill((pi2_pt-pi2r_pt)/pi2_pt);
         h1dpT->Fill((pi3_pt-pi3r_pt)/pi3_pt);
         h1deta->Fill(deta_max);
         h1dphi->Fill(dphi_max);
         h1dR->Fill(dR_max);

         h1ratiosig5->Fill(pi_pt[0]/pi_pt[1]);
         h1ratiosig6->Fill(pi_pt[1]/pi_pt[2]);

         h2resoptratio->Fill(pi1_pt, (pi1_pt-pi1r_pt)/pi1_pt);
         h2resoptratio->Fill(pi2_pt, (pi2_pt-pi2r_pt)/pi2_pt);
         h2resoptratio->Fill(pi3_pt, (pi3_pt-pi3r_pt)/pi3_pt);

         h2drvsptsig->Fill(Pt_tot, b_dr_max);
         h2detavsptsig->Fill(Pt_tot, b_eta_max);
         h2dphivsptsig->Fill(Pt_tot, b_phi_max);
         // h2dr12vsptsig->Fill(Pt_tot, dr12);
         // h2dr23vsptsig->Fill(Pt_tot, dr23);
         // h2dr31vsptsig->Fill(Pt_tot, dr31);
         h1weightpTdRsig->Fill(weighteddR);
         h1mrhosig->Fill(mrho);
         h1mpipisig->Fill(mpipi12);
         h1drsig->Fill(dr12);
         h1detasig->Fill(deta12);
         h1dphisig->Fill(dphi12);
         h1dzpi1sig->Fill(pi1r_dz);
         h1doca3dpi1sig->Fill(pi1r_doca3d);
      }
      else
      {
         h1ratiobg1->Fill(pi_pt[0]/Pt_tot);
         h1ratiobg2->Fill(pi_pt[1]/Pt_tot);
         h1ratiobg3->Fill(pi_pt[2]/Pt_tot);
         h1ratiobg4->Fill(Pt_tot/tau_pt);
         h1vismassbg->Fill(vis_mass);


         h1ratiobg5->Fill(pi_pt[0]/pi_pt[1]);
         h1ratiobg6->Fill(pi_pt[1]/pi_pt[2]);

         h2drvsptbg->Fill(Pt_tot, b_dr_max);
         h2detavsptbg->Fill(Pt_tot, b_eta_max);
         h2dphivsptbg->Fill(Pt_tot, b_phi_max);
         // h2dr12vsptbg->Fill(Pt_tot, dr12);
         // h2dr23vsptbg->Fill(Pt_tot, dr23);
         // h2dr31vsptbg->Fill(Pt_tot, dr31);
         h1weightpTdRbg->Fill(weighteddR);
         h1mrhobg->Fill(mrho);
         if(pi1r_q == pi2r_q && pi2r_q == pi3r_q)
            h1mpipisbg->Fill(mpipi12); // Useless as the above condition never statisfies
         else
            h1mpipidbg->Fill(mpipi12);
         h1drbg->Fill(b_dr_max);
         h1detabg->Fill(b_eta_max);
         h1dphibg->Fill(b_phi_max);
         h1dzpi1bg->Fill(pi1r_dz);
         h1doca3dpi1bg->Fill(pi1r_doca3d);
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

   TCanvas *can7 = new TCanvas("can7", "dr_vs_pT", 300,20,1000,750);
   can7->Divide(1,2);
   can7->cd(1);
   h2drvsptsig->Draw("colz");
   can7->cd(2);
   h2drvsptbg->Draw("colz");
   can7->SaveAs("dr_vs_pT.pdf");

   // TCanvas *can8 = new TCanvas("can8", "dr12_vs_pT", 300,20,1000,750);
   // can8->Divide(1,2);
   // can8->cd(1);
   // h2dr12vsptsig->Draw("colz");
   // can8->cd(2);
   // h2dr12vsptbg->Draw("colz");
   // // can8->SaveAs("dr12_vs_pT.pdf");

   // TCanvas *can9 = new TCanvas("can9", "dr23_vs_pT", 300,20,1000,750);
   // can9->Divide(1,2);
   // can9->cd(1);
   // h2dr23vsptsig->Draw("colz");
   // can9->cd(2);
   // h2dr23vsptbg->Draw("colz");
   // // can9->SaveAs("dr23_vs_pT.pdf");

   // TCanvas *can10 = new TCanvas("can10", "dr31_vs_pT", 300,20,1000,750);
   // can10->Divide(1,2);
   // can10->cd(1);
   // h2dr31vsptsig->Draw("colz");
   // can10->cd(2);
   // h2dr31vsptbg->Draw("colz");
   // // can10->SaveAs("dr31_vs_pT.pdf");

   TCanvas *can11 = new TCanvas("can11", "dptresvspt", 300,20,1000,750);
   h2resoptratio->Draw("colz");
   can11->SaveAs("dptresovspt.pdf");

   TCanvas *can12 = new TCanvas("can12", "pt_pi1/pt_pi2", 300,20,1000,750);
   can12->Divide(1,2);
   can12->cd(1);
   h1ratiosig5->Draw("colz");
   can12->cd(2);
   h1ratiobg5->Draw("colz");
   can12->SaveAs("1by2.pdf");  

   TCanvas *can13 = new TCanvas("can13", "pt_pi2/pt_pi3", 300,20,1000,750);
   can13->Divide(1,2);
   can13->cd(1);
   h1ratiosig6->Draw("colz");
   can13->cd(2);
   h1ratiobg6->Draw("colz");
   can13->SaveAs("2by3.pdf"); 

   TCanvas *can14 = new TCanvas("can14", "deta_vs_pT", 300,20,1000,750);
   can14->Divide(1,2);
   can14->cd(1);
   h2detavsptsig->Draw("colz");
   can14->cd(2);
   h2detavsptbg->Draw("colz");
   can14->SaveAs("deta_vs_pT.pdf"); 

   TCanvas *can15 = new TCanvas("can15", "dphi_vs_pT", 300,20,1000,750);
   can15->Divide(1,2);
   can15->cd(1);
   h2dphivsptsig->Draw("colz");
   can15->cd(2);
   h2dphivsptbg->Draw("colz");
   can15->SaveAs("dphi_vs_pT.pdf"); 

   TCanvas *can16 = new TCanvas("can16", "weighteddR", 300,20,1000,750);
   can16->Divide(1,2);
   can16->cd(1);
   h1weightpTdRsig->Draw();
   can16->cd(2);
   h1weightpTdRbg->Draw();
   can16->SaveAs("weighteddR.pdf"); 

   TCanvas *can17 = new TCanvas("can17", "mrho", 300,20,1000,750);
   can17->Divide(1,2);
   can17->cd(1);
   h1mrhosig->Draw();
   can17->cd(2);
   h1mrhobg->Draw();
   can17->SaveAs("mrho.pdf");

   TCanvas *can18 = new TCanvas("can18", "mpipi", 300,20,1000,750);
   can18->Divide(1,2);
   can18->cd(1);
   h1mpipisig->Draw();
   can18->cd(2);
   h1mpipidbg->Draw();
   can18->SaveAs("mpipi.pdf");

   TCanvas *can19 = new TCanvas("can19", "dr", 300,20,1000,750);
   can19->Divide(1,2);
   can19->cd(1);
   h1drsig->Draw();
   can19->cd(2);
   h1drbg->Draw();
   can19->SaveAs("dr.pdf");
   
   TCanvas *can20 = new TCanvas("can20", "deta", 300,20,1000,750);
   can20->Divide(1,2);
   can20->cd(1);
   h1detasig->Draw();
   can20->cd(2);
   h1detabg->Draw();
   can20->SaveAs("deta.pdf");

   TCanvas *can21 = new TCanvas("can21", "dphi", 300,20,1000,750);
   can21->Divide(1,2);
   can21->cd(1);
   h1dphisig->Draw();
   can21->cd(2);
   h1dphibg->Draw();
   can21->SaveAs("dphi.pdf");

   TCanvas *can22 = new TCanvas("can22", "dz", 300,20,1000,750);
   can22->Divide(1,2);
   can22->cd(1);
   h1dzpi1sig->Draw();
   can22->cd(2);
   h1dzpi1bg->Draw();
   can22->SaveAs("dz.pdf");

   TCanvas *can23 = new TCanvas("can23", "doca3d", 300,20,1000,750);
   can23->Divide(1,2);
   can23->cd(1);
   h1doca3dpi1sig->Draw();
   can23->cd(2);
   h1doca3dpi1bg->Draw();
   can23->SaveAs("doca3d.pdf");
}

float TauMatch::deltaPhi(float phi1, float phi2) {
  float result = phi1 - phi2;
  while (result > TMath::Pi()) result -= float(2.*TMath::Pi());
  while (result <= -TMath::Pi()) result += float(2.*TMath::Pi());
  float absresult=TMath::Abs(result);
  return absresult;
}

Float_t TauMatch::getDeltaR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2)
{
  Float_t dphi = deltaPhi(phi1, phi2);
  return sqrt(dphi*dphi + (eta1-eta2)*(eta1-eta2));
}
