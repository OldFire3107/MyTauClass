#define MyTauClass_cxx
#include "MyTauClass/interface/MyTauClass.h"
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

void MyTauClass::Loop()
{
  // FileIn->cd();
  if (fChain == 0) return;
  Int_t ievent=0;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;

  //f->Close();

  // TH1F *h1nprong = new TH1F("h1nprong","nprong of the tau",4,0,4);
// TH1F *h1nprongpi = new TH1F("h1nprongtau","nprong of the gen pi",4,0,4);
  //  TH1I *h1npfCand = new TH1I("h1npfCand","number of pfCand",100,0.,100.);
  //  TH1F *h1pfCandpt = new TH1F("h1pfCandpt","pt of the PF candidate",100,0.,40.);

  // start of the loop
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    ievent++;
    if(ievent%50000==0) cout << "processed " << ievent << endl;

    // h1pttau->Fill(JpsiTau_tau_pt->at(0));
    if(JpsiTau_gen_tau_pt->size()!=0){
      unsigned int sum = 0;
      if(JpsiTau_gen_tau_pt->at(0)>1.){ 
        tau_true_pt=JpsiTau_gen_tau_pt->at(0);
        h1pttau->Fill(tau_true_pt);
        h1etatau->Fill(JpsiTau_gen_tau_eta->at(0));
        h1phitau->Fill(JpsiTau_gen_tau_phi->at(0));
        h1nprongtau->Fill(JpsiTau_gen_tau_nprong->at(0));
        h1nmatchedtau->Fill(JpsiTau_gen_tau_nmatched->at(0));

        h2etapt->Fill(JpsiTau_gen_tau_pt->at(0), JpsiTau_gen_tau_eta->at(0));
        h2phipt->Fill(JpsiTau_gen_tau_pt->at(0), JpsiTau_gen_tau_phi->at(0));
        h2etaphi->Fill(JpsiTau_gen_tau_phi->at(0), JpsiTau_gen_tau_eta->at(0)); 

        float max_dr = 0;
        float eta[3];
        float phi[3];

        if (int(JpsiTau_gen_tau_nprong->at(0) == 2))
        {
          max_dr = getDeltaR(JpsiTau_gen_pion_eta->at(sum), JpsiTau_gen_pion_phi->at(sum), JpsiTau_gen_pion_eta->at(sum+1), JpsiTau_gen_pion_phi->at(sum+1));
          h2dRpt->Fill(JpsiTau_gen_tau_pt->at(0), max_dr);
        }

        for(int j=0; j < int(JpsiTau_gen_tau_nprong->at(0)); j++, sum++){
          
          if (int(JpsiTau_gen_tau_nprong->at(0) == 3))
          {
            pi_pt[j] = JpsiTau_gen_pion_pt->at(sum);
            pi_eta[j] = JpsiTau_gen_pion_eta->at(sum);
            pi_phi[j] = JpsiTau_gen_pion_phi->at(sum);
            eta[j] = pi_eta[j];
            phi[j] = pi_phi[j];
          }
          h1ptpi->Fill(JpsiTau_gen_pion_pt->at(sum));
          h1etapi->Fill(JpsiTau_gen_pion_eta->at(sum));
          h1phipi->Fill(JpsiTau_gen_pion_phi->at(sum));
        }

        if (int(JpsiTau_gen_tau_nprong->at(0) == 3))
        {
          float dr = getDeltaR(eta[0], phi[0], eta[1], phi[1]);
          max_dr = max_dr < dr ? dr : max_dr;
          dr = getDeltaR(eta[1], phi[1], eta[2], phi[2]);
          max_dr = max_dr < dr ? dr : max_dr;
          dr = getDeltaR(eta[2], phi[2], eta[0], phi[0]);
          max_dr = max_dr < dr ? dr : max_dr;
          max_dr = TMath::Sqrt(max_dr);
          h2dRpt->Fill(JpsiTau_gen_tau_pt->at(0), max_dr);
        }
        
      }
      else
      {
        sum+=int(JpsiTau_gen_tau_nprong->at(0));
      }

    }/*
    // count pfCand 
    if(pfCand_pt->size()>0){
      for(Int_t j=0; j<pfCand_pt->size(); j++){
  h1pfCandpt->Fill((*pfCand_pt)[j]);
        // add DeltaR function
      }
    };
    */
    //h1npfCand->Fill(pfCand_pt->size());
    
    // if (Cut(ientry) < 0) continue;
    
    // Comparing eta and phi since pt can be lost via Bremsstrahlung
    Float_t dR_min = 40;
    int min_pos = -1;

    vector<vector<Int_t>> num_comb;
    Int_t num_comb_min[3] = {-1, -1, -1};
    Int_t input[3] = {0, 1, 2};
    CombAdd(num_comb, input, 3);

    
    for(int i=0; i < JpsiTau_tau_pi1_pt->size(); i++){
      h1pt1->Fill(JpsiTau_tau_pi1_pt->at(i));
      h1pt2->Fill(JpsiTau_tau_pi2_pt->at(i));
      h1pt3->Fill(JpsiTau_tau_pi3_pt->at(i));
      pir_pt[0] = JpsiTau_tau_pi1_pt->at(i);
      pir_pt[1] = JpsiTau_tau_pi2_pt->at(i);
      pir_pt[2] = JpsiTau_tau_pi3_pt->at(i);
      pir_eta[0] = JpsiTau_tau_pi1_eta->at(i);
      pir_eta[1] = JpsiTau_tau_pi2_eta->at(i);
      pir_eta[2] = JpsiTau_tau_pi3_eta->at(i);
      pir_phi[0] = JpsiTau_tau_pi1_phi->at(i);
      pir_phi[1] = JpsiTau_tau_pi2_phi->at(i);
      pir_phi[2] = JpsiTau_tau_pi3_phi->at(i);
      pi_flag = false;

      for(const auto& in: num_comb){
        if(JpsiTau_gen_tau_nprong->size()!=0 && int(JpsiTau_gen_tau_nprong->at(0) == 3))
        {
          pi_flag = true;
          Float_t dR = 0;
          Float_t dR_sum = 0;
          for(int j=0; j < 3; j++)
          {
            dR = getDeltaR(pi_eta[in.at(j)], pi_phi[in.at(j)], pir_eta[j], pir_phi[j]);
            dR_sum += dR;
            if (dR > 0.05)
            {
              break;
            }
          }
          if(dR_sum < dR_min)
          {
            dR_min = dR_sum;
            min_pos = i;
            num_comb_min[0]=in.at(0);
            num_comb_min[1]=in.at(1);
            num_comb_min[2]=in.at(2);
          }
        }
      }
    }

    // Making it match
    Float_t temp_transfer[3];
    temp_transfer[0] = pi_pt[num_comb_min[0]];
    temp_transfer[1] = pi_pt[num_comb_min[1]];
    temp_transfer[2] = pi_pt[num_comb_min[2]];
    pi_pt[0] = temp_transfer[0];
    pi_pt[1] = temp_transfer[1];
    pi_pt[2] = temp_transfer[2];
    temp_transfer[0] = pi_eta[num_comb_min[0]];
    temp_transfer[1] = pi_eta[num_comb_min[1]];
    temp_transfer[2] = pi_eta[num_comb_min[2]];
    pi_eta[0] = temp_transfer[0];
    pi_eta[1] = temp_transfer[1];
    pi_eta[2] = temp_transfer[2];
    temp_transfer[0] = pi_phi[num_comb_min[0]];
    temp_transfer[1] = pi_phi[num_comb_min[1]];
    temp_transfer[2] = pi_phi[num_comb_min[2]];
    pi_phi[0] = temp_transfer[0];
    pi_phi[1] = temp_transfer[1];
    pi_phi[2] = temp_transfer[2];


    if(JpsiTau_gen_tau_nprong->size()!=0 && int(JpsiTau_gen_tau_nprong->at(0) == 3))
    {
      for(int i=0; i < JpsiTau_tau_pi1_pt->size(); i++){
        pir_pt[0] = JpsiTau_tau_pi1_pt->at(i);
        pir_pt[1] = JpsiTau_tau_pi2_pt->at(i);
        pir_pt[2] = JpsiTau_tau_pi3_pt->at(i);
        pir_eta[0] = JpsiTau_tau_pi1_eta->at(i);
        pir_eta[1] = JpsiTau_tau_pi2_eta->at(i);
        pir_eta[2] = JpsiTau_tau_pi3_eta->at(i);
        pir_phi[0] = JpsiTau_tau_pi1_phi->at(i);
        pir_phi[1] = JpsiTau_tau_pi2_phi->at(i);
        pir_phi[2] = JpsiTau_tau_pi3_phi->at(i);
        pir_q[0] = JpsiTau_tau_pi1_q->at(i);
        pir_q[1] = JpsiTau_tau_pi2_q->at(i);
        pir_q[2] = JpsiTau_tau_pi3_q->at(i);
        pir_dce3d[0] = JpsiTau_tau_pi1_doca3d->at(i);
        pir_dce3d[1] = JpsiTau_tau_pi2_doca3d->at(i);
        pir_dce3d[2] = JpsiTau_tau_pi3_doca3d->at(i);
        pir_dz[0] = JpsiTau_tau_pi1_dz->at(i);
        pir_dz[1] = JpsiTau_tau_pi1_dz->at(i);
        pir_dz[2] = JpsiTau_tau_pi1_dz->at(i);

        // To Do to sort all var in order of pirpt
        if (pir_pt[0] < pir_pt[1])
        {
          SwapValue(pi_pt[0], pi_pt[1]);
          SwapValue(pi_eta[0], pi_eta[1]);
          SwapValue(pi_phi[0], pi_phi[1]);
          SwapValue(pir_pt[0], pir_pt[1]);
          SwapValue(pir_eta[0], pir_eta[1]);
          SwapValue(pir_phi[0], pir_phi[1]);
          SwapValue(pir_q[0], pir_q[1]);
          SwapValue(pir_dce3d[0], pir_dce3d[1]);
          SwapValue(pir_dz[0], pir_dz[1]);
        }
        if (pir_pt[0] < pir_pt[2])
        {
          SwapValue(pi_pt[0], pi_pt[2]);
          SwapValue(pi_eta[0], pi_eta[2]);
          SwapValue(pi_phi[0], pi_phi[2]);
          SwapValue(pir_pt[0], pir_pt[2]);
          SwapValue(pir_eta[0], pir_eta[2]);
          SwapValue(pir_phi[0], pir_phi[2]);
          SwapValue(pir_q[0], pir_q[2]);
          SwapValue(pir_dce3d[0], pir_dce3d[2]);
          SwapValue(pir_dz[0], pir_dz[2]);
        }
        if (pir_pt[1] < pir_pt[2])
        {
          SwapValue(pi_pt[1], pi_pt[2]);
          SwapValue(pi_eta[1], pi_eta[2]);
          SwapValue(pi_phi[1], pi_phi[2]);
          SwapValue(pir_pt[1], pir_pt[2]);
          SwapValue(pir_eta[1], pir_eta[2]);
          SwapValue(pir_phi[1], pir_phi[2]);
          SwapValue(pir_q[1], pir_q[2]);
          SwapValue(pir_dce3d[1], pir_dce3d[2]);
          SwapValue(pir_dz[1], pir_dz[2]);
        }

        if(i == min_pos) 
          pi_flag = true; 
        else 
          pi_flag = false;

        tree1->Fill();
      }
    
    }

    // sinec above 1 prong and 2 prong messeages were never printed it is
    // it is taken to be same

    for(int i=0; i < JpsiTau_tau_pi1_pt->size(); i++){
      h1ratio->Fill(JpsiTau_tau_pi1_pt->at(i)/JpsiTau_tau_pt->at(i));
    }

    for(int i=0; i < JpsiTau_tau_pi1_pt->size(); i++){
      if(JpsiTau_tau_pt->at(i) < 5.)
      {
        h15pt->Fill(JpsiTau_tau_max_dr_3prong->at(i));
      }
      else if(JpsiTau_tau_pt->at(i) < 10.)
      {
        h110pt->Fill(JpsiTau_tau_max_dr_3prong->at(i));
      }
      else if(JpsiTau_tau_pt->at(i) < 20.)
      {
        h120pt->Fill(JpsiTau_tau_max_dr_3prong->at(i));
      }
      else
      {
        h1expt->Fill(JpsiTau_tau_max_dr_3prong->at(i));
      }
    }

  } // end loop on events

  FileOut->cd();
  h1pttau->Write();
  h1etatau->Write();
  h1phitau->Write(); 
  h1nprongtau->Write();
  h1nmatchedtau->Write();
  h1ptpi->Write();
  h1etapi->Write();
  h1phipi->Write();
  h1pt1->Write();
  h1pt2->Write();
  h1pt3->Write();
  h1ratio->Write();
  h15pt->Write();
  h110pt->Write();
  h120pt->Write();
  h1expt->Write();
  h2etapt->Write();
  h2phipt->Write();
  h2etaphi->Write();
  h2dRpt->Write();
  tree1->Write();
}

float MyTauClass::deltaPhi(float phi1, float phi2) {
  float result = phi1 - phi2;
  while (result > TMath::Pi()) result -= float(2.*TMath::Pi());
  while (result <= -TMath::Pi()) result += float(2.*TMath::Pi());
  float absresult=TMath::Abs(result);
  return absresult;
}

Float_t MyTauClass::getDeltaR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2)
{
  Float_t dphi = deltaPhi(phi1, phi2);
  return sqrt(dphi*dphi + (eta1-eta2)*(eta1-eta2));
}


void MyTauClass::MegaLoop(const char *str)
{
  gSystem->Exec("mkdir files");
  TSystemDirectory dir(str, str); 
  TList *files = dir.GetListOfFiles(); 
  if (files) { 
    TSystemFile *file; 
    TString fname; 
    TIter next(files); 
    while ((file=(TSystemFile*)next())) 
    { 
      fname = file->GetName(); 
      if (!file->IsDirectory() && fname.EndsWith(".root")) {
        TString filename = str;
        filename += fname;

        Int_t dot = fname.First('.');
        Int_t len = fname.Length();
        fname.Remove(dot,len-dot);
        Int_t under = fname.Last('_');
        fname.Remove(0,under+1);

        FileNameIn = filename;

        const char* fileo="files/histos_";
        TString fileout = fileo;
        fileout += fname;
        fileout += ".root";
        FileNameOut = fileout;

        FileOut = new TFile(fileout, "RECREATE");
        tree1=NULL;
        InitOut();

        FileIn = new TFile(filename);
        TTree *tree=NULL;

        TDirectory *dir = (TDirectory*)FileIn->Get("ntuplizer");
        dir->GetObject("tree",tree);

        Init(tree);
        Loop();
        delete tree;
        // delete tree1;
        // delete FileOut;
        FileIn->Close();
      }
    }
  }
}

void MyTauClass::show(const char * file="final.root")
{
// Book here my histograms
  TFile *f = new TFile(file);
  TH1F *h1pttau = (TH1F*)f->Get("h1pttau");
  TH1F *h1etatau = (TH1F*)f->Get("h1etatau");
  TH1F *h1phitau = (TH1F*)f->Get("h1phitau");
  TH1F *h1nprongtau = (TH1F*)f->Get("h1nprongtau");

  TH1F *h1ptpi = (TH1F*)f->Get("h1ptpi");
  TH1F *h1etapi = (TH1F*)f->Get("h1etapi");
  TH1F *h1phipi = (TH1F*)f->Get("h1phipi");
  

  TH1F *h1pt1 = (TH1F*)f->Get("h1pt1");
  TH1F *h1pt2 = (TH1F*)f->Get("h1pt2");
  TH1F *h1pt3 = (TH1F*)f->Get("h1pt3");

  TH1F *h1ratio = (TH1F*)f->Get("h1ratio");
  TH1F *h15pt = (TH1F*)f->Get("h15pt");
  TH1F *h110pt = (TH1F*)f->Get("h110pt");
  TH1F *h120pt =  (TH1F*)f->Get("h120pt");
  TH1F *h1expt = (TH1F*)f->Get("h1expt");

  // PART 2
  TH2F *h2etapt = (TH2F*)f->Get("h2etapt");
  TH2F *h2phipt = (TH2F*)f->Get("h2phipt");
  TH2F *h2etaphi = (TH2F*)f->Get("h2etaphi");
  TH2F *h2dRpt = (TH2F*)f->Get("h2dRpt");


  // TH1F *h1nprong = new TH1F("h1nprong","nprong of the tau",4,0,4);
// TH1F *h1nprongpi = new TH1F("h1nprongtau","nprong of the gen pi",4,0,4);
  //  TH1I *h1npfCand = new TH1I("h1npfCand","number of pfCand",100,0.,100.);
  //  TH1F *h1pfCandpt = new TH1F("h1pfCandpt","pt of the PF candidate",100,0.,40.);
  h1pttau->GetXaxis()->SetTitle("pt[GeV]");
  h1etatau->GetXaxis()->SetTitle("#eta");
  h1phitau->GetXaxis()->SetTitle("#phi");
  h1nprongtau->GetXaxis()->SetTitle("nprongs");

  h1ptpi->GetXaxis()->SetTitle("pt[GeV]");
  h1etapi->GetXaxis()->SetTitle("#eta");
  h1phipi->GetXaxis()->SetTitle("#phi");
  h1pt3->GetXaxis()->SetTitle("pt[GeV]");

  h2etaphi->GetXaxis()->SetTitle("#phi");
  h2etaphi->GetYaxis()->SetTitle("#eta");
  h2etapt->GetXaxis()->SetTitle("pt[GeV]");
  h2etapt->GetYaxis()->SetTitle("#eta");
  h2phipt->GetXaxis()->SetTitle("pt[GeV]");
  h2phipt->GetYaxis()->SetTitle("#phi");
  // start of the loop
  //Plot
  TCanvas *c1=new TCanvas("c1","GeneralCanvas",200,10,1400,1000);
  c1->Divide(2,2);
  c1->cd(1);
  h1pttau->Draw();
  c1->cd(2);
  h1etatau->Draw();
  c1->cd(3);
  h1phitau->Draw();
  c1->cd(4);
  h1nprongtau->Draw();
  c1->SaveAs("firstplot.pdf");
  c1->SaveAs("firstplot.png");

  TCanvas *c2=new TCanvas("c2","SecondCanvas",300,20,1400,1000);
  c2->Divide(2,2);
  c2->cd(1);
  h1ptpi->Draw();
  c2->cd(2);
  h1etapi->Draw();
  c2->cd(3);
  h1phipi->Draw();
  c2->SaveAs("secondplot.pdf");
  c2->SaveAs("secondplot.png");

  gStyle->SetOptStat(kFALSE);
  TCanvas *c3=new TCanvas("c3","ThirdCanvas",300,20,1000,750);
  // h1pt1->SetFillColor(8);
  // h1pt2->SetFillColor(7);
  // h1pt3->SetFillColor(2);
  h1pt1->SetFillColorAlpha(8, 0.35);
  h1pt2->SetFillColorAlpha(7, 0.35);
  h1pt3->SetFillColorAlpha(2, 0.35);
  h1pt3->Draw("");
  h1pt2->Draw("SAME");
  h1pt1->Draw("SAME");
  h1pt3->SetTitle("Pions");
  TLegend *mylegend= new TLegend(.6,.6,.8,.8,"Pions");
  mylegend->AddEntry(h1pt1,"pt of pion 1", "f");
  mylegend->AddEntry(h1pt2,"pt of pion 2", "f");
  mylegend->AddEntry(h1pt3,"pt of pion 3", "f");
  mylegend->Draw("SAME");
  c3->SaveAs("thirdplot.pdf");
  c3->SaveAs("thirdplot.png");

  gStyle->SetOptStat(kTRUE);
  TCanvas *c4=new TCanvas("c4","FourthCanvas",300,20,1000,750);
  h1ratio->Draw();
  c4->SaveAs("fourthplot.pdf");
  c4->SaveAs("fourthplot.png");

  TCanvas *c5=new TCanvas("c5","FifthCanvas",300,20,1400,1000);
  c5->Divide(2,2);
  c5->cd(1);
  h15pt->Draw();
  c5->cd(2);
  h110pt->Draw();
  c5->cd(3);
  h120pt->Draw();
  c5->cd(4);
  h1expt->Draw();
  c5->SaveAs("fifthplot.pdf");
  c5->SaveAs("fifthplot.png");


  TCanvas *c6 = new TCanvas("c6", "SixthCanvas", 300,20,1000,750);
  h2etapt->Draw("COLZ");
  c6->SaveAs("sixthplot.pdf");
  c6->SaveAs("sixthplot.png");

  TCanvas *c7 = new TCanvas("c7", "SeventhCanvas", 300,20,1000,750);
  h2phipt->Draw("COLZ");
  c7->SaveAs("seventhplot.pdf");
  c7->SaveAs("seventhplot.png");

  TCanvas *c8 = new TCanvas("c8", "EightCanvas", 300,20,1000,750);
  h2etaphi->Draw("COLZ");
  c8->SaveAs("eightplot.pdf");
  c8->SaveAs("eightplot.png");

  TCanvas *c9 = new TCanvas("c9", "NinthCanvas", 300,20,1000,750);
  h2dRpt->Draw("COLZ");
  c9->SaveAs("ninthplot.pdf");
  c9->SaveAs("ninthplot.png");

}

void MyTauClass::SwapValue(Float_t &a, Float_t &b) 
{
  Float_t t = a;
  a = b;
  b = t;
}

void MyTauClass::SwapValue(Int_t &a, Int_t &b) 
{
  Int_t t = a;
  a = b;
  b = t;
}

void MyTauClass::CombAdd(vector<vector<Int_t>>& comb, Int_t arr[3], int n)
{
  do {
    vector<Int_t> temp;
    for(int i=0; i<n; i++) temp.push_back(arr[i]);
    comb.push_back(temp);
  } while (next_permutation(arr, arr + n));
}
