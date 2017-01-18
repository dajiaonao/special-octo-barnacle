#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <iostream>
#include <math.h>       /* pow */
// #include "Math/GSLMinimizer.h"
#ifndef __CINT__
#include "Minuit2/Minuit2Minimizer.h"
#endif
#include "Math/Functor.h"
#include <TLatex.h>
#include <vector>

class effFitter{
   TH1F* hP1;
   TH1F* hF1;
   TH1F* hP2;
   TH1F* hF2;
public:
   float eff;
   float effErr;
   float TF;
   float TFErr;

   double getChi(const double *xx){
    float e = xx[0];
    float s = xx[1];

    float eP = (1.-e)/e;

    TH1F hA(*hP1);
    hA.Scale(eP);
    hA.Add(hF2, s);

    TH1F hB(*hP2);
    hB.Add(hF1, s*eP);

//     std::cout << "e=" << e <<"; s=" << s <<  " val=" << hA.Chi2Test(&hB, "WW") << std::endl; 

    double r=0.;
    for(int i=1; i<hA.GetNbinsX()+1;i++){
      r += pow(hA.GetBinContent(i)-hB.GetBinContent(i),2)/(pow(hA.GetBinError(i),2)+pow(hB.GetBinError(i),2));
    }

    return r;
//     return hA.Chi2Test(&hB, "WWCHI2");
    }

   double getEff(TH1F* t_hP1, TH1F* t_hF1, TH1F* t_hP2, TH1F* t_hF2){
     hP1 = t_hP1;
     hF1 = t_hF1;
     hP2 = t_hP2;
     hF2 = t_hF2;

//      hP1->Sumw2();
//      hP2->Sumw2();
//      hF1->Sumw2();
//      hF2->Sumw2();

//      ROOT::Math::GSLSimAnMinimizer min;
     ROOT::Minuit2::Minuit2Minimizer min( ROOT::Minuit2::kMigrad );
     min.SetMaxFunctionCalls(1000000);
     min.SetMaxIterations(100000);
     min.SetTolerance(0.1);
   
     ROOT::Math::Functor f(this, &effFitter::getChi,2); 

     double step[2] = {0.01,0.01};
     double variable[2] = {0.98,1.};
   
     min.SetFunction(f);
   
     // Set the free variables to be minimized!
     min.SetVariable(0,"x",variable[0], step[0]);
     min.SetVariable(1,"y",variable[1], step[1]);

     min.SetVariableLimits(0, 0.6, 1.);
     min.SetVariableLimits(1, 0.5, 10);
   
     min.Minimize(); 
   
     const double *xs = min.X();
//      std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " 
//           << getChi(xs) << std::endl;

     eff = xs[0];
     TF = xs[1];
     effErr = sqrt(min.CovMatrix(0,0));
     TFErr = sqrt(min.CovMatrix(1,1));
   }
};


int run_HistFitting(TH1F* hP1, TH1F* hF1, TH1F* hP2, TH1F* hF2){
  std::cout << "testing" << std::endl;

//   for(int i=0; i<hP1->GetNbinsX()+2; i++){
//     std::cout << hP1->GetBinContent(i)
//               << " " << hP2->GetBinContent(i)
//               << " " << hF1->GetBinContent(i)
//               << " " << hF2->GetBinContent(i)
//               << std::endl;
//    }

  float s = 0.8;
  float e = 0;

  for(int i=0; i<20; i++){
    e = 0.81+0.01*i;
    float eP = (1.-e)/e;

    TH1F hA(*hP1);
    hA.Scale(eP);
    hA.Add(hF2, s);

    TH1F hB(*hP2);
    hB.Add(hF1, s*eP);

    float chi2 = hA.Chi2Test(&hB, "WWCHI2");
    std::cout << i << " e=" << e << " chi2=" << chi2 << std::endl;
  }


  return 0;
}


int HistMore(){
  TFile fPass("Dist_Bkg_Mll_passIso.root","read");
  TFile fFail("Dist_Bkg_Mll_failIso.root", "read");

  TFile fout1("fout1.root","recreate");
  TGraphErrors gEff;
  TGraphErrors gTF;
  float pt0 = 4.5;

  effFitter j;
  std::vector< std::string > ptBins{"4_5","5_6","6_7","7_8","8_9","9_10","10_11","111_12","12_13","13_14","14_15"};
  for(auto x: ptBins){
    TCanvas* cav1 = (TCanvas*)fPass.Get(("Dist_Bkg_Data_ZMass_"+x).c_str());
    TH1F* hP1 = (TH1F*) cav1->GetPrimitive("eff_Data_Draw");
    TH1F* hF1 = (TH1F*) cav1->GetPrimitive("eff_bkg_Draw");

    TCanvas* cav2 = (TCanvas*)fFail.Get(("Dist_Bkg_Data_ZMass_"+x).c_str());
    TH1F* hP2 = (TH1F*) cav2->GetPrimitive("eff_Data_Draw");
    TH1F* hF2 = (TH1F*) cav2->GetPrimitive("eff_bkg_Draw");

    j.getEff(hP1, hF1, hP2, hF2);
    std::cout << x << " e=" << j.eff << "+/-" << j.effErr << " s=" << j.TF << "+/-" << j.TFErr << std::endl;

    int n = gEff.GetN();
    gEff.SetPoint(n, pt0, j.eff);
    gEff.SetPointError(n, 0.5, j.effErr);
    gTF.SetPoint(n, pt0, j.TF);
    gTF.SetPointError(n, 0.5, j.TFErr);

    pt0 += 1;
   }

  gEff.Write("eff");
  gTF.Write("TF");

  return 0;
}


int HistFitting(){

//   TFile fPass("/home/dzhang/links/download/Dist_Bkg_ptCut56_NoMllCut_passIso.root","read");
//   TFile fFail("/home/dzhang/links/download/Dist_Bkg_ptCut56_NoMllCut_failIso.root", "read");
  TFile fPass("Dist_Bkg_Mll_passIso.root","read");
  TFile fFail("Dist_Bkg_Mll_failIso.root", "read");
  std::string x("_13_14");
  TString lx("13<p_{T}<14 GeV");

  TCanvas* cav1 = (TCanvas*)fPass.Get(("Dist_Bkg_Data_ZMass"+x).c_str());
  TH1F* hP1 = (TH1F*) cav1->GetPrimitive("eff_Data_Draw");
  TH1F* hF1 = (TH1F*) cav1->GetPrimitive("eff_bkg_Draw");

  TCanvas* cav2 = (TCanvas*)fFail.Get(("Dist_Bkg_Data_ZMass"+x).c_str());
  TH1F* hP2 = (TH1F*) cav2->GetPrimitive("eff_Data_Draw");
  TH1F* hF2 = (TH1F*) cav2->GetPrimitive("eff_bkg_Draw");

//   return run_HistFitting(hP1, hF1, hP2, hF2);
//   effFitter j;
//   j.getEff(hP1, hF1, hP2, hF2);


//   0.910141,4.6401
//   0.967004,1.40231
// 4_5 e=0.959668+/-0.00641386 s=1.36733+/-0.00890238
// 5_6 e=0.967004+/-0.00491406 s=1.40231+/-0.0276626
// 6_7 e=0.985617+/-0.0020087 s=1.63938+/-0.0327642
// 7_8 e=0.957873+/-0.00280571 s=1.70809+/-0.0481239
// 8_9 e=0.964739+/-0.00231842 s=1.87594+/-0.0654838
// 9_10 e=0.969996+/-0.0020627 s=2.11101+/-0.0710226
// 10_11 e=0.995463+/-0.000865506 s=2.56746+/-0.179051
// 111_12 e=0.997166+/-0.000542415 s=2.87417+/-0.257211
// 12_13 e=0.99685+/-0.000325247 s=3.0257+/-0.273914
// 13_14 e=0.998426+/-0.000414012 s=4.08638+/-0.441891
// 14_15 e=0.998169+/-0.000422332 s=7.51555+/-1.2246//
  float e = 0.998426;
  float s = 4.08638;
  float eP = (1.-e)/e;

  TLatex lt;

  TCanvas* cxx1 = new TCanvas();
  cxx1->cd();

  TH1F* hA = (TH1F*)hP1->Clone("hA");
  hA->Scale(eP);
  hA->Add(hF2, s);

  TH1F* hB = (TH1F*)hP2->Clone("hB");
  hB->Add(hF1, s*eP);

  hA->SetLineColor(2);
  hA->SetMarkerColor(2);
  hA->SetMarkerStyle(24);
  hA->DrawCopy();
  hB->DrawCopy("same");
  lt.DrawLatexNDC(0.7,0.8, lx);

  cxx1->Update();

  TCanvas* cxx2 = new TCanvas();
  cxx2->cd();
  hP1->DrawCopy();

  TH1F* hC = (TH1F*)hF1->Clone("hC");
  hC->Scale(s);
  hC->SetLineColor(2);
  hC->SetMarkerColor(2);
  hC->SetMarkerStyle(24);
  hC->DrawCopy("same");

  TH1F* hP1c=(TH1F*)hP1->Clone("hP1c");
  hP1c->Add(hC,-1);

  hP1c->SetLineColor(4);
  hP1c->SetMarkerColor(4);
  hP1c->SetMarkerStyle(25);
  hP1c->DrawCopy("same");

  lt.DrawLatexNDC(0.7,0.8, lx);
  cxx2->Update();

  TCanvas* cxx3 = new TCanvas();
  cxx3->cd();
  hP2->DrawCopy();

  TH1F* hD = (TH1F*)hF2->Clone("hD");
  hD->Scale(s);
  hD->SetLineColor(2);
  hD->SetMarkerColor(2);
  hD->SetMarkerStyle(24);
  hD->DrawCopy("same");

  lt.DrawLatexNDC(0.7,0.8, lx);
  cxx3->Update();

  TCanvas* cxx4 = new TCanvas();
  cxx4->cd();
//   hP2->DrawCopy();
  hF1->SetLineColor(2);
  hF1->SetMarkerColor(2);
  hF1->SetMarkerStyle(24);
  hF2->DrawNormalized();
  hF1->DrawNormalized("same");

//   TH1F* hFx = (TH1F*)hF2->Clone("hFx");
//   hFx->Divide(hF1);
//   hFx->SetLineColor(4);
//   hFx->SetMarkerColor(4);
//   hFx->SetMarkerStyle(26);
//   hFx->DrawCopy("same");

  lt.DrawLatexNDC(0.7,0.8, lx);
  cxx4->Update();

//   TH1F h1("h1","h1",10,0,10);
//   TH1F* h2 = (TH1F*)h1.Clone("h2");
// 
//   h1.Fill(3);

//   for(int i=0; i<10;i++) h2->Fill(3);
//   for(int i=0; i<100;i++){h1.Fill(i); h2->Fill(i);}
//   for(int i=0; i<10000;i++){h2->Fill(i);}

//   std::cout << "chi2 = " << h1.Chi2Test(h2,"UU") << std::endl;

  return 0;
}

int main(){
  return HistFitting();
//   return HistMore();
}
