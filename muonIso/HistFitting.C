/// Use locations:
// /home/dzhang/work/muons/isolation/iso207/testing

// compile: g++ `root-config --cflags` `root-config --glibs` -lMinuit2 HistFitting.C -o HistFitting
// Compile to a library

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
#include <string>

class effFitter{
public:
   TH1F* hOS1;
   TH1F* hSS1;
   TH1F* hOS2;
   TH1F* hSS2;
   TH1F* hMC1;
   TH1F* hMC2;
   TH1F* hMCr;
   float eff;
   float effErr;
   float TF;
   float TFErr;
   TString fitMessage;

   void showValues(){
     std::cout << "eff = " << eff << "+/-" << effErr << "; TF = " << TF << "+/-" << TFErr << std::endl;
   }

   void setMC(TH1F* h1, TH1F* h2){hMC1 = h1; hMC2 = h2;}


   void showHists(TString savename="test"){
    float e = eff;
    float s = TF;

    TH1F e1(*hMC1);
    TH1F hMC0(*hMC1);
    hMC0.Add(hMC2);
    e1.Divide(&hMC0);

    /// multiply scalar efficiency
    std::cout << e1.GetMaximumBin() << "->" << e1.GetBinContent(e1.GetMaximumBin()) << " " << e << std::endl;
    e1.Scale(e/e1.GetBinContent(e1.GetMaximumBin()));

    TH1F e1m(e1); /// 1-eff
    for(int i=0;i<e1m.GetNbinsX()+2;i++){
      e1m.SetBinContent(i, 1.-e1m.GetBinContent(i));
     }

    /* The method:  ssPass*TF*(1-e) + osFail*e = osPass*(1-e) + ssFail*T*e is the formular*/
    TH1F hA1(*hSS1);
    hA1.Multiply(&e1m);
    hA1.Scale(s);

    TH1F hA2(*hOS2);
    hA2.Multiply(&e1);

    /// LHS
    TH1F hA(hA1);
    hA.Add(&hA2);

    TH1F hB1(*hOS1);
    hB1.Multiply(&e1m);

    TH1F hB2(*hSS2);
    hB2.Multiply(&e1);
    hB2.Scale(s);

    /// RHS
    TH1F hB(hB1);
    hB.Add(&hB2);

    TCanvas cv;
    cv.Divide(3,2);
    cv.cd(1);
    hA1.Draw();

    cv.cd(2);
    hA2.Draw();

    cv.cd(3);
    hA.Draw();

    cv.cd(4);
    hB1.Draw();

    cv.cd(5);
    hB2.Draw();

    cv.cd(6);
    hB.Draw();

    cv.SaveAs(savename+".png");
    
    return;
    }

   double getChi(const double *xx){
    float e = xx[0];
    float s = xx[1];
    float eP = (1.-e)/e;

    TH1F e1(*hMC1);
    TH1F hMC0(*hMC1);
    hMC0.Add(hMC2);
    e1.Divide(&hMC0);

    /// multiply scalar efficiency
    e1.Scale(e/e1.GetBinContent(e1.GetMaximumBin()));

    TH1F e1m(e1);
    for(int i=0;i<e1m.GetNbinsX()+2;i++){
      e1m.SetBinContent(i, 1.-e1m.GetBinContent(i));
     }

    /* The method:  ssPass*TF*(1-e) + osFail*e = osPass*(1-e) + ssFail*T*e is the formular*/
    TH1F hA(*hSS1);
    hA.Multiply(&e1m);
    hA.Scale(s);

    TH1F hA2(*hOS2);
    hA2.Multiply(&e1);

    /// LHS
    hA.Add(&hA2);

    TH1F hB(*hOS1);
    hB.Multiply(&e1m);

    TH1F hB2(*hSS2);
    hB2.Multiply(&e1);
    hB2.Scale(s);

    /// RHS
    hB.Add(&hB2);

    /// fit
    double r=0.;
    for(int i=1; i<hA.GetNbinsX()+1;i++){
      r += pow(hA.GetBinContent(i)-hB.GetBinContent(i),2)/(pow(hA.GetBinError(i),2)+pow(hB.GetBinError(i),2));
    }

    return r;
    }


   double getChiS(const double *xx){
    float e = xx[0];
    float s = xx[1];

    float eP = (1.-e)/e;

    TH1F hA(*hOS1);
    hA.Scale(eP);
    hA.Add(hSS2, s);

    TH1F hB(*hOS2);
    hB.Add(hSS1, s*eP);

//     std::cout << "e=" << e <<"; s=" << s <<  " val=" << hA.Chi2Test(&hB, "WW") << std::endl; 

    double r=0.;
    for(int i=1; i<hA.GetNbinsX()+1;i++){
      r += pow(hA.GetBinContent(i)-hB.GetBinContent(i),2)/(pow(hA.GetBinError(i),2)+pow(hB.GetBinError(i),2));
    }

    return r;
//     return hA.Chi2Test(&hB, "WWCHI2");
    }

   double getEff(TH1F* t_hOS1, TH1F* t_hSS1, TH1F* t_hOS2, TH1F* t_hSS2, std::string checkMode=""){
     hOS1 = t_hOS1;
     hSS1 = t_hSS1;
     hOS2 = t_hOS2;
     hSS2 = t_hSS2;

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

     min.SetVariableLimits(0, 0.6, 1);
     min.SetVariableLimits(1, 0.5, 10);
   
     min.Minimize(); 
   
     const double *xs = min.X();
//      std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " 
//           << getChi(xs) << std::endl;

     eff = xs[0];
     TF = xs[1];
     effErr = sqrt(min.CovMatrix(0,0));
     TFErr = sqrt(min.CovMatrix(1,1));

     if(checkMode != ""){
       checkEff1(checkMode);
      }

     return 0;
   }

  double checkEff1(TString savename="noSave"){

    float e = eff;
    float s = TF;

    TH1F e1(*hMC1);
    TH1F hMC0(*hMC1);
    hMC0.Add(hMC2);
    e1.Divide(&hMC0);

    /// multiply scalar efficiency
    e1.Scale(e/e1.GetBinContent(e1.GetMaximumBin()));

    TH1F e1m(e1); /// 1-eff
    for(int i=0;i<e1m.GetNbinsX()+2;i++){
      e1m.SetBinContent(i, 1.-e1m.GetBinContent(i));
     }

    /* The method:  ssPass*TF*(1-e) + osFail*e = osPass*(1-e) + ssFail*T*e is the formular*/
    TH1F hA(*hSS1);
    hA.Multiply(&e1m);
    hA.Scale(s);

    TH1F hA2(*hOS2);
    hA2.Multiply(&e1);

    /// LHS
    hA.Add(&hA2);

    TH1F hB(*hOS1);
    hB.Multiply(&e1m);

    TH1F hB2(*hSS2);
    hB2.Multiply(&e1);
    hB2.Scale(s);

    /// RHS
    hB.Add(&hB2);

/////////////////////////
    TLatex lt;

    TCanvas* cxx1 = new TCanvas();
    cxx1->Divide(3,2);

    /// check plot 1: the passed ones, The fit histograms
    cxx1->cd(1);

    hA.SetLineColor(2);
    hA.SetMarkerColor(2);
    hA.SetMarkerStyle(24);
    hA.DrawCopy();
    hB.DrawCopy("same");
    lt.DrawLatexNDC(0.5,0.9, fitMessage);
    lt.DrawLatexNDC(0.5,0.8, "Fit histograms");

    /// check plot 2: show signal, with estimated background subtracted 
    auto pad2 = cxx1->cd(2);
    hOS1->DrawCopy();

    TH1F* hC = (TH1F*)hSS1->Clone("hC");
    hC->Scale(s);
    hC->SetLineColor(2);
    hC->SetMarkerColor(2);
    hC->SetMarkerStyle(24);
    hC->DrawCopy("same");

    TH1F* hOS1c=(TH1F*)hOS1->Clone("hOS1c");
    hOS1c->Add(hC,-1);

    hOS1c->SetLineColor(4);
    hOS1c->SetMarkerColor(4);
    hOS1c->SetMarkerStyle(25);
    hOS1c->DrawCopy("same");

    if(hMC1){
      TH1F* hMC_t = (TH1F*)hMC1->Clone("hMC_temp");
      hMC_t->Scale(hOS1c->Integral()/hMC1->Integral());
      hMC_t->SetLineColor(6);
      hMC_t->Draw("same");

      cxx1->cd(5);
      TH1F* hOS1c_t = (TH1F*)hOS1c->Clone("hOS1c_t");
      hOS1c_t->Divide(hMC_t);
      hOS1c_t->Draw();
      hOS1c_t->GetYaxis()->SetRangeUser(0.5,1);

      if(hMC2){
        cxx1->cd(6);
        hMC1->Draw();
        hMC2->SetLineStyle(2);
        hMC2->Draw("same");

        cxx1->cd(5);
        TH1F* hMC_t3 = (TH1F*)hMC1->Clone("hMC_temp3");
        TH1F* hMC_t4 = (TH1F*)hMC1->Clone("hMC_temp4");
        hMC_t3->Add(hMC2);
        hMC_t4->Divide(hMC_t3);
        hMC_t4->Draw("same");
//         gPad->Update();
      }

      cxx1->cd(2);
     }

    pad2->SetLogy();
    lt.DrawLatexNDC(0.6,0.8, "Bkg subtraction");

    /// check plot 3
    cxx1->cd(3);
    hOS2->DrawCopy();

    TH1F* hD = (TH1F*)hSS2->Clone("hD");
    hD->Scale(s);
    hD->SetLineColor(2);
    hD->SetMarkerColor(2);
    hD->SetMarkerStyle(24);
    hD->DrawCopy("same");
    lt.DrawLatexNDC(0.6,0.8, "SS pass/fail");

    /// check plot 4
    cxx1->cd(4);
    hSS1->SetLineColor(2);
    hSS1->SetMarkerColor(2);
    hSS1->SetMarkerStyle(24);
    hSS2->DrawNormalized();
    hSS1->DrawNormalized("same");

    lt.DrawLatexNDC(0.6,0.8, "Bkg shapes");
    cxx1->Update();

    /// save the plots
    if(savename!="noSave"){
      cxx1->SaveAs(savename+".eps");
      cxx1->SaveAs(savename+".pdf");
      cxx1->SaveAs(savename+".png");
     }

    return 0;
   }

//   double checkEff(std::string savename="noSave"){
  double checkEff(TString savename="noSave"){
    float e = eff;
    float s = TF;
    float eP = (1.-e)/e;

    TLatex lt;

    TCanvas* cxx1 = new TCanvas();
    cxx1->Divide(3,2);

    /// check plot 1: the passed ones, The fit histograms
    cxx1->cd(1);

    TH1F* hA = (TH1F*)hOS1->Clone("hA");
    hA->Scale(eP);
    hA->Add(hSS2, s);

    TH1F* hB = (TH1F*)hOS2->Clone("hB");
    hB->Add(hSS1, s*eP);

    hA->SetLineColor(2);
    hA->SetMarkerColor(2);
    hA->SetMarkerStyle(24);
    hA->DrawCopy();
    hB->DrawCopy("same");
    lt.DrawLatexNDC(0.6,0.9, fitMessage);
    lt.DrawLatexNDC(0.6,0.8, "Fit histograms");

    /// check plot 2: show signal, with estimated background subtracted 
    auto pad2 = cxx1->cd(2);
    hOS1->DrawCopy();

    TH1F* hC = (TH1F*)hSS1->Clone("hC");
    hC->Scale(s);
    hC->SetLineColor(2);
    hC->SetMarkerColor(2);
    hC->SetMarkerStyle(24);
    hC->DrawCopy("same");

    TH1F* hOS1c=(TH1F*)hOS1->Clone("hOS1c");
    hOS1c->Add(hC,-1);

    hOS1c->SetLineColor(4);
    hOS1c->SetMarkerColor(4);
    hOS1c->SetMarkerStyle(25);
    hOS1c->DrawCopy("same");

    if(hMC1){
      TH1F* hMC_t = (TH1F*)hMC1->Clone("hMC_temp");
      hMC_t->Scale(hOS1c->Integral()/hMC1->Integral());
      hMC_t->SetLineColor(6);
      hMC_t->Draw("same");

      cxx1->cd(5);
      TH1F* hOS1c_t = (TH1F*)hOS1c->Clone("hOS1c_t");
      hOS1c_t->Divide(hMC_t);
      hOS1c_t->Draw();
      hOS1c_t->GetYaxis()->SetRangeUser(0.5,1);

      if(hMC2){
        cxx1->cd(6);
        hMC1->Draw();
        hMC2->SetLineStyle(2);
        hMC2->Draw("same");

        cxx1->cd(5);
        TH1F* hMC_t3 = (TH1F*)hMC1->Clone("hMC_temp3");
        TH1F* hMC_t4 = (TH1F*)hMC1->Clone("hMC_temp4");
        hMC_t3->Add(hMC2);
        hMC_t4->Divide(hMC_t3);
        hMC_t4->Draw("same");
//         gPad->Update();
      }

      cxx1->cd(2);
     }

    pad2->SetLogy();
    lt.DrawLatexNDC(0.6,0.8, "Bkg subtraction");

    /// check plot 3
    cxx1->cd(3);
    hOS2->DrawCopy();

    TH1F* hD = (TH1F*)hSS2->Clone("hD");
    hD->Scale(s);
    hD->SetLineColor(2);
    hD->SetMarkerColor(2);
    hD->SetMarkerStyle(24);
    hD->DrawCopy("same");
    lt.DrawLatexNDC(0.6,0.8, "SS pass/fail");

    /// check plot 4
    cxx1->cd(4);
    hSS1->SetLineColor(2);
    hSS1->SetMarkerColor(2);
    hSS1->SetMarkerStyle(24);
    hSS2->DrawNormalized();
    hSS1->DrawNormalized("same");

    lt.DrawLatexNDC(0.6,0.8, "Bkg shapes");
    cxx1->Update();

    /// save the plots
    if(savename!="noSave"){
      cxx1->SaveAs(savename+".eps");
      cxx1->SaveAs(savename+".pdf");
      cxx1->SaveAs(savename+".png");
     }

    return 0;
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

    j.fitMessage = x;
    j.getEff(hP1, hF1, hP2, hF2, "figs_2017Jan18/fit_"+x);
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
  effFitter j;
  j.getEff(hP1, hF1, hP2, hF2, "testing");
  return 0;


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
  cxx1->Divide(2,2);
  cxx1->cd(1);

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

//   TCanvas* cxx2 = new TCanvas();
//   cxx2->cd();
  cxx1->cd(2);
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
//   cxx2->Update();

  cxx1->cd(3);
//   TCanvas* cxx3 = new TCanvas();
//   cxx3->cd();
  hP2->DrawCopy();

  TH1F* hD = (TH1F*)hF2->Clone("hD");
  hD->Scale(s);
  hD->SetLineColor(2);
  hD->SetMarkerColor(2);
  hD->SetMarkerStyle(24);
  hD->DrawCopy("same");

  lt.DrawLatexNDC(0.7,0.8, lx);
//   cxx3->Update();

  cxx1->cd(4);
//   TCanvas* cxx4 = new TCanvas();
//   cxx4->cd();
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
//   cxx4->Update();

//   TH1F h1("h1","h1",10,0,10);
//   TH1F* h2 = (TH1F*)h1.Clone("h2");
// 
//   h1.Fill(3);

//   for(int i=0; i<10;i++) h2->Fill(3);
//   for(int i=0; i<100;i++){h1.Fill(i); h2->Fill(i);}
//   for(int i=0; i<10000;i++){h2->Fill(i);}

//   std::cout << "chi2 = " << h1.Chi2Test(h2,"UU") << std::endl;

  cxx1->Update();

  return 0;
}

int main(){
//   return HistFitting();
  return HistMore();
}
