#include "Riostream.h"
#include "TGraph.h"
#include "TObjArray.h"
#include "TAxis.h"
#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TVirtualFitter.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFrame.h"
#include "Fit/Fitter.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TRandom3.h"

using namespace  std;

Int_t NGraphs = 0;
Int_t Colors[10] = {kMagenta,kOrange-3,kGreen,kBlue,kRed,kBlack,kCyan-3,kRed-7};

TCanvas *MyCanvas = nullptr;
TCanvas *MyCanvas2 = nullptr;

Int_t DateMin=0;
Int_t DateMax=0;

TLegend *legend = nullptr;

TH1 *firstHist = nullptr;
TH1 *firstHist_tot = nullptr;

Double_t MaxY = 0;
Double_t MaxY_Tot = 0;

vector<TString> vDates;
vector<Double_t> vDeaths;
vector<Double_t> vDeaths_e;
vector<Double_t> vDeaths_Tot;
map<TString,Int_t> mapofdates;

Bool_t DoEst = false;
Bool_t DoParis = false;

TF1 *fTotalD = nullptr;
TF1 *fTotalD2 = nullptr;
TF1 *fTotalD2Full = nullptr;

void AnaPerCountry(TString theCountry, Int_t Smooth = 3, TString FitUpTo = "", Bool_t DoD = true, Bool_t DoD2 = false, Bool_t DoD2Full = false, Bool_t DoPlotSub = false);

Double_t FuncD(Double_t*xx,Double_t*pp);
Double_t FuncD2(Double_t*xx,Double_t*pp);
Double_t FuncD2Full(Double_t*xx,Double_t*pp);

Float_t fChi2D = 0.;
Float_t fChi2D2 = 0.;
Float_t fChi2D2Full = 0.;

void AnaPerCountry(TString theCountry, Int_t Smooth, TString FitUpTo, Bool_t DoD, Bool_t DoD2, Bool_t DoD2Full, Bool_t DoPlotSub) {

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    if(theCountry.EqualTo("SA",TString::kIgnoreCase)) theCountry = "South_Africa";

    TString Folder = "./worldometers/";

//    TString Folder = "/Users/dudouet/Documents/Perso/Divers/Worldometers";

    ifstream file(Form("%s/%s.csv",Folder.Data(),theCountry.Data()));

    if(!file) {
        cout<<theCountry<<" not found in "<<Folder<<endl;
        return;
    }

    MaxY = 0;
    MaxY_Tot = 0;

    DateMin=0;
    DateMax=0;

    NGraphs = 0;

    MyCanvas = new TCanvas("daily","daily",1600,1200);
    MyCanvas->SetLeftMargin(0.107769);
    MyCanvas->SetRightMargin(0.00125313);
    MyCanvas->SetTopMargin(0.00173913);
    MyCanvas->SetBottomMargin(0.135652);

    TString Buffer;
    string line;

    vDates.clear();
    vDeaths.clear();
    vDeaths_e.clear();
    vDeaths_Tot.clear();
    mapofdates.clear();

    Int_t Mounth=1;
    Int_t Day=1;
    TString Mounth_str[12] = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
    Int_t NDaysPerMounth[12] = {31,29,31,30,31,30,31,31,30,31,30,31};
    Int_t NDaysInYear=0;
    for(int i=0 ; i<12 ; i++) NDaysInYear += NDaysPerMounth[i];

    getline(file,line);

    int NToRemove = 0;

    Int_t Current_Year = 20;

    while(file) {
        getline(file,line);
        Buffer = line;

        TObjArray *arr = nullptr;

        if(Buffer.Contains(";")) {
            Buffer.Append(";");
            Buffer.ReplaceAll(";;","; ;");
            Buffer.ReplaceAll(";;","; ;");

            arr = Buffer.Tokenize(";");
        }
        else if(Buffer.Contains(",")) {
            Buffer.Append(",");
            Buffer.ReplaceAll(",,",", ,");
            Buffer.ReplaceAll(",,",", ,");

            arr = Buffer.Tokenize(",");
        }
        else continue;

        // The 2nd value of the array (index 1), corresponds to the date
        TString Date = (TString)arr->At(1)->GetName();
        TObjArray *arr2 = Date.Tokenize("$");
        Date = arr2->First()->GetName();
        delete arr2;
        Date.ReplaceAll(" ","-");

        TObjArray *temp = Date.Tokenize("-");
        TString Mounth_tmp = (TString)temp->At(0)->GetName();
        Day = ((TString)temp->At(1)->GetName()).Atoi();
        Date = Form("%d-%s-%d",Day,Mounth_tmp.Data(),Current_Year);
        delete temp;

        if(Date.BeginsWith("31-Dec")) Current_Year++;

        Int_t Deaths = ((TString)arr->At(3)->GetName()).Atoi()-NToRemove;
        delete arr;
        if(Deaths) {
            vDates.push_back(Date);
            vDeaths_Tot.push_back(Deaths);
        }
    }

    if(theCountry == "US") theCountry = "USA";
    if(theCountry == "UK") theCountry = "United Kingdom";
    theCountry.ReplaceAll("_"," ");

    bool stoptakingdata=false;

    for(size_t i=0 ; i<vDeaths_Tot.size() ; i++) {
        if(stoptakingdata) {
            vDeaths_Tot.erase(vDeaths_Tot.begin()+i,vDeaths_Tot.end());
            vDates.erase(vDates.begin()+i,vDates.end());
        }
        else if(vDates.at(i) == FitUpTo) stoptakingdata = true;
    }

    if(vDeaths_Tot.empty()) return;

    Int_t DeathsMin = 10;

    while(vDeaths_Tot.front()<DeathsMin) {
        vDates.erase(vDates.begin());
        vDeaths_Tot.erase(vDeaths_Tot.begin());
    }

    vDeaths.push_back(vDeaths_Tot.front());
    for(size_t i=1 ; i<vDates.size() ; i++) {
        if(vDeaths_Tot.size()>i && vDeaths_Tot.at(i)>0) vDeaths.push_back(vDeaths_Tot.at(i)-vDeaths_Tot.at(i-1));
    }

    vector<Double_t> vDeaths_Smooth;
    for(size_t i=0 ; i<vDeaths_Tot.size() ; i++) {
        Double_t NPoints = 0;
        Double_t Tot = 0.;
        Double_t Err2 = 0.;

        bool stop = false;

        if(Smooth>=7){
            if(i>3) {
                Tot += vDeaths_Tot.at(i-3);
                Err2 += TMath::Power(2*sqrt(vDeaths_Tot.at(i-3)),2);
                NPoints ++;
            }
            else stop = true;
        }
        if(Smooth>=5){
            if(i>2) {
                Tot += vDeaths_Tot.at(i-2);
                Err2 += TMath::Power(2*sqrt(vDeaths_Tot.at(i-2)),2);
                NPoints ++;
            }
            else stop = true;
        }
        if(Smooth>=3){
            if(i>1) {
                Tot += vDeaths_Tot.at(i-1);
                Err2 += TMath::Power(2*sqrt(vDeaths_Tot.at(i-1)),2);
                NPoints ++;
            }
            else stop = true;
        }
        if(true) {
            Tot += vDeaths_Tot.at(i);
            Err2 += TMath::Power(2*sqrt(vDeaths_Tot.at(i)),2);
            NPoints ++;
        }
        if(Smooth>=3){
            if(i<(vDeaths_Tot.size()-1)) {
                Tot += vDeaths_Tot.at(i+1);
                Err2 += TMath::Power(2*sqrt(vDeaths_Tot.at(i+1)),2);
                NPoints ++;
            }
            else stop = true;
        }
        if(Smooth>=5){
            if(i<(vDeaths_Tot.size()-2)) {
                Tot += vDeaths_Tot.at(i+2);
                Err2 += TMath::Power(2*sqrt(vDeaths_Tot.at(i+2)),2);
                NPoints ++;
            }
            else stop = true;
        }
        if(Smooth>=7){
            if(i<(vDeaths_Tot.size()-3)) {
                Tot += vDeaths_Tot.at(i+3);
                Err2 += TMath::Power(2*sqrt(vDeaths_Tot.at(i+3)),2);
                NPoints ++;
            }
            else stop = true;
        }

        Tot = Tot/NPoints;

        if(Tot>0. && !stop) {
            vDeaths_Smooth.push_back(Tot);
            vDeaths_e.push_back(sqrt(Err2)/NPoints);
        }
        else {
            vDeaths_Smooth.push_back(0);
            vDeaths_e.push_back(0);
        }
    }

    for(size_t i=0 ; i<vDeaths_Tot.size() ; i++) vDeaths_Tot.at(i) = vDeaths_Smooth.at(i);

    TH1D *hDeces_Tot = new TH1D(Form("DecesTot_%s",theCountry.Data()),Form("DecesTot_%s",theCountry.Data()),NDaysInYear*2,0,NDaysInYear*2);
    hDeces_Tot->GetYaxis()->SetTitle("TOTAL DEATHS");
    hDeces_Tot->GetYaxis()->CenterTitle();
    hDeces_Tot->GetXaxis()->SetLabelSize(0.06);
    hDeces_Tot->GetXaxis()->SetTitleOffset(1.);
    hDeces_Tot->GetXaxis()->SetTitleFont(132);
    hDeces_Tot->GetXaxis()->SetLabelFont(132);

    hDeces_Tot->GetYaxis()->SetLabelSize(0.05);
    hDeces_Tot->GetYaxis()->SetTitleSize(0.05);
    hDeces_Tot->GetYaxis()->SetTitleOffset(1.15);
    hDeces_Tot->GetYaxis()->SetTickSize(0.01);
    hDeces_Tot->GetXaxis()->SetTickSize(0.01);
    hDeces_Tot->GetYaxis()->SetTitleFont(132);
    hDeces_Tot->GetYaxis()->SetLabelFont(132);

    hDeces_Tot->SetDirectory(nullptr);
    hDeces_Tot->SetMarkerStyle(20);
    hDeces_Tot->SetMarkerColor(kBlack);
    hDeces_Tot->SetLineColor(kBlack);
    hDeces_Tot->SetBinContent(1,0.001);
    hDeces_Tot->SetBinContent(NDaysInYear,0.001);

    Int_t ibin=1;
    for(int year=20 ; year<=21 ; year++) {
        for(int i=0 ; i<12 ; i++) {
            if(year==21 && i==1) NDaysPerMounth[i] = 28;
            for(int j=0 ; j<NDaysPerMounth[i] ; j++) {
                TString Label = Form("%d-%s-%d",j+1,Mounth_str[i].Data(),year);
                hDeces_Tot->GetXaxis()->SetBinLabel(ibin,Label);
                ibin++;
            }
        }
    }

    // if not defined in the parameters, we fix the range from March 1rst to December 15th
    DateMin = hDeces_Tot->GetXaxis()->FindFixBin("1-Mar-20");
    DateMax = hDeces_Tot->GetXaxis()->FindFixBin("15-Feb-21");


    firstHist = hDeces_Tot;
    firstHist_tot = hDeces_Tot;

    TString LastDate = vDates.back();

    cout<<vDeaths.size()<<endl;
    cout<<left<<setw(10)<<"Date"<<setw(10)<<"Daily"<<setw(10)<<"Total"<<endl;
    for(size_t i=0 ; i<vDates.size() ; i++) {
        if(i<vDeaths_Tot.size() && vDeaths_Tot.at(i)) {
            Int_t Bin = firstHist->GetXaxis()->FindFixBin(vDates.at(i));
            if(Bin>0) {
                hDeces_Tot->SetBinContent(Bin,vDeaths_Tot.at(i));
                //                hDeces_Tot->SetBinError(Bin,2*sqrt(vDeaths_Tot.at(i)));
                hDeces_Tot->SetBinError(Bin,vDeaths_e.at(i));
//                LastDate = vDates.at(i);
            }
        }
        if(i<vDeaths.size() && vDeaths.at(i)) {
            Int_t Bin = firstHist->GetXaxis()->FindFixBin(vDates.at(i));
            if(Bin>0) {
                //                hDeces->SetBinContent(Bin,vDeaths.at(i));
                //                hDeces->SetBinError(Bin,vDeaths_e.at(i));
            }

            cout<<left<<setw(10)<<vDates.at(i)<<setw(10)<<vDeaths.at(i)<<setw(10)<<vDeaths_Tot.at(i)<<endl;
        }
    }

    MyCanvas->cd();
    //    gPad->SetGrid(1,1);

    hDeces_Tot->Draw("p");

    Double_t XMin2 = hDeces_Tot->GetXaxis()->GetBinLowEdge(hDeces_Tot->FindFirstBinAbove(1));
    Double_t XMax2 = hDeces_Tot->GetXaxis()->GetBinUpEdge(hDeces_Tot->FindLastBinAbove(1));

    // Fit params

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2","Migrad");
    //    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1e9);
    ROOT::Math::MinimizerOptions::SetDefaultErrorDef(2);
    //    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1);
    //        ROOT::Math::MinimizerOptions::SetDefaultPrecision(1e-9);

    //D Model
    if(DoD) {
        fTotalD = new TF1(Form("DModel_%s",hDeces_Tot->GetName()),FuncD,hDeces_Tot->GetXaxis()->GetXmin(),hDeces_Tot->GetXaxis()->GetXmax(),4);
        fTotalD->SetLineColor(kRed);
        fTotalD->SetNpx(1000);

        const int NPars = fTotalD->GetNpar();
        Double_t *Pars = new Double_t[NPars];
        //    Double_t ParsErr[NPars];

        Pars[0] = 10;
        Pars[1] = 15.;
        Pars[2] = 1e-3;
        Pars[3] = hDeces_Tot->GetXaxis()->GetBinLowEdge(hDeces_Tot->FindFirstBinAbove(1));

        fTotalD->SetParameter(0,Pars[0]);
        fTotalD->SetParameter(1,Pars[1]);
        fTotalD->SetParameter(2,Pars[2]);
        fTotalD->FixParameter(3,Pars[3]);

        fTotalD->SetParLimits(0,1.,1000);
        fTotalD->SetParLimits(1,0.1,50.);
        fTotalD->SetParLimits(2,1e-6,1);

        TFitResultPtr r = hDeces_Tot->Fit(fTotalD,"S0","",XMin2,XMax2);
        fChi2D = r->Chi2()/r->Ndf();

        r->Print("V");
        fTotalD->Draw("same");

        /*Create a histogram to hold the confidence intervals*/
        auto *herror = (TH1*)hDeces_Tot->Clone();
        herror->Reset();
        herror->SetName(((TString)hDeces_Tot->GetName()).Append("_error"));
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(herror);

        //Now the "hint" histogram has the fitted function values as the
        //bin contents and the confidence intervals as bin errors
        herror->SetStats(kFALSE);
        herror->SetFillColor(fTotalD->GetLineColor());
        herror->SetFillStyle(3002);
        herror->SetFillColorAlpha(fTotalD->GetLineColor(),0.5);
        herror->SetMarkerSize(0);
        herror->Draw("e3 same");
    }

    // D2 Model
    if(DoD2) {

        fTotalD2 = new TF1(Form("D2_%s",hDeces_Tot->GetName()),FuncD2,hDeces_Tot->GetXaxis()->GetXmin(),hDeces_Tot->GetXaxis()->GetXmax(),5);

        fTotalD2->SetLineColor(kBlue);
        fTotalD2->SetNpx(1000);

        const int NPars2 = fTotalD2->GetNpar();
        Double_t *Pars2 = new Double_t[NPars2];

        Pars2[0] = 50;
        Pars2[1] = 4.;
        Pars2[2] = 1e-3;
        Pars2[3] = 7;
        Pars2[4] = hDeces_Tot->GetXaxis()->GetBinLowEdge(hDeces_Tot->FindFirstBinAbove(1));

        fTotalD2->SetParameter(0,Pars2[0]);
        fTotalD2->SetParameter(1,Pars2[1]);
        fTotalD2->SetParameter(2,Pars2[2]);
        fTotalD2->SetParameter(3,Pars2[3]);
        fTotalD2->FixParameter(4,Pars2[4]);

        fTotalD2->SetParLimits(0,1,1000);
        fTotalD2->SetParLimits(1,0.1,50.);
        fTotalD2->SetParLimits(2,1e-6,1);
        fTotalD2->SetParLimits(3,0.1,50);

        TFitResultPtr r = hDeces_Tot->Fit(fTotalD2,"S0","",XMin2,XMax2);
        fChi2D2 = r->Chi2()/r->Ndf();

        r->Print("V");
        fTotalD2->Draw("same");

        /*Create a histogram to hold the confidence intervals*/
        auto *herror = (TH1*)hDeces_Tot->Clone();
        herror->Reset();
        herror->SetName(((TString)hDeces_Tot->GetName()).Append("_error"));
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(herror);

        //Now the "hint" histogram has the fitted function values as the
        //bin contents and the confidence intervals as bin errors
        herror->SetStats(kFALSE);
        herror->SetFillColor(fTotalD2->GetLineColor());
        herror->SetFillStyle(3002);
        herror->SetFillColorAlpha(fTotalD2->GetLineColor(),0.5);
        herror->SetMarkerSize(0);
        herror->Draw("e3 same");

        if(DoPlotSub) {
            TF1 *f1 = new TF1(Form("DModel_%s_1",hDeces_Tot->GetName()),FuncD,hDeces_Tot->GetXaxis()->GetXmin(),hDeces_Tot->GetXaxis()->GetXmax(),4);
            f1->SetParameters(fTotalD2->GetParameter(0),fTotalD2->GetParameter(1),fTotalD2->GetParameter(2),fTotalD2->GetParameter(4));
            f1->SetLineColor(fTotalD2->GetLineColor());
            f1->SetLineStyle(kDashed);
            f1->Draw("same");
            TF1 *f2 = new TF1(Form("DModel_%s_2",hDeces_Tot->GetName()),FuncD,hDeces_Tot->GetXaxis()->GetXmin(),hDeces_Tot->GetXaxis()->GetXmax(),4);
            f2->SetParameters(fTotalD2->GetParameter(0),fTotalD2->GetParameter(3),fTotalD2->GetParameter(2),fTotalD2->GetParameter(4));
            f2->SetLineColor(fTotalD2->GetLineColor());
            f2->SetLineStyle(kDashed);
            f2->Draw("same");
        }
    }

    // D2 Model
    if(DoD2Full) {
        fTotalD2Full = new TF1(Form("D2_%s",hDeces_Tot->GetName()),FuncD2Full,hDeces_Tot->GetXaxis()->GetXmin(),hDeces_Tot->GetXaxis()->GetXmax(),7);

        fTotalD2Full->SetLineColor(kGreen);
        fTotalD2Full->SetNpx(1000);

        const int NPars2 = fTotalD2Full->GetNpar();
        Double_t *Pars2 = new Double_t[NPars2];

        Pars2[0] = 50;
        Pars2[1] = 4.;
        Pars2[2] = 1e-3;

        Pars2[3] = 50;
        Pars2[4] = 10.;
        Pars2[5] = 1e-3;

        Pars2[6] = hDeces_Tot->GetXaxis()->GetBinLowEdge(hDeces_Tot->FindFirstBinAbove(1));

        fTotalD2Full->SetParameter(0,Pars2[0]);
        fTotalD2Full->SetParameter(1,Pars2[1]);
        fTotalD2Full->SetParameter(2,Pars2[2]);
        fTotalD2Full->SetParameter(3,Pars2[3]);
        fTotalD2Full->SetParameter(4,Pars2[4]);
        fTotalD2Full->SetParameter(5,Pars2[5]);
        fTotalD2Full->FixParameter(6,Pars2[6]);

        fTotalD2Full->SetParLimits(0,1,1000);
        fTotalD2Full->SetParLimits(1,0.1,50.);
        fTotalD2Full->SetParLimits(2,1e-6,1);
        fTotalD2Full->SetParLimits(3,1.,1000);
        fTotalD2Full->SetParLimits(4,0.1,50.);
        fTotalD2Full->SetParLimits(5,1e-6,1);

        TFitResultPtr r = hDeces_Tot->Fit(fTotalD2Full,"S0","",XMin2,XMax2);
        fChi2D2Full = r->Chi2()/r->Ndf();

        r->Print("V");
        fTotalD2Full->Draw("same");

        /*Create a histogram to hold the confidence intervals*/
        auto *herror = (TH1*)hDeces_Tot->Clone();
        herror->Reset();
        herror->SetName(((TString)hDeces_Tot->GetName()).Append("_error"));
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(herror);

        //Now the "hint" histogram has the fitted function values as the
        //bin contents and the confidence intervals as bin errors
        herror->SetStats(kFALSE);
        herror->SetFillColor(fTotalD2Full->GetLineColor());
        herror->SetFillStyle(3002);
        herror->SetFillColorAlpha(fTotalD2Full->GetLineColor(),0.5);
        herror->SetMarkerSize(0);
        herror->Draw("e3 same");

        if(DoPlotSub) {
            TF1 *f1 = new TF1(Form("DModel_%s_1",hDeces_Tot->GetName()),FuncD,hDeces_Tot->GetXaxis()->GetXmin(),hDeces_Tot->GetXaxis()->GetXmax(),4);
            f1->SetParameters(fTotalD2Full->GetParameter(0),fTotalD2Full->GetParameter(1),fTotalD2Full->GetParameter(2),fTotalD2Full->GetParameter(6));
            f1->SetLineColor(fTotalD2Full->GetLineColor());
            f1->SetLineStyle(kDashed);
            f1->Draw("same");
            TF1 *f2 = new TF1(Form("DModel_%s_2",hDeces_Tot->GetName()),FuncD,hDeces_Tot->GetXaxis()->GetXmin(),hDeces_Tot->GetXaxis()->GetXmax(),4);
            f2->SetParameters(fTotalD2Full->GetParameter(3),fTotalD2Full->GetParameter(4),fTotalD2Full->GetParameter(5),fTotalD2Full->GetParameter(6));
            f2->SetLineColor(fTotalD2Full->GetLineColor());
            f2->SetLineStyle(kDashed);
            f2->Draw("same");
        }
    }

    MaxY = firstHist->GetMaximum() * 1.2;
    firstHist->GetYaxis()->SetRangeUser(0,MaxY);
    firstHist->GetXaxis()->SetRange(DateMin,DateMax);

    for(int i=1 ; i<=firstHist->GetNbinsX() ; i+=7) {
        firstHist->GetXaxis()->SetBinLabel(i+1,"");
        firstHist->GetXaxis()->SetBinLabel(i+2,"");
        firstHist->GetXaxis()->SetBinLabel(i+3,"");
        firstHist->GetXaxis()->SetBinLabel(i+4,"");
        firstHist->GetXaxis()->SetBinLabel(i+5,"");
        firstHist->GetXaxis()->SetBinLabel(i+6,"");
    }

    MyCanvas->cd();

    MyCanvas->Modified();
    MyCanvas->Update();

    Float_t XVal = gPad->GetFrame()->GetX1()*1.02;
    Float_t YVal = gPad->GetFrame()->GetY2()*0.96;

    TLatex *text = new TLatex(XVal,YVal,Form("%s: %s",theCountry.Data(),LastDate.Data()));
    text->SetTextColor(kBlack);text->Draw();
    text->SetTextSize(0.05);
    text->SetTextFont(132);
    text->Draw();

    YVal = gPad->GetFrame()->GetY2()*0.9;
    XVal = gPad->GetFrame()->GetX1() * 1.05;
    Int_t NDY=0;
    Float_t DY = gPad->GetFrame()->GetY2()*0.04;
    Float_t TextSize = 0.03;

    // Print D
    if(DoD) {
        TLatex *text = new TLatex(XVal,YVal-DY*NDY,"D model");
        text->SetTextColor(fTotalD->GetLineColor());text->Draw();
        text->SetTextFont(132);
        text->SetTextSize(TextSize+0.01);
        NDY++;

        text = new TLatex(XVal,YVal-DY*NDY,Form("a = %.2f #pm %.2f",fTotalD->GetParameter(0),fTotalD->GetParError(0)));
        text->SetTextColor(fTotalD->GetLineColor());text->Draw();
        text->SetTextFont(132);
        text->SetTextSize(TextSize);
        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("b = %.2f #pm %.2f",fTotalD->GetParameter(1),fTotalD->GetParError(1)));
        text->SetTextColor(fTotalD->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("c = %.2e #pm %.2e",fTotalD->GetParameter(2),fTotalD->GetParError(2)));
        text->SetTextColor(fTotalD->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;

        text = new TLatex(XVal,YVal-DY*NDY,Form("Chi2/ndf = %.2f",fChi2D));
        text->SetTextColor(fTotalD->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;
        NDY++;
    }
    // Print D2
    if(DoD2) {
        TLatex *text = new TLatex(XVal,YVal-DY*NDY,"D2 model");
        text->SetTextColor(fTotalD2->GetLineColor());text->Draw();
        text->SetTextFont(132);
        text->SetTextSize(TextSize+0.01);
        NDY++;

        text = new TLatex(XVal,YVal-DY*NDY,Form("a = %.2f #pm %.2f",fTotalD2->GetParameter(0),fTotalD2->GetParError(0)));
        text->SetTextColor(fTotalD2->GetLineColor());text->Draw();
        text->SetTextFont(132);
        text->SetTextSize(TextSize);
        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("b1 = %.2f #pm %.2f",fTotalD2->GetParameter(1),fTotalD2->GetParError(1)));
        text->SetTextColor(fTotalD2->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("b2 = %.2f #pm %.2f",fTotalD2->GetParameter(3),fTotalD2->GetParError(3)));
        text->SetTextColor(fTotalD2->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("c = %.2e #pm %.2e",fTotalD2->GetParameter(2),fTotalD2->GetParError(2)));
        text->SetTextColor(fTotalD2->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("Chi2/ndf = %.2f",fChi2D2));
        text->SetTextColor(fTotalD2->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;
        NDY++;
    }
    // Print D2
    if(DoD2Full) {
        TLatex *text = new TLatex(XVal,YVal-DY*NDY,"D2 full model");
        text->SetTextColor(fTotalD2Full->GetLineColor());text->Draw();
        text->SetTextFont(132);
        text->SetTextSize(TextSize+0.01);
        NDY++;

        text = new TLatex(XVal,YVal-DY*NDY,Form("a1 = %.2f #pm %.2f",fTotalD2Full->GetParameter(0),fTotalD2Full->GetParError(0)));
        text->SetTextColor(fTotalD2Full->GetLineColor());text->Draw();
        text->SetTextFont(132);
        text->SetTextSize(TextSize);
        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("b1 = %.2f #pm %.2f",fTotalD2Full->GetParameter(1),fTotalD2Full->GetParError(1)));
        text->SetTextColor(fTotalD2Full->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("c1 = %.2e #pm %.2e",fTotalD2Full->GetParameter(2),fTotalD2Full->GetParError(2)));
        text->SetTextColor(fTotalD2Full->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;

        text = new TLatex(XVal,YVal-DY*NDY,Form("a2 = %.2f #pm %.2f",fTotalD2Full->GetParameter(3),fTotalD2Full->GetParError(3)));
        text->SetTextColor(fTotalD2Full->GetLineColor());text->Draw();
        text->SetTextFont(132);
        text->SetTextSize(TextSize);
        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("b2 = %.2f #pm %.2f",fTotalD2Full->GetParameter(4),fTotalD2Full->GetParError(4)));
        text->SetTextColor(fTotalD2Full->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("c2 = %.2e #pm %.2e",fTotalD2Full->GetParameter(5),fTotalD2Full->GetParError(5)));
        text->SetTextColor(fTotalD2Full->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;

        text = new TLatex(XVal,YVal-DY*NDY,Form("Chi2/ndf = %.2f",fChi2D2Full));
        text->SetTextColor(fTotalD2Full->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;
        NDY++;
    }

    TString FileName = Form("COVID_Total_Death_Ana_%s",theCountry.Data());
    (Smooth>1) ? FileName.Append(Form("_%dDaysSmooth",Smooth)) : FileName ;
    FileName.Append(Form("_%s.png",vDates.at(vDeaths_Tot.size()-1).Data()));
    gPad->GetCanvas()->SaveAs(FileName);
}

Double_t FuncD(Double_t*xx,Double_t*pp) {

    Double_t a1  = pp[0];
    Double_t b1  = pp[1];
    Double_t c1  = pp[2];

    Double_t t0 = pp[3];
    Double_t x   = xx[0] - t0;

    Double_t D1 = a1*TMath::Exp(x/b1)/(1+c1*TMath::Exp(x/b1));

    Double_t Daily = D1;

    return Daily;
}

Double_t FuncD2(Double_t*xx,Double_t*pp) {

    Double_t a1  = pp[0];
    Double_t b1  = pp[1];
    Double_t c1  = pp[2];
    Double_t b2  = pp[3];

    Double_t t0 = pp[4];
    Double_t x   = xx[0] - t0;

    Double_t D1 = a1*TMath::Exp(x/b1)/(1+c1*TMath::Exp(x/b1));
    Double_t D2 = a1*TMath::Exp(x/b2)/(1+c1*TMath::Exp(x/b2));

    Double_t Daily = D1+D2;

    return Daily;
}

Double_t FuncD2Full(Double_t*xx,Double_t*pp) {

    Double_t a1  = pp[0];
    Double_t b1  = pp[1];
    Double_t c1  = pp[2];

    Double_t a2  = pp[3];
    Double_t b2  = pp[4];
    Double_t c2  = pp[5];

    Double_t t0 = pp[6];
    Double_t x   = xx[0] - t0;

    Double_t D1 = a1*TMath::Exp(x/b1)/(1+c1*TMath::Exp(x/b1));
    Double_t D2 = a2*TMath::Exp(x/b2)/(1+c2*TMath::Exp(x/b2));

    Double_t Daily = D1+D2;

    return Daily;
}
