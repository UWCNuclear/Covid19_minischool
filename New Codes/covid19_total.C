#include "covid19_total.h"

///****************************************************************************************************************
///                                             User Guide
///****************************************************************************************************************
/// Main function:
///           Analyse(TString CountryName);
///             => country name need to correspond to a csv file where the data are downloaded from worldometers
/// Avalailble Options:
///           SetModels(Bool_t DoD, Bool_t DoD2, Bool_t DoESIR, Bool_t DoESIR2, Bool_t FullModel);
///             => Define the models that will be fitted on the data, default is D'2 and ESIR2 in full mode
///
///           SetSmoothing(Int_t Ndays);
///             => Number of average days in the sliding window. Default: 3
///
///           ReadDataRange(TString DateFrom,TString DateTo);
///             => Define the range of dates to read, default is all.
///             => ReadDataRange("","1-Mar-2021") means all up to 1-Mar-2021
///             => ReadDataRange("1-Aug-2020","") means all from to 1-Aug-2020
///
///           SetAxisRange(TString DateFrom,TString DateTo);
///             => Define the range of the histogram axis, default is adapted to the data.
///
///           SetFitRange(TString DateFrom,TString DateTo);
///             => Define the range of the histogram axis, default is adapted to the axis range
///
///****************************************************************************************************************

// Main fonction that plots the data and process the fits
void
Analyse(TString theCountry) {

    // histogram initialization
    InitHistograms();

    // to print the program's configuration in the terminal
    PrintParameters(theCountry);

    // style parameters, to remove unsed default titles and stat
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    // Path where the files form worldometers have been downloaded
    TString Folder = "./worldometers/";

    TString FileName = Form("%s/%s.csv",Folder.Data(),theCountry.Data());

    // now, the data file is read using the ReadData function
    bool data_ok = ReadData(FileName);
    if(data_ok == false) return;

    // if no data has been read, we exit
    if(vTotal_Deaths.empty()) {
        cout<<"OUPS, empty data"<<endl;
        return;
    }

    // we remove the first possible data points that are bellow the defined threshold
    while(vTotal_Deaths.front()<DeathsMin) {
        vDates.erase(vDates.begin());
        vTotal_Deaths.erase(vTotal_Deaths.begin());
    }

    // The fonction SmoothVector is then used to smooth the data on Smooth successive days
    SmoothVector(fNSmoothing,vTotal_Deaths,vTotal_Deaths_error);

    // for better printouts in the plots, we change the coutries names of US and UK
    if(theCountry.EqualTo("US",TString::kIgnoreCase)) theCountry = "USA";
    if(theCountry.EqualTo("UK",TString::kIgnoreCase)) theCountry = "United Kingdom";
    theCountry.ReplaceAll("_"," ");

    hTotal_Deaths->SetNameTitle(Form("TotalD_%s",theCountry.Data()),Form("TotalD_%s",theCountry.Data()));

    // get the bins corresponding to the defined range
    Int_t DateMin,DateMax;
    if(fAxisRangeFrom=="") DateMin = max(1,hDummyHist->GetXaxis()->FindFixBin(vDates.front())-5);
    else DateMin = hDummyHist->GetXaxis()->FindFixBin(fAxisRangeFrom);
    if(fAxisRangeTo=="") DateMax = min(hDummyHist->GetNbinsX(),hDummyHist->GetXaxis()->FindFixBin(vDates.back())+5);
    else DateMax = hDummyHist->GetXaxis()->FindFixBin(fAxisRangeTo);

    // define the fit range
    Int_t XMin, XMax;
    if(fFitRangeFrom=="") XMin = DateMin;
    else XMin = hDummyHist->GetXaxis()->FindFixBin(fFitRangeFrom);
    if(fFitRangeTo=="") XMax = DateMax;
    else XMax = hDummyHist->GetXaxis()->FindFixBin(fFitRangeTo);

    // the LastDate string is used to plot the last date of the data
    TString LastDate = vDates.back();
    for(size_t i=0 ; i<vDates.size() ; i++) {
        if(i<vTotal_Deaths.size() && vTotal_Deaths.at(i)) {
            Int_t Bin = hTotal_Deaths->GetXaxis()->FindFixBin(vDates.at(i));
            if(Bin>0) {
                hTotal_Deaths->SetBinContent(Bin,vTotal_Deaths.at(i));
                hTotal_Deaths->SetBinError(Bin,vTotal_Deaths_error.at(i));
            }
            LastDate = vDates.at(i);
        }
    }

    // We create the Canvas and margins in which all will be ploted
    TCanvas *MyCanvas = new TCanvas("Total","Total",1600,1200);
    MyCanvas->SetLeftMargin(0.107635);
    MyCanvas->SetRightMargin(0.00125156);
    MyCanvas->SetBottomMargin(0.13619);
    MyCanvas->SetTopMargin(0.00190476);

    // The Total deaths histogram is ploted
    hTotal_Deaths->Draw("p");

    // Chi2 definition
    Double_t fChi2D, fChi2D2;
    // functions definition
    TF1 *fTotal_D,*fTotal_D2;

    // Minimizer definition
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2","Migrad");
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(kMaxInt);
    ROOT::Math::MinimizerOptions::SetDefaultErrorDef(2);
    //    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-3);
    //    ROOT::Math::MinimizerOptions::SetDefaultPrecision(1e-9);

    // In the following, we define the different models, as a function of what has been asked in the Main fonctiuon parameters
    // DModel
    if(fDoD) {
        // The function is defined using the FuncD function, defined at the end of the file
        fTotal_D = new TF1(Form("D_%s",hTotal_Deaths->GetName()),FuncD,hTotal_Deaths->GetXaxis()->GetXmin(),hTotal_Deaths->GetXaxis()->GetXmax(),5);
        fTotal_D->SetLineColor(fColors[0]);
        fTotal_D->SetNpx(1000);

        const int NPars = fTotal_D->GetNpar();
        Double_t *Pars = new Double_t[NPars];

        if(fUseOffset) {
            fTotal_D->SetParameter(0,hTotal_Deaths->GetBinContent(hDummyHist->GetBinLowEdge(XMin)));
            fTotal_D->SetParLimits(0,0,hTotal_Deaths->GetBinContent(hDummyHist->GetBinLowEdge(XMax)));
        }
        else
            fTotal_D->FixParameter(0,0);

        // parameters initialization
        Pars[1] = 50;
        Pars[2] = 4.;
        Pars[3] = 1e-3;
        Pars[4] = XMin;//hTotal_Deaths->GetXaxis()->GetBinLowEdge(hTotal_Deaths->FindFirstBinAbove(1));

        fTotal_D->SetParameter(1,Pars[1]);
        fTotal_D->SetParameter(2,Pars[2]);
        fTotal_D->SetParameter(3,Pars[3]);
        fTotal_D->FixParameter(4,Pars[4]);

        fTotal_D->SetParLimits(1,0,1000);
        fTotal_D->SetParLimits(2,1.,20.);
        fTotal_D->SetParLimits(3,1e-6,1);

        // Fit of the histogram
        TFitResultPtr r = hTotal_Deaths->Fit(fTotal_D,"S0","",XMin,XMax);
        fChi2D = r->Chi2()/r->Ndf();

        r->Print("V");
        fTotal_D->Draw("same");

        /*Create a histogram to hold the confidence intervals*/
        auto *herrorD2 = (TH1*)hTotal_Deaths->Clone();
        herrorD2->Reset();
        herrorD2->SetName(((TString)hTotal_Deaths->GetName()).Append("_errorD"));
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(herrorD2);

        //Now the "hint" histogram has the fitted function values as the
        //bin contents and the confidence intervals as bin errors
        herrorD2->SetStats(kFALSE);
        herrorD2->SetFillColor(fTotal_D->GetLineColor());
        herrorD2->SetFillStyle(3002);
        herrorD2->SetFillColorAlpha(fTotal_D->GetLineColor(),0.5);
        herrorD2->SetMarkerSize(0);
        herrorD2->Draw("e3 same");
    }

    // D2Model
    if(fDoD2){
        if(fDoFullModel) {
            fTotal_D2 = new TF1(Form("D2_%s",hTotal_Deaths->GetName()),FuncD2Full,hTotal_Deaths->GetXaxis()->GetXmin(),hTotal_Deaths->GetXaxis()->GetXmax(),8);
        }
        else {
            fTotal_D2 = new TF1(Form("D2_%s",hTotal_Deaths->GetName()),FuncD2,hTotal_Deaths->GetXaxis()->GetXmin(),hTotal_Deaths->GetXaxis()->GetXmax(),6);
        }
        fTotal_D2->SetLineColor(fColors[1]);
        fTotal_D2->SetNpx(1000);

        const int NPars = fTotal_D2->GetNpar();
        Double_t *Pars = new Double_t[NPars];

        if(fUseOffset) {
            fTotal_D2->SetParameter(0,hTotal_Deaths->GetBinContent(hDummyHist->GetBinLowEdge(XMin)));
            fTotal_D2->SetParLimits(0,0,hTotal_Deaths->GetBinContent(hDummyHist->GetBinLowEdge(XMax)));
        }
        else
            fTotal_D2->FixParameter(0,0);

        if(fDoFullModel) {
            Pars[1] = 50;
            Pars[2] = 4.;
            Pars[3] = 1e-3;

            Pars[4] = 50;
            Pars[5] = 10.;
            Pars[6] = 1e-3;

            Pars[7] = XMin;//hTotal_Deaths->GetXaxis()->GetBinLowEdge(hTotal_Deaths->FindFirstBinAbove(1));

            fTotal_D2->SetParameter(1,Pars[1]);
            fTotal_D2->SetParameter(2,Pars[2]);
            fTotal_D2->SetParameter(3,Pars[3]);
            fTotal_D2->SetParameter(4,Pars[4]);
            fTotal_D2->SetParameter(5,Pars[5]);
            fTotal_D2->SetParameter(6,Pars[6]);
            fTotal_D2->FixParameter(7,Pars[7]);

            fTotal_D2->SetParLimits(1,1,1000);
            fTotal_D2->SetParLimits(2,1.,50.);
            fTotal_D2->SetParLimits(3,1e-6,1);
            fTotal_D2->SetParLimits(4,0,1000);
            fTotal_D2->SetParLimits(5,1.,50.);
            fTotal_D2->SetParLimits(6,1e-6,1);
        }
        else {
            Pars[1] = 50;
            Pars[2] = 4.;
            Pars[3] = 1e-3;
            Pars[4] = 7;
            Pars[5] = XMin;//hTotal_Deaths->GetXaxis()->GetBinLowEdge(hTotal_Deaths->FindFirstBinAbove(1));

            fTotal_D2->SetParameter(1,Pars[1]);
            fTotal_D2->SetParameter(2,Pars[2]);
            fTotal_D2->SetParameter(3,Pars[3]);
            fTotal_D2->SetParameter(4,Pars[4]);
            fTotal_D2->FixParameter(5,Pars[5]);

            fTotal_D2->SetParLimits(1,0,1000);
            fTotal_D2->SetParLimits(2,1.,50.);
            fTotal_D2->SetParLimits(3,1e-6,1);
            fTotal_D2->SetParLimits(4,3,50);
        }

        TFitResultPtr r = hTotal_Deaths->Fit(fTotal_D2,"S0","",XMin,XMax);
        fChi2D2 = r->Chi2()/r->Ndf();

        r->Print("V");
        fTotal_D2->Draw("same");

        /*Create a histogram to hold the confidence intervals*/
        auto *herrorD2 = (TH1*)hTotal_Deaths->Clone();
        herrorD2->Reset();
        herrorD2->SetName(((TString)hTotal_Deaths->GetName()).Append("_errorD2"));
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(herrorD2);

        //Now the "hint" histogram has the fitted function values as the
        //bin contents and the confidence intervals as bin errors
        herrorD2->SetStats(kFALSE);
        herrorD2->SetFillColor(fTotal_D2->GetLineColor());
        herrorD2->SetFillStyle(3002);
        herrorD2->SetFillColorAlpha(fTotal_D2->GetLineColor(),0.5);
        herrorD2->SetMarkerSize(0);
        herrorD2->Draw("e3 same");

        if(fDoFullModel) {
            TF1 *f1 = new TF1(Form("D2_%s_1",hTotal_Deaths->GetName()),FuncD,hTotal_Deaths->GetXaxis()->GetXmin(),hTotal_Deaths->GetXaxis()->GetXmax(),5);
            f1->SetParameters(fTotal_D2->GetParameter(0),fTotal_D2->GetParameter(1),fTotal_D2->GetParameter(2),fTotal_D2->GetParameter(3),fTotal_D2->GetParameter(7));
            f1->SetLineColor(fTotal_D2->GetLineColor());
            f1->SetLineStyle(kDashed);
            f1->Draw("same");
            TF1 *f2 = new TF1(Form("D2_%s_2",hTotal_Deaths->GetName()),FuncD,hTotal_Deaths->GetXaxis()->GetXmin(),hTotal_Deaths->GetXaxis()->GetXmax(),5);
            f2->SetParameters(fTotal_D2->GetParameter(0),fTotal_D2->GetParameter(4),fTotal_D2->GetParameter(5),fTotal_D2->GetParameter(6),fTotal_D2->GetParameter(7));
            f2->SetLineColor(fTotal_D2->GetLineColor());
            f2->SetLineStyle(kDashed);
            f2->Draw("same");
        }
        else {
            TF1 *f1 = new TF1(Form("D2_%s_1",hTotal_Deaths->GetName()),FuncD,hTotal_Deaths->GetXaxis()->GetXmin(),hTotal_Deaths->GetXaxis()->GetXmax(),5);
            f1->SetParameters(fTotal_D2->GetParameter(0),fTotal_D2->GetParameter(1),fTotal_D2->GetParameter(2),fTotal_D2->GetParameter(3),fTotal_D2->GetParameter(5));
            f1->SetLineColor(fTotal_D2->GetLineColor());
            f1->SetLineStyle(kDashed);
            f1->Draw("same");
            TF1 *f2 = new TF1(Form("D2_%s_2",hTotal_Deaths->GetName()),FuncD,hTotal_Deaths->GetXaxis()->GetXmin(),hTotal_Deaths->GetXaxis()->GetXmax(),5);
            f2->SetParameters(fTotal_D2->GetParameter(0),fTotal_D2->GetParameter(1),fTotal_D2->GetParameter(4),fTotal_D2->GetParameter(3),fTotal_D2->GetParameter(5));
            f2->SetLineColor(fTotal_D2->GetLineColor());
            f2->SetLineStyle(kDashed);
            f2->Draw("same");
        }
    }

    Double_t MaxY = hTotal_Deaths->GetMaximum() * 1.2;
    hTotal_Deaths->GetYaxis()->SetRangeUser(0,MaxY);
    hTotal_Deaths->GetXaxis()->SetRange(DateMin,DateMax);

    // Here, we remove some of the labels, in order to not have more than 40 labels on the graph

    Int_t NBinsInRange = DateMax-DateMin;
    Int_t Step = NBinsInRange/40;
    for(int i=1 ; i<=hTotal_Deaths->GetNbinsX() ; i+=Step) {
        for(int ii=1 ; ii<Step ; ii++) {
            if((i+ii) <= hTotal_Deaths->GetNbinsX()) {
                hTotal_Deaths->GetXaxis()->SetBinLabel(i+ii,"");
            }
        }
    }

    // We update the plot and we print all the fit parameters and Chi2 values
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

    YVal = gPad->GetFrame()->GetY1() + (gPad->GetFrame()->GetY2()-gPad->GetFrame()->GetY1()) * 0.9;
    XVal = gPad->GetFrame()->GetX1() + (gPad->GetFrame()->GetX2()-gPad->GetFrame()->GetX1()) * 0.01;
    Int_t NDY=0;

    Int_t NFuncs = (fDoD+fDoD2);
    Float_t DY = gPad->GetFrame()->GetY2()*0.05;
    Float_t TextSize = 0.04;

    if(NFuncs>2) {
        TextSize = 0.02;
        DY = gPad->GetFrame()->GetY2()*0.03;
    }
    else {
        TextSize = 0.03;
        DY = gPad->GetFrame()->GetY2()*0.04;
    }

    // Print D
    if(fDoD) {
        TLatex *text = new TLatex(XVal,YVal-DY*NDY,"D model");
        text->SetTextColor(fTotal_D->GetLineColor());text->Draw();
        text->SetTextFont(132);
        text->SetTextSize(TextSize+0.01);
        NDY++;

        if(fUseOffset) {
            text = new TLatex(XVal,YVal-DY*NDY,Form("Off = %.2f #pm %.2f",fTotal_D->GetParameter(0),fTotal_D->GetParError(0)));
            text->SetTextColor(fTotal_D->GetLineColor());text->Draw();
            text->SetTextFont(132);
            text->SetTextSize(TextSize);
            NDY++;
        }
        text = new TLatex(XVal,YVal-DY*NDY,Form("a = %.2f #pm %.2f",fTotal_D->GetParameter(1),fTotal_D->GetParError(1)));
        text->SetTextColor(fTotal_D->GetLineColor());text->Draw();
        text->SetTextFont(132);
        text->SetTextSize(TextSize);
        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("b = %.2f #pm %.2f",fTotal_D->GetParameter(2),fTotal_D->GetParError(2)));
        text->SetTextColor(fTotal_D->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("c = %.2e #pm %.2e",fTotal_D->GetParameter(3),fTotal_D->GetParError(3)));
        text->SetTextColor(fTotal_D->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;

        text = new TLatex(XVal,YVal-DY*NDY,Form("Chi2/ndf = %.2f",fChi2D));
        text->SetTextColor(fTotal_D->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;
        NDY++;
    }

    // Print D2
    if(fDoD2) {
        if(fDoFullModel) {

            TLatex *text = new TLatex(XVal,YVal-DY*NDY,"D2 full model");
            text->SetTextColor(fTotal_D2->GetLineColor());text->Draw();
            text->SetTextFont(132);
            text->SetTextSize(TextSize+0.01);
            NDY++;

            if(fUseOffset) {
                text = new TLatex(XVal,YVal-DY*NDY,Form("Off = %.2f #pm %.2f",fTotal_D2->GetParameter(0),fTotal_D2->GetParError(0)));
                text->SetTextColor(fTotal_D2->GetLineColor());text->Draw();
                text->SetTextFont(132);
                text->SetTextSize(TextSize);
                NDY++;
            }
            text = new TLatex(XVal,YVal-DY*NDY,Form("a1 = %.2f #pm %.2f",fTotal_D2->GetParameter(1),fTotal_D2->GetParError(1)));
            text->SetTextColor(fTotal_D2->GetLineColor());text->Draw();
            text->SetTextFont(132);
            text->SetTextSize(TextSize);
            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("b1 = %.2f #pm %.2f",fTotal_D2->GetParameter(2),fTotal_D2->GetParError(2)));
            text->SetTextColor(fTotal_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("c1 = %.2e #pm %.2e",fTotal_D2->GetParameter(3),fTotal_D2->GetParError(3)));
            text->SetTextColor(fTotal_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;

            text = new TLatex(XVal,YVal-DY*NDY,Form("a2 = %.2f #pm %.2f",fTotal_D2->GetParameter(4),fTotal_D2->GetParError(4)));
            text->SetTextColor(fTotal_D2->GetLineColor());text->Draw();
            text->SetTextFont(132);
            text->SetTextSize(TextSize);
            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("b2 = %.2f #pm %.2f",fTotal_D2->GetParameter(5),fTotal_D2->GetParError(5)));
            text->SetTextColor(fTotal_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("c2 = %.2e #pm %.2e",fTotal_D2->GetParameter(6),fTotal_D2->GetParError(6)));
            text->SetTextColor(fTotal_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;

            text = new TLatex(XVal,YVal-DY*NDY,Form("Chi2/ndf = %.2f",fChi2D2));
            text->SetTextColor(fTotal_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;
            NDY++;
        }
        else {
            TLatex *text = new TLatex(XVal,YVal-DY*NDY,"D2 model");
            text->SetTextColor(fTotal_D2->GetLineColor());text->Draw();
            text->SetTextFont(132);
            text->SetTextSize(TextSize+0.01);
            NDY++;

            if(fUseOffset) {
            text = new TLatex(XVal,YVal-DY*NDY,Form("Off = %.2f #pm %.2f",fTotal_D2->GetParameter(0),fTotal_D2->GetParError(0)));
            text->SetTextColor(fTotal_D2->GetLineColor());text->Draw();
            text->SetTextFont(132);
            text->SetTextSize(TextSize);
            NDY++;
            }
            text = new TLatex(XVal,YVal-DY*NDY,Form("a = %.2f #pm %.2f",fTotal_D2->GetParameter(1),fTotal_D2->GetParError(1)));
            text->SetTextColor(fTotal_D2->GetLineColor());text->Draw();
            text->SetTextFont(132);
            text->SetTextSize(TextSize);
            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("b1 = %.2f #pm %.2f",fTotal_D2->GetParameter(2),fTotal_D2->GetParError(3)));
            text->SetTextColor(fTotal_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("b2 = %.2f #pm %.2f",fTotal_D2->GetParameter(4),fTotal_D2->GetParError(4)));
            text->SetTextColor(fTotal_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("c = %.2e #pm %.2e",fTotal_D2->GetParameter(3),fTotal_D2->GetParError(3)));
            text->SetTextColor(fTotal_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("Chi2/ndf = %.2f",fChi2D2));
            text->SetTextColor(fTotal_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;
            NDY++;
        }
    }

    gSystem->mkdir("Pictures");
    TString OutputFileName = Form("Pictures/covid19_Total_deaths_%s",theCountry.Data());
    if(fNSmoothing>1) OutputFileName.Append(Form("_%dDaysSmooth",fNSmoothing));
    OutputFileName.Append(Form("_%s.png",vDates.at(vTotal_Deaths.size()-1).Data()));
    gPad->GetCanvas()->SaveAs(OutputFileName);
}

void InitHistograms() {

    // definitions of the used histograms, for the years 2020 and 2021
    vector<TString> vdates;

    // definition of the string dates, to be used for gaphical plot of the dates
    TString Mounth_str[12] = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
    Int_t NDaysPerMounth[12] = {31,29,31,30,31,30,31,31,30,31,30,31};

    // Here, we define the names of the bins of the histogram to use the dates
    for(int year=20 ; year<=21 ; year++) {
        for(int i=0 ; i<12 ; i++) {
            if(year==21 && i==1) NDaysPerMounth[i] = 28;
            for(int j=0 ; j<NDaysPerMounth[i] ; j++) {
                TString Label = Form("%d-%s-%d",j+1,Mounth_str[i].Data(),year);
                vdates.push_back(Label);
            }
        }
    }

    delete hTotal_Deaths;
    delete hDummyHist;

    // now we define the Total deaths histogram and we define the axis properties
    hTotal_Deaths = new TH1D("hTotal_Deaths","hTotal_Deaths",vdates.size(),0,vdates.size());
    hDummyHist = new TH1D("dummy_hist","dummy_hist",vdates.size(),0,vdates.size());

    for(size_t ibin=0 ; ibin<vdates.size() ; ibin++) {
        hTotal_Deaths->GetXaxis()->SetBinLabel(ibin+1,vdates.at(ibin));
        hDummyHist->GetXaxis()->SetBinLabel(ibin+1,vdates.at(ibin));
    }

    hTotal_Deaths->GetYaxis()->SetTitle("DEATHS / DAY");
    hTotal_Deaths->GetYaxis()->CenterTitle();
    hTotal_Deaths->GetXaxis()->SetLabelSize(0.04);
    hTotal_Deaths->GetXaxis()->SetTitleOffset(1.);
    hTotal_Deaths->GetXaxis()->SetTitleFont(132);
    hTotal_Deaths->GetXaxis()->SetLabelFont(132);

    hTotal_Deaths->GetYaxis()->SetLabelSize(0.05);
    hTotal_Deaths->GetYaxis()->SetTitleSize(0.05);
    hTotal_Deaths->GetYaxis()->SetTitleOffset(1.15);
    hTotal_Deaths->GetYaxis()->SetTickSize(0.01);
    hTotal_Deaths->GetXaxis()->SetTickSize(0.01);
    hTotal_Deaths->GetYaxis()->SetTitleFont(132);
    hTotal_Deaths->GetYaxis()->SetLabelFont(132);

    hTotal_Deaths->SetDirectory(nullptr);
    hTotal_Deaths->SetMarkerStyle(20);
    hTotal_Deaths->SetMarkerColor(kBlack);
    hTotal_Deaths->SetLineColor(kBlack);

    // These points is used force the function to start and ends at 0.
    hTotal_Deaths->SetBinContent(1,0.001);
    hTotal_Deaths->SetBinContent(hTotal_Deaths->GetNbinsX(),0.001);


}

void SetSmoothing(Int_t Ndays) {
    fNSmoothing = Ndays;

    INFO_MESS << "Smoothing set to " << fNSmoothing << "days" << ENDL;
}

void ReadDataRange(TString DateFrom,TString DateTo) {

    if(hDummyHist == nullptr) InitHistograms();
    if(DateFrom!="" && hDummyHist->GetXaxis()->FindFixBin(DateFrom)==-1) {
        WARN_MESS << DateFrom << "not found in the histogram date range, ignored" << ENDL;
        DateFrom = "";
    }
    if(DateTo!="" && hDummyHist->GetXaxis()->FindFixBin(DateTo)==-1) {
        WARN_MESS << DateTo << "not found in the histogram date range, ignored" << ENDL;
        DateTo = "";
    }

    fReadDataFrom = DateFrom;
    fReadDataTo = DateTo;

    if(fReadDataFrom=="" && fReadDataTo=="") INFO_MESS << "Read all the available data" << ENDL;
    else if(fReadDataFrom=="") INFO_MESS << "Read all data up to " << fReadDataTo << ENDL;
    else if(fReadDataTo=="") INFO_MESS << "Read all data from " << fReadDataFrom << ENDL;
    else INFO_MESS << "Read data from " << fReadDataFrom << " to " << fReadDataTo << ENDL;
}

void SetAxisRange(TString DateFrom,TString DateTo) {

    if(hDummyHist == nullptr) InitHistograms();
    if(DateFrom!="" && hDummyHist->GetXaxis()->FindFixBin(DateFrom)==-1) {
        WARN_MESS << DateFrom << "not found in the histogram date range, ignored" << ENDL;
        DateFrom = "";
    }
    if(DateTo!="" && hDummyHist->GetXaxis()->FindFixBin(DateTo)==-1) {
        WARN_MESS << DateTo << "not found in the histogram date range, ignored" << ENDL;
        DateTo = "";
    }

    fAxisRangeFrom = DateFrom;
    fAxisRangeTo = DateTo;

    if(fAxisRangeFrom=="" && fAxisRangeTo=="") INFO_MESS << "Axis range adapted to the data" << ENDL;
    else if(fAxisRangeFrom=="") INFO_MESS << "Axis range up to " << fAxisRangeTo << ENDL;
    else if(fAxisRangeTo=="") INFO_MESS << "Axis range up from " << fAxisRangeFrom << ENDL;
    else INFO_MESS << "Axis range from " << fAxisRangeFrom << " to " << fAxisRangeTo << ENDL;
}

void SetFitRange(TString DateFrom,TString DateTo) {

    if(hDummyHist == nullptr) InitHistograms();
    if(DateFrom!="" && hDummyHist->GetXaxis()->FindFixBin(DateFrom)==-1) {
        WARN_MESS << DateFrom << "not found in the histogram date range, ignored" << ENDL;
        DateFrom = "";
    }
    if(DateTo!="" && hDummyHist->GetXaxis()->FindFixBin(DateTo)==-1) {
        WARN_MESS << DateTo << "not found in the histogram date range, ignored" << ENDL;
        DateTo = "";
    }

    fFitRangeFrom = DateFrom;
    fFitRangeTo = DateTo;

    if(fFitRangeFrom=="" && fFitRangeTo=="") INFO_MESS << "Fit range adapted to axis range" << ENDL;
    else if(fFitRangeFrom=="") INFO_MESS << "Fit range up to " << fFitRangeTo << ENDL;
    else if(fFitRangeTo=="") INFO_MESS << "Fit range up from " << fFitRangeFrom << ENDL;
    else INFO_MESS << "Fit range from " << fFitRangeFrom << " to " << fFitRangeTo << ENDL;
}

void SetModels(Bool_t DoD, Bool_t DoD2, Bool_t FullModel, Bool_t UseOffset) {
    fDoFullModel = FullModel;
    fDoD = DoD;
    fDoD2 = DoD2;
    fUseOffset = UseOffset;

    INFO_MESS << "Models parameters: ";
    if(fDoD) cout << " D': On ";
    if(fDoD2) cout << " D2': On ";
    if(fDoFullModel) cout << " ==> Full parameters mode activated";
    if(fUseOffset) cout << " ==> With offset";
    cout << ENDL;
}

void PrintParameters(TString country_name) {

    INFO_MESS << "Analyse data from: " << country_name << ENDL << ENDL;

    INFO_MESS << "Models parameters: ";
    if(fDoD) cout << " D': On ";
    if(fDoD2) cout << " D2': On ";
    if(fDoFullModel) cout << " ==> Full parameters mode activated";
    cout << ENDL;

    INFO_MESS << "Smoothing set to " << fNSmoothing << "days" << ENDL;

    if(fReadDataFrom=="" && fReadDataTo=="") INFO_MESS << "Read all the available data" << ENDL;
    else if(fReadDataFrom=="") INFO_MESS << "Read all data up to " << fReadDataTo << ENDL;
    else if(fReadDataTo=="") INFO_MESS << "Read all data from " << fReadDataFrom << ENDL;
    else INFO_MESS << "Read data from " << fReadDataFrom << " to " << fReadDataTo << ENDL;

    if(fAxisRangeFrom=="" && fAxisRangeTo=="") INFO_MESS << "Axis range adapted to the data" << ENDL;
    else if(fAxisRangeFrom=="") INFO_MESS << "Axis range up to " << fAxisRangeTo << ENDL;
    else if(fAxisRangeTo=="") INFO_MESS << "Axis range up from " << fAxisRangeFrom << ENDL;
    else INFO_MESS << "Axis range from " << fAxisRangeFrom << " to " << fAxisRangeTo << ENDL;

    if(fFitRangeFrom=="" && fFitRangeTo=="") INFO_MESS << "Fit range adapted to axis range" << ENDL;
    else if(fFitRangeFrom=="") INFO_MESS << "Fit range up to " << fFitRangeTo << ENDL;
    else if(fFitRangeTo=="") INFO_MESS << "Fit range up from " << fFitRangeFrom << ENDL;
    else INFO_MESS << "Fit range from " << fFitRangeFrom << " to " << fFitRangeTo << ENDL;

    INFO_MESS << "Press a key to continue"<< ENDL;
    cin.get();
}

bool ReadData(TString filename)
{
    // The selected file is opended, and if not found, return with an error message
    ifstream file(filename);
    if(!file) {
        cout<<filename<<" not found"<<endl;
        return false;
    }

    // The vectors containing data are cleared from previous use
    vDates.clear();
    vTotal_Deaths.clear();
    vTotal_Deaths.clear();
    vTotal_Deaths_error.clear();

    TString Buffer;
    string line;

    // the first line is not used, we read it first to skip this line
    getline(file,line);

    Int_t Current_Year = 20;

    bool init = true;
    if(fReadDataFrom != "") init = false;

    // then we look on all the lines of the file, puting each line in the string: line
    while(file) {
        getline(file,line);
        // The string line, is then copied in a ROOT string (TString), on which specific methods can be used to easily play with the string
        Buffer = line;

        // The array arr will contain all the elements all the line separated by the separator
        // as a function of the sofware used to write the data, the separator can be either ; or ,
        // To not ignore posible empty lines, we separate ;; (or ,,) with a space to obtain an empty value in the array
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

        // To extact the day and mounth, we again cut the date in a new array, "temp"
        TObjArray *temp = Date.Tokenize("-");
        // the mounth number is stored
        TString Mounth = (TString)temp->At(0)->GetName();
        // and then the Day number
        Int_t Day = ((TString)temp->At(1)->GetName()).Atoi();
        // The date, in string format, is then stored in our prefered format
        Date = Form("%d-%s-%d",Day,Mounth.Data(),Current_Year);
        // The array temp is no more necessary, we delete it to free
        delete temp;

        if(init==false && Date == fReadDataFrom) init = true;
        if(init == false) continue;

        if(Date.BeginsWith("31-Dec")) Current_Year++;

        // Then, the total number of deaths is stored
        Int_t Deaths = ((TString)arr->At(3)->GetName()).Atoi();

        // The array arr is no more necessary, we delete it to free
        delete arr;

        // if the death number is well defined, we push this info (date + deaths) in the associated vectors
        if(Deaths) {
            vDates.push_back(Date);
            vTotal_Deaths.push_back(Deaths);
        }
        // if the date is the last that has been asked to be taken into acount, we stop reading the file
        if(fReadDataTo !="" && Date == fReadDataTo) return true;
    }

    return true;
}

void SmoothVector(Int_t N, vector<double> &data, vector<double> &data_err)
{
    vector<Double_t> vSmooth;

    for(size_t i=0 ; i<data.size() ; i++) {
        Double_t NPoints = 0;
        Double_t Tot = 0.;
        Double_t Err2 = 0.;

        if(i>=N) {
            for(int ii=0 ; ii<N ; ii++) {
                if((data.at(i-ii)>0)) {
                    Tot += data.at(i-ii);
                    Err2 += 2*TMath::Power(sqrt(data.at(i-ii)),2);
                    NPoints ++;
                }
            }
        }

        Tot = Tot/NPoints;

        if(Tot>0.) {
            vSmooth.push_back(Tot);
            data_err.push_back(sqrt(Err2)/NPoints);
        }
        else {
            vSmooth.push_back(0);
            data_err.push_back(0);
        }
    }
    for(size_t i=0 ; i<data.size() ; i++) data.at(i) = vSmooth.at(i);
}


Double_t FuncD(Double_t*xx,Double_t*pp) {

    Double_t offset  = pp[0];

    Double_t a1  = pp[1];
    Double_t b1  = pp[2];
    Double_t c1  = pp[3];

    Double_t t0 = pp[4];
    Double_t x   = xx[0] - t0;

    Double_t D1 = a1*TMath::Exp(x/b1)/(1+c1*TMath::Exp(x/b1));

    Double_t Total = offset+D1;

    return Total;
}

Double_t FuncD2(Double_t*xx,Double_t*pp) {

    Double_t offset  = pp[0];

    Double_t a1  = pp[1];
    Double_t b1  = pp[2];
    Double_t c1  = pp[3];
    Double_t b2  = pp[4];

    Double_t t0 = pp[5];
    Double_t x   = xx[0] - t0;

    Double_t D1 = a1*TMath::Exp(x/b1)/(1+c1*TMath::Exp(x/b1));
    Double_t D2 = a1*TMath::Exp(x/b2)/(1+c1*TMath::Exp(x/b2));

    Double_t Total = offset+D1+D2;

    return Total;
}

Double_t FuncD2Full(Double_t*xx,Double_t*pp) {

    Double_t offset  = pp[0];

    Double_t a1  = pp[1];
    Double_t b1  = pp[2];
    Double_t c1  = pp[3];

    Double_t a2  = pp[4];
    Double_t b2  = pp[5];
    Double_t c2  = pp[6];

    Double_t t0 = pp[7];
    Double_t x   = xx[0] - t0;

    Double_t D1 = a1*TMath::Exp(x/b1)/(1+c1*TMath::Exp(x/b1));
    Double_t D2 = a2*TMath::Exp(x/b2)/(1+c2*TMath::Exp(x/b2));

    Double_t Total = offset+D1+D2;

    return Total;
}
