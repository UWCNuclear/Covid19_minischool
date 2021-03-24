#include "COVID_Deaths_Ana.h"

// Main fonction that plots the data and process the fits
void AnaPerCountry(TString theCountry, Int_t Smooth, TString ReadUpTo, Bool_t FullModel, Bool_t DoD, Bool_t DoD2, Bool_t DoESIR, Bool_t DoESIR2, TString theDateMax, TString theDateMin,TString FitDateMin, TString FitDateMax) {

    // style parameters, to remove unsed default titles and stat
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    // Path where the files form worldometers have been downloaded
    TString Folder = "./worldometers/";
//    TString Folder = "/Users/dudouet/Documents/Perso/Divers/Worldometers";

    TString FileName = Form("%s/%s.csv",Folder.Data(),theCountry.Data());

    // now, the data file is read using the ReadData function
    bool data_ok = ReadData(FileName, ReadUpTo);
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

    // now, we calculate the daily data as the difference between two successive days
    vDaily_Deaths.push_back(vTotal_Deaths.front());
    for(size_t i=1 ; i<vDates.size() ; i++) {
        if(vTotal_Deaths.size()>i && vTotal_Deaths.at(i)>0) vDaily_Deaths.push_back(vTotal_Deaths.at(i)-vTotal_Deaths.at(i-1));
    }

    // The fonction SmoothVector is then used to smooth the data on Smooth successive days
    SmoothVector(Smooth,vDaily_Deaths,vDaily_Deaths_error);

    // definition of the string dates, to be used for gaphical plot of the dates
    TString Mounth_str[12] = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
    Int_t NDaysPerMounth[12] = {31,29,31,30,31,30,31,31,30,31,30,31};
    Int_t NDaysInYear=0;
    for(int i=0 ; i<12 ; i++) NDaysInYear += NDaysPerMounth[i];

    // for better printouts in the plots, we change the coutries names of US and UK
    if(theCountry.EqualTo("US",TString::kIgnoreCase)) theCountry = "USA";
    if(theCountry.EqualTo("UK",TString::kIgnoreCase)) theCountry = "United Kingdom";
    theCountry.ReplaceAll("_"," ");

    // now we define the daily deaths histogram and we define the axis properties
    TH1D *hDaily_Deaths = new TH1D(Form("Deces_%s",theCountry.Data()),Form("Deces_%s",theCountry.Data()),NDaysInYear*2,0,NDaysInYear*2);
    hDaily_Deaths->GetYaxis()->SetTitle("DEATHS / DAY");
    hDaily_Deaths->GetYaxis()->CenterTitle();
    hDaily_Deaths->GetXaxis()->SetLabelSize(0.04);
    hDaily_Deaths->GetXaxis()->SetTitleOffset(1.);
    hDaily_Deaths->GetXaxis()->SetTitleFont(132);
    hDaily_Deaths->GetXaxis()->SetLabelFont(132);

    hDaily_Deaths->GetYaxis()->SetLabelSize(0.05);
    hDaily_Deaths->GetYaxis()->SetTitleSize(0.05);
    hDaily_Deaths->GetYaxis()->SetTitleOffset(1.15);
    hDaily_Deaths->GetYaxis()->SetTickSize(0.01);
    hDaily_Deaths->GetXaxis()->SetTickSize(0.01);
    hDaily_Deaths->GetYaxis()->SetTitleFont(132);
    hDaily_Deaths->GetYaxis()->SetLabelFont(132);

    hDaily_Deaths->SetDirectory(nullptr);
    hDaily_Deaths->SetMarkerStyle(20);
    hDaily_Deaths->SetMarkerColor(kBlack);
    hDaily_Deaths->SetLineColor(kBlack);

    // These points is used force the function to start and ends at 0.
    hDaily_Deaths->SetBinContent(1,0.001);
    hDaily_Deaths->SetBinContent(NDaysInYear*2,0.001);

    // Here, we define the names of the bins of the histogram to use the dates
    Int_t ibin=1;
    for(int year=20 ; year<=21 ; year++) {
        for(int i=0 ; i<12 ; i++) {
            if(year==21 && i==1) NDaysPerMounth[i] = 28;
            for(int j=0 ; j<NDaysPerMounth[i] ; j++) {
                TString Label = Form("%d-%s-%d",j+1,Mounth_str[i].Data(),year);
                hDaily_Deaths->GetXaxis()->SetBinLabel(ibin,Label);
                ibin++;
            }
        }
    }

    // if not defined in the parameters, we fix the range from March 1rst to December 15th
    DateMin = hDaily_Deaths->GetXaxis()->FindFixBin("1-Mar-20");
    DateMax = hDaily_Deaths->GetXaxis()->FindFixBin("15-Feb-21");

    // if defined, we use the dates passed as parameters for the range
    if(theDateMin != "") DateMin = hDaily_Deaths->GetXaxis()->FindFixBin(theDateMin);
    if(theDateMax != "") DateMax = hDaily_Deaths->GetXaxis()->FindFixBin(theDateMax);

    // the LastDate string is used to plot the last date of the data
    TString LastDate = vDates.back();

    for(size_t i=0 ; i<vDates.size() ; i++) {
        if(i<vDaily_Deaths.size() && vDaily_Deaths.at(i)) {
            Int_t Bin = hDaily_Deaths->GetXaxis()->FindFixBin(vDates.at(i));
            if(Bin>0) {
                hDaily_Deaths->SetBinContent(Bin,vDaily_Deaths.at(i));
                hDaily_Deaths->SetBinError(Bin,vDaily_Deaths_error.at(i));
            }
            LastDate = vDates.at(i);
        }
    }

    // We create the Canvas and margins in which all will be ploted
    TCanvas *MyCanvas = new TCanvas("daily","daily",1600,1200);
    MyCanvas->SetLeftMargin(0.107635);
    MyCanvas->SetRightMargin(0.00125156);
    MyCanvas->SetBottomMargin(0.13619);
    MyCanvas->SetTopMargin(0.00190476);

    // The daily deaths histogram is ploted
    hDaily_Deaths->Draw("p");

    // The first and last data points are stored for defining then the fit range
    Double_t XMin = hDaily_Deaths->GetXaxis()->GetBinLowEdge(hDaily_Deaths->FindFirstBinAbove(0));
    Double_t XMax = hDaily_Deaths->GetXaxis()->GetBinUpEdge(hDaily_Deaths->FindLastBinAbove(0));

    if(FitDateMin!="") XMin = hDaily_Deaths->GetXaxis()->FindFixBin(FitDateMin);
    if(FitDateMax!="") XMax = hDaily_Deaths->GetXaxis()->FindFixBin(FitDateMax);

    cout<<XMin<<" "<<XMax<<endl;

    // Minimizer definition
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2","Migrad");
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(2147483647);
    ROOT::Math::MinimizerOptions::SetDefaultErrorDef(2);
    //    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-3);
    //    ROOT::Math::MinimizerOptions::SetDefaultPrecision(1e-9);

    // In the following, we define the different models, as a function of what has been asked in the Main fonctiuon parameters
    // DModel
    if(DoD) {
        // The function is defined using the FuncD function, defined at the end of the file
        fDaily_D = new TF1(Form("D'_%s",hDaily_Deaths->GetName()),FuncD,XMin,hDaily_Deaths->GetXaxis()->GetXmax(),4);
        fDaily_D->SetLineColor(Colors[0]);
        fDaily_D->SetNpx(1000);

        const int NPars2 = fDaily_D->GetNpar();
        Double_t *Pars2 = new Double_t[NPars2];

        // parameters initialization
        Pars2[0] = 50;
        Pars2[1] = 4.;
        Pars2[2] = 1e-3;
        Pars2[3] = XMin;//hDaily_Deaths->GetXaxis()->GetBinLowEdge(hDaily_Deaths->FindFirstBinAbove(1));

        fDaily_D->SetParameter(0,Pars2[0]);
        fDaily_D->SetParameter(1,Pars2[1]);
        fDaily_D->SetParameter(2,Pars2[2]);
        fDaily_D->FixParameter(3,Pars2[3]);

        fDaily_D->SetParLimits(0,0,1000);
        fDaily_D->SetParLimits(1,1.,20.);
        fDaily_D->SetParLimits(2,1e-6,1);

        // Fit of the histogram
        TFitResultPtr r = hDaily_Deaths->Fit(fDaily_D,"S0","",XMin,XMax);
        fChi2D = r->Chi2()/r->Ndf();

        r->Print("V");
        fDaily_D->Draw("same");

        /*Create a histogram to hold the confidence intervals*/
        auto *herrorD2 = (TH1*)hDaily_Deaths->Clone();
        herrorD2->Reset();
        herrorD2->SetName(((TString)hDaily_Deaths->GetName()).Append("_errorD"));
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(herrorD2);

        //Now the "hint" histogram has the fitted function values as the
        //bin contents and the confidence intervals as bin errors
        herrorD2->SetStats(kFALSE);
        herrorD2->SetFillColor(fDaily_D->GetLineColor());
        herrorD2->SetFillStyle(3002);
        herrorD2->SetFillColorAlpha(fDaily_D->GetLineColor(),0.5);
        herrorD2->SetMarkerSize(0);
        herrorD2->Draw("e3 same");
    }

    // D2Model
    if(DoD2){
        if(FullModel) {
            fDaily_D2 = new TF1(Form("D'2_%s",hDaily_Deaths->GetName()),FuncD2Full,XMin,hDaily_Deaths->GetXaxis()->GetXmax(),7);
        }
        else {
            fDaily_D2 = new TF1(Form("D'2_%s",hDaily_Deaths->GetName()),FuncD2,XMin,hDaily_Deaths->GetXaxis()->GetXmax(),5);
        }
        fDaily_D2->SetLineColor(Colors[1]);
        fDaily_D2->SetNpx(1000);

        const int NPars2 = fDaily_D2->GetNpar();
        Double_t *Pars2 = new Double_t[NPars2];

        if(FullModel) {
            Pars2[0] = 50;
            Pars2[1] = 4.;
            Pars2[2] = 1e-3;

            Pars2[3] = 50;
            Pars2[4] = 10.;
            Pars2[5] = 1e-3;

            Pars2[6] = XMin;//hDaily_Deaths->GetXaxis()->GetBinLowEdge(hDaily_Deaths->FindFirstBinAbove(1));

            fDaily_D2->SetParameter(0,Pars2[0]);
            fDaily_D2->SetParameter(1,Pars2[1]);
            fDaily_D2->SetParameter(2,Pars2[2]);
            fDaily_D2->SetParameter(3,Pars2[3]);
            fDaily_D2->SetParameter(4,Pars2[4]);
            fDaily_D2->SetParameter(5,Pars2[5]);
            fDaily_D2->FixParameter(6,Pars2[6]);

            fDaily_D2->SetParLimits(0,1,1000);
            fDaily_D2->SetParLimits(1,1.,50.);
            //            fDaily_D2->SetParLimits(2,1e-6,1);
            //            fDaily_D2->SetParLimits(3,0,1000);
            fDaily_D2->SetParLimits(4,1.,50.);
            //            fDaily_D2->SetParLimits(5,1e-6,1);
        }
        else {
            Pars2[0] = 50;
            Pars2[1] = 4.;
            Pars2[2] = 1e-3;
            Pars2[3] = 7;
            Pars2[4] = XMin;//hDaily_Deaths->GetXaxis()->GetBinLowEdge(hDaily_Deaths->FindFirstBinAbove(1));

            fDaily_D2->SetParameter(0,Pars2[0]);
            fDaily_D2->SetParameter(1,Pars2[1]);
            fDaily_D2->SetParameter(2,Pars2[2]);
            fDaily_D2->SetParameter(3,Pars2[3]);
            fDaily_D2->FixParameter(4,Pars2[4]);

            fDaily_D2->SetParLimits(0,0,1000);
            fDaily_D2->SetParLimits(1,1.,50.);
            fDaily_D2->SetParLimits(2,1e-6,1);
            fDaily_D2->SetParLimits(3,3,50);
        }

        TFitResultPtr r = hDaily_Deaths->Fit(fDaily_D2,"S0","",XMin,XMax);
        fChi2D2 = r->Chi2()/r->Ndf();

        r->Print("V");
        fDaily_D2->Draw("same");

        /*Create a histogram to hold the confidence intervals*/
        auto *herrorD2 = (TH1*)hDaily_Deaths->Clone();
        herrorD2->Reset();
        herrorD2->SetName(((TString)hDaily_Deaths->GetName()).Append("_errorD2"));
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(herrorD2);

        //Now the "hint" histogram has the fitted function values as the
        //bin contents and the confidence intervals as bin errors
        herrorD2->SetStats(kFALSE);
        herrorD2->SetFillColor(fDaily_D2->GetLineColor());
        herrorD2->SetFillStyle(3002);
        herrorD2->SetFillColorAlpha(fDaily_D2->GetLineColor(),0.5);
        herrorD2->SetMarkerSize(0);
        herrorD2->Draw("e3 same");

        if(FullModel) {
            TF1 *f1 = new TF1(Form("D2_%s_1",hDaily_Deaths->GetName()),FuncD,hDaily_Deaths->GetXaxis()->GetXmin(),hDaily_Deaths->GetXaxis()->GetXmax(),4);
            f1->SetParameters(fDaily_D2->GetParameter(0),fDaily_D2->GetParameter(1),fDaily_D2->GetParameter(2),fDaily_D2->GetParameter(6));
            f1->SetLineColor(fDaily_D2->GetLineColor());
            f1->SetLineStyle(kDashed);
            f1->Draw("same");
            TF1 *f2 = new TF1(Form("DModel_%s_2",hDaily_Deaths->GetName()),FuncD,hDaily_Deaths->GetXaxis()->GetXmin(),hDaily_Deaths->GetXaxis()->GetXmax(),4);
            f2->SetParameters(fDaily_D2->GetParameter(3),fDaily_D2->GetParameter(4),fDaily_D2->GetParameter(5),fDaily_D2->GetParameter(6));
            f2->SetLineColor(fDaily_D2->GetLineColor());
            f2->SetLineStyle(kDashed);
            f2->Draw("same");
        }
        else {
            TF1 *f1 = new TF1(Form("D2_%s_1",hDaily_Deaths->GetName()),FuncD,hDaily_Deaths->GetXaxis()->GetXmin(),hDaily_Deaths->GetXaxis()->GetXmax(),4);
            f1->SetParameters(fDaily_D2->GetParameter(0),fDaily_D2->GetParameter(1),fDaily_D2->GetParameter(2),fDaily_D2->GetParameter(4));
            f1->SetLineColor(fDaily_D2->GetLineColor());
            f1->SetLineStyle(kDashed);
            f1->Draw("same");
            TF1 *f2 = new TF1(Form("DModel_%s_2",hDaily_Deaths->GetName()),FuncD,hDaily_Deaths->GetXaxis()->GetXmin(),hDaily_Deaths->GetXaxis()->GetXmax(),4);
            f2->SetParameters(fDaily_D2->GetParameter(0),fDaily_D2->GetParameter(3),fDaily_D2->GetParameter(2),fDaily_D2->GetParameter(4));
            f2->SetLineColor(fDaily_D2->GetLineColor());
            f2->SetLineStyle(kDashed);
            f2->Draw("same");
        }
    }

    //ESIR
    if(DoESIR) {
        fDaily_ESIR = new TF1(Form("ESIR_%s",hDaily_Deaths->GetName()),FuncESIR,XMin,hDaily_Deaths->GetXaxis()->GetXmax(),6);
        fDaily_ESIR->SetLineColor(Colors[2]);
        fDaily_ESIR->SetNpx(1000);

        const int NPars = fDaily_ESIR->GetNpar();
        Double_t *Pars = new Double_t[NPars];

        Pars[0] = 5e-6;
        Pars[1] = 10.;
        Pars[2] = 5e-6;

        Pars[3] = 500.;
        Pars[4] = 1e-4;
        Pars[5] = XMin;//hDaily_Deaths->GetXaxis()->GetBinLowEdge(hDaily_Deaths->FindFirstBinAbove(1));

        fDaily_ESIR->SetParameter(0,Pars[0]);
        fDaily_ESIR->SetParameter(1,Pars[1]);
        fDaily_ESIR->SetParameter(2,Pars[2]);
        fDaily_ESIR->SetParameter(3,Pars[3]);
        fDaily_ESIR->SetParameter(4,Pars[4]);
        fDaily_ESIR->FixParameter(5,Pars[5]);

        fDaily_ESIR->SetParLimits(0,1e-15,1e-5);
        fDaily_ESIR->SetParLimits(1,1.,50.);
        fDaily_ESIR->SetParLimits(2,1e-15,1e-5);
        fDaily_ESIR->SetParLimits(3,1e1,1e7);
        fDaily_ESIR->SetParLimits(4,1e-8,0.1);

        TFitResultPtr r = hDaily_Deaths->Fit(fDaily_ESIR,"S0","",XMin,XMax);
        fChi2ESIR = r->Chi2()/r->Ndf();

        r->Print("V");
        fDaily_ESIR->Draw("same");

        /*Create a histogram to hold the confidence intervals*/
        auto *herror = (TH1*)hDaily_Deaths->Clone();
        herror->Reset();
        herror->SetName(((TString)hDaily_Deaths->GetName()).Append("_errorESIR"));
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(herror);

        //Now the "hint" histogram has the fitted function values as the
        //bin contents and the confidence intervals as bin errors
        herror->SetStats(kFALSE);
        herror->SetFillColor(fDaily_ESIR->GetLineColor());
        herror->SetFillStyle(3002);
        herror->SetFillColorAlpha(fDaily_ESIR->GetLineColor(),0.5);
        herror->SetMarkerSize(0);
        herror->Draw("e3 same");
    }

    //ESIR2
    if(DoESIR2) {
        if(FullModel) {
            fDaily_ESIR2 = new TF1(Form("ESIR2_%s",hDaily_Deaths->GetName()),FuncESIR2Full,XMin,hDaily_Deaths->GetXaxis()->GetXmax(),9);
        }
        else {
            fDaily_ESIR2 = new TF1(Form("ESIR2_%s",hDaily_Deaths->GetName()),FuncESIR2,XMin,hDaily_Deaths->GetXaxis()->GetXmax(),6);
        }
        fDaily_ESIR2->SetLineColor(Colors[3]);
        fDaily_ESIR2->SetNpx(1000);

        const int NPars = fDaily_ESIR2->GetNpar();
        Double_t *Pars = new Double_t[NPars];

        if(FullModel) {

            Pars[0] = 5e-6;
            Pars[1] = 5.;
            Pars[2] = 5e-6;

            Pars[3] = 5e-6;
            Pars[4] = 15.;
            Pars[5] = 5e-6;

            Pars[6] = 500.;
            Pars[7] = 1e-4;
            Pars[8] = XMin;//hDaily_Deaths->GetXaxis()->GetBinLowEdge(hDaily_Deaths->FindFirstBinAbove(1));

            fDaily_ESIR2->SetParameter(0,Pars[0]);
            fDaily_ESIR2->SetParameter(1,Pars[1]);
            fDaily_ESIR2->SetParameter(2,Pars[2]);
            fDaily_ESIR2->SetParameter(3,Pars[3]);
            fDaily_ESIR2->SetParameter(4,Pars[4]);
            fDaily_ESIR2->SetParameter(5,Pars[5]);
            fDaily_ESIR2->SetParameter(6,Pars[6]);
            fDaily_ESIR2->SetParameter(7,Pars[7]);
            fDaily_ESIR2->FixParameter(8,Pars[8]);

            fDaily_ESIR2->SetParLimits(0,1e-15,1e-2);
            fDaily_ESIR2->SetParLimits(1,1.,50.);
            fDaily_ESIR2->SetParLimits(2,1e-15,1e-2);
            fDaily_ESIR2->SetParLimits(3,1e-15,1e-2);
            fDaily_ESIR2->SetParLimits(4,1.,50.);
            fDaily_ESIR2->SetParLimits(5,1e-15,1e-2);
            fDaily_ESIR2->SetParLimits(6,1e1,1e7);
            fDaily_ESIR2->SetParLimits(7,1e-15,1e-1);
        }
        else {
            Pars[0] = 5e-6;
            Pars[1] = 5.;
            Pars[2] = 15.;
            Pars[3] = 500.;
            Pars[4] = 1e-4;
            Pars[5] = XMin;//hDaily_Deaths->GetXaxis()->GetBinLowEdge(hDaily_Deaths->FindFirstBinAbove(1));

            fDaily_ESIR2->SetParameter(0,Pars[0]);
            fDaily_ESIR2->SetParameter(1,Pars[1]);
            fDaily_ESIR2->SetParameter(2,Pars[2]);
            fDaily_ESIR2->SetParameter(3,Pars[3]);
            fDaily_ESIR2->SetParameter(4,Pars[4]);
            fDaily_ESIR2->FixParameter(5,Pars[5]);

            fDaily_ESIR2->SetParLimits(0,1e-15,1e-2);
            fDaily_ESIR2->SetParLimits(1,1.,50.);
            fDaily_ESIR2->SetParLimits(2,1.,50.);
            fDaily_ESIR2->SetParLimits(3,1e1,1e7);
            fDaily_ESIR2->SetParLimits(4,1e-8,0.1);
        }

        TFitResultPtr r = hDaily_Deaths->Fit(fDaily_ESIR2,"S0","",XMin,XMax);
        fChi2ESIR2 = r->Chi2()/r->Ndf();

        r->Print("V");
        fDaily_ESIR2->Draw("same");

        /*Create a histogram to hold the confidence intervals*/
        auto *herror = (TH1*)hDaily_Deaths->Clone();
        herror->Reset();
        herror->SetName(((TString)hDaily_Deaths->GetName()).Append("_errorESIR2"));
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(herror);

        //Now the "hint" histogram has the fitted function values as the
        //bin contents and the confidence intervals as bin errors
        herror->SetStats(kFALSE);
        herror->SetFillColor(fDaily_ESIR2->GetLineColor());
        herror->SetFillStyle(3002);
        herror->SetFillColorAlpha(fDaily_ESIR2->GetLineColor(),0.5);
        herror->SetMarkerSize(0);
        herror->Draw("e3 same");
    }

    MaxY = hDaily_Deaths->GetMaximum() * 1.2;
    hDaily_Deaths->GetYaxis()->SetRangeUser(0,MaxY);
    hDaily_Deaths->GetXaxis()->SetRange(DateMin,DateMax);

    // Here, we remove some of the labels, in order to not have more than 40 labels on the graph
    Int_t NBinsInRange = DateMax-DateMin;
    Int_t Step = NBinsInRange/40;
    for(int i=1 ; i<=hDaily_Deaths->GetNbinsX() ; i+=Step) {
        for(int ii=1 ; ii<Step ; ii++) {
            if((i+ii) <= hDaily_Deaths->GetNbinsX()) {
                hDaily_Deaths->GetXaxis()->SetBinLabel(i+ii,"");
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

    YVal = gPad->GetFrame()->GetY2() * 0.95;
    XVal = gPad->GetFrame()->GetX2() * 0.78;
    Int_t NDY=0;

    Int_t NFuncs = (DoD+DoD2+DoESIR+DoESIR2);
    Float_t DY = gPad->GetFrame()->GetY2()*0.05;
    Float_t TextSize = 0.04;

    if(NFuncs>2) {
        TextSize = 0.02;
        DY = gPad->GetFrame()->GetY2()*0.03;
        XVal = gPad->GetFrame()->GetX2() * 0.88;
        YVal = gPad->GetFrame()->GetY2() * 0.97;
    }
    else {
        TextSize = 0.03;
        DY = gPad->GetFrame()->GetY2()*0.04;
        XVal = gPad->GetFrame()->GetX2() * 0.83;
        YVal = gPad->GetFrame()->GetY2() * 0.96;
    }

    // Print D
    if(DoD) {
        TLatex *text = new TLatex(XVal,YVal-DY*NDY,"D' model");
        text->SetTextColor(fDaily_D->GetLineColor());text->Draw();
        text->SetTextFont(132);
        text->SetTextSize(TextSize+0.01);
        NDY++;

        text = new TLatex(XVal,YVal-DY*NDY,Form("a = %.2f #pm %.2f",fDaily_D->GetParameter(0),fDaily_D->GetParError(0)));
        text->SetTextColor(fDaily_D->GetLineColor());text->Draw();
        text->SetTextFont(132);
        text->SetTextSize(TextSize);
        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("b = %.2f #pm %.2f",fDaily_D->GetParameter(1),fDaily_D->GetParError(1)));
        text->SetTextColor(fDaily_D->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("c = %.2e #pm %.2e",fDaily_D->GetParameter(2),fDaily_D->GetParError(2)));
        text->SetTextColor(fDaily_D->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;

        text = new TLatex(XVal,YVal-DY*NDY,Form("Chi2/ndf = %.2f",fChi2D));
        text->SetTextColor(fDaily_D->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);
        NDY++;
        NDY++;
    }

    // Print D2
    if(DoD2) {
        if(FullModel) {

            TLatex *text = new TLatex(XVal,YVal-DY*NDY,"D'2 full model");
            text->SetTextColor(fDaily_D2->GetLineColor());text->Draw();
            text->SetTextFont(132);
            text->SetTextSize(TextSize+0.01);
            NDY++;

            text = new TLatex(XVal,YVal-DY*NDY,Form("a1 = %.2f #pm %.2f",fDaily_D2->GetParameter(0),fDaily_D2->GetParError(0)));
            text->SetTextColor(fDaily_D2->GetLineColor());text->Draw();
            text->SetTextFont(132);
            text->SetTextSize(TextSize);
            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("b1 = %.2f #pm %.2f",fDaily_D2->GetParameter(1),fDaily_D2->GetParError(1)));
            text->SetTextColor(fDaily_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("c1 = %.2e #pm %.2e",fDaily_D2->GetParameter(2),fDaily_D2->GetParError(2)));
            text->SetTextColor(fDaily_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;

            text = new TLatex(XVal,YVal-DY*NDY,Form("a2 = %.2f #pm %.2f",fDaily_D2->GetParameter(3),fDaily_D2->GetParError(3)));
            text->SetTextColor(fDaily_D2->GetLineColor());text->Draw();
            text->SetTextFont(132);
            text->SetTextSize(TextSize);
            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("b2 = %.2f #pm %.2f",fDaily_D2->GetParameter(4),fDaily_D2->GetParError(4)));
            text->SetTextColor(fDaily_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("c2 = %.2e #pm %.2e",fDaily_D2->GetParameter(5),fDaily_D2->GetParError(5)));
            text->SetTextColor(fDaily_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;

            text = new TLatex(XVal,YVal-DY*NDY,Form("Chi2/ndf = %.2f",fChi2D2));
            text->SetTextColor(fDaily_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;
            NDY++;
        }
        else {
            TLatex *text = new TLatex(XVal,YVal-DY*NDY,"D'2 model");
            text->SetTextColor(fDaily_D2->GetLineColor());text->Draw();
            text->SetTextFont(132);
            text->SetTextSize(TextSize+0.01);
            NDY++;

            text = new TLatex(XVal,YVal-DY*NDY,Form("a = %.2f #pm %.2f",fDaily_D2->GetParameter(0),fDaily_D2->GetParError(0)));
            text->SetTextColor(fDaily_D2->GetLineColor());text->Draw();
            text->SetTextFont(132);
            text->SetTextSize(TextSize);
            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("b1 = %.2f #pm %.2f",fDaily_D2->GetParameter(1),fDaily_D2->GetParError(1)));
            text->SetTextColor(fDaily_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("b2 = %.2f #pm %.2f",fDaily_D2->GetParameter(3),fDaily_D2->GetParError(3)));
            text->SetTextColor(fDaily_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("c = %.2e #pm %.2e",fDaily_D2->GetParameter(2),fDaily_D2->GetParError(2)));
            text->SetTextColor(fDaily_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("Chi2/ndf = %.2f",fChi2D2));
            text->SetTextColor(fDaily_D2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);
            NDY++;
            NDY++;
        }
    }

    // Print ESIR
    if(DoESIR) {
        text = new TLatex(XVal,YVal-DY*NDY,"ESIR model");
        text->SetTextColor(fDaily_ESIR->GetLineColor());text->Draw();
        text->SetTextFont(132);
        text->SetTextSize(TextSize+0.01);
        NDY++;

        text = new TLatex(XVal,YVal-DY*NDY,Form("a = %.2e #pm %.2e",fDaily_ESIR->GetParameter(0),fDaily_ESIR->GetParError(0)));
        text->SetTextColor(fDaily_ESIR->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);

        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("b = %.2f #pm %.2f",fDaily_ESIR->GetParameter(1),fDaily_ESIR->GetParError(1)));
        text->SetTextColor(fDaily_ESIR->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);

        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("c = %.2e #pm %.2e",fDaily_ESIR->GetParameter(2),fDaily_ESIR->GetParError(2)));
        text->SetTextColor(fDaily_ESIR->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);

        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("a2 = %.3g #pm %.3g",fDaily_ESIR->GetParameter(3),fDaily_ESIR->GetParError(3)));
        text->SetTextColor(fDaily_ESIR->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);

        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("b2 = %.2e #pm %.2e",fDaily_ESIR->GetParameter(4),fDaily_ESIR->GetParError(4)));
        text->SetTextColor(fDaily_ESIR->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);

        NDY++;
        text = new TLatex(XVal,YVal-DY*NDY,Form("Chi2/ndf = %.2f",fChi2ESIR));
        text->SetTextColor(fDaily_ESIR->GetLineColor());text->Draw();
        text->SetTextSize(TextSize);
        text->SetTextFont(132);

        NDY++;
        NDY++;
    }

    // Print ESIR2
    if(DoESIR2) {
        if(FullModel) {
            text = new TLatex(XVal,YVal-DY*NDY,"ESIR2 full model");
            text->SetTextColor(fDaily_ESIR2->GetLineColor());text->Draw();
            text->SetTextFont(132);
            text->SetTextSize(TextSize+0.01);
            NDY++;

            text = new TLatex(XVal,YVal-DY*NDY,Form("a = %.2e #pm %.2e",fDaily_ESIR2->GetParameter(0),fDaily_ESIR2->GetParError(0)));
            text->SetTextColor(fDaily_ESIR2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);

            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("b = %.2f #pm %.2f",fDaily_ESIR2->GetParameter(1),fDaily_ESIR2->GetParError(1)));
            text->SetTextColor(fDaily_ESIR2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);

            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("c = %.2e #pm %.2e",fDaily_ESIR2->GetParameter(2),fDaily_ESIR2->GetParError(2)));
            text->SetTextColor(fDaily_ESIR2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);

            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("a' = %.2e #pm %.2e",fDaily_ESIR2->GetParameter(3),fDaily_ESIR2->GetParError(3)));
            text->SetTextColor(fDaily_ESIR2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);

            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("b' = %.2f #pm %.2f",fDaily_ESIR2->GetParameter(4),fDaily_ESIR2->GetParError(4)));
            text->SetTextColor(fDaily_ESIR2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);

            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("c' = %.2e #pm %.2e",fDaily_ESIR2->GetParameter(5),fDaily_ESIR2->GetParError(5)));
            text->SetTextColor(fDaily_ESIR2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);

            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("a2 = %.3g #pm %.3g",fDaily_ESIR2->GetParameter(6),fDaily_ESIR2->GetParError(6)));
            text->SetTextColor(fDaily_ESIR2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);

            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("b2 = %.2e #pm %.2e",fDaily_ESIR2->GetParameter(7),fDaily_ESIR2->GetParError(7)));
            text->SetTextColor(fDaily_ESIR2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);

            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("Chi2/ndf = %.2f",fChi2ESIR2));
            text->SetTextColor(fDaily_ESIR2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);

            NDY++;
            NDY++;
        }
        else {
            text = new TLatex(XVal,YVal-DY*NDY,"ESIR2 model");
            text->SetTextColor(fDaily_ESIR2->GetLineColor());text->Draw();
            text->SetTextFont(132);
            text->SetTextSize(TextSize+0.01);
            NDY++;

            text = new TLatex(XVal,YVal-DY*NDY,Form("a = %.2e #pm %.2e",fDaily_ESIR2->GetParameter(0),fDaily_ESIR2->GetParError(0)));
            text->SetTextColor(fDaily_ESIR2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);

            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("b = %.2f #pm %.2f",fDaily_ESIR2->GetParameter(1),fDaily_ESIR2->GetParError(1)));
            text->SetTextColor(fDaily_ESIR2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);

            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("b' = %.2f #pm %.2f",fDaily_ESIR2->GetParameter(2),fDaily_ESIR2->GetParError(2)));
            text->SetTextColor(fDaily_ESIR2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);

            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("a2 = %.3g #pm %.3g",fDaily_ESIR2->GetParameter(3),fDaily_ESIR2->GetParError(3)));
            text->SetTextColor(fDaily_ESIR2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);

            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("b2 = %.2e #pm %.2e",fDaily_ESIR2->GetParameter(4),fDaily_ESIR2->GetParError(4)));
            text->SetTextColor(fDaily_ESIR2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);

            NDY++;
            text = new TLatex(XVal,YVal-DY*NDY,Form("Chi2/ndf = %.2f",fChi2ESIR2));
            text->SetTextColor(fDaily_ESIR2->GetLineColor());text->Draw();
            text->SetTextSize(TextSize);
            text->SetTextFont(132);

            NDY++;
            NDY++;
        }
    }

    gSystem->mkdir("Pictures");
    TString OutputFileName = Form("Pictures/COVID_Deaths_Ana_%s",theCountry.Data());
    if(Smooth>1) OutputFileName.Append(Form("_%dDaysSmooth",Smooth));
    OutputFileName.Append(Form("_%s.png",vDates.at(vTotal_Deaths.size()-1).Data()));
    gPad->GetCanvas()->SaveAs(OutputFileName);
}

bool ReadData(TString filename, TString ReadUpTo)
{
    // The selected file is opended, and if not found, return with an error message
    ifstream file(filename);
    if(!file) {
        cout<<filename<<" not found"<<endl;
        return false;
    }

    // The vectors containing data are cleared from previous use
    vDates.clear();
    vDaily_Deaths.clear();
    vTotal_Deaths.clear();
    vDaily_Deaths_error.clear();

    TString Buffer;
    string line;

    // the first line is not used, we read it first to skip this line
    getline(file,line);

    Int_t Current_Year = 20;

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
        if(ReadUpTo !="" && Date == ReadUpTo) return true;
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


Double_t FuncESIR2(Double_t*xx,Double_t*pp) {

    Double_t a  = pp[0];
    Double_t b  = pp[1];
    Double_t c  = pp[2];
    Double_t a2  = pp[3];
    Double_t b2  = pp[4];

    Double_t t0 = pp[5];
    Double_t x   = xx[0] - t0;

    Double_t r = a/(2*a + TMath::Exp(-x/b)) + a/(2*a + TMath::Exp(-x/c));

    Double_t Daily = a2 * (1 - TMath::Exp(-r/b2) - r);

    return Daily;
}

Double_t FuncESIR(Double_t*xx,Double_t*pp) {

    Double_t a  = pp[0];
    Double_t b  = pp[1];
    Double_t c  = pp[2];

    Double_t a2  = pp[3];
    Double_t b2  = pp[4];

    Double_t t0 = pp[5];
    Double_t x   = xx[0] - t0;

    Double_t r = a/(2*c + TMath::Exp(-x/b));

    Double_t Daily = a2 * (1 - TMath::Exp(-r/b2) - r);

    return Daily;
}

Double_t FuncESIR2Full(Double_t*xx,Double_t*pp) {

    Double_t a  = pp[0];
    Double_t b  = pp[1];
    Double_t c  = pp[2];

    Double_t ap  = pp[3];
    Double_t bp  = pp[4];
    Double_t cp  = pp[5];

    Double_t a2  = pp[6];
    Double_t b2  = pp[7];

    Double_t t0 = pp[8];
    Double_t x   = xx[0] - t0;

    Double_t r = a/(2*c + TMath::Exp(-x/b)) + ap/(2*cp + TMath::Exp(-x/bp));

    Double_t Daily = a2 * (1 - TMath::Exp(-r/b2) - r);

    return Daily;
}

Double_t FuncD2(Double_t*xx,Double_t*pp) {

    Double_t a1  = pp[0];
    Double_t b1  = pp[1];
    Double_t c1  = pp[2];
    Double_t b2  = pp[3];

    Double_t t0 = pp[4];
    Double_t x   = xx[0] - t0;

    Double_t D1 = a1*TMath::Exp(x/b1)/(b1*TMath::Power(1+c1*TMath::Exp(x/b1),2));
    Double_t D2 = a1*TMath::Exp(x/b2)/(b2*TMath::Power(1+c1*TMath::Exp(x/b2),2));

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

    Double_t D1 = a1*TMath::Exp(x/b1)/(b1*TMath::Power(1+c1*TMath::Exp(x/b1),2));
    Double_t D2 = a2*TMath::Exp(x/b2)/(b2*TMath::Power(1+c2*TMath::Exp(x/b2),2));

    Double_t Daily = D1+D2;

    return Daily;
}

Double_t FuncD(Double_t*xx,Double_t*pp) {

    Double_t a1  = pp[0];
    Double_t b1  = pp[1];
    Double_t c1  = pp[2];

    Double_t t0 = pp[3];
    Double_t x   = xx[0] - t0;

    Double_t D1 = a1*TMath::Exp(x/b1)/(b1*TMath::Power(1+c1*TMath::Exp(x/b1),2));

    Double_t Daily = D1;

    return Daily;
}
