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
#include "TSystem.h"
#include "TGraphErrors.h"

using namespace  std;

////////////////////////////////////
/// Global parameters definition ///
////////////////////////////////////

// number of average days in the sliding window
Int_t fNSmoothing = 7;

// Models parameters
Bool_t fDoFullModel = true;
Bool_t fDoD = false;
Bool_t fDoD2 = true;
Bool_t fDoESIR = false;
Bool_t fDoESIR2 = true;

// Range of dates to be read from the input files
TString fReadDataFrom = "";
TString fReadDataTo = "";

// Range of dates for the X axis of the histogram
TString fAxisRangeFrom = "";
TString fAxisRangeTo = "";

// Range of dates for the fit
TString fFitRangeFrom = "";
TString fFitRangeTo = "";

///////////////////////////////////
/// Global variables definition ///
///////////////////////////////////

// dummy histogram to make some tests on the defined dates
TH1D *hDummyHist = nullptr;
TH1D *hDaily_Deaths = nullptr;

// declaration of global variables used in the code for D, D2, ESIR, ESIR2
Int_t fColors[4] = {kMagenta,kGreen,kBlue,kRed};

// Minimal number of deaths to start to be taken into acount
Int_t DeathsMin = 10;

// vectors containing the data
vector<TString> vDates;
vector<Double_t> vTotal_Deaths;
vector<Double_t> vDaily_Deaths;
vector<Double_t> vDaily_Deaths_error;

////////////////////////////
/// Functions definition ///
////////////////////////////

// to print all the parameters in the terminal
void PrintParameters(TString country_name);

// to define the models we want to fit
void SetModels(Bool_t DoD=false, Bool_t DoD2=true, Bool_t DoESIR=false, Bool_t DoESIR2=true, Bool_t FullModel=true);

// to change fNSmoothing
void SetSmoothing(Int_t Ndays=7);

// to change the range of dates to be read
void ReadDataRange(TString DateFrom="",TString DateTo="");

// to change axis range
void SetAxisRange(TString DateFrom="",TString DateTo="");

// to change the fit range
void SetFitRange(TString DateFrom="",TString DateTo="");

// Init histograms
void InitHistograms();

// Fonction used to read the data files
bool ReadData(TString filename);

// fonction to smooth the data on N sucessive days
void SmoothVector(Int_t N, vector<double> &data, vector<double> &data_err);

// Fit Functions definition
Double_t FuncD(Double_t*xx,Double_t*pp);
Double_t FuncD2(Double_t*xx,Double_t*pp);
Double_t FuncD2Full(Double_t*xx,Double_t*pp);

Double_t FuncESIR(Double_t*xx,Double_t*pp);
Double_t FuncESIR2(Double_t*xx,Double_t*pp);
Double_t FuncESIR2Full(Double_t*xx,Double_t*pp);

////////////////////////////
/// Printouts coloration ///
////////////////////////////

#define ERR_MESS  std::cout<<"\e[0;3;91m -- ERROR   : "
#define WARN_MESS std::cout<<"\e[0;3;93m -- WARNNING: "
#define INFO_MESS std::cout<<"\e[0;3;92m -- INFO    : \e[0;3;94m"
#define TITLE_MESS std::cout<<" \e[0;4;92m"
#define END_MESS  "\e[0;3m"
#define ENDL END_MESS<<std::endl
