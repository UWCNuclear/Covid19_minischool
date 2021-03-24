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

using namespace  std;

// declaration of global variables used in the code

Int_t NGraphs = 0;
Int_t Colors[4] = {kGreen,kRed,kSpring+10,kAzure+1};

// Minimal number of deaths to start to be taken into acount
Int_t DeathsMin = 10;

Int_t DateMin=0;
Int_t DateMax=0;

TLegend *legend = nullptr;

Double_t MaxY = 0;
Double_t MaxY_Tot = 0;

vector<TString> vDates;
vector<Double_t> vTotal_Deaths;
vector<Double_t> vDaily_Deaths;
vector<Double_t> vDaily_Deaths_error;

Bool_t DoEst = false;
Bool_t DoParis = false;

TF1 *fDaily_ESIR = nullptr;
TF1 *fDaily_ESIR2 = nullptr;

TF1 *fDaily_D = nullptr;
TF1 *fDaily_D2 = nullptr;

Double_t fChi2D=0;
Double_t fChi2D2=0;
Double_t fChi2ESIR=0;
Double_t fChi2ESIR2=0;

// Definition of the functions used in the code, and initialization of the default parameters

// Main function used to plot and fit the data with our prefered models
void AnaPerCountry(TString theCountry, Int_t Smooth = 3, TString ReadUpTo = "", Bool_t FullModel = false, Bool_t DoD=0, Bool_t DoD2=0, Bool_t DoESIR=0, Bool_t DoESIR2=1, TString theDateMax="",TString theDateMin="",TString FitDateMin="", TString FitDateMax="");

// Fonction used to read the data files
bool ReadData(TString filename, TString ReadUpTo);

// fonction to smooth the data on N sucessive days
void SmoothVector(Int_t N, vector<double> &data, vector<double> &data_err);

Double_t FuncESIR(Double_t*xx,Double_t*pp);
Double_t FuncESIR2(Double_t*xx,Double_t*pp);
Double_t FuncESIR2Full(Double_t*xx,Double_t*pp);

Double_t FuncD(Double_t*xx,Double_t*pp);
Double_t FuncD2(Double_t*xx,Double_t*pp);
Double_t FuncD2Full(Double_t*xx,Double_t*pp);
