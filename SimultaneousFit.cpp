
#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooGamma.h>
#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooFit.h>
#include <RooCategory.h>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>
#include "TAxis.h"
#include "RooSimultaneous.h"
#include "RooConstVar.h"
#include "RooPolyVar.h"
#include "RooFormulaVar.h"
#include "RooArgList.h"

using namespace RooFit;
using namespace std;

void SimultaneousFit() {

    vector<const char*> file_paths = {
        "/home/comp-lab3/root/simFit/M166p5.root",
        "/home/comp-lab3/root/simFit/M169p5.root"};

    vector<string> fileNames = {"M166p5", "M169p5" };

    TH1F *hist_M166p5 = new TH1F("M166p5", "M166p5", 30, 0, 400);
    TH1F *hist_M169p5 = new TH1F("M169p5", "M169p5", 30, 0, 400);

 
    RooCategory sample("sample", "Sample Category");
    sample.defineType("M166p5");
    sample.defineType("M169p5");

    TFile* file_M166p5 = TFile::Open(file_paths[0]);
    TFile* file_M169p5 = TFile::Open(file_paths[1]);

    TTree* tree_file_M166p5 = (TTree*)file_M166p5->Get("Events");
    TTree* tree_file_M169p5 = (TTree*)file_M169p5->Get("Events");

    Float_t ToMaAn_MlMu_value; 

    tree_file_M166p5->SetBranchAddress("ToMaAn_MlMu", &ToMaAn_MlMu_value);
    tree_file_M169p5->SetBranchAddress("ToMaAn_MlMu", &ToMaAn_MlMu_value);

    Long64_t nentries_M166p5 = tree_file_M166p5->GetEntries();
    cout << "Number of entries for M166p5: " << nentries_M166p5 << endl;

    Long64_t nentries_M169p5 = tree_file_M169p5->GetEntries();
    cout << "Number of entries for M169p5: " << nentries_M169p5 << endl;

    for (Long64_t j = 0; j < nentries_M166p5; j++) {
        tree_file_M166p5->GetEntry(j);
        hist_M166p5->Fill(ToMaAn_MlMu_value); 
    }

    for (Long64_t t = 0; t < nentries_M169p5; t++) {
        tree_file_M169p5->GetEntry(t);
        hist_M169p5->Fill(ToMaAn_MlMu_value); 
    }
    TCanvas* cM166 = new TCanvas("cM166", "M166p5", 800, 600);
    hist_M166p5->Draw("Hist");
    
    TCanvas* cM169 = new TCanvas("cM169", "M169p5", 800, 600);
    hist_M169p5->Draw("Hist");


    double weightValue166 =833.9*38.25*1000*4/(9*nentries_M166p5);
    double weightValue169 =833.9*38.25*1000*4/(9*nentries_M169p5);

    RooRealVar weight_M166p5("weight_M166p5", "Weight for M166p5", weightValue166);
    RooRealVar weight_M169p5("weight_M169p5", "Weight for M169p5", weightValue169);

    cout << "--------------------------------------------------"<<endl;
    cout << "weight_M166p5 = "<<weight_M166p5<<endl;
    cout << "weight_M169p5 = "<<weight_M169p5<<endl;
    cout << "--------------------------------------------------"<<endl;

    RooRealVar x("x", "x", 0, 400);
    // RooRealVar para1Gauss166("para1Gauss166", "first coefficient166",  0,5);
    // RooRealVar para2Gauss166("para2Gauss166", "second coefficient166", 0,5); 

    // RooRealVar para1Gauss169("para1Gauss169", "first coefficient169",  0, 5);
    // RooRealVar para2Gauss169("para2Gauss169", "second coefficient169", 0, 5); 


    RooRealVar gamma_gamma166("gamma_gamma166", "shape parameter166",5, 5, 50);
    RooRealVar gamma_beta166("gamma_beta166", "scale parameter166" ,5, 5, 50);
    RooRealVar gamma_mu166("gamma_mu166","shift parameter166", 4, 0, 50);
    RooFormulaVar gammamu166("mu166p5", "x[0] < x[1] ? x[0] : x[1]", RooArgList(x, gamma_mu166));
    RooGamma gamma166("gamma166", "gamma PDF166", x, gamma_gamma166, gamma_beta166, gammamu166);

    RooRealVar gamma_gamma169("gamma_gamma169", "shape parameter169", 5 ,5, 20);
    RooRealVar gamma_beta169("gamma_beta169", "scale parameter169" ,5 , 5, 20);
    RooRealVar gamma_mu169("gamma_mu169","shift parameter169",4, 0, 50);
    RooFormulaVar gammamu169("mu169p5", "x[0] < x[1] ? x[0] : x[1]", RooArgList(x, gamma_mu169));
    RooGamma gamma169("gamma169", "gamma PDF169", x, gamma_gamma169, gamma_beta169, gammamu169);


    // RooRealVar M_top166("M_top166", "Top quark mass 166", 166.5);
    // RooRealVar M_top169("M_top169", "Top quark mass 169", 169.5);

    RooRealVar muGauss166p5("muGauss166p5", "Mean of Gaussian as a function of M_top", 0,80);
    // RooPolyVar muGauss166p5("muGauss166p5", "Mean of Gaussian as a function of M_top", M_top166, RooArgList(para1Gauss166, para2Gauss166));
    RooRealVar sigmaGauss166p5("sigmaGauss166p5", "width of gaussian for M166.5", 5,2, 30);
    RooGaussian gauss166p5("gauss166p5", "Gaussian PDF for M166.5", x, muGauss166p5, sigmaGauss166p5);

    RooRealVar muGauss169p5("muGauss169p5", "Mean of Gaussian as a function of M_top",0,80);
    // RooPolyVar muGauss169p5("muGauss169p5", "Mean of Gaussian as a function of M_top", M_top169, RooArgList( para1Gauss169, para2Gauss169));
    RooRealVar sigmaGauss169p5("sigmaGauss169p5", "width of gaussian for M169.5", 5,2, 30);
    RooGaussian gauss169p5("gauss169p5", "Gaussian PDF for M169.5", x, muGauss169p5, sigmaGauss169p5);


    RooRealVar alpha166("alpha166", "fraction of gaussian",0,1);
    RooRealVar alpha169("alpha169", "fraction of gaussian",0,1);


    RooAddPdf model166p5("model166", "gauss166+gamma", RooArgList(gauss166p5, gamma166), RooArgList(alpha166));
    RooAddPdf model169p5("model169", "gauss169+gamma", RooArgList(gauss169p5, gamma169), RooArgList(alpha169));

    RooDataHist* data_M166p5 = new RooDataHist("data_M166p5", "Dataset_M166p5", x, Import(*hist_M166p5));
    RooDataHist* data_M169p5 = new RooDataHist("data_M169p5", "Dataset_M169p5", x, Import(*hist_M169p5));


    RooDataSet* weighted_data_M166p5 = new RooDataSet("weighted_data_M166p5", "Weighted Dataset M166p5",
                                                    x, Import(*data_M166p5), WeightVar("weight_M166p5"));

    RooDataSet* weighted_data_M169p5 = new RooDataSet("weighted_data_M169p5", "Weighted Dataset M166p5",
                                                    x, Import(*data_M169p5), WeightVar("weight_M169p5"));

    RooDataSet* combData = new RooDataSet("combData", "combined data", x, Index(sample),Import("M166p5",*weighted_data_M166p5),Import("M169p5",*weighted_data_M169p5));
    // RooDataSet* combData = new RooDataSet("combData", "combined data", x, Index(sample),Import("M166p5",*data_M166p5),Import("M169p5",*data_M169p5));

    RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample) ;

    simPdf.addPdf(model166p5,"M166p5") ;
    simPdf.addPdf(model169p5,"M169p5") ;
 
    simPdf.fitTo(*combData, SumW2Error(true)); 

    // RooFitResult* fitResult = simPdf.fitTo(*combData, Save());
    // fitResult->Print("v");

    RooPlot *frame_166 = x.frame(Bins(30), Title("M166p5"));
    RooPlot *frame_169 = x.frame(Bins(30), Title("M169p5"));

    combData->plotOn(frame_166, Cut("sample==sample::M166p5"));
    simPdf.plotOn(frame_166, Slice(sample, "M166p5"), ProjWData(sample, *combData));
    simPdf.plotOn(frame_166, Slice(sample, "M166p5"), ProjWData(sample, *combData));

    combData->plotOn(frame_169, Cut("sample==sample::M169p5"));
    simPdf.plotOn(frame_169, Slice(sample, "M169p5"), ProjWData(sample, *combData));
    simPdf.plotOn(frame_169, Slice(sample, "M169p5"), ProjWData(sample, *combData));


    TCanvas* c = new TCanvas("TopMass","TopMass",800,400) ;
    c->Divide(2) ;
    c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame_166->GetYaxis()->SetTitleOffset(1.4) ; frame_166->Draw() ;
    c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame_169->GetYaxis()->SetTitleOffset(1.4) ; frame_169->Draw() ;
             c->SaveAs("/home/comp-lab3/root/simFit/twoRootSimultaneousFit.pdf");
    // file_M166p5->Close();
    // file_M169p5->Close();

 
}


