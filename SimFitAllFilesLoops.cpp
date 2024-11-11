#include <RooGaussian.h>
#include <RooGamma.h>
#include <RooAddPdf.h>
#include <RooFit.h>
#include <TH1I.h>
#include <TH1F.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <ios>
#include <fstream>
#include "TH2.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooConstVar.h"
#include "RooFormulaVar.h"
#include "RooCategory.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TColorWheel.h"
#include "TColorGradient.h"
#include "RooPlot.h"
#include "RooRandom.h"

using namespace RooFit;
using namespace std;

void SimFitAllFilesLoops() {
//gDebug=2;
    RooRealVar x("ToMaAn_MlMu", "ToMaAn_MlMu", 0, 250);
    RooCategory sample("sample", "Sample Category");
    vector<RooFormulaVar*> weight_semilepttbar;
    vector<RooRealVar*> wVars;
    vector<RooDataSet*> datas;	
    vector<RooDataSet*> Wdatas;

    vector<RooAddPdf*> models;
    stringstream fname;

    vector<double> weightValue;

    vector<const char*> file_paths = {
	"/home/comp-lab3/root/simFit/M166p5.root",
	"/home/comp-lab3/root/simFit/M169p5.root",
	"/home/comp-lab3/root/simFit/M171p5.root",
	"/home/comp-lab3/root/simFit/M172p5.root",
	"/home/comp-lab3/root/simFit/M173p5.root",
	"/home/comp-lab3/root/simFit/M175p5.root",
	"/home/comp-lab3/root/simFit/M178p5.root" };


    vector<string> fileNames = {"M166p5", "M169p5", "M171p5", "M172p5", "M173p5", "M175p5", "M178p5"};
    vector<double> TopMasses = {166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5};

    vector<double> allentries;
    vector<TH1F*> allHists;


    Float_t ToMaAn_MlMu_value;

    for (int h = 0; h < fileNames.size(); h++) {

        TH1F* hist = new TH1F(fileNames[h].c_str(), fileNames[h].c_str(), 20, 0, 250);

        TFile* file_Top = TFile::Open(file_paths[h]);
        TTree* tree_Top = (TTree*)file_Top->Get("Events");
        tree_Top->SetBranchAddress("ToMaAn_MlMu", &ToMaAn_MlMu_value);

        Long64_t nentries_Top = tree_Top->GetEntries();
        cout << "Number of entries for Tree---> " << fileNames[h] << ": " << nentries_Top << endl;

        TH1I * his =(TH1I*)file_Top->Get("nEles");
        double nEntries = his->GetEntries();



        RooRealVar * nGen = new RooRealVar("nGen","nGen",nEntries);

        fname.str("");
		fname << "data" << fileNames[h];
        RooFormulaVar *weight_Top= new RooFormulaVar(fname.str().c_str(), "(833.9*38.25*1000*4/(9*@0))", RooArgSet(*nGen)); 
        cout<<"________ "<<fname.str().c_str()<<" sample: "<<nGen->getVal()<<" and nEntries: "<<(int)nEntries<<",   weight should be: "<<833.9*38.25*1000*4/(9*nEntries)<<endl;
        weight_semilepttbar.push_back(weight_Top);
        sample.defineType(fileNames[h].c_str());

///////////////////////////////////////////////////////
        // TCanvas* c = new TCanvas(( "c" + fileNames[h]).c_str(), fileNames[h].c_str(), 800, 600);
        // hist->Draw("Hist");
        // c->SaveAs(("/home/comp-lab3/root/simFit/"+fileNames[h] + "_hist.pdf").c_str());
///////////////////////////////////////////////////////


	    datas.push_back(new RooDataSet(fname.str().c_str(), fname.str().c_str(),tree_Top, x, ""));
		datas[datas.size()-1]->Print();
		wVars.push_back((RooRealVar*)datas[datas.size()-1]->addColumn(*weight_semilepttbar[h]));
		cout<<" weight_semilepttbar:  "<<weight_semilepttbar[h]->evaluate()<<endl;
		Wdatas.push_back(new RooDataSet(fname.str().c_str(), fname.str().c_str(),datas[datas.size()-1],*datas[datas.size()-1]->get(),0, wVars[wVars.size()-1]->GetName()));
		cout<<" should be the same with: "<<wVars[wVars.size()-1]->getVal()<<endl;
		cout<<" , name:  "<<wVars[wVars.size()-1]->GetName()<<endl;
		Wdatas[Wdatas.size()-1]->Print();		
        }

	 RooDataSet combDataWeighted("combDataWeighted", "combined data",  RooArgSet(x,*wVars[0]), Index(sample),
            Import("M166p5",*Wdatas[0]),WeightVar("dataM166p5"));
	 RooDataSet combData169p5("combDataWeighted169p51", "combined data", RooArgSet(x,*wVars[1]), Index(sample),
            Import("M169p5",*Wdatas[1]),WeightVar("dataM169p5"));
	 RooDataSet combData171p5("combDataWeighted171p5", "combined data" ,  RooArgSet(x,*wVars[2]), Index(sample),
             Import("M171p5",*Wdatas[2]),WeightVar("dataM171p5"));
	 RooDataSet combData172p5("combDataWeighted172p51", "combined data" ,  RooArgSet(x,*wVars[3]), Index(sample),
             Import("M172p5",*Wdatas[3]),WeightVar("dataM172p5"));
	 RooDataSet combData173p5("combDataWeighted173p5", "combined data" ,  RooArgSet(x,*wVars[4]), Index(sample),
             Import("M173p5",*Wdatas[4]),WeightVar("dataM173p5"));
	 RooDataSet combData175p5("combDataWeighted175p51", "combined data" , RooArgSet(x,*wVars[5]), Index(sample),
             Import("M175p5",*Wdatas[5]),WeightVar("dataM175p5"));
	 RooDataSet combData178p5("combDataWeighted178p51", "combined data" , RooArgSet(x,*wVars[6]), Index(sample),
             Import("M178p5",*Wdatas[6]),WeightVar("dataM178p5"));

        combDataWeighted.append(combData169p5);
        combDataWeighted.append(combData171p5);
        combDataWeighted.append(combData172p5);
        combDataWeighted.append(combData173p5);
        combDataWeighted.append(combData175p5);
        combDataWeighted.append(combData178p5);

// To fix the fraction of Gausssin function
        RooRealVar* alpha = new RooRealVar("alpha", "fraction of gaussian", 0, 1); 

    for (size_t w = 0; w < fileNames.size(); w++) {
// To calculate different alphas for each model
        // RooRealVar* alpha = new RooRealVar(("alpha_" + fileNames[w]).c_str(), "fraction of gaussian", 0, 1); 

//To examine the dependency of the Gaussian mean value on the Top mass 
        // RooRealVar* para1Gauss = new RooRealVar(("para1Gauss" + fileNames[w]).c_str(), ("first coefficient" + fileNames[w]).c_str(), 0, 1);
        // RooRealVar* para2Gauss = new RooRealVar(("para2Gauss" + fileNames[w]).c_str(), ("second coefficient" + fileNames[w]).c_str(), 1, 50);

        RooRealVar* M_top = new RooRealVar(("M_top" + fileNames[w]).c_str(), ("Top quark mass" + fileNames[w]).c_str(), TopMasses[w]);
        // RooPolyVar* muGauss = new RooPolyVar(("muGauss" + fileNames[w]).c_str(), ("Mean of Gaussian as a function of M_top" + fileNames[w]).c_str(), *mass, RooArgList(*para1Gauss, *para2Gauss));
        RooRealVar* muGauss = new RooRealVar(("mean_" + fileNames[w]).c_str(), ("mean_" + fileNames[w]).c_str() ,10, 80);
        RooRealVar* sigmaGauss = new RooRealVar(("sigmaGauss" + fileNames[w]).c_str(), ("width of gaussian for" + fileNames[w]).c_str() ,5,0, 50);
        RooGaussian* gauss = new RooGaussian(("gauss166" + fileNames[w]).c_str(), ("Gaussian PDF for" + fileNames[w]).c_str(), x, *muGauss, *sigmaGauss);

        RooRealVar* gamma_gamma = new RooRealVar(("gamma_" + fileNames[w]).c_str(), ("shape parameter" + fileNames[w]).c_str(),5, 2, 60);
        RooRealVar* gamma_beta = new RooRealVar(("beta_" + fileNames[w]).c_str(), ("scale parameter" + fileNames[w]).c_str(),5, 3, 30);
        RooRealVar* gamma_mu = new RooRealVar(("mu_" + fileNames[w]).c_str(), ("shift parameter" + fileNames[w]).c_str(), 4, 0, 30);

        RooFormulaVar* gammamu = new RooFormulaVar((fileNames[w] + "_formula").c_str(), "x[0] < x[1] ? x[0] : x[1]", RooArgList(x, *gamma_mu));
        RooGamma* gamma = new RooGamma(("gamma " + fileNames[w]).c_str(), ("gamma PDF for " + fileNames[w]).c_str(), x, *gamma_gamma, *gamma_beta, *gammamu);

        RooAddPdf* model = new RooAddPdf(fileNames[w].c_str(), "Gaussian + Gamma PDF", RooArgList(*gauss, *gamma), *alpha); 

        models.push_back(model);

    }



        RooSimultaneous simPdf("simPdf", "simultaneous pdf", sample);

        for (int m = 0; m < models.size(); m++) {
            simPdf.addPdf(*models[m], fileNames[m].c_str());
        }

        simPdf.fitTo(combDataWeighted, SumW2Error(true));
        
        for (int p = 0; p < fileNames.size(); p++) {

            TCanvas* c = new TCanvas(("fit_" + fileNames[p]).c_str(), ("Fit " + fileNames[p]).c_str(), 800, 600);
            RooPlot* frame = x.frame();
            combDataWeighted.plotOn(frame, Cut(Form("sample==sample::%s", fileNames[p].c_str())));
            simPdf.plotOn(frame, Slice(sample, fileNames[p].c_str()), ProjWData(sample,combDataWeighted));
            frame->Draw();
            frame->SetTitle(Form("Fit Result %s", fileNames[p].c_str()));
            c->SaveAs(("/home/comp-lab3/root/simFit/SimultaneousFit_" + fileNames[p] + ".pdf").c_str());
        }
        
    }

