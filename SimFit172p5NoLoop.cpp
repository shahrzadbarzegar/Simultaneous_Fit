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
#include <TLegend.h>
#include "THStack.h"
#include <RooMinimizer.h>

using namespace RooFit;
using namespace std;
//int debug = 1;


void SimFit172p5NoLoop() {
    RooRealVar  mass("ToMaAn_MlMu", "ToMaAn_MlMu", 0, 300); 

	vector<RooRealVar*> wVars;
	vector<RooDataSet*> datas;	
	vector<RooDataSet*> Wdatas;

	vector<RooRealVar*> wVarsBkg;
	vector<RooDataSet*> datasBkg;	
	vector<RooDataSet*> WdatasBkg;

	vector<RooRealVar*> wVarsSTanti;
	vector<RooDataSet*> datasSTanti;	
	vector<RooDataSet*> WdatasSTanti;

	vector<RooRealVar*> wVarsST;
	vector<RooDataSet*> datasST;	
	vector<RooDataSet*> WdatasST;


	stringstream fname;
	stringstream fnameBkg;
	stringstream fnameSTanti;
	stringstream fnameST;


    const char* file_paths ={"/home/comp-lab3/root/simFit/M172p5.root"};
    string fileNames = {"M172p5"};

    const char* BkgFile ={"/home/comp-lab3/root/simFit/bkgSamples/TTTo2L2Nu_172p5.root"};
    string Bkg_Names={"TTTo2L2Nu_172p5"};

    const char* Santitop_File ={"/home/comp-lab3/root/simFit/bkgSamples/ST_t-channel_antitop_172p5.root"};
    string Santitop_Names ={"ST_t-channel_antitop_172p5"};

    const char* Stop_File ={"/home/comp-lab3/root/simFit/bkgSamples/ST_t-channel_top_172p5.root"};
    string Stop_Names ={"ST_t-channel_top_172p5"};

    vector<string> BkgName= {"TTTo2L2Nu_","ST_t-channel_top_","ST_t-channel_antitop_"};

    Float_t ToMaAn_MlMu_value;

	// std::vector<double> bkg_expected;
	// bkg_expected.push_back(833.9*16.8*1000/9); //dilepton  https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
	// bkg_expected.push_back(134.2*16.8*1000); //single top  https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopNNLORef
	// bkg_expected.push_back(80.*16.8*1000); // anti top
	// bkg_expected.push_back(79.3*16.8*1000/2);
	// bkg_expected.push_back(79.3*16.8*1000/2);


	
	// --- reading dileptonic bkg root files ---	
        TFile* file_TopBkg = TFile::Open(BkgFile);
        TTree* tree_TopBkg = (TTree*)file_TopBkg->Get("Events");
        tree_TopBkg->SetBranchAddress("ToMaAn_MlMu", &ToMaAn_MlMu_value);

        Long64_t nentries_TopBkg = tree_TopBkg->GetEntries();
        cout << "Number of entries for Tree---> " << BkgFile << ": " << nentries_TopBkg << endl;

        TH1I * hisBkg =(TH1I*)file_TopBkg->Get("nEles");
        double nEntriesBkg = hisBkg->GetEntries();

        RooRealVar *nGenBkg = new RooRealVar("nGenBkg","nGenBkg",nEntriesBkg);

        fnameBkg.str("");
		fnameBkg << "data_" << Bkg_Names;
        RooFormulaVar *weight_TopBkg= new RooFormulaVar(fnameBkg.str().c_str(), "(833.9*16.8*1000/(9*@0))", RooArgSet(*nGenBkg)); 
        cout<<"________ "<<"Bkg sample: "<<nGenBkg->getVal()<<" and nEntriesBkg: "<<(int)nEntriesBkg<<",   weight should be: "<< 80.*16.8*1000/nEntriesBkg <<endl;

	    datasBkg.push_back(new RooDataSet(fnameBkg.str().c_str(), fnameBkg.str().c_str(),tree_TopBkg, mass, ""));
		datasBkg[datasBkg.size()-1]->Print();
		wVarsBkg.push_back((RooRealVar*)datasBkg[datasBkg.size()-1]->addColumn(*weight_TopBkg));
		cout<<" weight Bkg Dileptonic:  "<<weight_TopBkg->evaluate()<<endl;
		WdatasBkg.push_back(new RooDataSet(fnameBkg.str().c_str(), fnameBkg.str().c_str(),datasBkg[datasBkg.size()-1],*datasBkg[datasBkg.size()-1]->get(),0, wVarsBkg[wVarsBkg.size()-1]->GetName()));
		cout<<" should be the same with: "<<wVarsBkg[wVarsBkg.size()-1]->getVal()<<endl;
		cout<<" , name:  "<<wVarsBkg[wVarsBkg.size()-1]->GetName()<<endl;
		WdatasBkg[WdatasBkg.size()-1]->Print();		
        


	// --- reading Single anti_top root files ---	

        TFile* file_TopSTanti = TFile::Open(Santitop_File);
        TTree* tree_TopSTanti = (TTree*)file_TopSTanti->Get("Events");
        tree_TopSTanti->SetBranchAddress("ToMaAn_MlMu", &ToMaAn_MlMu_value);

        Long64_t nentries_TopSTanti = tree_TopSTanti->GetEntries();
        cout << "Number of entries for Tree---> " << Santitop_File << ": " << nentries_TopSTanti << endl;

        TH1I * hisSTanti =(TH1I*)file_TopSTanti->Get("nEles");
        double nEntriesSTanti = hisSTanti->GetEntries();

        RooRealVar *nGenSTanti = new RooRealVar("nGenSTanti","nGenSTanti",nEntriesSTanti);
    
        fnameSTanti.str("");
		fnameSTanti << "data_" << Santitop_Names;
        RooFormulaVar *weight_TopSTanti= new RooFormulaVar(fnameSTanti.str().c_str(), "(80.*16.8*1000/@0)", RooArgSet(*nGenSTanti)); 
        cout<<"________ "<<"ST_anti sample: "<<nGenSTanti->getVal()<<" and nEntriesSTanti: "<<(int)nEntriesSTanti<<",   weight should be: "<< 79.3*16.8*1000/(2*nEntriesSTanti) <<endl;

	    datasSTanti.push_back(new RooDataSet(fnameSTanti.str().c_str(), fnameSTanti.str().c_str(),tree_TopSTanti, mass, ""));
		datasSTanti[datasSTanti.size()-1]->Print();
		wVarsSTanti.push_back((RooRealVar*)datasSTanti[datasSTanti.size()-1]->addColumn(*weight_TopSTanti));
		cout<<" weight Single Anti top:  "<<weight_TopSTanti->evaluate()<<endl;
		WdatasSTanti.push_back(new RooDataSet(fnameSTanti.str().c_str(), fnameSTanti.str().c_str(),datasSTanti[datasSTanti.size()-1],*datasSTanti[datasSTanti.size()-1]->get(),0, wVarsSTanti[wVarsSTanti.size()-1]->GetName()));
		cout<<" should be the same with: "<<wVarsSTanti[wVarsSTanti.size()-1]->getVal()<<endl;
		cout<<" , name:  "<<wVarsSTanti[wVarsSTanti.size()-1]->GetName()<<endl;
		WdatasSTanti[WdatasSTanti.size()-1]->Print();		
        


	// --- reading Single top root files ---	
        TFile* file_TopST = TFile::Open(Stop_File);
        TTree* tree_TopST = (TTree*)file_TopST->Get("Events");
        tree_TopST->SetBranchAddress("ToMaAn_MlMu", &ToMaAn_MlMu_value);

        Long64_t nentries_TopST = tree_TopST->GetEntries();
        cout << "Number of entries for Tree---> " << Stop_File << ": " << nentries_TopST << endl;

        TH1I * hisST =(TH1I*)file_TopST->Get("nEles");
        double nEntriesST = hisST->GetEntries();

        RooRealVar *nGenST = new RooRealVar("nGenST","nGenST",nEntriesST);
    
        fnameST.str("");
		fnameST << "data_" << Stop_Names;
        RooFormulaVar *weight_TopST= new RooFormulaVar(fnameST.str().c_str(), "(134.2*16.8*1000/@0)", RooArgSet(*nGenST)); 
        cout<<"________ "<<"ST sample: "<<nGenST->getVal()<<" and nEntriesST: "<<(int)nEntriesST<<",   weight should be: "<< 79.3*16.8*1000/(2*nEntriesST) <<endl;

	    datasST.push_back(new RooDataSet(fnameST.str().c_str(), fnameST.str().c_str(),tree_TopST, mass, ""));
		datasST[datasST.size()-1]->Print();
		wVarsST.push_back((RooRealVar*)datasST[datasST.size()-1]->addColumn(*weight_TopST));
		cout<<" weight Single top :  "<<weight_TopST->evaluate()<<endl;
		WdatasST.push_back(new RooDataSet(fnameST.str().c_str(), fnameST.str().c_str(),datasST[datasST.size()-1],*datasST[datasST.size()-1]->get(),0, wVarsST[wVarsST.size()-1]->GetName()));
		cout<<" should be the same with: "<<wVarsST[wVarsST.size()-1]->getVal()<<endl;
		cout<<" , name:  "<<wVarsST[wVarsST.size()-1]->GetName()<<endl;
		WdatasST[WdatasST.size()-1]->Print();		
        


// --- reading Signal root files ---	
        TH1F* hist = new TH1F((fileNames).c_str(), (fileNames).c_str(), 100, 0, 250);

        TFile* file_Top = TFile::Open(file_paths);
        TTree* tree_Top = (TTree*)file_Top->Get("Events");
        tree_Top->SetBranchAddress("ToMaAn_MlMu", &ToMaAn_MlMu_value);


        Long64_t nentries_Top = tree_Top->GetEntries();
        cout << "Number of entries for Tree---> " << fileNames << ": " << nentries_Top << endl;

        TH1I * his =(TH1I*)file_Top->Get("nEles");
        double nEntries = his->GetEntries();

        RooRealVar * nGen = new RooRealVar("nGen","nGen",nEntries);

        fname.str("");
		fname << "data_" << fileNames;
        RooFormulaVar *weight_Top= new RooFormulaVar(fname.str().c_str(), "(833.9*16.8*1000*4/(9*@0))", RooArgSet(*nGen)); 
        cout<<"________ "<<fname.str().c_str()<<" sample: "<<nGen->getVal()<<" and nEntries: "<<(int)nEntries<<",   weight should be: "<<833.9*16.8*1000*4/(9*nEntries)<<endl;

	    datas.push_back(new RooDataSet(fname.str().c_str(), fname.str().c_str(),tree_Top, mass, ""));
		datas[datas.size()-1]->Print();
		wVars.push_back((RooRealVar*)datas[datas.size()-1]->addColumn(*weight_Top));
		cout<<" weight Signal:  "<<weight_Top ->evaluate()<<endl;
		Wdatas.push_back(new RooDataSet(fname.str().c_str(), fname.str().c_str(),datas[datas.size()-1],*datas[datas.size()-1]->get(),0, wVars[wVars.size()-1]->GetName()));
		cout<<" should be the same with: "<<wVars[wVars.size()-1]->getVal()<<endl;
		cout<<" , name:  "<<wVars[wVars.size()-1]->GetName()<<endl;
		Wdatas[Wdatas.size()-1]->Print();		



        RooDataSet combDataWeighted(*Wdatas[0]);
        combDataWeighted.append(*WdatasBkg[0]);
        combDataWeighted.append(*WdatasSTanti[0]);
        combDataWeighted.append(*WdatasST[0]);
            


    //Without Errors for Parameters
        RooRealVar* alpha = new RooRealVar("alpha", "fraction of Gaussian", 0.463723);
        alpha->setConstant(kTRUE);

        RooRealVar* gamma_gamma = new RooRealVar("gamma_Shape", "shape parameter",2.29167);
        gamma_gamma->setConstant(kTRUE);

        RooRealVar* gamma_beta = new RooRealVar("gamma_scale", "scale parameter" ,33.8293);
        gamma_beta->setConstant(kTRUE);

        RooRealVar* gamma_mu = new RooRealVar("gamma_shift" , "shift parameter" , 9.29586);
        gamma_mu->setConstant(kTRUE);

    //With Errors for Parameters
        // RooRealVar* alpha = new RooRealVar("alpha", "fraction of Gaussian", 0.463723 , 0.4582638 , 0.4691822 );
        // RooRealVar* gamma_gamma = new RooRealVar("gamma_Shape", "shape parameter",2.29167, 2.2695881 , 2.3137519);
        // RooRealVar* gamma_beta = new RooRealVar("gamma_scale", "scale parameter" ,33.8293, 33.497177, 34.161423);
        // RooRealVar* gamma_mu = new RooRealVar("gamma_shift" , "shift parameter" , 9.29586, 9.114378, 9.477342);


        // double slope_mean = 0.47396;
        // double intercept_mean = -19.522;
        // double slope_sigma = 0.0620741;
        // double intercept_sigma = 6.98514;
        // RooFormulaVar linear_mean("linearFunc_mean", "Linear Function", Form("%f * @0 + %f", slope_mean, intercept_mean), RooArgList(mass));
        // RooFormulaVar linear_sigma("linearFunc_sigma", "Linear Function", Form("%f * @0 + %f", slope_sigma, intercept_sigma), RooArgList(mass));

        RooRealVar* M_top = new RooRealVar("M_top", "Top quark mass" , 172, 150, 180);

        RooFormulaVar mean_gauss("mean_gauss", "Mean of Gaussian","(0.47396) * @0 + (-19.522)", RooArgList(*M_top));
        RooFormulaVar sigma_gauss("sigma_gauss", "Sigma of Gaussian","(0.0620741) * @0 + (6.98514)", RooArgList(*M_top));

        RooGaussian* gauss = new RooGaussian(("gauss" + fileNames).c_str(), ("Gaussian PDF for" + fileNames).c_str(), mass, mean_gauss, sigma_gauss);

        RooFormulaVar* gammamu = new RooFormulaVar((fileNames + "_formula").c_str(), "x[0] < x[1] ? x[0] : x[1]", RooArgList(mass, *gamma_mu));
        RooGamma* gamma = new RooGamma(("gamma " + fileNames).c_str(), ("gamma PDF for " + fileNames).c_str(), mass, *gamma_gamma, *gamma_beta, *gammamu);

        RooAddPdf* model = new RooAddPdf(fileNames.c_str(), "Gaussian + Gamma PDF", RooArgList(*gauss, *gamma), *alpha); 



        model->fitTo(combDataWeighted, SumW2Error(true), PrintLevel(3),Binning(100));


        TCanvas* c = new TCanvas("fit_172p5", "Fit for Top Quark Mass 172.5 GeV", 800, 600);

        RooPlot* frame = mass.frame();
        combDataWeighted.plotOn(frame, RooFit::SumW2Error(true)); 
        model->plotOn(frame);           

        frame->SetTitle("Fit Result for Top Quark Mass 172.5 GeV");
        frame->Draw();

        c->SaveAs("/home/comp-lab3/root/simFit/Fit_172p5.pdf");
        // delete frame;
        // delete c;





        // Pull distribution
        RooHist* pullHist = frame->pullHist();
        RooPlot* pullFrame = mass.frame(Title("Pull Distribution"));
        pullFrame->addPlotable(pullHist, "P");

        TCanvas* cPull = new TCanvas("cPull", "Pull Distribution", 800, 600);
        pullFrame->GetYaxis()->SetTitle("Pull");
        pullFrame->GetYaxis()->SetRangeUser(-5, 5); 
        pullFrame->Draw();
        cPull->SaveAs("/home/comp-lab3/root/simFit/Pull_172p5.pdf");




        // Create the NLL function
        RooAbsReal* nll = model->createNLL(combDataWeighted, RooFit::Extended());

        double minMTop = 167;
        double maxMTop = 179;
        int nPoints = 100;
        double stepSize = (maxMTop - minMTop) / nPoints;

        std::vector<double> mTopValues;
        std::vector<double> nllValues;

        for (int i = 0; i <= nPoints; ++i) {
            double mTop = minMTop + i * stepSize;
            M_top->setVal(mTop);      
            M_top->setConstant(true); 
            nllValues.push_back(nll->getVal());
            mTopValues.push_back(mTop);
            M_top->setConstant(false); 
        }

        double minNLL = *std::min_element(nllValues.begin(), nllValues.end());
        for (auto& val : nllValues) {
            val -= minNLL; 
        }

        TGraph* nllGraph = new TGraph(nPoints + 1, mTopValues.data(), nllValues.data());
        nllGraph->SetTitle("Negative Log-Likelihood Curve;M_{top} (GeV);-Log(Likelihood)");
        nllGraph->SetLineColor(kBlue);
        nllGraph->SetLineWidth(2);

        TCanvas* cNLL = new TCanvas("cNLL", "Negative Log-Likelihood", 800, 600);
        nllGraph->Draw("AL");
        cNLL->SaveAs("/home/comp-lab3/root/simFit/NLL_Curve_172p5.pdf");

        // Clean up
        // delete nll;





    }

