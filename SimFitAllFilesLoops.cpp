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
// int debug = 1;
void SimFitAllFilesLoops() {
//gDebug=2;
    RooRealVar x("ToMaAn_MlMu", "ToMaAn_MlMu", 0, 250);
    // RooRealVar * mass = new RooRealVar("ToMaAn_MlMu", "ToMaAn_MlMu", 0, 300); 
    RooCategory sample("sample", "Sample Category");
	vector<RooFormulaVar*> weight_semilepttbar;
	vector<RooRealVar*> wVars;
	vector<RooDataSet*> datas;	
	vector<RooDataSet*> Wdatas;
    // vector<RooDataHist*> datas;	
	// vector<RooDataHist*> Wdatas;
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

// if(debug==1)cout << "debug  1 :"<<endl;

    // for (const auto& name : fileNames) {
    //     sample.defineType(name.c_str());
    // }
// if(debug==1)cout << "debug  2 :"<<endl;
    Float_t ToMaAn_MlMu_value;

    for (int h = 0; h < fileNames.size(); h++) {

        TH1F* hist = new TH1F(fileNames[h].c_str(), fileNames[h].c_str(), 20, 0, 250);

        TFile* file_Top = TFile::Open(file_paths[h]);
        TTree* tree_Top = (TTree*)file_Top->Get("Events");
        tree_Top->SetBranchAddress("ToMaAn_MlMu", &ToMaAn_MlMu_value);

        // double nEntries = hist->GetEntries();

        Long64_t nentries_Top = tree_Top->GetEntries();
        cout << "Number of entries for Tree---> " << fileNames[h] << ": " << nentries_Top << endl;

        TH1I * his =(TH1I*)file_Top->Get("nEles");
        double nEntries = his->GetEntries();

        //////////////////////////////////////////////////////

        // for (Long64_t j = 0; j < nentries_Top; j++) {
        //     tree_Top->GetEntry(j);
        //     // hist->GetEntry(j);
        //     hist->Fill(ToMaAn_MlMu_value);

        // }
        //////////////////////////////////////////////////////


        // vector<double> weightValue(nEntries);
        // RooRealVar * nGen = new RooRealVar("nGen","nGen",nentries_Top);
        RooRealVar * nGen = new RooRealVar("nGen","nGen",nEntries);

        fname.str("");
		fname << "data" << fileNames[h];
        RooFormulaVar *weight_Top= new RooFormulaVar(fname.str().c_str(), "(833.9*38.25*1000*4/(9*@0))", RooArgSet(*nGen)); 
        ////////////////////////////////////////////
        cout<<"________ "<<fname.str().c_str()<<" sample: "<<nGen->getVal()<<" and nEntries: "<<(int)nEntries<<",   weight should be: "<<833.9*38.25*1000*4/(9*nEntries)<<endl;
        // cout<<"________ "<<fname.str().c_str()<<" sample: "<<nGen->getVal()<<" and nentries_Top: "<<(int)nentries_Top<<",   weight should be: "<<833.9*38.25*1000*4/(9*nentries_Top)<<endl;
        weight_semilepttbar.push_back(weight_Top);

        sample.defineType(fileNames[h].c_str());


///////////////////////////////////////////////////////
        // TCanvas* c = new TCanvas(( "c" + fileNames[h]).c_str(), fileNames[h].c_str(), 800, 600);
        // hist->Draw("Hist");
        // c->SaveAs(("/home/comp-lab3/root/simFit/"+fileNames[h] + "_hist.pdf").c_str());
///////////////////////////////////////////////////////


        // RooRealVar weight(fname.str().c_str(), fname.str().c_str(), nEntries);
        // for (int i=0; i<weight_semilepttbar.size();i++){
	    datas.push_back(new RooDataSet(fname.str().c_str(), fname.str().c_str(),tree_Top, x, ""));
		datas[datas.size()-1]->Print();
		wVars.push_back((RooRealVar*)datas[datas.size()-1]->addColumn(*weight_semilepttbar[h]));
		cout<<" weight_semilepttbar:  "<<weight_semilepttbar[h]->evaluate()<<endl;
		Wdatas.push_back(new RooDataSet(fname.str().c_str(), fname.str().c_str(),datas[datas.size()-1],*datas[datas.size()-1]->get(),0, wVars[wVars.size()-1]->GetName()));
		cout<<" should be the same with: "<<wVars[wVars.size()-1]->getVal()<<endl;
		cout<<" , name:  "<<wVars[wVars.size()-1]->GetName()<<endl;
		Wdatas[Wdatas.size()-1]->Print();		
        }
		// for (unsigned int b = 0; b < fileNames.size(); b++) {
        //                	cout<<"  ... appending "<<bkgFiles[b]<<endl;
        //                	Wdatas[i]->append(*bkgdatas[b]);
        //                	Wdatas[Wdatas.size()-1]->Print();
		// }

        // RooDataSet* weightedData = new RooDataSet(fileNames[h].c_str(), ("Weighted Dataset " + fileNames[h]).c_str(),
        //                                       x, Import(*hist), WeightVar("weight"));

        //  for (Long64_t j = 0; j < nEntries; j++) {
        //     weightedData->add(*dataTopMass->get(j),/*nullptr,*/ weightValue[j]); // Add each entry with its corresponding weight
        // }
            // Wdatas.push_back(weightedData);
        

        /////////////////////////////////      
        // RooDataHist* dataTopMass = new RooDataHist(fileNames[h].c_str(), ("Dataset " + fileNames[h]).c_str(), x, Import(*hist));
        /////////////////////////////////

        // RooDataSet* dataTopMass = new RooDataSet(fileNames[h].c_str(), ("Dataset " + fileNames[h]).c_str(), x, Import(*hist));
        
        /////////////////////////////////////////   
        // dataTopMass->add(*dataTopMass, nullptr, weight_Top.getVal());
        ////////////////////////////////////////

        // datas.push_back(dataTopMass);

		// wVars.push_back((RooRealVar*)datas[datas.size()-1]->addColumn(*weight_semilepttbar[h]));
		// Wdatas.push_back(new RooDataHist(fileNames[h].c_str(), fileNames[h].c_str(),datas[datas.size()-1],*datas[datas.size()-1]->get(),0, wVars[wVars.size()-1]->GetName()));


    

// if(debug==1)cout << "debug  3 :"<<endl;

        // for (int f = 0; f < allentries.size(); f++) {
        //             weightValue[f] =833.9*38.25*1000*4/(9*allentries[f]);
        // }
// if(debug==1)cout << "debug  4 :"<<endl;

        // RooRealVar weight_M166p5("weight_M166p5", "Weight for M166p5", weightValue[0]);
        // RooRealVar weight_M169p5("weight_M169p5", "Weight for M169p5", weightValue[1]);
        // RooRealVar weight_M171p5("weight_M171p5", "Weight for M171p5", weightValue[2]);
        // RooRealVar weight_M172p5("weight_M172p5", "Weight for M172p5", weightValue[3]);
        // RooRealVar weight_M173p5("weight_M173p5", "Weight for M173p5", weightValue[4]); 
        // RooRealVar weight_M175p5("weight_M175p5", "Weight for M175p5", weightValue[5]);
        // RooRealVar weight_M178p5("weight_M178p5", "Weight for M178p5", weightValue[6]); 

        // RooDataHist* dataTopM166p5 = new RooDataHist("M166p5", "M166p5", x, Import(*allHists[0]));
        // // dataTopM166p5->add(*dataTopM166p5, nullptr, weight_M166p5.getVal());

        // RooDataHist* dataTopM169p5 = new RooDataHist("M169p5", "M169p5", x, Import(*allHists[1]));
        // // dataTopM169p5->add(*dataTopM169p5, nullptr, weight_M169p5.getVal());

        // RooDataHist* dataTopM171p5 = new RooDataHist("M171p5", "M171p5", x, Import(*allHists[2]));
        // // dataTopM171p5->add(*dataTopM171p5, nullptr, weight_M171p5.getVal());

        // RooDataHist* dataTopM172p5 = new RooDataHist("M172p5", "M172p5", x, Import(*allHists[3]));
        // // dataTopM172p5->add(*dataTopM172p5, nullptr, weight_M172p5.getVal());

        // RooDataHist* dataTopM173p5 = new RooDataHist("M173p5", "M1731p5", x, Import(*allHists[4]));
        // // dataTopM173p5->add(*dataTopM173p5, nullptr, weight_M173p5.getVal());

        // RooDataHist* dataTopM175p5 = new RooDataHist("M175p5", "M175p5", x, Import(*allHists[5]));
        // // dataTopM175p5->add(*dataTopM175p5, nullptr, weight_M175p5.getVal());

        // RooDataHist* dataTopM178p5 = new RooDataHist("M178p5", "M178p5", x, Import(*allHists[6]));
        // // dataTopM178p5->add(*dataTopM178p5, nullptr, weight_M178p5.getVal());




        // RooDataSet* weighted_data_M166p5 = new RooDataSet("weighted_data_M166p5", "Weighted Dataset M166p5",
        //                                                 x, Import(*dataTopM166p5), WeightVar("weight_M166p5"));

        // RooDataSet* weighted_data_M169p5 = new RooDataSet("weighted_data_M169p5", "Weighted Dataset M166p5",
        //                                                 x, Import(*dataTopM169p5), WeightVar("weight_M169p5"));
        
        // RooDataSet* weighted_data_M171p5 = new RooDataSet("weighted_data_M171p5", "Weighted Dataset M171p5",
        //                                                 x, Import(*dataTopM171p5), WeightVar("weight_M171p5"));
                                                                                                        
        // RooDataSet* weighted_data_M172p5 = new RooDataSet("weighted_data_M172p5", "Weighted Dataset M172p5",
        //                                                 x, Import(*dataTopM172p5), WeightVar("weight_M172p5"));

        // RooDataSet* weighted_data_M173p5 = new RooDataSet("weighted_data_M173p5", "Weighted Dataset M173p5",
        //                                                 x, Import(*dataTopM173p5), WeightVar("weight_M173p5"));

        // RooDataSet* weighted_data_M175p5 = new RooDataSet("weighted_data_M175p5", "Weighted Dataset M175p5",
        //                                                 x, Import(*dataTopM175p5), WeightVar("weight_M175p5"));

        // RooDataSet* weighted_data_M178p5 = new RooDataSet("weighted_data_M178p5", "Weighted Dataset M178p5",
        //                                            x, Import(*dataTopM178p5), WeightVar("weight_M178p5"));        




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
    
// // RooDataSet *combDataWeighted = new RooDataSet("combData", "combined data", x, Index(sample));
// RooDataSet *combDataWeighted = new RooDataSet("combData", "combined data", x, Index(sample),Import("M166p5", *weighted_data_M166p5),
// Import("M169p5", *weighted_data_M169p5),Import("M171p5", *weighted_data_M171p5),Import("M172p5", *weighted_data_M172p5),
// Import("M173p5", *weighted_data_M173p5), Import("M175p5", *weighted_data_M175p5), Import("M178p5", *weighted_data_M178p5));

// /*
//                                       Import("M166p5", *dataTopM166p5, WeightVar(weight_M166p5.GetName())),
//                                       Import("M169p5", *dataTopM169p5, WeightVar(weight_M169p5.GetName())),
//                                       Import("M171p5", *dataTopM171p5, WeightVar(weight_M171p5.GetName())),
//                                       Import("M172p5", *dataTopM172p5, WeightVar(weight_M172p5.GetName())),
//                                       Import("M173p5", *dataTopM173p5, WeightVar(weight_M173p5.GetName())),
//                                       Import("M175p5", *dataTopM175p5, WeightVar(weight_M175p5.GetName())),
//                                       Import("M178p5", *dataTopM178p5, WeightVar(weight_M178p5.GetName()))
//                                       );
// */
        // combDataWeighted.append(weighted_data_M166p5);
        // combDataWeighted.append(weighted_data_M169p5);
        // combDataWeighted.append(weighted_data_M171p5);
        // combDataWeighted.append(weighted_data_M172p5);
        // combDataWeighted.append(weighted_data_M173p5);
        // combDataWeighted.append(weighted_data_M175p5);
        // combDataWeighted.append(weighted_data_M178p5);

        RooRealVar* alpha = new RooRealVar("alpha", "fraction of gaussian", 0, 1); 

    for (size_t w = 0; w < fileNames.size(); w++) {

        // RooRealVar* alpha = new RooRealVar(("alpha_" + fileNames[w]).c_str(), "fraction of gaussian", 0, 1); 
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

        RooAddPdf* model = new RooAddPdf(fileNames[w].c_str(), "Gaussian + Gamma PDF",
                                        RooArgList(*gauss, *gamma), 
                                        *alpha); 

        models.push_back(model);

    }



        RooSimultaneous simPdf("simPdf","simultaneous pdf", sample);

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

