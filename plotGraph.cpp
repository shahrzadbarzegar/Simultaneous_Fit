#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>

void plotGraph() {
    const int n = 7;


    string fileNames[n] = {"M166p5", "M169p5", "M171p5", "M172p5", "M173p5", "M175p5", "M178p5"};
    double mass_top[n] = {166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5};
    
    double alpha_GaussFrac[n] ={ 0.725065	,0.726588, 0.730813	 , 0.729184	 ,0.735256	,0.733477, 0.738999};	  	
    double alphaErrors[n] = {0.000291574, 0.000260022, 0.000254181, 0.000270751, 0.000228737, 0.000202304, 0.000263292};

    TGraphErrors *gr = new TGraphErrors(n, mass_top, alpha_GaussFrac, nullptr, alphaErrors);

    TCanvas *c1 = new TCanvas("c1", "alpha_GaussFrac and Top Mass", 800, 600);
    
    gr->Draw("AP");

    gr->SetMarkerStyle(21); 
    gr->SetMarkerSize(1);   

    gr->SetLineColor(kBlue);
    
    TF1 *fitFunction = new TF1("fitFunction", "pol1", mass_top[0], mass_top[n-1]);
    
    gr->Fit(fitFunction);
    fitFunction->SetLineColor(kRed); 
    
    c1->Update();
    c1->SaveAs("/home/comp-lab3/root/simFit/alpha_GaussFrac.pdf");
//scale parameter of Gamma function
    double gamma_beta[n] = {30,30, 30,30,30,30,30};
    double gamma_beta_errors[n] = {0.0544049, 0.047422, 0.0498006, 0.0420236, 0.0427316, 0.0441441, 0.0441572};

    TGraphErrors *g_gamma_beta = new TGraphErrors(n, mass_top, gamma_beta, nullptr, gamma_beta_errors);

    TCanvas *c_gamma_beta = new TCanvas("cgamma_beta", "gamma_beta and Top Mass", 800, 600);
    
    g_gamma_beta->Draw("AP");

    g_gamma_beta->SetMarkerStyle(21); 
    g_gamma_beta->SetMarkerSize(1);   

    g_gamma_beta->SetLineColor(kBlue);
    
    // TF1 *fitFunction_gamma_beta = new TF1("fitFunction", "pol1", mass_top[0], mass_top[n-1]);
    
    g_gamma_beta->Fit(fitFunction);
    fitFunction->SetLineColor(kBlue); 
    
    c_gamma_beta->Update();
    c_gamma_beta->SaveAs("/home/comp-lab3/root/simFit/gamma_beta.pdf");

//shape parameter

    double gamma_gamma[n] = {2.44644, 2.48229, 2.49058, 2.51321, 2.52707, 2.54333, 2.57203};
    double gamma_gamma_errors[n] = {0.0178804, 0.017378, 0.0174462, 0.017519, 0.0174969, 0.0169983, 0.0171131};

    TGraphErrors *g_gamma_gamma = new TGraphErrors(n, mass_top, gamma_gamma, nullptr, gamma_gamma_errors);

    TCanvas *c_gamma_gamma = new TCanvas("cgamma_gamma", "gamma_gamma and Top Mass", 800, 600);
    
    g_gamma_gamma->Draw("AP");

    g_gamma_gamma->SetMarkerStyle(21); 
    g_gamma_gamma->SetMarkerSize(1);   

    g_gamma_gamma->SetLineColor(kBlue);
    
    // TF1 *fitFunction_gamma_beta = new TF1("fitFunction", "pol1", mass_top[0], mass_top[n-1]);
    
    g_gamma_gamma->Fit(fitFunction);
    fitFunction->SetLineColor(kBlue); 
    
    c_gamma_gamma->Update();
    c_gamma_gamma->SaveAs("/home/comp-lab3/root/simFit/gamma_gamma.pdf");

    double mean_gauss[n] = { 59.5717 , 60.7702, 61.6627, 61.974, 62.5632, 63.1756 , 64.3926 };
    double mean_gauss_errors[n] = {0.216097, 0.218871, 0.219972, 0.216401, 0.217587, 0.218671, 0.217722};
    TGraphErrors *g_mean_gauss = new TGraphErrors(n, mass_top, mean_gauss, nullptr, mean_gauss_errors);

    TCanvas *c_mean_gauss = new TCanvas("cmean_gauss", "mean_gauss and Top Mass", 800, 600);
    
    g_mean_gauss->Draw("AP");

    g_mean_gauss->SetMarkerStyle(21); 
    g_mean_gauss->SetMarkerSize(1);   

    g_mean_gauss->SetLineColor(kBlue);
    
    // TF1 *fitFunction_gamma_beta = new TF1("fitFunction", "pol1", mass_top[0], mass_top[n-1]);
    
    g_mean_gauss->Fit(fitFunction);
    fitFunction->SetLineColor(kBlue); 
    
    c_mean_gauss->Update();
    c1->SaveAs("/home/comp-lab3/root/simFit/gamma_mean.pdf");

//shift parameter 
    double mu_gamma[n] = {8.33757, 8.40653, 8.45317	, 8.29288, 8.16278 , 8.33167,  8.252};
    double mu_gaamma_errors[n] = {0.28403, 0.27647, 0.273213, 0.272109,  0.27658, 0.260442, 0.268632};
    TGraphErrors *g_mu_gamma = new TGraphErrors(n, mass_top, mu_gamma, nullptr, mu_gaamma_errors);

    TCanvas *c_mu_gamma = new TCanvas("cmu_gamma", "mu_gamma and Top Mass", 800, 600);
    
    g_mu_gamma->Draw("AP");

    g_mu_gamma->SetMarkerStyle(21); 
    g_mu_gamma->SetMarkerSize(1);   

    g_mu_gamma->SetLineColor(kBlue);
    
    // TF1 *fitFunction_gamma_beta = new TF1("fitFunction", "pol1", mass_top[0], mass_top[n-1]);
    
    g_mu_gamma->Fit(fitFunction);
    fitFunction->SetLineColor(kBlue); 
    
    c_mu_gamma->Update();
    c_mu_gamma->SaveAs("/home/comp-lab3/root/simFit/mu_gamma.pdf");




    double sigma_gauss[n]={16.5561, 16.9337, 17.1608, 17.0376, 17.2309, 17.414, 17.6062};
    double sigma_gauss_errors[n] = {0.163794, 0.165361 ,0.165082, 0.1639 , 0.16413, 0.164414, 0.165027 };
    TGraphErrors *g_sigma_gauss = new TGraphErrors(n, mass_top, sigma_gauss, nullptr, sigma_gauss_errors);

    TCanvas *c_sigma_gauss = new TCanvas("csigma_gauss", "sigma_gauss and Top Mass", 800, 600);
    
    g_sigma_gauss->Draw("AP");

    g_sigma_gauss->SetMarkerStyle(21); 
    g_sigma_gauss->SetMarkerSize(1);   

    g_sigma_gauss->SetLineColor(kBlue);
    
    // TF1 *fitFunction_gamma_beta = new TF1("fitFunction", "pol1", mass_top[0], mass_top[n-1]);
    
    g_sigma_gauss->Fit(fitFunction);
    fitFunction->SetLineColor(kBlue); 
    
    c_sigma_gauss->Update();
    c_sigma_gauss->SaveAs("/home/comp-lab3/root/simFit/sigma_gauss.pdf");
}
/*
root [0] 
Processing plotGraph.cpp...
****************************************
Minimizer is Linear / Migrad
Chi2                      =      4.50175
NDf                       =            5
p0                        =     -11.6872   +/-   1.92813     
p1                        =     0.432722   +/-   0.011161    
root [1] .q

*/