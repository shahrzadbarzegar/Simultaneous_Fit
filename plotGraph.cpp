#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <string>
#include <vector>

void plotGraph() {
    const int n = 7;
    const int numGraphs = 6;

    std::string fileNames[n] = {"M166p5", "M169p5", "M171p5", "M172p5", "M173p5", "M175p5", "M178p5"};
    double mass_top[n] = {166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5};

    double data[numGraphs][n] = {
        {0.725065, 0.726588, 0.730813, 0.729184, 0.735256, 0.733477, 0.738999},   // alpha_GaussFrac
        {30, 30, 30, 30, 30, 30, 30},                                             // gamma_beta
        {2.44644, 2.48229, 2.49058, 2.51321, 2.52707, 2.54333, 2.57203},          // gamma_gamma
        {59.5717, 60.7702, 61.6627, 61.974, 62.5632, 63.1756, 64.3926},           // mean_gauss
        {8.33757, 8.40653, 8.45317, 8.29288, 8.16278, 8.33167, 8.252},            // mu_gamma
        {16.5561, 16.9337, 17.1608, 17.0376, 17.2309, 17.414, 17.6062}            // sigma_gauss
    };
    
    double errors[numGraphs][n] = {
        {0.000291574, 0.000260022, 0.000254181, 0.000270751, 0.000228737, 0.000202304, 0.000263292}, // alphaErrors
        {0.0544049, 0.047422, 0.0498006, 0.0420236, 0.0427316, 0.0441441, 0.0441572},                // gamma_beta_errors
        {0.0178804, 0.017378, 0.0174462, 0.017519, 0.0174969, 0.0169983, 0.0171131},                 // gamma_gamma_errors
        {0.216097, 0.218871, 0.219972, 0.216401, 0.217587, 0.218671, 0.217722},                     // mean_gauss_errors
        {0.28403, 0.27647, 0.273213, 0.272109, 0.27658, 0.260442, 0.268632},                        // mu_gaamma_errors
        {0.163794, 0.165361, 0.165082, 0.1639, 0.16413, 0.164414, 0.165027}                         // sigma_gauss_errors
    };

    std::string graphNames[numGraphs] = {"alpha_GaussFrac", "gamma_beta", "gamma_gamma", "mean_gauss", "mu_gamma", "sigma_gauss"};

    for (int i = 0; i < numGraphs; i++) {
        TGraphErrors *graph = new TGraphErrors(n, mass_top, data[i], nullptr, errors[i]);
        TCanvas *canvas = new TCanvas(("c_" + graphNames[i]).c_str(), (graphNames[i] + " and Top Mass").c_str(), 800, 600);

        graph->Draw("AP");
        graph->SetMarkerStyle(21);
        graph->SetMarkerSize(1);
        graph->SetLineColor(kBlue);

        TF1 *fitFunction = new TF1("fitFunction", "pol1", mass_top[0], mass_top[n-1]);
        graph->Fit(fitFunction);
        fitFunction->SetLineColor(kRed);
        canvas->Update();
        canvas->SaveAs(("/home/comp-lab3/root/simFit/" + graphNames[i] + ".pdf").c_str());
    }
}
