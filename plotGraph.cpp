#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <string>
#include <vector>

void plotGraph() {
    const int n = 6;
    const int numGraphs = 5;

    // Names of files and variables
    std::string fileNames[n] = {"M166p5", "M169p5", "M171p5", "M173p5", "M175p5", "M178p5"};
    double mass_top[n] = {166.5, 169.5, 171.5, 173.5, 175.5, 178.5};

    // Data arrays for each parameter and their errors
    double data[numGraphs][n] = {
        {30, 30, 30, 30, 30, 30},                                             // gamma_beta
        {2.44805, 2.48395, 2.49202, 2.52864, 2.54485, 2.5736},          // gamma_gamma
        {59.5612, 60.7595, 61.6529, 62.5521, 63.1634, 64.3832},           // mean_gauss
        {8.33001, 8.39546, 8.44387, 8.15394, 8.32168, 8.24161},            // mu_gamma
        {16.5713, 16.9489, 17.1758, 17.2469, 17.43, 17.6222}            // sigma_gauss
    };

    double errors[numGraphs][n] = {
        {0.0543043, 0.0473491, 0.0492577, 0.0425215, 0.0438611, 0.0446834},                // gamma_beta_errors
        {0.0178725, 0.0173959, 0.0175328, 0.0176342, 0.0171388, 0.0172597},                 // gamma_gamma_errors
        {0.216601, 0.219246, 0.220383, 0.218051, 0.219187, 0.218247},                     // mean_gauss_errors
        {0.281536, 0.275236, 0.274313, 0.278992, 0.261227, 0.269804},                        // mu_gaamma_errors
        {0.164747, 0.16648, 0.166316, 0.165473, 0.165613, 0.166598}                         // sigma_gauss_errors
    };

    std::string graphNames[numGraphs] = { "scale parameter", "shape parameter", "mean_gauss", "shift parameter", "sigma_gauss"};

    // Loop through each graph to create and plot
    for (int i = 0; i < numGraphs; i++) {
        TGraphErrors *graph = new TGraphErrors(n, mass_top, data[i], nullptr, errors[i]);
        TCanvas *canvas = new TCanvas(("c_" + graphNames[i]).c_str(), (graphNames[i] + " and Top Mass").c_str(), 800, 600);

        // Draw graph
        graph->Draw("AP");
        graph->SetMarkerStyle(21);
        graph->SetMarkerSize(1);
        graph->SetLineColor(kBlue);
        graph->SetTitle((graphNames[i]).c_str());
        graph->GetXaxis()->SetTitle("Top Quark Mass");
	    graph->GetYaxis()->SetTitle((graphNames[i]).c_str());
        // Fit function
        TF1 *fitFunction = new TF1("fitFunction", "pol1", mass_top[0], mass_top[n-1]);
        graph->Fit(fitFunction);
        fitFunction->SetLineColor(kRed);

        // Update and save canvas
        canvas->Update();
        canvas->SaveAs(("/home/comp-lab3/root/simFit/" + graphNames[i] + ".pdf").c_str());
    }
}
