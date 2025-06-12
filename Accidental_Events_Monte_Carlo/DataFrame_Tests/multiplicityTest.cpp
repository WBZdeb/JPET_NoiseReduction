#include <ROOT/RDataFrame.hxx>
#include <TROOT.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLine.h>
#include <TStyle.h>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include "../testDTW.cpp"


static void multiplicityTest(int window_count, double activity) {
    std::string fileName = generate_DataFrame(window_count, activity);

    std::string treeName = "FlatTree";	
    TChain chain(treeName.c_str());
    chain.Add(fileName.c_str());

    ROOT::RDataFrame df(chain);
	
    // Zliczamy wystąpienia poszczególnych timeWindowNumber
    std::map<int, int> windowCounts;

    auto timeWindowNumbers = df.Take<int>("timeWindowNumber");
    for (auto num : *timeWindowNumbers) {
        windowCounts[num]++;
    }

    // Histogram
    int minWin = windowCounts.begin()->first;
    int maxWin = windowCounts.rbegin()->first;

    int nbins = maxWin - minWin + 1;
    TH1F* h = new TH1F("hTimeWindows", "Histogram krotnosci", nbins, minWin - 0.5, maxWin + 0.5);

    double sum = 0;
    for (const auto& [windowNum, count] : windowCounts) {
        h->SetBinContent(h->FindBin(windowNum), count);
        sum += count;
    }

    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c", "Histogram", 800, 600);
    h->SetFillColor(kAzure + 1);
    h->Draw();

    // Wyrysowanie średniej
    double mean = sum / windowCounts.size();
    TLine* line = new TLine(minWin - 0.5, mean, maxWin + 0.5, mean);
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->SetLineStyle(2);
    line->Draw("same");

	c->SaveAs("multiplicity_hist.png");
    c->Draw();
}
