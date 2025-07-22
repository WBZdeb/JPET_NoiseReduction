#include <ROOT/RDataFrame.hxx>
#include <TROOT.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <string>
#include "../testDTW.cpp"

static void plotTrueIntervalsFromDF(int window_count, double activity){
	std::string fileName = generate_DataFrame(window_count, activity);

    std::string treeName = "FlatTree";	
    TChain chain(treeName.c_str());
    chain.Add(fileName.c_str());

    ROOT::RDataFrame df(chain);

    // Histogram dla różnic czasów prompt - gamma
    TH1F* hDeltaT = new TH1F("hDeltaT", "Prompt-gamma time diff", 50, -1000, 1000);

    // Wypelnienie histogramu
    auto times = df.Take<std::vector<float>>("time");
    auto energies = df.Take<std::vector<float>>("energy");

    size_t nEvents = times->size();
    for (size_t i = 0; i < nEvents; ++i) {
        const auto& tVec = (*times)[i];
        const auto& eVec = (*energies)[i];

        if (tVec.size() != 3 || eVec.size() != 3) continue;

        int promptIdx = -1;
        std::vector<int> gammaIdxs;

        for (size_t j = 0; j < eVec.size(); ++j) {
            if (eVec[j] > 511.0) {
                promptIdx = j;
            } else {
                gammaIdxs.push_back(j);
            }
        }

        if (promptIdx >= 0 && gammaIdxs.size() == 2) {
            float promptTime = tVec[promptIdx];
            for (int gammaIdx : gammaIdxs) {
                float deltaT = tVec[gammaIdx] - promptTime;
                hDeltaT->Fill(deltaT);
            }
        }
    }

    // Rysowanie histogramu
    gStyle->SetOptStat(1110);
    TCanvas* c = new TCanvas("c", "DeltaT Histogram", 800, 600);
    hDeltaT->Draw();

	// Dodanie info do stat box'a
	c->Update();
	TPaveStats* stats = (TPaveStats*)hDeltaT->FindObject("stats");
	TPaveStats* newStats = (TPaveStats*)stats->Clone("newStats");
	TText* t = newStats->AddText(Form("Event count: %zu", nEvents));
	stats->SetBit(kCanDelete);
    stats->Delete();
    newStats->Draw("same");
	c->Modified();
	c->Update();

    c->SaveAs("deltaT_True_hist.png");
    
    delete hDeltaT;
    delete c;
}

static void plotRandomsIntervalsFromGen(int window_count, double activity){
	auto hitsVec = generateEvents(window_count, activity);
	std::vector<std::vector<gammaP>> pairs;
	
	int nEvents = findRandoms(hitsVec, &pairs);
	
	// Histogram dla różnic czasów prompt - gamma
    TH1F* hDeltaT = new TH1F("hDeltaT", "Prompt-gamma time diff", 50, -1000, 1000);
    
    // Wypelnienie histogramu
    for (int i = 0; i < nEvents; ++i) {
    	const auto& event = pairs[i];
    	
    	int promptIdx = -1;
        std::vector<int> gammaIdxs;

		if (event.size() != 3) continue;

        for (size_t j = 0; j < event.size(); ++j) {
            if (event[j].isPrompt()) {
                promptIdx = j;
            } else {
                gammaIdxs.push_back(j);
            }
        }

        float promptTime = event[promptIdx].getTime();
        for (int gammaIdx : gammaIdxs) {
        	float deltaT = event[gammaIdx].getTime() - promptTime;
			hDeltaT->Fill(deltaT);
		}
    }
    
    // Rysowanie histogramu
    gStyle->SetOptStat(1110);
    TCanvas* c = new TCanvas("c", "DeltaT Histogram", 800, 600);
    hDeltaT->Draw();

	// Dodanie info do stat box'a
	c->Update();
	TPaveStats* stats = (TPaveStats*)hDeltaT->FindObject("stats");
	TPaveStats* newStats = (TPaveStats*)stats->Clone("newStats");
	TText* t = newStats->AddText(Form("Event count: %zu", hitsVec.size()/3));
	stats->SetBit(kCanDelete);
    stats->Delete();
    newStats->Draw("same");
	c->Modified();
	c->Update();

    c->SaveAs("deltaT_Rand_Gen_hist.png");
    
    delete hDeltaT;
    delete c;
}

static void plotTruesIntervalsFromGen(int window_count, double activity){
	auto hitsVec = generateEvents(window_count, activity);
	std::vector<std::vector<gammaP>> pairs;
	
	int nEvents = findTrues(hitsVec, &pairs);
	
	// Histogram dla różnic czasów prompt - gamma
    TH1F* hDeltaT = new TH1F("hDeltaT", "Prompt-gamma time diff", 50, -1000, 1000);
    
    // Wypelnienie histogramu
    for (int i = 0; i < nEvents; ++i) {
    	const auto& event = pairs[i];
    	
    	int promptIdx = -1;
        std::vector<int> gammaIdxs;

		if (event.size() != 3) continue;

        for (size_t j = 0; j < event.size(); ++j) {
            if (event[j].isPrompt()) {
                promptIdx = j;
            } else {
                gammaIdxs.push_back(j);
            }
        }

        float promptTime = event[promptIdx].getTime();
        for (int gammaIdx : gammaIdxs) {
        	float deltaT = event[gammaIdx].getTime() - promptTime;
			hDeltaT->Fill(deltaT);
		}
    }
    
    // Rysowanie histogramu
    gStyle->SetOptStat(1110);
    TCanvas* c = new TCanvas("c", "DeltaT Histogram", 800, 600);
    hDeltaT->Draw();

	// Dodanie info do stat box'a
	c->Update();
	TPaveStats* stats = (TPaveStats*)hDeltaT->FindObject("stats");
	TPaveStats* newStats = (TPaveStats*)stats->Clone("newStats");
	TText* t = newStats->AddText(Form("Event count: %zu", hitsVec.size()/3));
	stats->SetBit(kCanDelete);
    stats->Delete();
    newStats->Draw("same");
	c->Modified();
	c->Update();

    c->SaveAs("deltaT_Trues_Gen_hist.png");
    
    delete hDeltaT;
    delete c;
}

static void plotPromptGammaTimeDiff(int window_count, double activity) {
    plotTrueIntervalsFromDF(window_count, activity);
    plotRandomsIntervalsFromGen(window_count, activity);
    plotTruesIntervalsFromGen(window_count, activity);
}

