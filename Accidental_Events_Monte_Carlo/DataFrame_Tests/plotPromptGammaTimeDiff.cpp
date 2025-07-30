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
#include "../../DTW_algorithm/calc_lifetime.h"

using DTW_func = void (*)(TH1F*, RNode, int);

//Helper funcs
static void drawHistogram(TH1F* hDeltaT, std::string histName, size_t eventCount){
	gStyle->SetOptStat(1110);
    TCanvas* c = new TCanvas("c", "DeltaT Histogram", 800, 600);
    hDeltaT->Draw();
	hDeltaT->SetMinimum(0);
	hDeltaT->SetMaximum(hDeltaT->GetMaximum() * 1.1);

	// Dodanie info do stat box'a
	c->Update();
	TPaveStats* stats = (TPaveStats*)hDeltaT->FindObject("stats");
	TPaveStats* newStats = (TPaveStats*)stats->Clone("newStats");
	TText* t = newStats->AddText(Form("Event count: %zu", eventCount));
	stats->SetBit(kCanDelete);
    stats->Delete();
    newStats->Draw("same");
	c->Modified();
	c->Update();

    c->SaveAs(histName.c_str());
    
    delete c;
}

//Second 511 delayed - random is prompt and second 511
static void fillHistForDTW_1(TH1F* hDeltaT, RNode df, int skips){
    std::vector<std::vector<float>> times;
    DTW_type1(df, skips, &times);

    size_t nEvents = times.size();
    for (size_t i = 0; i < nEvents; ++i) {
		auto promptTime = times[i][0];
		float deltaT = times[i][2] - promptTime;
		hDeltaT->Fill(deltaT);
    }    
    
    drawHistogram(hDeltaT, "deltaT_Rand_DTW1_hist.png", nEvents);
}

//Second and prompt delayed - random is prompt and first 511
static void fillHistForDTW_2(TH1F* hDeltaT, RNode df, int skips){
    std::vector<std::vector<float>> times;
    DTW_type2(df, skips, &times);

    size_t nEvents = times.size();
    for (size_t i = 0; i < nEvents; ++i) {
		auto promptTime = times[i][0];
		float deltaT = times[i][1] - promptTime;
		hDeltaT->Fill(deltaT);
    }
    
    drawHistogram(hDeltaT, "deltaT_Rand_DTW2_hist.png", nEvents);
}

//Second and prompt delayed (diff shifts) - all randoms
static void fillHistForDTW_3(TH1F* hDeltaT, RNode df, int skips){
    std::vector<std::vector<float>> times;
    DTW_type3(df, skips, &times);

    size_t nEvents = times.size();
    for (size_t i = 0; i < nEvents; ++i) {
		auto promptTime = times[i][0];
		for (int j = 1; j < 3; j++){
			float deltaT = times[i][j] - promptTime;
			hDeltaT->Fill(deltaT);
		}
    }    
    
    drawHistogram(hDeltaT, "deltaT_Rand_DTW3_hist.png", nEvents);
}

//Prompr delayed - all randoms
static void fillHistForDTW_4(TH1F* hDeltaT, RNode df, int skips){
    std::vector<std::vector<float>> times;
    DTW_type4(df, skips, &times);

    size_t nEvents = times.size();
    for (size_t i = 0; i < nEvents; ++i) {
		auto promptTime = times[i][0];
		for (int j = 1; j < 3; j++){
			float deltaT = times[i][j] - promptTime;
			hDeltaT->Fill(deltaT);
		}
    }    
    
    drawHistogram(hDeltaT, "deltaT_Rand_DTW4_hist.png", nEvents);
}

//Plotting funcs
static void plotTrueIntervalsFromDF(int window_count, double activity){
	std::string fileName = generate_DataFrame(window_count, activity);

    std::string treeName = "FlatTree";	
    TChain chain(treeName.c_str());
    chain.Add(fileName.c_str());

    ROOT::RDataFrame df(chain);

    // Histogram dla różnic czasów prompt - gamma
    TH1F* hDeltaT = new TH1F("hDeltaT", "Prompt-gamma time diff", 50, -400, 1000);

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
    
	drawHistogram(hDeltaT, "deltaT_True_hist.png", nEvents);
    delete hDeltaT;
}

static void plotRandomsIntervalsFromDTW(int window_count, double activity, DTW_func fillHistForDTW){
	std::string fileName = generate_DataFrame(window_count, activity);

    std::string treeName = "FlatTree";	
    TChain chain(treeName.c_str());
    chain.Add(fileName.c_str());

    ROOT::RDataFrame df(chain);

    // Histogram dla różnic czasów prompt - gamma
    TH1F* hDeltaT = new TH1F("hDeltaT", "Prompt-gamma time diff", 50, -2500, 2500);
    fillHistForDTW(hDeltaT, df, 2);
    delete hDeltaT;
}

static void plotRandomsIntervalsFromGen(int window_count, double activity){
	auto hitsVec = generateEvents(window_count, activity);
	std::vector<std::vector<gammaP>> pairs;
	
	int nEvents = findRandoms(hitsVec, &pairs);
	
	// Histogram dla różnic czasów prompt - gamma
    TH1F* hDeltaT = new TH1F("hDeltaT", "Prompt-gamma time diff", 50, -2500, 2500);
    
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
		    if ( event[gammaIdx].getEventNum() != event[promptIdx].getEventNum() ){
		    	float deltaT = event[gammaIdx].getTime() - promptTime;
				hDeltaT->Fill(deltaT);
			}
		}
    }
    
	drawHistogram(hDeltaT, "deltaT_Rand_Gen_hist.png", hitsVec.size()/3);
    delete hDeltaT;
}

static void plotTruesIntervalsFromGen(int window_count, double activity){
	auto hitsVec = generateEvents(window_count, activity);
	std::vector<std::vector<gammaP>> pairs;
	
	int nEvents = findTrues(hitsVec, &pairs);
	
	// Histogram dla różnic czasów prompt - gamma
    TH1F* hDeltaT = new TH1F("hDeltaT", "Prompt-gamma time diff", 50, -2500, 2500);
    
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
    
    drawHistogram(hDeltaT, "deltaT_Trues_Gen_hist.png", hitsVec.size()/3);
    delete hDeltaT;
}


static void plotTruesAndRandoms(int window_count, double activity){
	auto hitsVec = generateEvents(window_count, activity);
	std::vector<std::vector<gammaP>> pairs;
	
	int nEvents = findCoincidences(hitsVec, &pairs);
	int nEvents_True = 0, nEvents_Rand = 0;
	
	// Histogramy
    TH1F* hDeltaT_True = new TH1F("hDeltaT_True", "Prompt-gamma time diff", 50, -2500, 2500);
    TH1F* hDeltaT_Rand = new TH1F("hDeltaT_Rand", "Prompt-gamma time diff", 50, -2500, 2500);
    
    // Wypelnienie histogramow
    for (int i = 0; i < nEvents; ++i) {
    	const auto& event = pairs[i];
    	
    	int promptIdx = -1;
        std::vector<int> gammaIdxs;

		if (event.size() != 3) continue;
		
		//Czy event to random?
		bool isRandom = true;
		if( (event[0].getEventNum() == event[1].getEventNum()) && 
			(event[1].getEventNum() == event[2].getEventNum()) ) {
			
			isRandom = false;
		}

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
        	
        	if(isRandom){
        		if ( event[gammaIdx].getEventNum() != event[promptIdx].getEventNum() ){
					hDeltaT_Rand->Fill(deltaT);
					nEvents_Rand++;
				}
        	} else {
        		hDeltaT_True->Fill(deltaT);
        		nEvents_True++;
			}
		}
    }
    
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c", "DeltaT Comparison", 800, 600);
    
    // Trues hist
    hDeltaT_True->SetLineColor(kBlue);
    hDeltaT_True->SetLineWidth(2);
    hDeltaT_True->SetMinimum(0);
    hDeltaT_True->SetMaximum(std::max(hDeltaT_True->GetMaximum(), hDeltaT_Rand->GetMaximum()) * 1.1);
    
    hDeltaT_True->Draw();
    
    // Rand hist
    hDeltaT_Rand->SetLineColor(kRed);
    hDeltaT_Rand->SetLineWidth(2);
    hDeltaT_Rand->Draw("SAME");
    
    // Legenda
    TLegend* legend = new TLegend(0.65, 0.75, 0.88, 0.88);
    legend->AddEntry(hDeltaT_True, Form("True coinc. (n = %d)", nEvents_True), "l");
    legend->AddEntry(hDeltaT_Rand, Form("Random coinc. (n = %d)", nEvents_Rand), "l");
    legend->Draw();
    
    c->SaveAs("deltaT_True_and_Rand_hist.png");
    delete c;
    delete hDeltaT_True;
    delete hDeltaT_Rand;
}

//Main
static void plotPromptGammaTimeDiff(int window_count, double activity) {
    plotTrueIntervalsFromDF(window_count, activity);
    plotRandomsIntervalsFromGen(window_count, activity);
    plotTruesIntervalsFromGen(window_count, activity);
    plotRandomsIntervalsFromDTW(window_count, activity, fillHistForDTW_1);
    plotRandomsIntervalsFromDTW(window_count, activity, fillHistForDTW_2);
    plotRandomsIntervalsFromDTW(window_count, activity, fillHistForDTW_3);
    plotRandomsIntervalsFromDTW(window_count, activity, fillHistForDTW_4);
    plotTruesAndRandoms(window_count, activity);
}

