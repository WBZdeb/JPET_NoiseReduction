#include <ROOT/RDataFrame.hxx>
#include <TROOT.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TH1F.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <string>
#include "../testDTW.cpp"
#include "../../DTW_algorithm/calc_lifetime.h"

using DTW_func = void (*)(TH1F*, RNode, int);

enum class RandomType {
    I = 0,
    IIa = 1,
    IIb = 2,
    III = 3
};

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
	
	//Dopasowanie osi X

	int firstBin = hDeltaT->FindFirstBinAbove(0);
	int lastBin  = hDeltaT->FindLastBinAbove(0);

	double minX = hDeltaT->GetBinLowEdge(firstBin);
	double maxX = hDeltaT->GetBinLowEdge(lastBin + 1);
	double margin = (maxX - minX) * 0.1;
	
	hDeltaT->GetXaxis()->SetRangeUser(minX - margin, maxX + margin);

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
    
    TH1F* hDeltaT_ref = new TH1F("hDeltaT_ref", "Prompt-gamma time diff (Reference only)", 100, -2000, 3000);
    TH1F* hDeltaT_nonRef = new TH1F("hDeltaT_nonRef", "Prompt-gamma time diff (Non-reference only)", 100, -2000, 3000);

    size_t nEvents = times.size();
    for (size_t i = 0; i < nEvents; ++i) {
		auto promptTime = times[i][0];
		for (int j = 1; j < 3; j++){
			float deltaT = times[i][j] - promptTime;
			hDeltaT->Fill(deltaT);
			
			if(j == 1) {
				hDeltaT_ref->Fill(deltaT);
			} else {
				hDeltaT_nonRef->Fill(deltaT);
			}
		}
    }    
    
    drawHistogram(hDeltaT, "deltaT_Rand_DTW3_hist.png", nEvents);
    drawHistogram(hDeltaT_ref, "deltaT_Rand_ref_DTW3_hist.png", nEvents);
    drawHistogram(hDeltaT_nonRef, "deltaT_Rand_nonRef_DTW3_hist.png", nEvents);
    
    delete hDeltaT_ref;
    delete hDeltaT_nonRef;
}

//Prompr delayed - all randoms
static void fillHistForDTW_4(TH1F* hDeltaT, RNode df, int skips){
    std::vector<std::vector<float>> times;
    DTW_type4(df, skips, &times);
    
    TH1F* hDeltaT_ref = new TH1F("hDeltaT_ref", "Prompt-gamma time diff (Reference only)", 100, -2000, 3000);
    TH1F* hDeltaT_nonRef = new TH1F("hDeltaT_nonRef", "Prompt-gamma time diff (Non-reference only)", 100, -2000, 3000);

    size_t nEvents = times.size();
    for (size_t i = 0; i < nEvents; ++i) {
		auto promptTime = times[i][0];
		for (int j = 1; j < 3; j++){
			float deltaT = times[i][j] - promptTime;
			hDeltaT->Fill(deltaT);
			
			if(j == 1) {
				hDeltaT_ref->Fill(deltaT);
			} else {
				hDeltaT_nonRef->Fill(deltaT);
			}
		}
    }    
    
    drawHistogram(hDeltaT, "deltaT_Rand_DTW4_hist.png", nEvents);
    drawHistogram(hDeltaT_ref, "deltaT_Rand_ref_DTW4_hist.png", nEvents);
    drawHistogram(hDeltaT_nonRef, "deltaT_Rand_nonRef_DTW4_hist.png", nEvents);
    
    delete hDeltaT_ref;
    delete hDeltaT_nonRef;
}

//Assign random type, given event = {prompt, gamma_1, gamma_2}
static int assignRandomType(const std::vector<gammaP>& event){
	RandomType randomType = RandomType::III;
	if( event[1].getEventNum() == event[2].getEventNum() ) randomType = RandomType::I;
	if( event[1].getEventNum() == event[0].getEventNum() ) randomType = RandomType::IIa;
	if( event[0].getEventNum() == event[2].getEventNum() ) randomType = RandomType::IIb;
	
	return static_cast<int>(randomType);
}

//Plotting funcs
static void plotTrueIntervalsFromDF(int window_count, double activity){
	std::string fileName = generate_DataFrame(window_count, activity);

    std::string treeName = "FlatTree";	
    TChain chain(treeName.c_str());
    chain.Add(fileName.c_str());

    ROOT::RDataFrame df(chain);

    // Histogram dla różnic czasów prompt - gamma
    TH1F* hDeltaT = new TH1F("hDeltaT", "Prompt-gamma time diff", 50, -1000, 1500);

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
    TH1F* hDeltaT = new TH1F("hDeltaT", "Prompt-gamma time diff", 100, -2000, 3000);
    fillHistForDTW(hDeltaT, df, 2);
    delete hDeltaT;
}

static void plotRandomsIntervalsFromGen(int window_count, double activity){
	auto hitsVec = generateEvents(window_count, activity);
	std::vector<std::vector<gammaP>> pairs;
	
	int nEvents = findRandoms(hitsVec, &pairs);
	
	// Histogramy dla różnic czasów prompt - gamma
    TH1F* hDeltaT = new TH1F("hDeltaT", "Prompt-gamma time diff", 100, -2000, 3000);
    TH1F* hDeltaT_ref = new TH1F("hDeltaT_ref", "Prompt-gamma time diff (Reference only)", 100, -2000, 3000);
    TH1F* hDeltaT_nonRef = new TH1F("hDeltaT_nonRef", "Prompt-gamma time diff (Non-reference only)", 100, -2000, 3000);
    
    TH1F* hDeltaT_typeI = new TH1F("hDeltaT_typeI", "Prompt-gamma time diff (type I)", 100, -2000, 3000);
    TH1F* hDeltaT_typeIIa = new TH1F("hDeltaT_typeIIa", "Prompt-gamma time diff (type IIa)", 100, -2000, 3000);
    TH1F* hDeltaT_typeIIb = new TH1F("hDeltaT_typeIIb", "Prompt-gamma time diff (type IIb)", 100, -2000, 3000);
    TH1F* hDeltaT_typeIII = new TH1F("hDeltaT_typeIII", "Prompt-gamma time diff (type III)", 100, -2000, 3000);
    
    // Wypelnienie histogramu
    for (int i = 0; i < nEvents; ++i) {
    	const auto& event = pairs[i];

		if (event.size() != 3) continue;

        float promptTime = event[2].getTime();
        for (int j = 0; j < 2; j++) {
		    if ( event[j].getEventNum() != event[2].getEventNum() ){
		    	float deltaT = event[j].getTime() - promptTime;
				hDeltaT->Fill(deltaT);
				
				if(j == 0) hDeltaT_ref->Fill(deltaT);
				if(j == 1) hDeltaT_nonRef->Fill(deltaT);
				
				switch(assignRandomType(event)) {
					case 0:
						hDeltaT_typeI->Fill(deltaT);
						break;
					case 1:
						hDeltaT_typeIIa->Fill(deltaT);
						break;
					case 2:
						hDeltaT_typeIIb->Fill(deltaT);
						break;
					case 3:
						hDeltaT_typeIII->Fill(deltaT);
						break;
					
				}
			}
		}
    }
    
	drawHistogram(hDeltaT, "deltaT_Rand_Gen_hist.png", hitsVec.size()/3);
	drawHistogram(hDeltaT_ref, "deltaT_Rand_ref_Gen_hist.png", hitsVec.size()/3);
	drawHistogram(hDeltaT_nonRef, "deltaT_Rand_nonRef_Gen_hist.png", hitsVec.size()/3);
	drawHistogram(hDeltaT_typeI, "deltaT_Rand_typeI_Gen_hist.png", hitsVec.size()/3);
	drawHistogram(hDeltaT_typeIIa, "deltaT_Rand_typeIIa_Gen_hist.png", hitsVec.size()/3);
	drawHistogram(hDeltaT_typeIIb, "deltaT_Rand_typeIIb_Gen_hist.png", hitsVec.size()/3);
	drawHistogram(hDeltaT_typeIII, "deltaT_Rand_typeIII_Gen_hist.png", hitsVec.size()/3);
    delete hDeltaT;
    delete hDeltaT_ref;
    delete hDeltaT_nonRef;
    delete hDeltaT_typeI;
    delete hDeltaT_typeIIa;
    delete hDeltaT_typeIIb;
    delete hDeltaT_typeIII;
}

static void plotTruesIntervalsFromGen(int window_count, double activity){
	auto hitsVec = generateEvents(window_count, activity);
	std::vector<std::vector<gammaP>> pairs;
	
	int nEvents = findTrues(hitsVec, &pairs);
	
	// Histogramy dla różnic czasów prompt - gamma
    TH1F* hDeltaT = new TH1F("hDeltaT", "Prompt-gamma time diff", 100, -2000, 3000);
    TH1F* hDeltaT_Avg = new TH1F("hDeltaT_Avg", "Prompt-gamma time diff (Avg of anih. time)", 100, -2000, 3000);
    
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
        float avgTime = 0.0f;
        for (int gammaIdx : gammaIdxs) {
        	float deltaT = event[gammaIdx].getTime() - promptTime;
        	avgTime += event[gammaIdx].getTime();
			hDeltaT->Fill(deltaT);
		}
		
		hDeltaT_Avg->Fill(avgTime/2 - promptTime);
    }
    
    drawHistogram(hDeltaT, "deltaT_Trues_Gen_hist.png", hitsVec.size()/3);
    drawHistogram(hDeltaT_Avg, "deltaT_Trues_Gen_hist_avg.png", hitsVec.size()/3);
    delete hDeltaT;
    delete hDeltaT_Avg;
}


static void plotTruesAndRandoms(int window_count, double activity){
	auto hitsVec = generateEvents(window_count, activity);
	std::vector<std::vector<gammaP>> pairs;
	std::vector<int> DTWCounts = {0, 0, 0, 0};
	
	int nEvents = findCoincidences(hitsVec, &pairs);
	int nEvents_True = 0, nEvents_Rand = 0;
	
	// Histogramy
    TH1F* hDeltaT_True = new TH1F("hDeltaT_True", "Prompt-gamma time diff", 100, -2000, 3000);
    TH1F* hDeltaT_Rand = new TH1F("hDeltaT_Rand", "Prompt-gamma time diff", 100, -2000, 3000);
    
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
				DTWCounts[assignRandomType(event)]++;
        	} else {
        		hDeltaT_True->Fill(deltaT);
        		nEvents_True++;
			}
		}
    }
    
    // Zapis danych o randomach
    std::ofstream outFile("randomTypeData.txt");
    if (!outFile) {
        std::cerr << "Error opening file!\n";
    }

    outFile << "Random type: I	|	Count == " << DTWCounts[0] << std::endl;
    outFile << "Random type: IIa	|	Count == " << DTWCounts[1] << std::endl;
    outFile << "Random type: IIb	|	Count == " << DTWCounts[2] << std::endl;
    outFile << "Random type: III	|	Count == " << DTWCounts[3] << std::endl;

    outFile.close();
    
    //Hist
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c", "DeltaT Comparison", 800, 600);
	// Fix zero/negative bins
	for (int i = 0; i <= hDeltaT_True->GetNbinsX() + 1; ++i) {
		if (hDeltaT_True->GetBinContent(i) <= 0)
		    hDeltaT_True->SetBinContent(i, 1e-5);
	}
	for (int i = 0; i <= hDeltaT_Rand->GetNbinsX() + 1; ++i) {
		if (hDeltaT_Rand->GetBinContent(i) <= 0)
		    hDeltaT_Rand->SetBinContent(i, 1e-5);
	}
    c->SetLogy();
    // Set display minimum
	hDeltaT_True->SetMinimum(1e-5);
	hDeltaT_Rand->SetMinimum(1e-5);
    
    // Trues hist
    hDeltaT_True->SetLineColor(kBlue);
    hDeltaT_True->SetLineWidth(2);
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

