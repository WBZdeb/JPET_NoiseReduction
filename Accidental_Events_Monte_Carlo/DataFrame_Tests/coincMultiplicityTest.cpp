#include <iostream>
#include <vector>
#include <map>
#include <string>
#include "../testDTW.cpp"

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

static void coincMultiplicityTest(int window_count, double activity){
	auto hitsVec = generateEvents(window_count, activity);
	std::vector<std::vector<gammaP>> pairs;
	std::map<int, int> coincMultiplicity;
	int aboveOneCount = 0;
	int numOfMissedCoincidences = 0;
	int numOfPairs = 0;
	
	int nEvents = findAllPossibleCoincidences(hitsVec, &coincMultiplicity, &pairs);
	
	std::ofstream outFile("coincMult.txt");

    if (!outFile) {
        std::cerr << "Error opening file.\n";
        exit(0);
    }

    for (const auto& pair : coincMultiplicity) {
    	if( pair.second > 0){
        	outFile << "Gamma with index " << pair.first << "	|	Count: " << pair.second << "\n";
        	numOfPairs++;
        }
        if( pair.second > 1){
        	aboveOneCount++;
        	numOfMissedCoincidences += pair.second - 1;
        }
    }
    
    std::cout << "Non-propmt gammas count:	" << coincMultiplicity.size() * 2 << std::endl;
    std::cout << "Number of pairs:	" << numOfPairs << std::endl;
    std::cout << "Number of gammas with multi coincidences:	" << aboveOneCount << std::endl;
    std::cout << "Number of missed coincidences:	" << numOfMissedCoincidences << std::endl;
    
    TH1F* hDeltaT = new TH1F("hDeltaT", "Prompt-gamma time diff", 100, -10000, 6000);
    
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
    
    drawHistogram(hDeltaT, "deltaT_AllPossibleCoinc_hist.png", hitsVec.size()/3);
    delete hDeltaT;

    outFile.close();
}
