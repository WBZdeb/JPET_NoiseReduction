#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ROOT/RDataFrame.hxx>
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TRandom.h>

// [WIN_LEN] = ps
#define WIN_LEN -50000000.0

class gammaP {
private:
	float time;
	int timeWindowNum;
	int eventNum;
	bool bIsPrompt;
	bool paired;
public:
	gammaP(float time, int timeWindowNum, int eventNum, bool isPrompt = false)
    	: time(time), timeWindowNum(timeWindowNum), eventNum(eventNum), bIsPrompt(isPrompt), paired(false) {}
	
	float getTime() const{
		return this->time;
	}
	
	int getTimeWindowNum() const{
		return this->timeWindowNum;
	}
	
	int getEventNum() const{
		return this->eventNum;
	}
	
	bool isPrompt() const{
		return this->bIsPrompt;
	}
	
	bool isPaired() const{
		return this->paired;
	}
	
	void setPaired(bool val){
		paired = val;
	}
};

//Function for gammaP sorting
static bool compareGammaP(const gammaP& a, const gammaP& b){
	return a.getTime() < b.getTime();
}


//Function for finding randoms; returns number, stores pairs inside passed vector
int findRandoms(std::vector<gammaP>& hitsVec, std::vector<std::vector<gammaP>>* pairs = nullptr) {
    int randomsCount = 0;
    int windowStartIndex = 0;
    float gammaCoincWindowLen = 5000.0f;
    float promptCoincWindowLen = 13000.0f;	

    //Iterate through each gammaP in hitsVec
    for( int i = 0; i < hitsVec.size(); ++i ) {
    
    	//Adjust starting index of time window
    	if( hitsVec[i].getTimeWindowNum() > hitsVec[windowStartIndex].getTimeWindowNum() ){
    		windowStartIndex = i;
    	}
    
        //Start by finding not yet paired 511 gamma
        if( hitsVec[i].isPrompt() == false && !hitsVec[i].isPaired() ) {
        
            //Find another 511 gamma withing same time window, up to 0.5 ns away
            for( int j = i + 1; j < hitsVec.size(); ++j ) {
                if( hitsVec[j].isPrompt() == false && !hitsVec[j].isPaired() &&
                    hitsVec[j].getTimeWindowNum() == hitsVec[i].getTimeWindowNum() &&
                    hitsVec[j].getTime() <= hitsVec[i].getTime() + gammaCoincWindowLen ) {
                    
                    //Find a not-paired prompt gamma within 1.8 ns around reference 511
                    for( int k = windowStartIndex; k < hitsVec.size(); ++k ) {
                        if( hitsVec[k].isPrompt() && !hitsVec[k].isPaired() &&
                            hitsVec[k].getTime() >= hitsVec[i].getTime() - promptCoincWindowLen &&
                            hitsVec[k].getTime() <= hitsVec[i].getTime() + gammaCoincWindowLen ) {
                            
                            //Compare eventNum
                            if( hitsVec[i].getEventNum() != hitsVec[j].getEventNum() ||
                                hitsVec[i].getEventNum() != hitsVec[k].getEventNum() ||
                                hitsVec[j].getEventNum() != hitsVec[k].getEventNum() ) {
                                
                                //If any unequal, increment the random count
                                randomsCount++;
                            }

                            //Mark all three as paired
                            hitsVec[i].setPaired(true);
                            hitsVec[j].setPaired(true);
                            hitsVec[k].setPaired(true);

                            //Store the paired gammas
                            if( pairs != nullptr ) {
                                pairs->push_back({ hitsVec[i], hitsVec[j], hitsVec[k] });
                            }
                            break;
                        }
                    }
                    break;
                }
            }
        }
    }
    return randomsCount;
}


//Function for finding true pairs; returns number, stores pairs inside passed vector
int findTrues(std::vector<gammaP>& hitsVec, std::vector<std::vector<gammaP>>* pairs = nullptr) {
    int truesCount = 0;
    int windowStartIndex = 0;
    float gammaCoincWindowLen = 500.0f;
    float promptCoincWindowLen = 1300.0f;	

    //Iterate through each gammaP in hitsVec
    for( int i = 0; i < hitsVec.size(); ++i ) {
    
    	//Adjust starting index of time window
    	if( hitsVec[i].getTimeWindowNum() > hitsVec[windowStartIndex].getTimeWindowNum() ){
    		windowStartIndex = i;
    	}
    
        //Start by finding not yet paired 511 gamma
        if( hitsVec[i].isPrompt() == false && !hitsVec[i].isPaired() ) {
        
            //Find another 511 gamma withing same time window, up to 0.5 ns away
            for( int j = i + 1; j < hitsVec.size(); ++j ) {
                if( hitsVec[j].isPrompt() == false && !hitsVec[j].isPaired() &&
                    hitsVec[j].getTimeWindowNum() == hitsVec[i].getTimeWindowNum() &&
                    hitsVec[j].getTime() <= hitsVec[i].getTime() + gammaCoincWindowLen ) {
                    
                    //Find a not-paired prompt gamma within 1.8 ns around reference 511
                    for( int k = windowStartIndex; k < hitsVec.size(); ++k ) {
                        if( hitsVec[k].isPrompt() && !hitsVec[k].isPaired() &&
                            hitsVec[k].getTime() >= hitsVec[i].getTime() - promptCoincWindowLen &&
                            hitsVec[k].getTime() <= hitsVec[i].getTime() + gammaCoincWindowLen ) {
                            
                            //Compare eventNum
                            if( hitsVec[i].getEventNum() == hitsVec[j].getEventNum() &&
                                hitsVec[i].getEventNum() == hitsVec[k].getEventNum() ) {
                                
                                //If all equal, increment the random count
                                truesCount++;
                            }

                            //Mark all three as paired
                            hitsVec[i].setPaired(true);
                            hitsVec[j].setPaired(true);
                            hitsVec[k].setPaired(true);

                            //Store the paired gammas
                            if( pairs != nullptr ) {
                                pairs->push_back({ hitsVec[i], hitsVec[j], hitsVec[k] });
                            }
                            break;
                        }
                    }
                    break;
                }
            }
        }
    }
    return truesCount;
}


// Calculate intervals between two 511
static std::vector<float> calcPromptIntervals(std::vector<gammaP>& hitsVec, int window) {
	float prevTime = 0.0f;
	std::vector<float> intervals;
	
	//Iterate through each gammaP in hitsVec
    for(const auto& hit : hitsVec) {	
    	//Consider only gammaP within given window
    	if( hit.getTimeWindowNum() != window ) continue;
    	
    	//Check if prompt
    	if( hit.isPrompt() ){
    		intervals.push_back(hit.getTime() - prevTime);
        	prevTime = hit.getTime();
    	}
    }
	
	//First element is invalid, need to be removed
	if( !intervals.empty() ){
		intervals.erase( intervals.begin() );
	}
	
	return intervals;
}


static std::vector<float> calcGammaIntervals(std::vector<gammaP>& hitsVec, int window) {
	int lastPair = -1;
	std::vector<float> intervals;
	
	//Iterate through each gammaP in hitsVec
    for(int i = 0; i < hitsVec.size(); ++i) {
    	
    	//Consider only gammaP within given window
    	if( hitsVec[i].getTimeWindowNum() != window ){
    		continue;
    	}
    	
    	//Check if 511
    	if( !hitsVec[i].isPrompt() && hitsVec[i].getEventNum() != lastPair ){
    		//look for 511 with identical eventNum
    		for(int j = i+1; j < hitsVec.size(); ++j) {
    		
    			if( hitsVec[i].getEventNum() == hitsVec[j].getEventNum() ){
    			
    				if( !hitsVec[j].isPrompt() ){
    					intervals.push_back( hitsVec[j].getTime() - hitsVec[i].getTime() );
    					lastPair = hitsVec[i].getEventNum();
    				}
    			} else {
    				break;
    			}
    		}
    	}
    }
	
	return intervals;
}


std::vector<gammaP> generateEvents(int windowCount, double activity) {
	const double kParaDecayTime = 125.0;
	const double kTimeResolution = 62.5;	// 250 ps

	double meanItPerWindow = activity * (-WIN_LEN) * 1e-12; //Number of events per window
	std::vector<gammaP> hitsVec, windowVec;
	int eventNum = 0;
	TRandom gRand(time(0));

	for(int window = 0; window < windowCount; window++){
        //itPerWindow -> to jest lambda = <N> w oknie         	
        //e{-n*lambda} * lambda^n/n!
        auto itPerWindow = gRand.Poisson(meanItPerWindow); //może być 0
		//generate hits for a given window
		for(int iter = 0; iter < itPerWindow; iter++){
		
			//generate prompt
			float time = gRand.Uniform(WIN_LEN, 0.0);
			//std::cout << "Pushing prompt:   " << time << "	" << window << std::endl;
			windowVec.push_back( gammaP(time, window, eventNum, true) );
			
			//generate decay time (in ps)
			float decayTime = gRand.Exp(kParaDecayTime);
			
			//Apply gauss twice to generate two 511-gammas
			for (int i = 0; i < 2; i++) {
				float t_511 = gRand.Gaus(time + decayTime, kTimeResolution);
				if (t_511 < 0) windowVec.push_back(gammaP(t_511, window, eventNum, false));
			}
			eventNum++;
		}
		//sort vector before appending
		std::sort(windowVec.begin(), windowVec.end(), compareGammaP);
		
		//append window vector to vector of all hits; clear after appending
		hitsVec.insert(hitsVec.end(), windowVec.begin(), windowVec.end());
		windowVec.clear();
	}
	return hitsVec;
}


std::string generate_DataFrame(int window_count, double activity){
	//If window_count or activity is too small, terminate macro
	if(window_count <= 0 || activity <= 0){
		std::cout << "Command line arguments should be greater than zero" << std::endl;
		exit(1);
	}
	
	//Populate vector
	auto hitsVec = generateEvents(window_count, activity);
	
	// Group hits by eventNum
	std::map<int, std::vector<gammaP>> groupedHits;
	for(const auto& hit : hitsVec){
		groupedHits[hit.getEventNum()].push_back(hit);
	}
	
	// Prepare vectors to fill dataframe columns
	std::vector<std::vector<float>> allTimes;
	std::vector<std::vector<float>> allEnergies;
	std::vector<int> allWindowNums;
	
	for(const auto& [eventNum, hits] : groupedHits){
		std::vector<float> times;
		std::vector<float> energies;
		
		for(const auto& hit : hits){
			times.push_back(hit.getTime());
			energies.push_back(hit.isPrompt() ? 900.0f : 300.0f);
		}
		
		allTimes.push_back(times);
		allEnergies.push_back(energies);
		allWindowNums.push_back(hits.front().getTimeWindowNum());
	}
	
	//Create empty dataframe
	std::vector<int> indices(allTimes.size());
	std::iota(indices.begin(), indices.end(), 0);
	ROOT::RDataFrame empty_df(allTimes.size());
	
	//Fill dataframe based on generated hits
	auto df = empty_df
		.Define("index", [&indices](ULong64_t i) { return indices[static_cast<int>(i)]; }, {"rdfentry_"})
		.Define("time", [&allTimes](int i) { return allTimes[i]; }, {"index"})
    	.Define("energy", [&allEnergies](int i) { return allEnergies[i]; }, {"index"})
    	.Define("timeWindowNumber", [&allWindowNums](int i) { return allWindowNums[i]; }, {"index"});
	
	//Save the dataframe (will overwrite)
	std::string outputFileName = "generatedData.root";
	auto file = TFile::Open(outputFileName.c_str(), "RECREATE");
	df.Snapshot("FlatTree", outputFileName.c_str());
	file->Close();
	
	return outputFileName;
}


// [activity] = Bq
// A = 700000.0 Bq = 0.7 MBq
// T_electronic_window = -50000000 = 50 *10^6 ps = 5000 ns = 5 us = 5 * 10^-6 s
// 10 ^{-12} = ps
// <N> = 5*10^-6 * 0.7*10^6 = 35
// T_anih = 3 lub 4 ns
// T_prompt_anihi = 10 ns
void testDTW(int window_count = 10, double activity = 700000.0){
	//If window_count or activity is too small, terminate macro
	if(window_count <= 0 || activity <= 0){
		std::cout << "Command line arguments should be greater than zero" << std::endl;
		exit(1);
	}
	
	//Populate vector
	auto hitsVec = generateEvents(window_count, activity);
	
	
	//Test random count
	std::vector<std::vector<gammaP>> pairs;
	
	std::cout << "Random count for:		window_count = " << window_count << "	activity = " << activity << std::endl;
	std::cout << "All randoms: " << findRandoms(hitsVec, &pairs)/ (double)window_count << std::endl;
	
	//Calc different types of randoms
	int rTypes[4] = {0, 0, 0, 0};
	
	for(int pair = 0; pair < pairs.size(); pair++){
		//Type I
		if( pairs[pair][0].getEventNum() == pairs[pair][1].getEventNum() &&
			pairs[pair][0].getEventNum() != pairs[pair][2].getEventNum() &&
			pairs[pair][1].getEventNum() != pairs[pair][2].getEventNum() ) rTypes[0]++;
			
		//Type IIa
		if( pairs[pair][0].getEventNum() != pairs[pair][1].getEventNum() &&
			pairs[pair][0].getEventNum() == pairs[pair][2].getEventNum() &&
			pairs[pair][1].getEventNum() != pairs[pair][2].getEventNum() ) rTypes[1]++;
					
		//Type IIb
		if( pairs[pair][0].getEventNum() != pairs[pair][1].getEventNum() &&
			pairs[pair][0].getEventNum() != pairs[pair][2].getEventNum() &&
			pairs[pair][1].getEventNum() == pairs[pair][2].getEventNum() ) rTypes[2]++;
					
		//Type III
		if( pairs[pair][0].getEventNum() != pairs[pair][1].getEventNum() &&
			pairs[pair][0].getEventNum() != pairs[pair][2].getEventNum() &&
			pairs[pair][1].getEventNum() != pairs[pair][2].getEventNum() ) rTypes[3]++;
					
	}
	
	std::cout << "Type I: " << rTypes[0]/ (double)window_count << ",  Type IIa: " << rTypes[1]/ (double)window_count; 
	std::cout << ",  Type IIb: " << rTypes[2]/ (double)window_count << ",  Type III: " << rTypes[3]/ (double)window_count << std::endl;
/*
	//Draw prompt interval histograms
	std::vector<const char*> hNames = {"Window1", "Window2", "Window3", "Window4"};
	
	for(int win = 1; win < 5; win++){
		std::vector<double> intervals = calcPromptIntervals(hitsVec, win);
	
		TH1D hist(hNames[win-1], hNames[win-1], 30, 0.0, 12000000.0);
		
		for(double val : intervals){
			hist.Fill(val);
		}
		
		TCanvas canvas;
		hist.Draw();
		
		std::string filePath = "promptIntervals/window" + std::to_string(win) + ".jpg";
		canvas.SaveAs(filePath.c_str());
	}

	//Draw 511 interval histograms
	std::vector<double> allIntervals;
	
	for(int win = 0; win < window_count; win++){
		std::vector<double> intervals = calcGammaIntervals(hitsVec, win);
		
		allIntervals.insert(allIntervals.end(), intervals.begin(), intervals.end());
	}
	
	//Save vector to txt file for reading in ROOT
	std::ofstream outFile("promptIntervals/Intervals_511.txt");
	for(double interval :  allIntervals){
		outFile << interval << std::endl;
	}
	outFile.close();
	
	
	TH1D hist("Intervals_511", "Intervals_511", 30, 0.0, 2000.0);
		
	for(double val : allIntervals){
		hist.Fill(val);
	}
		
	TCanvas canvas;
	hist.Draw();
		
	std::string filePath = "promptIntervals/Intervals_511.jpg";
	canvas.SaveAs(filePath.c_str());
*/
}




