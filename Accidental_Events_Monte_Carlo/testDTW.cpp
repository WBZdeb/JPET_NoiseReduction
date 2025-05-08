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
	double time;
	int timeWindowNum;
	int eventNum;
	bool bIsPrompt;
	bool paired;
public:
	gammaP(double time, int timeWindowNum, int eventNum, bool isPrompt = false)
    	: time(time), timeWindowNum(timeWindowNum), eventNum(eventNum), bIsPrompt(isPrompt), paired(false) {}
	
	double getTime() const{
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
bool compareGammaP(const gammaP& a, const gammaP& b){
	return a.getTime() < b.getTime();
}


//Function for finding randoms; returns number, stores pairs inside passed vector
int findRandoms(std::vector<gammaP>& hitsVec, std::vector<std::vector<gammaP>>* pairs = nullptr) {
    int randomsCount = 0;
    int windowStartIndex = 0;
    double gammaCoincWindowLen = 500.0;
    double promptCoincWindowLen = 1300.0;	

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

// Calculate intervals between two 511
std::vector<double> calcPromptIntervals(std::vector<gammaP>& hitsVec, int window) {
	double prevTime = 0.0;
	std::vector<double> intervals;
	
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


std::vector<double> calcGammaIntervals(std::vector<gammaP>& hitsVec, int window) {
	int lastPair = -1;
	std::vector<double> intervals;
	
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
	const double kTimeResolution = 250.0;

	double meanItPerWindow = activity * (-WIN_LEN) * 1e-12; //Number of events per window
	std::vector<gammaP> hitsVec, windowVec;
	int eventNum = 0;
	TRandom gRand;

	for(int window = 0; window < windowCount; window++){
        //itPerWindow -> to jest lambda = <N> w oknie         	
        //e{-n*lambda} * lambda^n/n!
        auto itPerWindow = gRand.Poisson(meanItPerWindow);
		//generate hits for a given window
		for(int iter = 0; iter < itPerWindow; iter++){
		
			//generate prompt
			double time = gRand.Uniform(WIN_LEN, 0.0);
			//std::cout << time << "	" << window << std::endl;
			windowVec.push_back( gammaP(time, window, eventNum, true) );
			
			//generate decay time (in ps)
			double decayTime = gRand.Exp(kParaDecayTime);
			
			//Apply gauss twice to generate two 511-gammas
			for (int i = 0; i < 2; i++) {
				double t_511 = gRand.Gaus(time + decayTime, kTimeResolution);
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


// [activity] = Bq
// A = 700000.0 Bq = 0.7 MBq
// T_electronic_window = -50000000 = 50 *10^6 ps = 5000 ns = 5 us = 5 * 10^-6 s
// 10 ^{-12} = ps
// <N> = 5*10^-6 * 0.7*10^6 = 3.5
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
	
/*	
	//Create empty dataframe
	ROOT::RDataFrame empty_df(hitsVec.size());
	
	//Fill dataframe based on generated hits
	auto df = empty_df.Define("time",[&hitsVec]() {
						std::vector<double> time;
						for(const auto& hit : hitsVec){
							time.push_back(hit.getTime());
						}
						return time;
					}).Define("energy",[&hitsVec]() {
						std::vector<double> energy;
						double En;
						for(const auto& hit : hitsVec){
							if(hit.isPrompt()){
								En = 900.0;
							} else {
								En = 300.0;
							}
							energy.push_back(En);
						}
						return energy;
					}).Define("timeWindowNumber",[&hitsVec]() {
						std::vector<double> twn;
						for(const auto& hit : hitsVec){
							twn.push_back(hit.getTimeWindowNum());
						}
						return twn;
					});
	
	//Save the dataframe (will overwrite)
	auto file = TFile::Open("testData.root", "RECREATE");
	df.Snapshot("testTree", "testData.root");
	file->Close();
*/	
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




