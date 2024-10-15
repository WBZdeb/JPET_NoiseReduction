#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <TRandom.h>
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
	gammaP(double time, int timeWindowNum, int eventNum, bool isPrompt = false){
		this->time = time;
		this->timeWindowNum = timeWindowNum;
		this->eventNum = eventNum;
		this->bIsPrompt = isPrompt;
		paired = false;
	}
	
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
bool compareGammaP(gammaP& a, gammaP& b){
	return a.getTime() < b.getTime();
}


//Function for finding randoms; returns number, stores pairs inside passed vector
int findRandoms(std::vector<gammaP>& hitsVec, std::vector<std::vector<gammaP>>* pairs = nullptr) {
    int randomsCount = 0;
    int windowStartIndex = 0;

    //Iterate through each gammaP in hitsVec
    for( int i = 0; i < hitsVec.size(); ++i ) {
    
    	//Adjust starting index of time window
    	if( hitsVec[i].getTimeWindowNum() > hitsVec[windowStartIndex].getTimeWindowNum() ){
    		windowStartIndex = i;
    	}
    
        //Start by finding not yet paired 511 gamma
        if( hitsVec[i].isPrompt() == false && !hitsVec[i].isPaired() ) {
        
            //Find another 511 gamma withing same time window, up to 5000 ns away
            for( int j = i + 1; j < hitsVec.size(); ++j ) {
                if( hitsVec[j].isPrompt() == false && !hitsVec[j].isPaired() &&
                    hitsVec[j].getTimeWindowNum() == hitsVec[i].getTimeWindowNum() &&
                    hitsVec[j].getTime() <= hitsVec[i].getTime() + 3000000 ) {
                    
                    //Find a not-paired prompt gamma within 23000 ns around reference 511
                    for( int k = windowStartIndex; k < hitsVec.size(); ++k ) {
                        if( hitsVec[k].isPrompt() && !hitsVec[k].isPaired() &&
                            hitsVec[k].getTime() >= hitsVec[i].getTime() - 9000000 &&
                            hitsVec[k].getTime() <= hitsVec[i].getTime() + 3000000 ) {
                            
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


std::vector<double> calcPromptIntervals(std::vector<gammaP>& hitsVec, int window) {
	double prevTime = 0.0;
	std::vector<double> intervals;
	
	//Iterate through each gammaP in hitsVec
    for(int i = 0; i < hitsVec.size(); ++i) {
    	
    	//Consider only gammaP within given window
    	if( hitsVec[i].getTimeWindowNum() != window ){
    		continue;
    	}
    	
    	//Check if prompt
    	if( hitsVec[i].isPrompt() ){
    		intervals.push_back( hitsVec[i].getTime() - prevTime );
    		prevTime = hitsVec[i].getTime();
    	}
    }
	
	//First element is invalid, need to be removed
	if( !intervals.empty() ){
		intervals.erase( intervals.begin() );
	}
	
	return intervals;
}

// [activity] = Bq
// A = 700000.0 Bq = 0.7 MBq
// T_electronic_window = -50000000 = 50 *10^6 ps = 5000 ns = 5 us = 5 * 10^-6 s
// 10 ^{-12} = ps
// <N> = 5*10^-6 * 0.7*10^6 = 3.5
// T_anih = 3 lub 4 ns
// T_prompt_anihi = 10 ns
void testDTW(int window_count = 20, double activity = 700000.0){
	//If window_count or activity is too small, terminate macro
	if(window_count <= 0 || activity <= 0){
		std::cout << "Command line arguments should be greater than zero" << std::endl;
		exit(1);
	}
	
	double meanItPerWindow = activity * (-WIN_LEN) * TMath::Power(10.0, -12.0); //Number of events per window
	std::vector<gammaP> hitsVec;	//Vector storing all generated hits of gamma particles (both 511 and prompt)
	std::vector<gammaP> windowVec;	//Vector storing particles generated for specific window
	int eventNum = 0;
	
	TRandom gRand;
    const double kParaDecayTime = 125;
    const double kTimeResolution = 250;
	
	//Populate vector
	for(int window = 0; window < window_count; window++){
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
			double time_511 = gRand.Gaus(time+decayTime, kTimeResolution);
			if(time_511 < 0){
				windowVec.push_back( gammaP(time_511, window, eventNum, false) );
			}
			
			time_511 = gRand.Gaus(time+decayTime, kTimeResolution);
			if(time_511 < 0){
				windowVec.push_back( gammaP(time_511, window, eventNum++, false) );
			}
		}
		//sort vector before appending
		std::sort(windowVec.begin(), windowVec.end(), compareGammaP);
		
		//append window vector to vector of all hits
		hitsVec.insert(hitsVec.end(), windowVec.begin(), windowVec.end());
	}
	
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
	
	//Test random count
	std::cout << "Random count for:		window_count = " << window_count << "	activity = " << activity << std::endl;
	std::cout << "All randoms: " << findRandoms(hitsVec) << std::endl;	
	

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

}


