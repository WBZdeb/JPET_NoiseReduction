#include <iostream>
#include <cstdlib>
#include <vector>
#include <TRandom.h>
#include <ROOT/RDataFrame.hxx>
#include <TTree.h>
#include <TFile.h>
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
public:
	gammaP(double time, int timeWindowNum, int eventNum, bool isPrompt = false){
		this->time = time;
		this->timeWindowNum = timeWindowNum;
		this->eventNum = eventNum;
		this->bIsPrompt = isPrompt;
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
};

//Function for gammaP sorting
bool compareGammaP(gammaP& a, gammaP& b){
	return a.getTime() < b.getTime();
}


// [activity] = Bq
void testDTW(int window_count = 20, double activity = 700000.0){
	//If window_count or activity is too small, terminate macro
	if(window_count <= 0 || activity <= 0){
		std::cout << "Command line arguments should be greater than zero" << std::endl;
		exit(1);
	}
	
	double itPerWindow = activity * (-WIN_LEN) * TMath::Power(10.0, -12.0); //Number of events per window
	std::vector<gammaP> hitsVec;	//Vector storing all generated hits of gamma particles (both 511 and prompt)
	std::vector<gammaP> windowVec;	//Vector storing particles generated for specific window
	int eventNum = 0;
	
	TRandom gRand;
	
	//Create empty dataframe
	ROOT::RDataFrame empty_df(0);
	
	//Populate vector
	for(int window = 0; window < window_count; window++){
		//generate hits for a given window
		for(int iter = 0; iter < itPerWindow; iter++){
			//generate prompt
			double time = gRand.Uniform(WIN_LEN, 0.0);
			windowVec.push_back( *(new gammaP(time, window, eventNum, true)) );
			
			//generate decay time (in ps)
			double decayTime = gRand.Exp(125.0);
			
			//Apply gauss twice to generate two 511-gammas
			double time_511 = gRand.Gaus(time+decayTime, 250);
			if(time_511 < 0){
				windowVec.push_back( *(new gammaP(time_511, window, eventNum, true)) );
			}
			
			time_511 = gRand.Gaus(time+decayTime, 250);
			if(time_511 < 0){
				windowVec.push_back( *(new gammaP(time_511, window, eventNum++, true)) );
			}
		}
		//sort vector before appending
		std::sort(windowVec.begin(), windowVec.end(), compareGammaP);
		
		//append window vector to vector of all hits
		hitsVec.insert(windowVec.end(), windowVec.begin(), windowVec.end());
	}
	
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
	
	return 0;
}


