/// guards
#ifndef CALC_LIFETIME
#define CALC_LIFETIME

#include <iostream>
#include <vector>
#include <TROOT.h>
#ifdef R__HAS_VDT
#undef R__HAS_VDT
#endif
#include "ROOT/RDataFrame.hxx"
#ifdef R__HAS_VDT
#undef R__HAS_VDT
#endif

//TAK NIE ROBIMY! -> namespace pollution
//using namespace std;

std::vector<float> calc_lifetime(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> df){
	float timeRescale = 1.0/1000.0;	//rescale time unit to nanoseconds
	const float energyTH = 350.0;
	float c = 30;		//speed of light [cm/ns]
	auto time_Vec = *(df.Take<std::vector<float>>("time"));
	auto x_Vec = *(df.Take<std::vector<float>>("x"));
	auto y_Vec = *(df.Take<std::vector<float>>("y"));
	auto z_Vec = *(df.Take<std::vector<float>>("z"));
	auto energy = *(df.Take<std::vector<float>>("energy"));
	std::vector<float> lifetimes;	//vector for calculated lifetimes
	
	float tFl1, tFl2, tEm1, tEm2, refTime, lenAB, origin[3] = {0.0, 0.0, 0.0};
	std::vector<float> vecAB, pDist, d1, d2;
	int iMap[3] = {0, 1, 2};
	
	for (int i = 0; i < time_Vec.size(); i++){
		//set index mapping based on prompt index
		if(energy[i][0] > energyTH){
			iMap[0] = 0; 
			if(time_Vec[i][1] > time_Vec[i][2]){ iMap[1] = 1; iMap[2] = 2;}
			else { iMap[1] = 2; iMap[2] = 1;}
		}
		else if(energy[i][1] > energyTH){
			iMap[0] = 1;
			if(time_Vec[i][0] > time_Vec[i][2]){ iMap[1] = 0; iMap[2] = 2;}
			else { iMap[1] = 2; iMap[2] = 0;}
		}
		else if(energy[i][2] > energyTH){
			iMap[0] = 2;
			if(time_Vec[i][0] > time_Vec[i][1]){ iMap[1] = 0; iMap[2] = 1;}
			else { iMap[1] = 1; iMap[2] = 0;}
		}
		else{
			continue;
		}
	
		//determine decay origin
		vecAB = { (x_Vec[i][iMap[1]] - x_Vec[i][iMap[2]]), (y_Vec[i][iMap[1]] - y_Vec[i][iMap[2]]), (z_Vec[i][iMap[1]] - z_Vec[i][iMap[2]]) };
		lenAB = TMath::Sqrt( vecAB[0]*vecAB[0] + vecAB[1]*vecAB[1] + vecAB[2]*vecAB[2]);
		
		origin[0] = x_Vec[i][iMap[2]] + (0.5 - c*timeRescale*(time_Vec[i][iMap[1]]-time_Vec[i][iMap[2]]) / (2*lenAB) ) * vecAB[0];
		origin[1] = y_Vec[i][iMap[2]] + (0.5 - c*timeRescale*(time_Vec[i][iMap[1]]-time_Vec[i][iMap[2]]) / (2*lenAB) ) * vecAB[1];
		origin[2] = z_Vec[i][iMap[2]] + (0.5 - c*timeRescale*(time_Vec[i][iMap[1]]-time_Vec[i][iMap[2]]) / (2*lenAB) ) * vecAB[2];
		//std::cout << origin[0] << "	" << origin[1] << "	" << origin[2] << std::endl;
	
		//calculate "machine time" T
		pDist = { x_Vec[i][iMap[0]] - origin[0], y_Vec[i][iMap[0]] - origin[1], z_Vec[i][iMap[0]] - origin[2] };
		refTime = timeRescale*time_Vec[i][iMap[0]] - TMath::Sqrt( pDist[0]*pDist[0] + pDist[1]*pDist[1] + pDist[2]*pDist[2]) / c;
		
		//calculate lifetime
		pDist = { x_Vec[i][iMap[1]] - origin[0], y_Vec[i][iMap[1]] - origin[1], z_Vec[i][iMap[1]] - origin[2] };
		tFl1 = TMath::Sqrt( pDist[0]*pDist[0] + pDist[1]*pDist[1] + pDist[2]*pDist[2]) / c;
		tEm1 = timeRescale*time_Vec[i][iMap[1]] - tFl1 - refTime;
		
		pDist = { x_Vec[i][iMap[2]] - origin[0], y_Vec[i][iMap[2]] - origin[1], z_Vec[i][iMap[2]] - origin[2] };
		tFl2 = TMath::Sqrt( pDist[0]*pDist[0] + pDist[1]*pDist[1] + pDist[2]*pDist[2]) / c;
		tEm2 = timeRescale*time_Vec[i][iMap[2]] - tFl2 - refTime;
		
		lifetimes.push_back( (tEm1 + tEm2)/2);
	}
	
	return lifetimes;
}


void DTW_type1(ROOT::RDataFrame df){
	const float lookahead = 5000, energyTH = 350.0;
	auto time_Vec = *(df.Take<std::vector<float>>("time"));
	auto energy = *(df.Take<std::vector<float>>("energy"));
	auto window_num = *(df.Take<int>("timeWindowNumber"));
	
	int event_num = 0, coinc_num = 0, lw = window_num.back();
	int cw, iStart, iEnd;
	
	//vector of vectors for tagging paired hits
	std::vector<std::vector<int>> paired = {};
	
	for(int i = 0; i < time_Vec.size(); i++){
		std::vector<int> v(time_Vec[i].size(), 0);
		paired.push_back(v);
	}
	
	for(int out = 0; out < time_Vec.size(); out++){
		for(int in = 0; in < time_Vec[out].size(); in++){
			event_num++;
			//Exit when at the second-to-last window; skip paired hits or prompts
			if(window_num[event_num-1] == (lw-1)) break;
			if(paired[out][in] == 1) break;
			if(energy[out][in] > energyTH) break;
			
			//set range of indexes to search for coincidence
			cw = window_num[event_num-1];
			for(int i = (out + 1); i < time_Vec.size(); i++){
				if(time_Vec[i][0] == (cw+2)){
					iStart = i;
					if(cw == lw){
						iEnd = time_Vec.size();
						break;
					}
				}
				if(time_Vec[i][0] == (cw+3)){
					iEnd = i;
					break;
				}
			}
			
			//look for match
			for(int j = iStart; j < iEnd; j++){
				for(int k = 0; k < time_Vec[j].size(); k++){
					if(paired[j][k] == 1) break;
					if(energy[j][k] > energyTH) break;
					if( (time_Vec[j][k] - time_Vec[out][in]) < lookahead){
						coinc_num+=2;
						paired[out][in] = 1;
						paired[j][k] = 1;
					}
				}
			}
		}
	}
	std::cout << "All hits: " << event_num << std::endl;
	std::cout << "Randoms (type 1): " << coinc_num << std::endl;
	std::cout << "Percentage: " << coinc_num/event_num << "%" << std::endl;
}


#endif
