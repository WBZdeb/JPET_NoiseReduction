/// guards
#ifndef CALC_LIFETIME
#define CALC_LIFETIME

#include <iostream>
#include <fstream>
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

using RDF = ROOT::RDataFrame;
using RNode = ROOT::RDF::RNode;

std::vector<float> calc_lifetime(RNode rnode){
	auto df = rnode;

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
		//std::ci << origin[0] << "	" << origin[1] << "	" << origin[2] << std::endl;
	
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

//Only second 511 is delayed --> gives randoms of types IIa and III
void DTW_type1(RNode rnode, int skips){
	auto df = rnode;
	if(skips < 1) skips = 1;

	const float lookahead = 5000, p_lookbehind = -18000, energyTH = 350.0, promptTH = 511;
	auto time_Vec = *(df.Take<std::vector<float>>("time"));
	auto energy = *(df.Take<std::vector<float>>("energy"));
	auto window_num = *(df.Take<int>("timeWindowNumber"));
	
	int event_num = 0, coinc_num = 0, lw;
	int cw, wPrev = window_num[0], wNum = 1, wStartInd = 0;
	int save_pos[] = {-1, -1};
	ofstream stats;
	
	//vector of vectors for tagging paired hits
	std::vector<std::vector<int>> paired = {};
	for(int ind = 0; ind < time_Vec.size(); ind++){
		std::vector<int> v(time_Vec[ind].size(), 0);
		paired.push_back(v);
	}
	
	if(skips == 2) {	
		stats.open("DTW_stats.txt");
		stats << "Pair   |   TWindow   |   HitTime   |   Energy\n";
	}
	
	//Re-mapping of window numbers to avoid gaps
	for(int wInd = 0; wInd < window_num.size(); wInd++){
		if(window_num[wInd] == wPrev) {
			window_num[wInd] = wNum;
		} else {
			wPrev = window_num[wInd];
			wNum++;
			window_num[wInd] = wNum;
		}
	}
	lw = window_num.back();
	
	for(int i = 0; i < time_Vec.size(); i++){
		save_pos[0] = -1;
		save_pos[1] = -1;
		
		//Tracking first index of current window, so we can check entire window for prompt matching
		if( window_num[i] > window_num[wStartInd] ){
			wStartInd = i;
		}
		
		for(int j = 0; j < time_Vec[i].size(); j++){
		//Exit when at the second-to-last window; skip paired hits or prompts
			event_num++;
			if(window_num[i] == (lw-skips+1)) continue;
			if(paired[i][j] == 1) continue;
			if(energy[i][j] > energyTH) continue;
				
			cw = window_num[i];
				
			//look for match
			for(int k = (i+1); k < time_Vec.size(); k++){
				//if(paired[i][j] == 1) break;
				if(window_num[k] > (cw+skips)) break;
				if(window_num[k] < (cw+skips)) continue;
				
				for(int l = 0; l < time_Vec[k].size(); l++){
					//if(paired[i][j] == 1) break;
					if(paired[k][l] == 1) continue;
					if(energy[k][l] > energyTH) continue;
					if( TMath::Abs(time_Vec[k][l] - time_Vec[i][j]) < lookahead){
						save_pos[0] = k;
						save_pos[1] = l;
						
						if(skips == 2){
							stats << coinc_num << "   |   " << cw << "   |   " << time_Vec[i][j] << "   |   " << energy[i][j] << "\n";
							stats << coinc_num << "   |   " << window_num[k] << "   |   " << time_Vec[k][l] << "   |   " << energy[k][l] << "\n";
						}
					}
				}
			}
			//if no match - skip to next iteration
			if( save_pos[0] == -1 ) continue;
			
			//look for prompt
			for(int k = wStartInd; k < time_Vec.size(); k++){
				//if(paired[i][j] == 1) break;
				if(window_num[k] > cw) break;
				
				for(int l = 0; l < time_Vec[k].size(); l++){
					//if(paired[i][j] == 1) break;
					if(paired[k][l] == 1) continue;
					if(energy[k][l] <= promptTH) continue;
					
					//Prompt can be 5 ns after or 18 ns before reference 511
					if( ((time_Vec[k][l] - time_Vec[i][j]) < lookahead)  &&
					     	((time_Vec[k][l] - time_Vec[i][j]) > p_lookbehind) ){
					    
						coinc_num+=1;
						paired[i][j] = 1;
						paired[k][l] = 1;
						paired[ save_pos[0] ][ save_pos[1] ] = 1;
					}
				}
			}
		}
	}
	std::cout << "All hits: " << event_num << std::endl;
	std::cout << "Randoms (type 1): " << coinc_num << std::endl;
	std::cout << "Percentage: " << 100.0*coinc_num/event_num << "%" << std::endl;
	
	if(skips == 2) stats.close();
}

//Both second 511 and prompt are delayed --> gives randoms of types IIb and III
void DTW_type2(RNode rnode, int skips){
	auto df = rnode;
	if(skips < 1) skips = 1;

	const float lookahead = 5000, p_lookbehind = -18000, energyTH = 350.0, promptTH = 511;
	auto time_Vec = *(df.Take<std::vector<float>>("time"));
	auto energy = *(df.Take<std::vector<float>>("energy"));
	auto window_num = *(df.Take<int>("timeWindowNumber"));
	
	int event_num = 0, coinc_num = 0, lw;
	int cw, wPrev = window_num[0], wNum = 1, wStartInd = 0;
	int save_pos[] = {-1, -1};
	ofstream stats;
	
	//vector of vectors for tagging paired hits
	std::vector<std::vector<int>> paired = {};
	for(int ind = 0; ind < time_Vec.size(); ind++){
		std::vector<int> v(time_Vec[ind].size(), 0);
		paired.push_back(v);
	}
	
	if(skips == 2) {	
		stats.open("DTW_stats.txt");
		stats << "Pair   |   TWindow   |   HitTime   |   Energy\n";
	}
	
	//Re-mapping of window numbers to avoid gaps
	for(int wInd = 0; wInd < window_num.size(); wInd++){
		if(window_num[wInd] == wPrev) {
			window_num[wInd] = wNum;
		} else {
			wPrev = window_num[wInd];
			wNum++;
			window_num[wInd] = wNum;
		}
	}
	lw = window_num.back();
	
	for(int i = 0; i < time_Vec.size(); i++){
		save_pos[0] = -1;
		save_pos[1] = -1;
		
		//Tracking first index of current window, so we can check entire window for prompt matching
		if( window_num[i] > window_num[wStartInd] ){
			wStartInd = i;
		}
		
		for(int j = 0; j < time_Vec[i].size(); j++){
		//Exit when at the second-to-last window; skip paired hits or prompts
			event_num++;
			if(window_num[i] == (lw-skips+1)) continue;
			if(paired[i][j] == 1) continue;
			if(energy[i][j] > energyTH) continue;
				
			cw = window_num[i];
				
			//look for match
			for(int k = (i+1); k < time_Vec.size(); k++){
				//if(paired[i][j] == 1) break;
				if(window_num[k] > (cw+skips)) break;
				if(window_num[k] < (cw+skips)) continue;
				
				for(int l = 0; l < time_Vec[k].size(); l++){
					//if(paired[i][j] == 1) break;
					if(paired[k][l] == 1) continue;
					if(energy[k][l] > energyTH) continue;
					if( TMath::Abs(time_Vec[k][l] - time_Vec[i][j]) < lookahead){
						save_pos[0] = k;
						save_pos[1] = l;
						
						if(skips == 2){
							stats << coinc_num << "   |   " << cw << "   |   " << time_Vec[i][j] << "   |   " << energy[i][j] << "\n";
							stats << coinc_num << "   |   " << window_num[k] << "   |   " << time_Vec[k][l] << "   |   " << energy[k][l] << "\n";
						}
					}
				}
			}
			//if no match - skip to next iteration
			if( save_pos[0] == -1 ) continue;
			
			//look for prompt
			for(int k = (i+1); k < time_Vec.size(); k++){
				//if(paired[i][j] == 1) break;
				if(window_num[k] > (cw+skips)) break;
				if(window_num[k] < (cw+skips)) continue;
				
				for(int l = 0; l < time_Vec[k].size(); l++){
					//if(paired[i][j] == 1) break;
					if(paired[k][l] == 1) continue;
					if(energy[k][l] <= promptTH) continue;
					
					//Prompt can be 5 ns after or 18 ns before reference 511
					if( ((time_Vec[k][l] - time_Vec[i][j]) < lookahead)  &&
					     	((time_Vec[k][l] - time_Vec[i][j]) > p_lookbehind) ){
					    
						coinc_num+=1;
						paired[i][j] = 1;
						paired[k][l] = 1;
						paired[ save_pos[0] ][ save_pos[1] ] = 1;
					}
				}
			}
		}
	}
	std::cout << "All hits: " << event_num << std::endl;
	std::cout << "Randoms (type 1): " << coinc_num << std::endl;
	std::cout << "Percentage: " << 100.0*coinc_num/event_num << "%" << std::endl;
	
	if(skips == 2) stats.close();
}

//Second 511 and prompt are delayed with different shifts --> gives randoms of type III
void DTW_type3(RNode rnode, int skips){
	auto df = rnode;
	if(skips < 1) skips = 1;
	int p_skips = skips+1;

	const float lookahead = 5000, p_lookbehind = -18000, energyTH = 350.0, promptTH = 511;
	auto time_Vec = *(df.Take<std::vector<float>>("time"));
	auto energy = *(df.Take<std::vector<float>>("energy"));
	auto window_num = *(df.Take<int>("timeWindowNumber"));
	
	int event_num = 0, coinc_num = 0, lw;
	int cw, wPrev = window_num[0], wNum = 1, wStartInd = 0;
	int save_pos[] = {-1, -1};
	ofstream stats;
	
	//vector of vectors for tagging paired hits
	std::vector<std::vector<int>> paired = {};
	for(int ind = 0; ind < time_Vec.size(); ind++){
		std::vector<int> v(time_Vec[ind].size(), 0);
		paired.push_back(v);
	}
	
	if(skips == 2) {	
		stats.open("DTW_stats.txt");
		stats << "Pair   |   TWindow   |   HitTime   |   Energy\n";
	}
	
	//Re-mapping of window numbers to avoid gaps
	for(int wInd = 0; wInd < window_num.size(); wInd++){
		if(window_num[wInd] == wPrev) {
			window_num[wInd] = wNum;
		} else {
			wPrev = window_num[wInd];
			wNum++;
			window_num[wInd] = wNum;
		}
	}
	lw = window_num.back();
	
	for(int i = 0; i < time_Vec.size(); i++){
		save_pos[0] = -1;
		save_pos[1] = -1;
		
		//Tracking first index of current window, so we can check entire window for prompt matching
		if( window_num[i] > window_num[wStartInd] ){
			wStartInd = i;
		}
		
		for(int j = 0; j < time_Vec[i].size(); j++){
		//Exit when at the second-to-last window; skip paired hits or prompts
			event_num++;
			if(window_num[i] == (lw-p_skips+1)) continue;
			if(paired[i][j] == 1) continue;
			if(energy[i][j] > energyTH) continue;
				
			cw = window_num[i];
				
			//look for match
			for(int k = (i+1); k < time_Vec.size(); k++){
				//if(paired[i][j] == 1) break;
				if(window_num[k] > (cw+skips)) break;
				if(window_num[k] < (cw+skips)) continue;
				
				for(int l = 0; l < time_Vec[k].size(); l++){
					//if(paired[i][j] == 1) break;
					if(paired[k][l] == 1) continue;
					if(energy[k][l] > energyTH) continue;
					if( TMath::Abs(time_Vec[k][l] - time_Vec[i][j]) < lookahead){
						save_pos[0] = k;
						save_pos[1] = l;
						
						if(skips == 2){
							stats << coinc_num << "   |   " << cw << "   |   " << time_Vec[i][j] << "   |   " << energy[i][j] << "\n";
							stats << coinc_num << "   |   " << window_num[k] << "   |   " << time_Vec[k][l] << "   |   " << energy[k][l] << "\n";
						}
					}
				}
			}
			//if no match - skip to next iteration
			if( save_pos[0] == -1 ) continue;
			
			//look for prompt
			for(int k = (i+1); k < time_Vec.size(); k++){
				//if(paired[i][j] == 1) break;
				if(window_num[k] > (cw+p_skips)) break;
				if(window_num[k] < (cw+p_skips)) continue;
				
				for(int l = 0; l < time_Vec[k].size(); l++){
					//if(paired[i][j] == 1) break;
					if(paired[k][l] == 1) continue;
					if(energy[k][l] <= promptTH) continue;
					
					//Prompt can be 5 ns after or 18 ns before reference 511
					if( ((time_Vec[k][l] - time_Vec[i][j]) < lookahead)  &&
					     	((time_Vec[k][l] - time_Vec[i][j]) > p_lookbehind) ){
					    
						coinc_num+=1;
						paired[i][j] = 1;
						paired[k][l] = 1;
						paired[ save_pos[0] ][ save_pos[1] ] = 1;
					}
				}
			}
		}
	}
	std::cout << "All hits: " << event_num << std::endl;
	std::cout << "Randoms (type 1): " << coinc_num << std::endl;
	std::cout << "Percentage: " << 100.0*coinc_num/event_num << "%" << std::endl;
	
	if(skips == 2) stats.close();
}


//Only prompt is delayed --> gives randoms of types I and III
void DTW_type4(RNode rnode, int skips){
	auto df = rnode;
	if(skips < 1) skips = 1;

	const float lookahead = 5000, p_lookbehind = -18000, energyTH = 350.0, promptTH = 511;
	auto time_Vec = *(df.Take<std::vector<float>>("time"));
	auto energy = *(df.Take<std::vector<float>>("energy"));
	auto window_num = *(df.Take<int>("timeWindowNumber"));
	
	int event_num = 0, coinc_num = 0, lw;
	int cw, wPrev = window_num[0], wNum = 1;
	int save_pos[] = {-1, -1};
	ofstream stats;
	
	//vector of vectors for tagging paired hits
	std::vector<std::vector<int>> paired = {};
	for(int ind = 0; ind < time_Vec.size(); ind++){
		std::vector<int> v(time_Vec[ind].size(), 0);
		paired.push_back(v);
	}
	
	if(skips == 2) {	
		stats.open("DTW_stats.txt");
		stats << "Pair   |   TWindow   |   HitTime   |   Energy\n";
	}
	
	//Re-mapping of window numbers to avoid gaps
	for(int wInd = 0; wInd < window_num.size(); wInd++){
		if(window_num[wInd] == wPrev) {
			window_num[wInd] = wNum;
		} else {
			wPrev = window_num[wInd];
			wNum++;
			window_num[wInd] = wNum;
		}
	}
	lw = window_num.back();
	
	for(int i = 0; i < time_Vec.size(); i++){
		save_pos[0] = -1;
		save_pos[1] = -1;
		
		for(int j = 0; j < time_Vec[i].size(); j++){
		//Exit when at the second-to-last window; skip paired hits or prompts
			event_num++;
			if(window_num[i] == (lw-skips+1)) continue;
			if(paired[i][j] == 1) continue;
			if(energy[i][j] > energyTH) continue;
				
			cw = window_num[i];
				
			//look for match
			for(int k = i; k < time_Vec.size(); k++){
				//if(paired[i][j] == 1) break;
				if(window_num[k] > cw) break;
				
				for(int l = 0; l < time_Vec[k].size(); l++){
					if( (k == i) && (l == j) ) continue;
					if(paired[k][l] == 1) continue;
					if(energy[k][l] > energyTH) continue;
					if( TMath::Abs(time_Vec[k][l] - time_Vec[i][j]) < lookahead){
						save_pos[0] = k;
						save_pos[1] = l;
						
						if(skips == 2){
							stats << coinc_num << "   |   " << cw << "   |   " << time_Vec[i][j] << "   |   " << energy[i][j] << "\n";
							stats << coinc_num << "   |   " << window_num[k] << "   |   " << time_Vec[k][l] << "   |   " << energy[k][l] << "\n";
						}
					}
				}
			}
			//if no match - skip to next iteration
			if( save_pos[0] == -1 ) continue;
			
			//look for prompt
			for(int k = (i+1); k < time_Vec.size(); k++){
				//if(paired[i][j] == 1) break;
				if(window_num[k] > (cw+skips)) break;
				if(window_num[k] < (cw+skips)) continue;
				
				for(int l = 0; l < time_Vec[k].size(); l++){
					//if(paired[i][j] == 1) break;
					if(paired[k][l] == 1) continue;
					if(energy[k][l] <= promptTH) continue;
					
					//Prompt can be 5 ns after or 18 ns before reference 511
					if( ((time_Vec[k][l] - time_Vec[i][j]) < lookahead)  &&
					     	((time_Vec[k][l] - time_Vec[i][j]) > p_lookbehind) ){
					    
						coinc_num+=1;
						paired[i][j] = 1;
						paired[k][l] = 1;
						paired[ save_pos[0] ][ save_pos[1] ] = 1;
					}
				}
			}
		}
	}
	std::cout << "All hits: " << event_num << std::endl;
	std::cout << "Randoms (type 1): " << coinc_num << std::endl;
	std::cout << "Percentage: " << 100.0*coinc_num/event_num << "%" << std::endl;
	
	if(skips == 2) stats.close();
}

#endif
