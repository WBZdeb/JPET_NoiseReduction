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


vector<float> calc_lifetime(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> df){
	float timeRescale = 1/1000;	//rescale time unit to seconds
	float posRescale = 1;	//rescale position to centimeters
	float c = 30;		//speed of light [cm/ns]
	auto time_Vec = *(df.Take<vector<float>>("time"));
	auto x_Vec = *(df.Take<vector<float>>("x"));
	auto y_Vec = *(df.Take<vector<float>>("y"));
	auto z_Vec = *(df.Take<vector<float>>("z"));
	vector<float> lifetimes;	//vector for calculated lifetimes
	
	float tFl1, tFl2, tEm1, tEm2, refTime, lenAB, origin[3] = {0.0, 0.0, 0.0};
	vector<float> vecAB, pDist, d1, d2;
	
	for (int i = 0; i < time_Vec.size(); i++){
		//determine decay origin
		vecAB = { (x_Vec[i][1] - x_Vec[i][2])*posRescale, (y_Vec[i][1] - y_Vec[i][2])*posRescale, (z_Vec[i][1] - z_Vec[i][2])*posRescale };
		lenAB = TMath::Sqrt( vecAB[0]*vecAB[0] + vecAB[1]*vecAB[1] + vecAB[2]*vecAB[2]);
		
		origin[0] = x_Vec[i][1]*posRescale + (0.5 + c*timeRescale*(time_Vec[i][1]-time_Vec[i][2]) / (2*lenAB) ) * vecAB[0];
		origin[1] = y_Vec[i][1]*posRescale + (0.5 + c*timeRescale*(time_Vec[i][1]-time_Vec[i][2]) / (2*lenAB) ) * vecAB[1];
		origin[2] = z_Vec[i][1]*posRescale + (0.5 + c*timeRescale*(time_Vec[i][1]-time_Vec[i][2]) / (2*lenAB) ) * vecAB[2];
	
		//calculate "machine time" T
		pDist = { x_Vec[i][0]*posRescale - origin[0], y_Vec[i][0]*posRescale - origin[1], z_Vec[i][0]*posRescale - origin[2] };
		refTime = timeRescale*time_Vec[i][0] - TMath::Sqrt( pDist[0]*pDist[0] + pDist[1]*pDist[1] + pDist[2]*pDist[2]) / c;
		
		//calculate lifetime
		pDist = { x_Vec[i][1]*posRescale - origin[0], y_Vec[i][1]*posRescale - origin[1], z_Vec[i][1]*posRescale - origin[2] };
		tFl1 = TMath::Sqrt( pDist[0]*pDist[0] + pDist[1]*pDist[1] + pDist[2]*pDist[2]) / c;
		tEm1 = timeRescale*time_Vec[i][1] - tFl1 - refTime;
		
		pDist = { x_Vec[i][2]*posRescale - origin[0], y_Vec[i][2]*posRescale - origin[1], z_Vec[i][2]*posRescale - origin[2] };
		tFl2 = TMath::Sqrt( pDist[0]*pDist[0] + pDist[1]*pDist[1] + pDist[2]*pDist[2]) / c;
		tEm2 = timeRescale*time_Vec[i][2] - tFl2 - refTime;
		
		lifetimes.push_back( (tEm1 + tEm2)/2 - timeRescale*time_Vec[i][0]);
	}
	
	return lifetimes;
}
