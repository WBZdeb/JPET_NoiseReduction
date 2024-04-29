#include <TROOT.h>
#include <TChain.h>
#ifdef R__HAS_VDT
#undef R__HAS_VDT
#endif
#include "ROOT/RDataFrame.hxx"
#ifdef R__HAS_VDT
#undef R__HAS_VDT
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <cassert>
#include <memory>

#include "calc_lifetime.h"
using namespace std;


void DTW(){
	gROOT->Reset();
	
	//open and read file with flatTree paths
	ifstream inFile("flatTrees.txt");
	
	if(!inFile.is_open() ){
		cerr << "Error opening the file! Are you sure file 'flatTrees.txt' exists?" << endl;
		return 1;
	}	
	
	std::string filePath;
	std::string treeName = "FlatTree";
	
	TChain chain(treeName.c_str());
	
	while(getline(inFile, filePath)){
		chain.Add(filePath.c_str());
	}
	
	ROOT::RDataFrame df(chain);
	
	//close file
	inFile.close();
	
	//3-hit Pickoff events, no scatters, has prompt
	auto df_fltr = df.Filter("numberOfHits == 3 && isPickOff && !isScattered && containsPrompt");
	
	//Split between true events and randoms
	auto df_true = df_fltr.Filter("!isAcc");
	auto df_acc = df_fltr.Filter("isAcc");
	
	
	//Calculate and plot lifetimes
	auto lf_True = calc_lifetime(df_true);
	auto lf_Acc = calc_lifetime(df_acc);
	
	//ture
	std::unique_ptr<TCanvas> canv(new TCanvas("canv", "canv", 1920, 1080));
	TH1F *histTrue = new TH1F("lifetime_true","lifetime_true",150,-10,10);
	
	for (int i = 0; i < lf_True.size(); i++){
		histTrue->Fill(lf_True[i]);
	}
	histTrue->GetXaxis()->SetTitle("Lifetime");
	histTrue->GetYaxis()->SetTitle("Count");
	histTrue->Draw();
	canv->SaveAs("lifetime_true.png");
	
	//acc
	std::unique_ptr<TCanvas> canv2(new TCanvas("canv2", "canv2", 1920, 1080));
	TH1F *histAcc = new TH1F("lifetime_acc","lifetime_acc",150,-10,10);
	
	for (int i = 0; i < 20; i++){
		cout << lf_Acc[i] << endl;
	}
	
	for (int i = 0; i < lf_Acc.size(); i++){
		histAcc->Fill(lf_Acc[i]);
	}
	histAcc->GetXaxis()->SetTitle("Lifetime");
	histAcc->GetYaxis()->SetTitle("Count");
	histAcc->Draw();
	canv2->SaveAs("lifetime_acc.png");
	
	
	delete gROOT->FindObject("lifetime_true");
	delete gROOT->FindObject("lifetime_acc");
}

/*
vector<float> calc_lifetime(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> df){
	int timeRescale = 1;	//rescale time unit to seconds
	int c = TMath::C();		//speed of light [m/s]
	auto time_Vec = *(df.Take<vector<float>>("time"));
	auto x_Vec = *(df.Take<vector<float>>("x"));
	auto y_Vec = *(df.Take<vector<float>>("y"));
	auto z_Vec = *(df.Take<vector<float>>("z"));
	vector<float> lifetimes;	//vector for calculated lifetimes
	
	float tFl1, tFl2, tEm1, tEm2, refTime, lenAB, origin[3] = {0.0, 0.0, 0.0};
	vector<float> vecAB, pDist, d1, d2;
	
	for (int i = 0; i < time_Vec.size(); i++){
		//determine decay origin
		vecAB = { x_Vec[i][1] - x_Vec[i][2], y_Vec[i][1] - y_Vec[i][2], z_Vec[i][1] - z_Vec[i][2] };
		lenAB = TMath::Sqrt( vecAB[0]*vecAB[0] + vecAB[1]*vecAB[1] + vecAB[2]*vecAB[2]);
		
		origin[0] = x_Vec[i][1] + (0.5 + c*timeRescale*(time_Vec[i][1]-time_Vec[i][2]) / (2*lenAB) ) * vecAB[0];
		origin[1] = y_Vec[i][1] + (0.5 + c*timeRescale*(time_Vec[i][1]-time_Vec[i][2]) / (2*lenAB) ) * vecAB[1];
		origin[2] = z_Vec[i][1] + (0.5 + c*timeRescale*(time_Vec[i][1]-time_Vec[i][2]) / (2*lenAB) ) * vecAB[2];
	
		//calculate "machine time" T
		pDist = { x_Vec[i][0] - origin[0], y_Vec[i][0] - origin[1], z_Vec[i][0] - origin[2] };
		refTime = timeRescale*time_Vec[i][0] - TMath::Sqrt( pDist[0]*pDist[0] + pDist[1]*pDist[1] + pDist[2]*pDist[2]) / c;
		
		//calculate lifetime
		pDist = { x_Vec[i][1] - origin[0], y_Vec[i][1] - origin[1], z_Vec[i][1] - origin[2] };
		tFl1 = TMath::Sqrt( pDist[0]*pDist[0] + pDist[1]*pDist[1] + pDist[2]*pDist[2]) / c;
		tEm1 = timeRescale*time_Vec[i][1] - tFl1 - refTime;
		
		pDist = { x_Vec[i][2] - origin[0], y_Vec[i][2] - origin[1], z_Vec[i][2] - origin[2] };
		tFl2 = TMath::Sqrt( pDist[0]*pDist[0] + pDist[1]*pDist[1] + pDist[2]*pDist[2]) / c;
		tEm2 = timeRescale*time_Vec[i][2] - tFl2 - refTime;
		
		lifetimes.push_back( (tEm1 + tEm2)/2 - timeRescale*time_Vec[i][0]);
	}
	
	return lifetimes;
}*/

