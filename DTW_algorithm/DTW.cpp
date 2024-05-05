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

float get_origin_x(const vector<float>& x, const vector<float>& y, const  vector<float>& z, const vector<float>& time)
{
	float lenAB, origin[3] = {0.0, 0.0, 0.0};

	float c = 30;		//speed of light [cm/ns]
	float timeRescale = 1/1000;
	vector<float> vecAB, pDist, d1, d2;
        //determine decay origin
        vecAB = { (x[1] - x[2]), (y[1] - y[2]), (z[1] - z[2]) };
        lenAB = TMath::Sqrt( vecAB[0]*vecAB[0] + vecAB[1]*vecAB[1] + vecAB[2]*vecAB[2]);
        origin[0] = x[1] + (0.5 + c*timeRescale*(time[1]-time[2]) / (2*lenAB) ) * vecAB[0];
        //origin[1] = y[1] + (0.5 + c*timeRescale*(time[1]-time[2]) / (2*lenAB) ) * vecAB[1];
        //origin[2] = z[1] + (0.5 + c*timeRescale*(time[1]-time[2]) / (2*lenAB) ) * vecAB[2];
        return origin[0];
}


void DTW(){
	gROOT->Reset();
	
	//open and read file with flatTree paths
	ifstream inFile("flatTrees.txt");
	
	if(!inFile.is_open() ){
		cerr << "Error opening the file! Are you sure file 'flatTrees.txt' exists?" << endl;
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
        auto df_vertex = df_fltr.Define("origin_x", get_origin_x, {"x", "y", "z", "time"});
	
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

