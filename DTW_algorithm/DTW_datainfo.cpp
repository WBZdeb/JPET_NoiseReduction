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

void DTW_datainfo(){
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
	
	//========================
	//	Info about data
	//========================
	
	std::cout << "Acc: " << (*(df.Filter("isAcc").Take<int>("eventNumber"))).size() << std::endl;
	std::cout << "All: " << (*(df.Take<int>("eventNumber"))).size() << std::endl;
	
	// Ratio of Acc/All
	std::cout << 100.0*(*(df.Filter("isAcc").Take<int>("eventNumber"))).size() / (*(df.Take<int>("eventNumber"))).size() << std::endl;
	
	// Ratio of Acc/All with only accidental pairings of 511
	std::cout << 100.0*(*(df.Filter("isAcc && (numberOfHits == 2) && (containsPrompt)").Take<int>("eventNumber"))).size() / (*(df.Take<int>("eventNumber"))).size() << std::endl;
	
	//Random coincidences (type 1)
	for (int i = 2; i < 6; i++){
		DTW_type1(df, i);
	}
	
	//Registered source activity
	auto window_num = *(df.Take<int>("timeWindowNumber"));
	auto window_count = window_num.back();
	double runtime = 0.00005 * window_count;	//[s]
	std::cout << "Registered source activity: "	<< 1022681.0/runtime << " Bq" << std::endl;
	
}

int main()
{
  DTW_datainfo();
  return 0;  
}

