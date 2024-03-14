#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;


void MC_coincidences(double activity = 1000000.0, double time = 0.001, double timeFrame = 0.0){

	TRandom2 *rand = new TRandom2(0);
	int sample_size = activity*time;	
	double tau = 142;	//measured in [ns] 
	
	vector<double> uniform, exponent, exp_sorted;
	
	//generate time (in [ms]) of decay into ortho-positronium using uniform distribution
	for(int i = 0; i < sample_size; i++){
		uniform.push_back(rand->Uniform(time*1000));	//convert time to [ms]
	}
	
	//sorting generated decays in ascending order
	sort(uniform.begin(), uniform.end());
	
	//for each decay, generate time of ortho-positronium's decay 
	//using exponential distribution (in [ms]) 
	for(int i = 0; i < sample_size; i++){
		exponent.push_back( uniform[i] + rand->Exp(tau/1000000) );	//converting tau to [ms]
	}
	
	exp_sorted = exponent;
	sort(exp_sorted.begin(), exp_sorted.end());
	
	//calculate percentage of mismatched creation-decay times
	int matched = 0, mismatched = 0, pos = 0;
	
	for(int i = 0; i < sample_size-1; i++){
		//move decay index to the one thats bigger than the creation time we are currently matching
		while(uniform[i] > exp_sorted[pos]) pos++;
	
		//skip if the next time-stamp is a time of creation - can't be matched
		if(uniform[i+1] < exp_sorted[pos]) continue;
		
		//skip if the match doesn't satisfy timeFrame parameter
		if( (timeFrame > 0.0) && (exp_sorted[pos] - uniform[i] > timeFrame) ) continue;
		
		// if we have a mismatch - increment counter by 1
		if( exp_sorted[i] != exponent[pos] ){ 
			mismatched++;
			cout << "Mismatch between t" << i << "=" << uniform[i] << " and t" << find(exponent.begin(), exponent.end(), exp_sorted[pos]) - exponent.begin();
			cout << "'=" << exp_sorted[pos]  << endl;
			cout << "Correct match: " << uniform[i] << " and " << exponent[i] << endl << endl;
		}
		//on both match and mismatch increment all matches and index position by 1
		pos++;
		matched++;
	}
	
	//handle last creation time matching
	while(uniform[sample_size-1] > exp_sorted[pos]) pos++;
	
	if( (timeFrame > 0.0) && (exp_sorted[pos] - uniform[sample_size-1] < timeFrame) ){
		matched++;
		if( exp_sorted[sample_size-1] != exponent[pos]) {
			mismatched++;		
			cout << "Mismatch between t" << sample_size-1 << "=" << uniform[sample_size-1] << " and t";
			cout << find(exponent.begin(), exponent.end(), exp_sorted[pos]) - exponent.begin() << "'=" << exp_sorted[pos] << endl;
			cout << "Correct match: " << uniform[sample_size-1] << " and " << exponent[sample_size-1] << endl << endl;	
		}
	}
	
	if( timeFrame <= 0.0 ){
		matched++;
		if( exp_sorted[sample_size-1] != exponent[pos]) {
			mismatched++;		
			cout << "Mismatch between t" << sample_size-1 << "=" << uniform[sample_size-1] << " and t";
			cout << find(exponent.begin(), exponent.end(), exp_sorted[pos]) - exponent.begin() << "'=" << exp_sorted[pos] << endl;
			cout << "Correct match: " << uniform[sample_size-1] << " and " << exponent[sample_size-1] << endl << endl;	
		}
	}
	
	
	cout << "All matches = " << matched << endl;
	cout << "Mismatched = " << mismatched << endl;
	cout << "Percentage = " << ((double)mismatched)*100 / matched << "%" << endl;
	
	delete rand;
}
