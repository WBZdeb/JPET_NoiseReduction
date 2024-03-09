#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;


void MC_coincidences(double activity = 1000000.0, double time = 0.001){

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
	int matched = 1, mismatched = 0, pos = 0;
	
	for(int i = 0; i < sample_size-1; i++){
		//move decay index to the one thats bigger than the creation time we are currently matching
		while(uniform[i] > exp_sorted[pos]) pos++;
	
		//skip if the next time-stamp is a time of creation - can't be matched
		if(uniform[i+1] < exp_sorted[pos]) continue;
		
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
	
	//Last element - we only need to check whether there is such k that:
	//		t_n < t_k' < t_n'
	//
	//where n = number of creations/decays. If it exists, we will have a mismatch; otherwise, a correct match.
	//Since we only test for existence, we can just check k = n-1
	if( exponent[sample_size-2] > uniform[sample_size-1] ){
		mismatched++;
	}
	
	cout << "All matches = " << matched << endl;
	cout << "Mismatched = " << mismatched << endl;
	cout << "Percentage = " << ((double)mismatched)*100 / matched << "%" << endl;
	
	delete rand;
}
