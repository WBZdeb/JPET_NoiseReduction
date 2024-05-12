#include <iostream>
#include <vector>
#include <cmath>


using namespace std;

vector<float> calculateAnnihilationPoint(const vector<float>& hitA, const vector<float>& hitB, float tof)
{
  const float kLightVelocity_cm_ps = 0.0299792458;
  vector<float> middleOfLOR = {0.5 *(hitB[0] + hitA[0]),0.5 *( hitB[1] + hitA[1]),0.5 *( hitB[2] + hitA[2])};
  float len = sqrt((hitB[0] - hitA[0])*(hitB[0] - hitA[0]) + (hitB[1] - hitA[1])*(hitB[1] - hitA[1]) + (hitB[2] - hitA[2])*(hitB[2] - hitA[2]));
  vector<float> versorOnLOR = {(hitB[0] - hitA[0])/len, (hitB[1] - hitA[1])/len, (hitB[2] - hitA[2])/len};

  float shift = 0.5 * tof  * kLightVelocity_cm_ps;
  vector<float> annihilationPoint = {middleOfLOR[0] + shift * versorOnLOR[0],
                             middleOfLOR[1] + shift * versorOnLOR[1],
                             middleOfLOR[2] + shift * versorOnLOR[2]};
  return annihilationPoint;
}

vector<float> get_origin_x(const vector<float>& hitA, const vector<float>& hitB, const vector<float>& time)
{
	float lenAB;

	float c = 30;		//speed of light [cm/ns]
	float timeRescale = 1/1000;
	vector<float> vecAB, origin = {0.0, 0.0, 0.0};
	
        //determine decay origin
        vecAB = { (hitA[0] - hitB[0]), (hitA[1] - hitB[1]), (hitA[2] - hitB[2]) };
        lenAB = sqrt( vecAB[0]*vecAB[0] + vecAB[1]*vecAB[1] + vecAB[2]*vecAB[2]);
        origin[0] = hitB[0] + (0.5 + c*timeRescale*(time[1]-time[2]) / (2*lenAB) ) * vecAB[0];
        origin[1] = hitB[1] + (0.5 + c*timeRescale*(time[1]-time[2]) / (2*lenAB) ) * vecAB[1];
        origin[2] = hitB[2] + (0.5 + c*timeRescale*(time[1]-time[2]) / (2*lenAB) ) * vecAB[2];
   return origin;
}

int main(){
	vector<float> hitA = {-27.10531, -50.71047, -13.61408};
	vector<float> hitB = {30.363901, -48.82912, 8.6609554};
	vector<float> time = {-30558104, -30557534, -30557468};
	
	vector<float> p1 = calculateAnnihilationPoint(hitA, hitB, (time[1]-time[2])/2.0);
	vector<float> p2 = get_origin_x(hitA, hitB, time);
	
	std::cout << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
	std::cout << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
	
	return 0;
}
