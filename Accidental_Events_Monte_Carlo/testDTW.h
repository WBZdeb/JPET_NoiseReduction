#ifndef TESTDTW_H
#define TESTDTW_H

#include <vector>
#include <string>
#include <ROOT/RDataFrame.hxx>

#define WIN_LEN -50000000.0

class gammaP {
private:
    double time;
    int timeWindowNum;
    int eventNum;
    bool bIsPrompt;
    bool paired;

public:
    gammaP(double time, int timeWindowNum, int eventNum, bool isPrompt = false);

    double getTime() const;
    int getTimeWindowNum() const;
    int getEventNum() const;
    bool isPrompt() const;
    bool isPaired() const;
    void setPaired(bool val);
};

bool compareGammaP(const gammaP& a, const gammaP& b);

int findRandoms(std::vector<gammaP>& hitsVec, std::vector<std::vector<gammaP>>* pairs = nullptr);

int findTrues(std::vector<gammaP>& hitsVec, std::vector<std::vector<gammaP>>* pairs = nullptr);

int findCoincidences(std::vector<gammaP>& hitsVec, std::vector<std::vector<gammaP>>* pairs = nullptr);

std::vector<double> calcPromptIntervals(std::vector<gammaP>& hitsVec, int window);

std::vector<double> calcGammaIntervals(std::vector<gammaP>& hitsVec, int window);

std::vector<gammaP> generateEvents(int windowCount, double activity);

std::string generate_DataFrame(int window_count, double activity);

void testDTW(int window_count = 10, double activity = 700000.0);

#endif // TESTDTW_H

