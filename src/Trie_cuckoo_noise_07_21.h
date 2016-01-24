/**
trie_cuckoo_noise_07_21.cpp
create by: Yunhong
create time: 07/21/2015
*/

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <string>
#include <vector>
#include <sys/time.h>
#include <arpa/inet.h>
#include "aggregation_add_cuckoo.h"
#include "RL.h"
# include <time.h>
#include <unistd.h>

using namespace std;

typedef vector<vector<int> > intss;

const int SWITCHNUM = 3;

struct BkInfo{
    vector<string> bks;
    vector<int> bkPres;
    vector<int> bkActions;
    vector<size_t> bkSizes;
};

struct UPrefix{
    vector<int> uPres;
};

typedef vector<UPrefix> VUPrefix;

typedef vector<BkInfo> VBkInfo;

void feedbackBlackkey(vector<string>& overBigKeys);

void feedbackBlackkey1(strings& overBigKeys);

bool readFile0(ifstream& infile, vector<string> &flow, vector<size_t> &flow_cnt, size_t readNum, bool& isEndFlag, size_t& i_seq);

void initRLearn(RLearn* rLearn);

void updateBlacklist(vector<string>& overBigKeys, vector<int>& overActions, RLearn* rLearn, int actionSeq
, vector<string>& blackkeyPres, vector<int>& blackActionPres, ofstream& blackKeyFileOut,  size_t& slotNum, CuckooTable& cuckooBlackKeyTable);

void feedbackBlackkeyRL(VBkInfo& vBkInfo, RLearn** rLearn[],
int actionSize, int switchNum, size_tss& slotNums, intss& flowActionUnique, size_t line, CuckooTable* cuckooBlackKeyTable);

void printQList(RLearn* rLearn);

void selectAction(RLearn** rLearn[], int actionSize, int switchNum, size_tss& slotNums);

void findMax(vector<QSum>& qSums, size_ts& slotNums, int actionSize, int switchNum );

void printVec(vector<size_t>& vec);

void loadKeys2Filter(string& inFileName, vector<size_t>& mask, VUPrefix& vuniquePrefix,
VUPrefix& vuniqueAggPrefix, char mL0[][ACTIONSIZE][20], char * argv[], int& finger, int& finger0,int switchNum,
int actionSize, intss& uniqueActs,
                      CuckooTable* cuckooTableKey, CuckooTable* cuckooBlackKeyTable, CuckooTable* cuckooAggrKeyTable,
                      CuckooFilter* cuckooFilter, CuckooFilter* cuckooFilterInit0);

double nextTime(double rateParameter);


