/**
otherFun.h

create by: Yunhong
create time: 05/22/2015
update time: 05/25/2015
*/
#include <stdio.h>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <vector>
#include <sys/time.h>
#include <stdint.h>
#include <arpa/inet.h>
#include <math.h>
#include <iterator>
# include <omp.h>
#include "hash_function.h"



using namespace std;

extern  double BIGNONKEYTHLD;// = 0.00000769230769230769;//0.00000169558403281403;//0.00000070483398121908;//0.00000140966796243816;//0.00117;100 packets
extern  double BIGKEYTHLD;// = 0.8;
extern  double KEYTHLD;// = 1;
extern  double GLOBAL_INITKEYTHLD;// = 1;
extern  int GLOBAL_BIGNONKEYNUM ;//= 3;
extern  double gFPR;// = 0.01;
extern  double nonkey_weight_sum; // nonkey packet ratio
extern  double key_weight_sum; // nonkey packet ratio
extern double para_k;    // control parameters
extern double diff_FPR; // FPR difference for the first iteration
extern double g_ebsr; // control parameters, if fabs(gFPR - expFPR)< g_ebsr, stop
extern double g_blackkey_weight_sum;
extern ofstream outfile;
extern size_t GLOBAL_KEY_SIZE;
extern size_t FLOW_EST_COUNT;
extern size_t FLOW_EST_SIZE;
extern size_t CUCKOO_BLACK_SIZE;
extern int AGGR_PREFIX;

extern string BLACKFILENAME;
extern string AGGRFILENAME;

extern vector<double> g_vweightThld;
extern vector<double> g_vnonkey_weight_sum; // nonkey packet ratio
extern vector<double> g_vkey_weight_sum; // nonkey packet ratio
extern vector<double> g_vnonkey_weight_sum0; // initial nonkey packet ratio
extern vector<double> g_vgFPR; // nonkey packet ratio
extern vector<double> g_vpara_k; // nonkey packet ratio
extern vector<long> g_vcountkey;
extern vector<long> g_vcountblackkey;
extern vector<long> g_vcountkey_init; //initial value
extern vector<double> g_vblackkey_weight_sum;

const int ACTIONSIZE = 2;

enum keyType
{
    iskey = 1,
    isnonkey = 2,
    isnonkeyprefix = 3,
    isaggregatekey = 4,
    isblackkey = 5,
    isbignonkey = 6,
    invalid = 0,
    invalidkey = 7,
    cannotaggr = 8
};

struct SubIPv
{

    vector<string> subIPUint;
};

struct ActionOrder
{

    vector<uint16_t> aTrieOder;
};


long print_elapsed(const char* desc, struct timeval* start, struct timeval* end, int niters);

bool IsPrime(int number);

std::string num2str(int num);

int str2num(string str);

void readCSV( FILE *fp, std::vector<std::string>& vls, std::vector<string>& vcnt );

void readCSVKey( FILE *fp, std::vector<std::string>& vls, std::vector<string>& vcnt );

uint32_t parseIPV4string(const char* ipAddress);

int char_to_ipv4(char *src);

string DecToBin(uint64_t dec);

string parseIPV42bin(const char* ipAddress);

uint32_t binConv(string binNum);

string parsebin2IPV4(const char* ipAddress);

string parsedec2IPV4(uint64_t dec);

void readCSVKeySlash( FILE *fp, std::vector<std::string>& vls);

void ctlPrefix(vector<string> &key, vector<string>& key_cnt,vector<string> &keyprefix,vector<keyType> &keytype,vector<string> &action,vector<string> &prefix,
               vector<string> &key_out,vector<string>& key_cnt_out,vector<string> &keyprefix_out,vector<keyType> &keytype_out,vector<string> &action_out );

void randKey(long& n, long& ip_num, vector<string>& flow, vector<string>& flow_cnt, vector<string> &key, vector<string> &key_cnt,vector<keyType> keytype);

void readCSVcnt( FILE *fp, std::vector<std::string>& vls, std::vector<std::string>& vpr, std::vector<string>& vcnt );

void readCSVaction( FILE *fp, std::vector<std::string>& vls, std::vector<std::string>& vpr, std::vector<string>& action, std::vector<string>& vcnt );

void readCSVactioncnt( FILE *fp, std::vector<std::string>& vls, std::vector<std::string>& vpr, std::vector<string>& action, std::vector<string>& keytype, std::vector<string>& vcnt );

float ctlPara(double givenFPR, double expFPR, float threshold_old, float threshold_old_old, int action);

void clearglobal_value();

void setGlobalValue(int actionSize);

double gaussrand(double miu, double sigma);

