/**
cuckoo_filter.h

create by: Yunhong
create time: 05/22/2015
update time: 05/25/2015
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <string>
#include <set>
#include <vector>
#include <math.h>
#include <sys/time.h>
#include "cuckoo_filter.h"
using namespace std;

class CuckooTable
{

// constructor
public:
    CuckooTable();
    void CuckooTableInit(long m,int f,int bc,long MaxNumKicks);
    ~CuckooTable();

// private variabeles
//private:
private:


    // Hash table storage
    int    mf;                   // fingerprint bits number
    int    mk;                  // current number of key-value pairs

    // Hash table parameters
    int    mm;                  // number of buckets
    int    mbc;                  // number of cells per bucket

    // List pointers
    vector<vector<uint32_t> > mL;
    vector<vector<int> > maction;
    vector<vector<int> > mprefix;
    vector<vector<size_t> > mcount;

    // Hash function parameters
    long    mp;                  // hash function parameter
    long    ma;                  // hash function parameter
    long    mb;                  // hash function parameter
    long    mc;                  // hash function parameter


    // Hash function parameters
    long    mpf;                  // hash function parameter
    long    maf;                  // hash function parameter
    long    mbf;                  // hash function parameter
    long    mcf;                  // hash function parameter

    // maximum number of kicking out an elt
    long    mMaxNumKicks;

    // high mask and low mask, used for check whether h2 exceeds range
    vector<vector<bool> > mhigh;


// public variables
private:

// private method
//private:
    long HashCode(const string& key);

    long Fingerprint(const string& key);

    bool InitHashFunction();
// public method
public:

    bool AddKey(const string& key,int action);

    bool AddKeyPrefix(const string& key,int prefix, int action);

    bool LookUpKey(const string& key,int &prefix);

    int LookUpKeyCol(const string& key,int & iaction);

    bool LookUpKeyAction(const string& key, int &prefix, int& action);

    bool LookUpKeyActionCount(const string& key, int &prefix, int& action, size_t& keyNo);

    bool RemoveKey(const string& key);

    void ClearTable();

    void clearCount();

    void returnKey(vector<string>& keys, vector<size_t>& keyNos,
                            vector<int>& keyPrefixes, vector<int>& keyActions);

    float occupyRate();
};
