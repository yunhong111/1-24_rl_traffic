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
#include "otherFun.h"
using namespace std;

class CuckooFilter
{

// constructor
public:
    CuckooFilter();
    void cuckooFilterInit(long m,int f,int bc,long MaxNumKicks);
    ~CuckooFilter();

// private variabeles
private:


    // Hash table storage
    int    mf;                   // fingerprint bits number
    int    mk;                  // current number of key-value pairs

    // Hash table parameters
    int    mm;                  // number of buckets
    int    mbc;                  // number of cells per bucket

    // List pointers
    //vector<vector<string> > mL0;
    //const char *mL0[100000][4];
    vector<vector<int> > mL;
    vector<vector<int> > maction;
    //vector<vector<size_t> > mcount;
    //vector<vector<size_t> > mcount0;
    //vector<vector<size_t> > mcountdiff;

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
public:



// private method
private:
    long HashCode(const string& key);
    long HashCodeXor(const string& key);

    long Fingerprint(const string& key);

    bool InitHashFunction();
// public method
public:

    bool AddKey(const string& key,int action);

    bool AddKeyCount(const string& key,int action, char mL0[][ACTIONSIZE][20]);

    bool LookUpKey(const string& key);

    int LookUpKeyCol(const string& key,int & action_out);

    int LookUpKeyActions(const string& key,vector<int>& action);

    int LookUpKeyActionsCount(const string& key,vector<int>& action, size_t& flowNo
                              ,char mL0[][ACTIONSIZE][20], vector<vector<size_t> >& mcount, vector<vector<size_t> >& mcount0,
                            vector<vector<size_t> >& mcountdiff, long& hBuck, int& slot_i);

    int LookUpKeyBack(const string& key, long& posx, int posy);

    bool RemoveKey(const string& key);

    bool RemovePos(long& xPos, int& yPos);

    void returnKey(vector<string>& keys, vector<int>& keyActions,
                              char mL0[][ACTIONSIZE][20], vector<vector<size_t> >& mcount, vector<vector<size_t> >& mcount0,
                            vector<vector<size_t> >& mcountdiff,vector<size_t>& keyCounts, vector<size_t>& keyCounts0,
                            vector<size_t>& keyCountDiffs);

    void ClearTable();
};
