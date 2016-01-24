/**
aggregation_add_cuckoo.cpp
create by: Yunhong
create time: 05/25/2015
*/


#include "aggregation_add_cuckoo.h"

/*CuckooFilter cuckooFilter2;  // cuckoo filter for aggregation
CuckooFilter cuckooFilterInit02; // cuckoo filter for non aggregation
CuckooTable cuckooBlackKeyTable2;    // cuckoo table for black key
CuckooTable cuckooTableKey2; // cuckoo table for non aggregation

CuckooFilter cuckooFilter3;  // cuckoo filter for aggregation
CuckooFilter cuckooFilterInit03; // cuckoo filter for non aggregation
CuckooTable cuckooBlackKeyTable3;    // cuckoo table for black key
CuckooTable cuckooTableKey3; // cuckoo table for non aggregation*/

// descend
struct CmpMass
{
    CmpMass(vector<size_t>& vec) : values(vec) {}
    bool operator() (const size_t& a, const size_t& b) const
    {
        return (values[a]>values[b]);
    }
    vector<size_t>& values;
};


// descend
struct CmpMassInt
{
    CmpMassInt(vector<int>& vec) : values(vec) {}
    bool operator() (const size_t& a, const size_t& b) const
    {
        return (values[a]<values[b]);
    }
    vector<int>& values;
};


int initCuckoo(vector<string> &keys,vector<int> &keyPrefixes,
                vector<int> &keyActions,float &storage, int& finger, char mL0[][ACTIONSIZE][20], CuckooFilter& cuck, CuckooFilter& cuck0)
{
    //-----------------------------------------------------------
    // Parameters for cuckoo filter
    float loadFactor = 0.90;
    int slotNo = 4;
    size_t keySize = keys.size();
    size_t bucketNo = size_t(keySize/(loadFactor*slotNo))+1;
    int fingerprint = storage*1024.0*8.0*loadFactor/(keys.size())-3;
    finger = fingerprint;
    long maxNumKicks = 1500;
    cout<<"* Fingerprint length: "<<fingerprint<<endl;
    //mL0 =vector<vector<string> > (bucketNo, vector<string>(slotNo, " "));

    // --------------------------------------------------------
    //init cuckoo filter
    cuck.ClearTable();
    cuck.cuckooFilterInit(bucketNo,fingerprint,slotNo,maxNumKicks);

    cuck0.ClearTable();
    cuck0.cuckooFilterInit(bucketNo,fingerprint,slotNo,
                                       maxNumKicks);

    // add keys to cuckoo filer
    addCuckooFilter(keys, keyPrefixes, keyActions, cuck, cuck0);

    return 1;

}
/**
The first aggregation
*/
size_t initAggregation(vector<string> &keyIns,vector<int> &keyPrefixIns,
                        vector<int> &keyActionIns, vector<size_t> &maskes,
                      int actionSize, float &storage, bool isInit,
                      int &fingerprintOld,vector<int> &uniqueAggKeyprefixes,char mL0[][ACTIONSIZE][20],
                      CuckooFilter& cuck, CuckooTable& cuckAggr, ints& flowactionunique)
{
    // ----------------------------------------
    // Get indexes
    vector<size_t> indexes;
    //vector<int> flowactionunique;
    clusterAction(keyIns, keyActionIns, maskes,indexes, flowactionunique);

    // ----------------------------------------------------------
    // Define tries
    uint16_t trieNum = indexes.size();   // trie number
    Trie *bTrie;            // define trees
    bTrie = new Trie[trieNum];            // define trees

    // ---------------------------------------------------------
    /* Insert keys to trie */
    insertWordTrieSimple(bTrie, trieNum, indexes, keyIns, keyPrefixIns,
                         keyActionIns, flowactionunique);

    // --------------------------------------------------------
    // classify tree according to action
    ActionOrder  actionOrder[actionSize];
    cout<<"* ActionOder clustering ..."<<endl;
    for(int ai = 0; ai < actionSize; ai++)
    {
        actionOrder[ai].aTrieOder.clear();
    }

    for(int ai = 0; ai < actionSize; ai++)
    {
        for (int ti = 0; ti < indexes.size(); ti++)
        {
            if (bTrie[ti].maction == flowactionunique[ai])
            {
                actionOrder[ai].aTrieOder.push_back(ti);
            }
        }

    }

    // --------------------------------------------------
    int cuckooBlackSize = CUCKOO_BLACK_SIZE;
    float loadFactor = 0.90;

    // -------------------------------------------
    size_t countKey = 0;
    size_t countAggregateKey = 0;
    size_t countBlackKey =0;
    size_t countOriKey =0;

    vector<string> keys;
    vector<int> keyActions;
    vector<string> blackKeys;
    vector<int> blackkeyPrefixes;
    vector<string> aggregateKeys;

    g_vcountkey.assign(actionSize,0);
    g_vcountblackkey.assign(actionSize,0);
    g_vweightThld.assign(actionSize,0);


    for(int ai = 0; ai < actionSize; ai++ )
    {
        g_vweightThld[ai] =  0.0001+ 0.0001/*ai*/;
    }

    countKey = keyIns.size();
    int fingerprint = (storage*1024.0f*8.0f-(cuckooBlackSize)*39/0.9-FLOW_EST_SIZE*17/0.9)*loadFactor/(countKey) -3;
    int fingerprintTarget = 15;

    // no need to compress
    if(fingerprint+2>=fingerprintTarget)
    {
        fingerprintTarget = fingerprint+2;
    }
    cout<<" Fingerprint target: "<<fingerprintTarget<<endl;
    if(fingerprint>14)
    {
        cout<<" No aggregation needed!"<<endl;
        return 0;
    }
       // return 0;
    int aggrprefixlength = 19;
    while(fingerprint<fingerprintTarget)
    {
        //fingerprintTarget --;
        aggrprefixlength --;

        cout<<"prefixlength "<<aggrprefixlength<<endl;
        int iteartion = 0;
        int maxIteration = 70;
        for(int ai = 0; ai < actionSize; ai++ )
        {
            g_vweightThld[ai] =  0.000004+0.000001;

            if(flowactionunique[ai] == 1)
            {
                g_vweightThld[ai] =  0.00001+0.000001;
            }
        }
        while(fingerprint!=fingerprintTarget && iteartion<maxIteration)
        {
            iteartion++;

            // -----------------------------------------
            // aggregate all other keys
            /* Aggregation */
            /* Init variables */
            countKey = 0;
            countAggregateKey = 0;
            countBlackKey =0;
            countOriKey =0;

            keys.clear();
            keyActions.clear();
            blackKeys.clear();
            blackkeyPrefixes.clear();
            aggregateKeys.clear();

            for(int ai = 0; ai < actionSize; ai++)
            {
                g_vcountkey[ai] = 0;
                g_vcountblackkey[ai] = 0;
                if(fingerprint>fingerprintTarget && iteartion!=1) // decompree, increase threshold
                    g_vweightThld[ai] +=  0.1*g_vweightThld[ai]*(fingerprint - fingerprintTarget)*(1) ;
                else if( g_vweightThld[ai] > 0&& iteartion!=1)
                {
                    g_vweightThld[ai] +=  0.01*g_vweightThld[ai]*(fingerprint - fingerprintTarget)*(1) ;
                    if(g_vweightThld[ai]<pow(10,-6.5))
                        g_vweightThld[ai] = pow(10,-6.5);

                }

                // ----------------------------------------------------------
                /* aggregate Trie */
                aggregateTrie(bTrie, ai,  actionOrder, countKey,countAggregateKey,
                              countBlackKey,countOriKey, keys, keyActions,blackKeys,
                              blackkeyPrefixes, aggregateKeys,aggrprefixlength, 1, isInit);
                cout<<"Action: "<<ai<<" threshold: "<<g_vweightThld[ai]<<" ";

            }

            cout<<endl;
            cout<<"* orikey  "<<countOriKey<<" agg  "<<
                countAggregateKey<<" blkey  "<<countBlackKey<<" countKey  "<<countKey<<endl;
            fingerprint = (storage*1024.0f*8.0f-(cuckooBlackSize)*39/0.9-FLOW_EST_SIZE*17/0.9)*loadFactor/(countKey) -3;
            cout<<"Fingerprint: "<<fingerprint<<endl;
            bool isbreakFlag = 0;
            for(int ai = 0; ai < actionSize; ai++)
            if(g_vweightThld[ai] <= pow(10,-6.5) && fingerprint < fingerprintTarget)
            {
                isbreakFlag = 1;
                break;
            }
            if(isbreakFlag)
            {
                break;
            }
        }

    }
    AGGR_PREFIX = aggrprefixlength;
    // --------------------------------------
    /* Insert to cuckoo filter */
    // parameters for cuckoo filter

    int slotNo = 4;
    size_t keySize = keys.size();
    long bucketNo = long(keySize/(loadFactor*slotNo))+1;

    fingerprint = (storage*1024.0f*8.0f-(cuckooBlackSize)*39/0.9-FLOW_EST_SIZE*17/0.9)*loadFactor/(countKey) -3;
    fingerprintOld = fingerprint;
    long maxNumKicks = 1000;
    cout<<"* Fingerprint length: "<<fingerprint<<endl;

    //init cuckoo filter
    cuck.ClearTable();
    cuck.cuckooFilterInit(bucketNo,fingerprint,slotNo,maxNumKicks);

    // add key to cuckoo filter
    addCuckooFiltermL(keys, keyActions,mL0, cuck);

    // ----------------------------------
    // test cuckoo filter
    /*long keyCount = 0;
    for(size_t i = 0; i < keyIns.size(); i++)
    {
        for(int mi = 0; mi <= 32 -8; mi++)
        {
                size_t ip = parseIPV4string(keyIns[i].c_str());
                size_t subIP = ip & maskes[mi];
                string flowstr = parsedec2IPV4(subIP)+"/"+num2str((mi+8));
                int flag_look = cuckooFilter.LookUpKey(flowstr);
                if (flag_look != 0)
                {
                    keyCount += 1;
                    break;
                }
        }
    }
    cout<<"* Test key in cuckoo filters: "<<keyCount<<endl;*/

    // ----------------------------
    // Compute unique prefix length
    vector<int> keyPrefixes;
    for(size_t i = 0; i< keys.size(); i++)
    {
        size_t found = keys[i].find('/');
        string str = keys[i].substr(0,found);
        string prefixStr = keys[i].substr(found+1, keys[i].size()-found);
        keyPrefixes.push_back(str2num(prefixStr));
    }

    prefixNum(keyPrefixes,uniqueAggKeyprefixes);

    // ----------------------------------------
    // Add aggregate keys to cuckooTable
    size_t aggregateKeySize = aggregateKeys.size();
    bucketNo = long(aggregateKeySize/(loadFactor*slotNo));
    cuckAggr.CuckooTableInit(bucketNo,fingerprint,slotNo,maxNumKicks);

    for(size_t i = 0; i < aggregateKeySize; i++)
    {
        size_t found = aggregateKeys[i].find('/');
        string str = aggregateKeys[i].substr(0,found);
        string prefixstr = aggregateKeys[i].substr(found+1,
                                                   aggregateKeys[i].size()-found);

        cuckAggr.AddKey(str,str2num(prefixstr));
    }

    keyPrefixes.clear();
    keys.clear();
    keyActions.clear();

    // --------------------------------------
    // Add blackkey to cuckooTable
    /*cout<<"* Add blackkey to cuckooTable!"<<endl;
    size_t blackKeySize = blackKeys.size();
    size_t bucketSize = int((blackKeySize)/(loadFactor*slotNo))+1;

    int MaxKickoutNum = 1000;
    cuckooBlackKeyTable.ClearTable();
    cuckooBlackKeyTable.CuckooTableInit(bucketSize,fingerprint,slotNo,
                                        MaxKickoutNum);
    for(size_t i = 0; i < blackKeySize; i++)
    {
        cuckooBlackKeyTable.AddKeyPrefix(blackKeys[i],blackkeyPrefixes[i], 4);
    }

    blackKeys.clear();
    blackkeyPrefixes.clear();*/

    // ---------------------------------------
    for(int i = 0; i < trieNum; i++)
    {
        bTrie[i].deleteChild(bTrie[i].root);
    }
    if(!bTrie)
        delete[] bTrie;

    return countKey;

}

/*size_t aggregation(vector<string> &keyIns,vector<int> &keyPrefixIns,
                        vector<int> &keyActionIns, vector<size_t> &maskes,
                      int actionSize, float &storage, bool isInit,
                      int &fingerprintOld,vector<int> &uniqueAggKeyprefixes, strings& overKeys, size_ts& overKeyNos,
                      floats& haoOvers, float target)
{
    // ---------------------------
    vector<size_t> indexes;
    clusterAction(keyIns, keyActionIns, maskes,indexes);

    // -----------------------------
    // find min value and max value for each subtrie
    vector<size_t> IPUpperBound;
    vector<size_t> IPLowerBound;
    IPLowerBound.push_back(parseIPV4string("0.0.0.0"));
    IPUpperBound.push_back(parseIPV4string("0.255.255.255"));
    for(uint16_t i = 0; i < indexes.size()-1; i ++)
    {
        size_t lowerBound = parseIPV4string(keyIns[indexes[i]].c_str()) & maskes[0];
        size_t upperBound = lowerBound + (1<<24)-1;
        IPLowerBound.push_back(lowerBound);
        IPUpperBound.push_back(upperBound);

    }

    // ----------------------------------------------------------
    // Define Tries
    uint16_t trieNum = indexes.size();   // trie number
    Trie *bTrie;            // define trees
    bTrie = new Trie[trieNum];            // define trees

    // ---------------------------------------------------------
    /* Insert keys to trie * /
    insertWordTrieSimple(bTrie, trieNum, indexes, keyIns, keyPrefixIns,
                         keyActionIns);

    // -------------------------------------------------------------
    // classify tries according to actions
    ActionOrder  actionOrder[actionSize];
    cout<<"* ActionOder clustering ..."<<endl;
    for(int ai = 0; ai < actionSize; ai++)
    {
        actionOrder[ai].aTrieOder.clear();
    }

    for(int ai = 0; ai < actionSize; ai++)
    {
        for (int ti = 0; ti < indexes.size(); ti++)
        {
            if (bTrie[ti].maction == ai)
            {
                actionOrder[ai].aTrieOder.push_back(ti);
            }
        }

    }

    // -------------------------------------------------------------------------
    /* Aggregation * /
    /* Init variables * /
    size_t countKey = 0;
    size_t countAggregateKey = 0;
    size_t countBlackKey =0;
    size_t countOriKey =0;

    vector<string> keys;
    vector<int> keyActions;
    vector<string> blackKeys;
    vector<int> blackkeyPrefixes;
    vector<string> aggregateKeys;

    int aggrPrefixlength = AGGR_PREFIX;

    for(int ai = 0; ai < actionSize; ai++)
    {
        g_vcountkey[ai] = 0;
        g_vcountblackkey[ai] = 0;

        // ----------------------------------------------------------
        /* aggregate Trie * /
        aggregateTrie(bTrie, ai,  actionOrder, countKey,countAggregateKey,
                      countBlackKey,countOriKey, keys,
        keyActions,blackKeys,blackkeyPrefixes, aggregateKeys, aggrPrefixlength, 1, isInit);
        cout<<"Action: "<<ai<<" threshold: "<<g_vweightThld[ai]<<" ";

    }
    cout<<endl;
    cout<<"* orikey  "<<countOriKey<<" agg  "<<
        countAggregateKey<<" blkey  "<<countBlackKey<<" countKey  "<<countKey<<endl;

    // ------------------------------------------------
    // Find the actions to compress and decompress
    ints compressActions;
    ints decompressActions;
    compressAct(haoOvers, target, actionSize, compressActions, decompressActions);

    // --------------------------------------------
    // Decompress aggregate keys
    strings aggrIPs;
    ints aggrPrefixes;
    size_t aggrCount = 0;

    cout<<"overbig Keys size: "<<overKeys.size()<<endl;
    for(size_t i = 0; i< overKeys.size(); i++)
    {
        size_t ipInt = parseIPV4string(overKeys[i].c_str());

        // lookup key in aggregate table, if inside, get the aggregate key
        bool isAggregatekey = 0;
        //if(cuckooAggrKeyTable.mm > 1)
        {
            for(int mi = 0; mi < uniqueAggKeyprefixes.size(); mi++)
            {
                size_t subIP = ipInt & maskes[uniqueAggKeyprefixes[mi]-8];
                string flowstr = parsedec2IPV4(subIP);
                int prefix = uniqueAggKeyprefixes[mi];
                isAggregatekey = cuckooAggrKeyTable.LookUpKey(flowstr,prefix);
                if (isAggregatekey)
                {

                    aggrIPs.push_back(flowstr);
                    aggrPrefixes.push_back(prefix);

                    // Get the aggregate keys and decompress them
                    for(uint16_t ti = 0; ti < trieNum; ti++)
                    {
                         if(ipInt >= IPLowerBound[ti] && ipInt < IPUpperBound[ti])
                        {
                            // If the action is to compress, no need to decompress the aggregate keys
                            for(int ai = 0; ai <compressActions.size(); ai++)
                            {
                                if(ti/64 == compressActions[ai])
                                {
                                    break;
                                }
                            }

                            bTrie[ti].searchAggrPrefix(DecToBin(subIP),prefix, aggrCount);
                            break;
                        }
                    }

                    // Find one match and break
                    break;

                }
            }
        }
    }

     // ------------------------------------------------
    // load aggregate previous keys and decompress them
    ifstream aggrFile(AGGRFILENAME.c_str());
    string aggrIPStr;
    int aggrPrefix;
    while(aggrFile>>aggrIPStr>>aggrPrefix && aggrIPs.size()<3500)
    {
        aggrIPs.push_back(aggrIPStr);
        aggrPrefixes.push_back(aggrPrefix);

        size_t aggrIPInt = parseIPV4string(aggrIPStr.c_str());
        // get the aggregate key
        for(uint16_t ti = 0; ti < trieNum; ti++)
        {
            if(aggrIPInt >= IPLowerBound[ti] && aggrIPInt < IPUpperBound[ti])
            {
                /*for(int ai = 0; ai <compressActions.size(); ai++)
                {
                    if(ti/64 == compressActions[ai])
                    {
                        break;
                    }
                }* /

                bTrie[ti].searchAggrPrefix(DecToBin(aggrIPInt),aggrPrefix, aggrCount);
                break;
            }
        }

    }
    aggrFile.clear();
    aggrFile.close();
    vector<size_t>().swap(IPLowerBound);
    vector<size_t>().swap(IPUpperBound);

    //---------------------------------------
    // write to aggr file;
    ofstream aggrFileOut("aggrFile");
    for(size_t i = 0; i < aggrIPs.size(); i++)
    {
        aggrFileOut<<aggrIPs[i]<<" "<<aggrPrefixes[i]<<endl;
    }
    aggrFileOut.clear();
    aggrFileOut.close();
    vector<string>().swap(aggrIPs);
    vector<int>().swap(aggrPrefixes);

    // -------------------------------------------
    // If decompress alters fingerprint length, compress other keys
    countKey += aggrCount;  // key number after decompress
    cout<<"Increase aggrcount: "<<aggrCount<<endl;

    uint16_t cuckooBlackSize = CUCKOO_BLACK_SIZE;
    float loadFactor = 0.90;
    int fingerprint = (storage*1024.0f*8.0f-(cuckooBlackSize)*39/0.9-FLOW_EST_SIZE*17/0.9)*loadFactor/(countKey) -3;

    int iteartion = 0;
    int maxIteration = 12;
    ints prefixlength;
    prefixlength.assign(actionSize,20);
    while(fingerprint!=fingerprintOld && iteartion<maxIteration)
    {
        iteartion++;

        // -----------------------------------------
        // aggregate all other keys
        /* Aggregation * /
        /* Init variables * /
        countKey = 0;
        countAggregateKey = 0;
        countBlackKey =0;
        countOriKey =0;

        keys.clear();
        keyActions.clear();
        blackKeys.clear();
        blackkeyPrefixes.clear();
        aggregateKeys.clear();

        for(int ai = 0; ai < actionSize; ai++)
        {
            g_vcountkey[ai] = 0;
            g_vcountblackkey[ai] = 0;
            if(fingerprint>fingerprintOld) // decompress, increase threshold
            {
                g_vweightThld[ai] +=  0.2*g_vweightThld[ai]*(fingerprint - fingerprintOld);
                prefixlength[ai] ++;
            }

            else if( g_vweightThld[ai] > 0)
            {
                g_vweightThld[ai] +=  0.1*g_vweightThld[ai]*(fingerprint - fingerprintOld)  ;
                prefixlength[ai] --;
            }

            //g_vweightThld[ai] = 0;
            // ----------------------------------------------------------
            /* aggregate Trie * /
            aggregateTrie(bTrie, ai,  actionOrder, countKey,countAggregateKey,
                          countBlackKey,countOriKey, keys,keyActions,blackKeys,blackkeyPrefixes, aggregateKeys, aggrPrefixlength, 1, isInit);
            cout<<"* threshold: "<<g_vweightThld[ai]<<" prefixlength: "<<prefixlength[ai]<<endl;

        }
        cout<<"  coutkey  "<<"agg  "<<"blkey  "<<"orikey  "<<countOriKey<<" "<<
        countAggregateKey<<" "<<countBlackKey<<" "<<countKey<<endl;
        fingerprint = (storage*1024.0f*8.0f*loadFactor-(cuckooBlackSize)*39/0.9-FLOW_EST_SIZE*17/0.9)/(countKey) -3;
        cout<<"Fingerprint: "<<fingerprint<<endl;
    }

    // -----------------------------------------
    // Get aggregate result
    countKey = 0;
    countAggregateKey = 0;
    countBlackKey =0;
    countOriKey =0;

    keys.clear();
    keyActions.clear();
    blackKeys.clear();
    blackkeyPrefixes.clear();
    aggregateKeys.clear();
    vector<char> word;
    size_t recoverCount = 0;

    for(uint16_t ti = 0; ti < trieNum; ti++)
    {
         bTrie[ti].nodeCount(bTrie[ti].root,countKey,countAggregateKey,countBlackKey,countOriKey);
         bTrie[ti].printNode(bTrie[ti].root,word,keys,keyActions,blackKeys,blackkeyPrefixes,aggregateKeys);
         bTrie[ti].recoverTrie(bTrie[ti].root8, recoverCount);
    }
    cout<<"  coutkey  "<<"agg  "<<"blkey  "<<"orikey  "<<countOriKey<<" "<<
    countAggregateKey<<" "<<countBlackKey<<" "<<countKey<<endl;

    // ----------------------------------------------
    // Compute Unique prefixes
    vector<int> keyPrefixes;
    for(size_t i = 0; i < keys.size(); i++)
    {
        size_t found = keys[i].find('/');
        string prefixstr = keys[i].substr(found+1,
                           keys[i].size()-found);

        keyPrefixes.push_back(str2num(prefixstr));
    }

    uniqueAggKeyprefixes.clear();
    prefixNum(keyPrefixes,uniqueAggKeyprefixes);

    // ----------------------------------
    /* Insert to cuckoo filter * /
    // parameters for cuckoo filter
    int slotNo = 4;
    size_t keySize = keys.size();
    size_t bucketNo = size_t(keySize/(loadFactor*slotNo))+1;
    fingerprint = (storage*1024.0f*8.0f-(cuckooBlackSize)*39/0.9-FLOW_EST_SIZE*17/0.9)*loadFactor/(countKey) -3;
    //fingerprintOld = fingerprint;
    long maxNumKicks = 1000;
    cout<<"* Fingerprint length: "<<fingerprint<<endl;

    //init cuckoo filter
    cuckooFilter.ClearTable();
    cuckooFilter.cuckooFilterInit(bucketNo,fingerprint,slotNo,maxNumKicks);

    // add key to cuckoo filter
    addCuckooFilter(keys, keyActions);

    // --------------------------------------
    // Add aggregate keys to cuckooTable
    size_t aggregateKeySize = aggregateKeys.size();
    bucketNo = long(aggregateKeySize/(loadFactor*slotNo));
    cuckooAggrKeyTable.CuckooTableInit(bucketNo,fingerprint,slotNo,maxNumKicks);

    for(size_t i = 0; i < aggregateKeySize; i++)
    {
        size_t found = aggregateKeys[i].find('/');
        string str = aggregateKeys[i].substr(0,found);
        string prefixstr = aggregateKeys[i].substr(found+1,
                                                   aggregateKeys[i].size()-found);

        cuckooAggrKeyTable.AddKey(str,str2num(prefixstr));
    }

    vector<string>().swap(keys);
    vector<int>().swap(keyActions);
    vector<string>().swap(blackKeys);
    vector<int>().swap(blackkeyPrefixes);
    vector<string>().swap(aggregateKeys);

    // --------------------------------------
    // Add blackkey to cuckooTable
    /*cout<<"* Add blackkey to cuckooTable!"<<endl;
    size_t blackKeySize = blackKeys.size();
    size_t bucketSize = int((blackKeySize)/(loadFactor*slotNo))+1;

    int MaxKickoutNum = 1000;
    cuckooBlackKeyTable.ClearTable();
    cuckooBlackKeyTable.CuckooTableInit(bucketSize,fingerprint,slotNo,
                                        MaxKickoutNum);
    for(size_t i = 0; i < blackKeySize; i++)
    {
        cuckooBlackKeyTable.AddKeyPrefix(blackKeys[i],blackkeyPrefixes[i], 4);
    }

    blackKeys.clear();
    blackkeyPrefixes.clear();* /

    // ---------------------------------------
    for(int i = 0; i < trieNum; i++)
    {
        bTrie[i].deleteChild(bTrie[i].root);
    }
    if(!bTrie)
        delete[] bTrie;

    cout<<"* Aggregation Return! "<<endl;
    return countKey;

}*/

bool clusterAction(vector<string> &flow, vector<int> &flowaction,
                   vector<size_t> &mask, vector<size_t> &index, vector<int>& flowactionunique)
{
    vector<string> subIpBinary;
    size_t ip_num = flow.size();
    const char *charIP;
    size_t ip, subIP;
    for (size_t i = 0; i < ip_num; i++)
    {
        charIP = flow[i].c_str(); // get IP string
        ip = parseIPV4string(charIP); //convert IP to int

        int mi = 0;
        subIP = ip & mask[mi]; // mask
        subIpBinary.push_back(DecToBin(subIP)); // store /8
    }

    //--------------------------------------------------------
    // unique items in subIPv8
    cout<<"* Assign flows to subtrees according to prefix/8 and actions ... ..."
    <<endl;
    vector<string> subIPv8;
    subIPv8 = subIpBinary; // copy prefix length /8 prefix and group
    {
        vector<string>::iterator it;
        it = unique (subIPv8.begin(), subIPv8.end());
        subIPv8.resize( distance(subIPv8.begin(),it));
    }

    // ------------------------------------------------------
    // clustering according to prefix length /8
    index.clear();
    {
        size_t  pos = 0;
        vector<string>::iterator posIt = subIpBinary.begin();
        for(size_t i = 1; i < subIPv8.size(); i++)
        {
            pos = find(posIt,subIpBinary.end(),subIPv8[i])-subIpBinary.begin();
            index.push_back(pos);
            posIt = subIpBinary.begin();
            advance (posIt,pos);
        }
        index.push_back(subIpBinary.end()-subIpBinary.begin());
    }
    cout<<"* /8 index size:       "<<index.size()<<endl;

    //---------------------------------------------------------------
    // Cluster according to actions

    flowactionunique = flowaction;
    {
        vector<int>::iterator it;
        it = unique (flowactionunique.begin(), flowactionunique.end());
        flowactionunique.resize( distance(flowactionunique.begin(),it));
    }
    cout<<"* Unique flowaction size: "<<flowactionunique.size()<<endl;

    vector<size_t> index1;
    for(size_t i = 1; i < flowactionunique.size(); i++)
    {
        size_t It = find(flowaction.begin(),flowaction.end(),flowactionunique[i])
        -flowaction.begin();//find(~L(:,h1),1); // the first 0 in bucket h1
        index1.push_back(It);
    }

    // ----------------------------------------------------------------
    // merge index and index1
    vector<size_t> index2;
    merge(index.begin(), index.end(), index1.begin(), index1.end(),
          back_inserter(index2));
    {
        vector<size_t>::iterator it;
        it = unique (index2.begin(), index2.end());
        index2.resize( distance(index2.begin(),it));
    }

    index.clear();
    index = index2;
    cout<<"* Index size: "<<index.size()<<endl;

    // Release space
    vector<size_t>().swap(index2);
    vector<string>().swap(subIPv8);
    //vector<int>().swap(flowactionunique);
    vector<string>().swap(subIpBinary);
    return 1;
}

bool addCuckooFilter(vector<string> &keys, vector<int> &keyPrefixes,
                      vector<int> &keyActions, CuckooFilter& cuck, CuckooFilter& cuck0)
{
    cout<<"* Cuckoo filter adding ..."<<endl;
    // define variables
    bool flagAdd,flagAdd0;
    int countFail = 0;
    string str;
    size_t keySize = keys.size();

    for(size_t i = 0; i < keySize; i++)
    {
        str = keys[i]+'/'+num2str(keyPrefixes[i]);

        //add keys to cuckoo filter
        flagAdd = cuck.AddKey(str,(keyActions[i]));
        if(flagAdd == 0)
        {
            countFail++;
            cout << "* Line 818 Cuckoo filter flag_add fail"<<endl;
            cout<<"* Order: "<<i<<"  ";
            break;
        }

        flagAdd0 =cuck0.AddKey(str,(keyActions[i]));
        if(flagAdd0 == 0)
        {
            countFail++;
            cout << "* Line 827 CuckooFilter0 flag_add fail"<<endl;
            cout<<"* Order: "<<i<<"  ";
            break;
        }


    }
    cout<<"* Line 834 Count fail num: "<<countFail<<endl;
    return 1;

}

bool addCuckooFilter(vector<string> &keys, vector<int> &keyActions, CuckooFilter& cuck)
{
    // define variables
    cout<<"* Cuckoo filter adding ..."<<endl;
    bool flagAdd;
    int countFail = 0;
    string str;
    size_t keySize = keys.size();

    for(size_t i = 0; i < keySize; i++)
    {
        str = keys[i];
        //add keys to cuckoo filter
        flagAdd = cuck.AddKey(str,(keyActions[i]));
        if(flagAdd == 0)
        {
            countFail++;
            cout << "* 664 Cuckoo filter flag_add fail"<<endl;
            cout<<"* Order: "<<i<<"  ";
            break;
        }

    }
    cout<<"* Line 670 Count fail num: "<<countFail<<endl;
    return 1;

}

bool addCuckooFiltermL(vector<string> &keys, vector<int> &keyActions,char mL0[][ACTIONSIZE][20] , CuckooFilter& cuck)
{
    // define variables
    cout<<"* Cuckoo filter adding ..."<<endl;
    bool flagAdd;
    int countFail = 0;
    string str;
    size_t keySize = keys.size();

    for(size_t i = 0; i < keySize; i++)
    {
        str = keys[i];
        //add keys to cuckoo filter
        flagAdd = cuck.AddKeyCount(str,(keyActions[i]),mL0);
        if(flagAdd == 0)
        {
            countFail++;
            cout << "* 664 Cuckoo filter flag_add fail"<<endl;
            cout<<"* Order: "<<i<<"  ";
            break;
        }

    }
    cout<<"* Line 670 Count fail num: "<<countFail<<endl;
    return 1;

}

bool addCuckooFilter0(vector<string> &keys, vector<int> &keyPrefixes,
                      vector<int> &keyActions, char mL0[][ACTIONSIZE][20], CuckooFilter& cuck, CuckooFilter& cuck0)
{
    cout<<"* Cuckoo filter adding ..."<<endl;

    // define variables
    bool flagAdd,flagAdd0;
    int countFail = 0;
    string str;
    size_t keySize = keys.size();

    for(size_t i = 0; i < keySize; i++)
    {
        str = keys[i]+'/'+num2str(keyPrefixes[i]);

        //add keys to cuckoo filter
        flagAdd = cuck.AddKeyCount(str,(keyActions[i]),mL0);
        if(flagAdd == 0)
        {
            countFail++;
            cout << "* Line 818 Cuckoo filter flag_add fail"<<endl;
            cout<<"* Order: "<<i<<"  ";
            break;
        }

        flagAdd0 =cuck0.AddKey(str,(keyActions[i]));
        if(flagAdd0 == 0)
        {
            countFail++;
            cout << "* Line 827 CuckooFilter0 flag_add fail"<<endl;
            cout<<"* Order: "<<i<<"  ";
            break;
        }


    }
    cout<<"* Line 834 Count fail num: "<<countFail<<endl;
    return 1;

}

bool insertWordTrieSimple(Trie *bTrie, int trieNum, vector<size_t> &index,
                          vector<string> &flow,vector<int> &keyprefixlength,
                          vector<int> &flowaction, vector<int>& flowactionunique)
{
    cout<<"* Trie insert word ... ..."<<endl;
    vector<string> subIPvtmp;
    vector<int>     prelentmp;
    vector<int>      actiontmp;
    string prefix;
    const char *charIP;

    for(uint16_t ti = 0; ti<trieNum; ti++)
    {
        if (ti == 0)
            bTrie[ti].maction = flowactionunique[0];//(ti)/64;
        else
            bTrie[ti].maction = (flowaction[index[ti-1]]);

        int mi = 32-8;              // only insert the original flow prefix
        size_t beginIndex = 0;     // the start index for the first subtree
        if (ti > 0)
            beginIndex = index[ti-1];   // start index for other subtrees

        /* initialize values*/
        subIPvtmp.clear();
        prelentmp.clear();
        actiontmp.clear();


        /* get the values for one subtree*/
        for (int i = beginIndex; i < index[ti]; i++)
        {
            subIPvtmp.push_back(parseIPV42bin(flow[i].c_str()));
            prelentmp.push_back((keyprefixlength[i]));
            actiontmp.push_back(flowaction[i]);


        }
        for(int i = 0; i <subIPvtmp.size(); i++)
        {

            charIP = subIPvtmp[i].c_str();
            prefix = string(charIP, 0, prelentmp[i]);
            bTrie[ti].addWord(prefix,0, iskey, prelentmp[i],actiontmp[i]);
        }


    }
    subIPvtmp.clear();
    prelentmp.clear();
    actiontmp.clear();
    return 1;
}
bool aggregateTrie(Trie *bTrie, int ai, ActionOrder *actionOrder, size_t &countkey,
                   size_t &countaggregatekey, size_t &countblackkey,
                   size_t &countorikey, vector<string> &key,
                   vector<int> &keyaction,vector<string> &blackkey,
                   vector<int> &blackkeyPrefixes, vector<string> &aggregatekey,
                   int& prefixlength,
                   bool isPrint, bool isInit)
{
    //cout<<"* Aggregate keys ... ..."<<endl;
    vector<char> word;
    for(uint16_t ti = 0; ti < actionOrder[ai].aTrieOder.size(); ti++)
    {
        // compute keynum for each node
        size_t recoverCount = 0;
        bTrie[actionOrder[ai].aTrieOder[ti]].recoverTrie(bTrie[actionOrder[ai].aTrieOder[ti]].root8, recoverCount);

        bTrie[actionOrder[ai].aTrieOder[ti]].computeKeyNum(bTrie[actionOrder[ai].aTrieOder[ti]].root8,iskey);

        bTrie[actionOrder[ai].aTrieOder[ti]].computeNodeNum(bTrie[actionOrder[ai].aTrieOder[ti]].root);

        bTrie[actionOrder[ai].aTrieOder[ti]].arregatePrefix8(bTrie[actionOrder[ai].aTrieOder[ti]].root8, g_vweightThld[ai], prefixlength,  isInit);

        bTrie[actionOrder[ai].aTrieOder[ti]].nodeCount(bTrie[actionOrder[ai].aTrieOder[ti]].root,countkey,countaggregatekey,countblackkey,countorikey);

        // print key and blackkey to file
        if(isPrint)
            bTrie[actionOrder[ai].aTrieOder[ti]].printNode(bTrie[actionOrder[ai].aTrieOder[ti]].root,word,key,keyaction,blackkey,blackkeyPrefixes,aggregatekey);


    }

    word.clear();
    return 1;

}

void assignAction(vector<string> &flow,vector<size_t> &flow_cnt,
                  vector<int> &flowaction,vector<size_t> &pkt_num_action,
                  int actionSize)
{
    cout<<"* Assign actions ... ..."<<endl<<endl;
    vector<string> actionBound;
    actionBound.push_back("64.0.0.0");// = {"25.0.0.0","32.0.0.0"};
    actionBound.push_back("128.0.0.0");
    actionBound.push_back("192.0.0.0");// = {"25.0.0.0","32.0.0.0"};
    actionBound.push_back("255.255.255.255");

    pkt_num_action.assign(actionSize,0);

    for(int i = 0; i < flow.size(); i++)
    {
        for (int ai = 0; ai < actionSize; ai ++)
        {
            if(parseIPV4string(flow[i].c_str())<
               parseIPV4string(actionBound[ai].c_str()))
            {
                flowaction.push_back((ai));
                pkt_num_action[ai] += (flow_cnt[i]);
                break;
            }
        }
    }
}


/*bool lookupKey(size_t &line, int actionSize, size_t &flow_int,size_t &flow_cnt_int, vector<size_t> &mask,
               double &key_sum, double &pkt_sum, double &aggrSum, floats& keySums, floats& countIps, double &countIP,
               size_t &countNum, double &countIP0, size_t &countNum0, size_t &countblack_tmp,
               vector<string> &overbigKey, vector<size_t> &overbigKey_cnt, Trie* trie, int2s& bigNonCounts, long2s& timeCounts,
               vector<int> &UniquePrefix,vector<int> &UniqueAggPrefix)
{

    int prefix;
    int flag_look,flag_look0, flag_lookkey;
    bool flag_lookblack = 0;
    uint32_t subIP;
    string flowstr;
    string flowIpv4 = parsedec2IPV4(flow_int);

    uint32_t ip = flow_int; //convert IP to int
    string keyType_cur = "0";
    string flow_action_str;

    //pkt_sum += flow_cnt_int;
    //cout<<"1"<<endl;

    // -------------------------------
    // look up blackkey
    int iflowaction;
    //if(CUCKOO_BLACK_SIZE > 0)
    /*{
        for(int mi = 0; mi <= 32-8; mi++)
        {

            subIP = ip & mask[mi]; // mask
            flowstr = parsedec2IPV4(ip);//s.str();
            prefix = mi+8;

            flag_lookblack = cuckooBlackKeyTable.LookUpKeyAction(flowstr,prefix,iflowaction);
            if (flag_lookblack)
            {
                countblack_tmp++;
                break;
            }
        }


    }* /
    //cout<<"2"<<endl;

    flag_look = 0;
    //if (flag_lookblack == 0) // not a blackkey
    {
        // -------------------------------
        // look up key
        for(int mi = 0; mi < UniquePrefix.size(); mi++)
        {
            subIP = ip & mask[UniquePrefix[mi]-8];
            flowstr = parsedec2IPV4(subIP);
            prefix = UniquePrefix[mi];

            flag_lookkey = cuckooTableKey.LookUpKeyActionCount(flowstr,prefix,iflowaction,flow_cnt_int);
            if (flag_lookkey)
            {

                keyType_cur = "1";
                flow_action_str = num2str(iflowaction);
                #pragma omp atomic
                key_sum += (flow_cnt_int);
                for(int ai = 0; ai <actionSize; ai++)
                {
                    if(iflowaction == ai)
                    {
                        #pragma omp atomic
                        keySums[ai] += (flow_cnt_int);
                    }

                }
                break;

            }
        }
        //cout<<"3"<<endl;
        // -------------------------------
        // look up aggregate key
        bool isAggregatekey = 0;
        //if(cuckooAggrKeyTable.mm > 1)
        {
            for(int mi = 0; mi < UniqueAggPrefix.size(); mi++)
            {
                subIP = ip & mask[UniqueAggPrefix[mi]-8];
                flowstr = parsedec2IPV4(subIP);
                prefix = UniqueAggPrefix[mi];
                isAggregatekey = cuckooAggrKeyTable.LookUpKey(flowstr,prefix);
                if (isAggregatekey && !flag_lookkey)
                {
                    #pragma omp atomic
                    aggrSum += (flow_cnt_int);
                    break;

                }
            }
        }
        //cout<<"4"<<endl;
        // -------------------------------
        // lookup key from cuckoo filter
        vector<int> iactions;
        bool isAKey = 0;

        for(int mi = 0; mi < UniqueAggPrefix.size(); mi++)
        {
            bool overbigFlag = 0;
            subIP = ip & mask[UniqueAggPrefix[mi]-8];
            flowstr = parsedec2IPV4(subIP)+"/"+num2str(UniqueAggPrefix[mi]);
            iactions.clear();
            flag_look = cuckooFilter.LookUpKeyActions(flowstr,iactions);
            if (flag_look != 0)
            {
                isAKey = 1;

                #pragma omp atomic
                countNum += flag_look;
                #pragma omp atomic
                countIP += flag_look*(flow_cnt_int);

                // each action lookups

                for(int ai = 0; ai <actionSize; ai++)
                {
                    for(int li = 0; li < iactions.size(); li++)
                    {
                        if(iactions[li] == ai)
                        {
                            #pragma omp atomic
                            countIps[ai] += (flow_cnt_int);
                        }

                    }
                }
                //cout<<"41"<<endl;
                // blackkeys
                if(!flag_lookkey && !overbigFlag)
                {
                    overbigFlag = 1;
                    //#pragma omp critical
                    {
                        //trie->addWordCount(DecToBin(flow_int),32, isnonkey, 0);
                    }

                }
            }
        }
        // ---------------------------------------------------
        // add to flow est filter
        if(isAKey == 0)  // not inside cuckooFilter
        {
            // add to filter
            //if(FLOW_EST_SIZE>0)
            //addFlowEst(flowIpv4, bigNonCounts, timeCounts, overbigKey );

        }

    }

    //cout<<"5"<<endl;
    // -------------------------------
    // lookup key from cuckoo filterInit0

    vector<int> iactions;
    for(int mi = 0; mi < UniquePrefix.size(); mi++)
    {
            subIP = ip & mask[UniquePrefix[mi]-8];
            flowstr = parsedec2IPV4(subIP)+"/"+num2str(UniquePrefix[mi]);
            iactions.clear();
            flag_look0 =cuckooFilterInit0.LookUpKeyActions(flowstr,iactions);
            if (flag_look0)
            {
                #pragma omp atomic
                countNum0 += flag_look0;
                #pragma omp atomic
                countIP0 += flag_look0*(flow_cnt_int);

            }

    }

    //cout<<"6"<<endl;

    return 1;

}*/

bool assignAction(vector<string> &key,vector<int> &keyaction,int &actionSize)
{
    cout<<"* Assign action ... ..."<<endl;
    vector<string> actionBound;
    actionBound.push_back("64.0.0.0");
    actionBound.push_back("128.0.0.0");
    actionBound.push_back("192.0.0.0");
    actionBound.push_back("255.255.255.255");
    keyaction.clear();

    for(int i = 0; i < key.size(); i++)
    {
        for (int ai = 0; ai < actionSize; ai ++)
        {
            if(parseIPV4string(key[i].c_str())<
               parseIPV4string(actionBound[ai].c_str()))
            {
                keyaction.push_back((ai));
                break;
            }
        }

    }
    return 1;
}


void prefixNum(vector<int> &keyprefix, vector<int> &UniquePrefix)
{
    vector<int > keyprefixSort = keyprefix;
    sort(keyprefixSort.begin(),keyprefixSort.end());
    vector<int>::iterator it;
    it = unique (keyprefixSort.begin(), keyprefixSort.end());
    keyprefixSort.resize( distance(keyprefixSort.begin(),it) );
    UniquePrefix.clear();
    UniquePrefix = keyprefixSort;

    for(int i = 0; i < UniquePrefix.size(); i++)
        cout<<"prefix: "<<UniquePrefix[i]<<"    ";
    cout<<endl;
}

void keySort(vector<string>& keys, vector<size_t>& keyNos)
{
    vector<size_t> indexes;
    size_t keySize = keys.size();
    for(size_t i = 0; i < keys.size(); i++)
    {
        indexes.push_back(i);
    }
    sort(indexes.begin(), indexes.end(), CmpMass(keyNos));

    size_t curip_num = 0;
    vector<string> keySorts;
    vector<size_t> keyNoSorts;
    while(curip_num < keySize)
    {
        keySorts.push_back(keys[indexes[curip_num]]);
        keyNoSorts.push_back(keyNos[indexes[curip_num]]);
        curip_num ++;
    }
    keys.clear();
    keyNos.clear();
    keys = keySorts;
    keyNos = keyNoSorts;
}

void keySort3(vector<string>& keys, vector<size_t>& keyNos, vector<int>& keyActions)
{
    vector<size_t> indexes;
    size_t keySize = keys.size();
    for(size_t i = 0; i < keys.size(); i++)
    {
        indexes.push_back(i);
    }
    sort(indexes.begin(), indexes.end(), CmpMass(keyNos));

    size_t curip_num = 0;
    vector<string> keySorts;
    vector<size_t> keyNoSorts;
    vector<int> keyActSorts;

    while(curip_num < keySize)
    {
        keySorts.push_back(keys[indexes[curip_num]]);
        keyNoSorts.push_back(keyNos[indexes[curip_num]]);
        keyActSorts.push_back(keyActions[indexes[curip_num]]);

        curip_num ++;
    }
    keys.clear();
    keyNos.clear();
    keyActions.clear();

    keys = keySorts;
    keyNos = keyNoSorts;
    keyActions = keyActSorts;
}


void keySortInt(vector<string>& keys, vector<int>& keyNos)
{
    vector<size_t> indexes;
    size_t keySize = keys.size();
    for(size_t i = 0; i < keys.size(); i++)
    {
        indexes.push_back(i);
    }
    sort(indexes.begin(), indexes.end(), CmpMassInt(keyNos));

    size_t curip_num = 0;
    vector<string> keySorts;
    vector<int> keyNoSorts;
    while(curip_num < keySize)
    {
        keySorts.push_back(keys[indexes[curip_num]]);
        keyNoSorts.push_back(keyNos[indexes[curip_num]]);
        curip_num ++;
    }
    keys.clear();
    keyNos.clear();
    keys = keySorts;
    keyNos = keyNoSorts;
}


/*void addFlowEst(string& key, int2s& bigNonCounts, long2s& timeCounts, strings& overBigKeys )
{
    // generate a random number, if >= 1,  do nothing
    int v1 = rand() % 1000;         // v1 in the range 0 to 99
    if(v1 >29)
        return;

    // time++ and find maximum value
    long nRow = timeCounts.size();
    int nCol = timeCounts[0].size();
    long xMax = 0;
    int yMax = 0;
    long tMax = timeCounts[0][0];
    long cMax = bigNonCounts[0][0];
    for(int i = 0; i < nRow; i++)
    {
        for(int j = 0; j < nCol; j++)
        {
            if(timeCounts[i][j] != 0)
            {
                timeCounts[i][j] ++;
            }
            if(timeCounts[i][j]>tMax)
            {
                   xMax = i;
                   yMax = j;
                   cMax = bigNonCounts[i][j];
            }
            if(timeCounts[i][j]==tMax && bigNonCounts[i][j]>cMax)
            {
                   xMax = i;
                   yMax = j;
                   cMax = bigNonCounts[i][j];
            }
        }

    }

    // if inside cuckooFIlterFlowEst, count++, return
    long posx;
    int posy;
    if(cuckooFilterFlowEst.LookUpKeyBack(key,posx,posy))
    {
        bigNonCounts[posx][posy] ++;
        //cout<<"* FLOW_EST_COUNT: "<<FLOW_EST_COUNT<<endl;
        // if count==7, feedback and delete the item
        if(bigNonCounts[posx][posy] == 2)
        {
            // feedback to blackkey
            overBigKeys.push_back(key);
            // delete the item
            cuckooFilterFlowEst.RemoveKey(key);
            bigNonCounts[posx][posy] = 0;
            timeCounts[posx][posy] = 0;
            FLOW_EST_COUNT --;

        }
        return;
    }


    // if not overflow, insert
    if(FLOW_EST_COUNT<FLOW_EST_SIZE)
    {
        FLOW_EST_COUNT++;
        cuckooFilterFlowEst.AddKey(key,0);
        if(cuckooFilterFlowEst.LookUpKeyBack(key,posx,posy))
        {
            bigNonCounts[posx][posy] = 1;
            timeCounts[posx][posy] = 1;
            FLOW_EST_COUNT ++;
        }
        return;
    }

    // if overflow, delete one, and insert
    if(FLOW_EST_COUNT>=FLOW_EST_SIZE)
    {
        //cout<<"* FLOW_EST_COUNT>=FLOW_EST_SIZE!"<<FLOW_EST_COUNT<<endl<<endl;
        // delete one, find the one with maximum timeCounts and minimum value
        cuckooFilterFlowEst.RemovePos(xMax,yMax);
        bigNonCounts[xMax][yMax] = 0;
        timeCounts[xMax][yMax] = 0;

        // insert one
        cuckooFilterFlowEst.AddKey(key,0);
        if(cuckooFilterFlowEst.LookUpKeyBack(key,posx,posy))
        {
            bigNonCounts[posx][posy] = 1;
            timeCounts[posx][posy] = 1;
        }
        return;
    }

}*/

void compressAct(floats& haoOvers, float target, int actionSize, ints& compressActions, ints& decompressActions)
{
    // --------------------------------------------
    // determine which to remain unchanged
    for(int ai = 0; ai < actionSize; ai++)
    {
        // compute aggregate overrate
        float increament = (haoOvers[ai] - target);

            // compress, lower the threshold
            if(increament > 0.0f)
            {
                decompressActions.push_back(ai);
                cout<<"Decompress actions: "<<ai<<" ";
            }
            else
            {
                compressActions.push_back(ai);
                cout<<"Compress actions: "<<ai<<" ";
            }

    }
    cout<<endl;

}





