/**
trie_cuckoo_noise_07_21.cpp
create by: Yunhong
create time: 07/21/2015
*/


#include "Trie_cuckoo_noise_07_21.h"

ofstream outfileR[3];

// time
double runTime = 0.0;

// packet generation rate: pps
const float rateParameter0 = 297298.0;

int main(int argc, char * argv[])
{
    cout << "\n";
    cout << "PRIME_OPENMP\n";
    cout << "  C++/OpenMP version\n";

    cout << "\n";
    cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
    cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";

    CuckooFilter cuckooFilter[3];  // cuckoo filter for aggregation
    CuckooFilter cuckooFilterInit0[3]; // cuckoo filter for non aggregation
    CuckooTable cuckooBlackKeyTable[3];    // cuckoo table for black key
    CuckooTable cuckooTableKey[3]; // cuckoo table for non aggregation
    CuckooTable cuckooAggrKeyTable[3]; // cuckoo table for aggregated key
    CuckooFilter cuckooFilterFlowEst[3];

    // -----------------------------------------------------------------------80
    // Global variables
    CUCKOO_BLACK_SIZE = strtof(argv[6],NULL); // black table size
    FLOW_EST_SIZE = 0;       // flow estimator size
    float overTarget = 0.01;    // target overselection rate

    // -----------------------------------------------------------------------80
    /* initialize random seed: */
    srand (time(NULL));

    // time v.r.
    struct timeval gen_start, gen_end; /* gettimeofday stuff */
    struct timezone tzp;
    long timeInv = 0;

    // -----------------------------------------------------------------------80
    int switchNum = 3;
    const int actionSize = ACTIONSIZE;

    // -----------------------------------------------------------------------80
    /*define mask*/
    vector<size_t> mask;
    size_t maskIP;
    mask.clear();
    cout<<"* Compute prefixes main ... ..."<<endl<<endl;;
    for(int i = 8; i <= 32; i++)
    {
        maskIP = (size_t(pow(2,i))-1)<<(32-i);
        mask.push_back(maskIP);
    }

    // -----------------------------------------------------------------------80
    // load key files and add keys to cuckoo filter
    std::string inFileName = argv[1];
    VUPrefix vuniquePrefix;
    VUPrefix vuniqueAggPrefix;

    vuniquePrefix.clear();
    vuniqueAggPrefix.clear();

    long but = 180000/3;
    char mL0[but][actionSize][20];//
    bzero(&mL0,sizeof(mL0));

    int finger = 0;
    int finger0 = 0;

    intss flowactionunique;
    flowactionunique.resize(switchNum);

    ints flowActions0; // action 0 1 2
    for(int i = 0; i < actionSize+1; i++)
    {
        flowActions0.push_back(i);
    }

    loadKeys2Filter(inFileName, mask, vuniquePrefix, vuniqueAggPrefix, mL0, argv, finger, finger0, switchNum, actionSize,
    flowactionunique, cuckooTableKey,  cuckooBlackKeyTable,  cuckooAggrKeyTable,
                       cuckooFilter,  cuckooFilterInit0);

    for(int si = 0; si < switchNum ; si ++)
    {
        for(int ai = 0; ai < actionSize; ai ++)
        {
            cout<<"* A: "<<flowactionunique[si][ai]<<" ";
        }
        cout<<endl;
    }

    // ----------------------------------------------------------------------80
    vector<vector<size_t> > keyCountcs;// = vector<vector<size_t> > (but, vector<size_t>(actionSize, 0));
    vector<vector<size_t> > keyCountcs0;//= vector<vector<size_t> > (but, vector<size_t>(actionSize, 0));
    vector<vector<size_t> > keyCountDiffs;//= vector<vector<size_t> > (but, vector<size_t>(actionSize, 0));

    uint16_t blackBackSize = 0;
    float feedSumPortion[switchNum];

    // -----------------------------------------------------------------------80
    // File name for caida trace
    const char * fp[] = {"mapped-caida-1","mapped-caida-6", "mapped-caida-11",
                         "mapped-caida-16","mapped-caida-21","mapped-caida-26",
                         "mapped-caida-31","mapped-caida-36","mapped-caida-41",
                         "mapped-caida-46","mapped-caida-51","mapped-caida-56"
                        };

    // -----------------------------------------------------------------------80
    // Open out file, and write result into it
    std::ofstream outfile0[switchNum];
    for(int si = 0; si < switchNum; si++)
    {
        string outfileName = ("../test/outfile0_simple_"+ string(argv[2]) + '_' +
                              string(argv[4]))+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6])+"_s"+num2str(si)+"_c"+string(argv[7])+".csv";
        outfile0[si].open(outfileName.c_str());

        // wirte resource allocation
        //std::ofstream outfileR;
        outfileName = ("../test/resouce_assign_"+ string(argv[2]) + '_' +
                       string(argv[4]))+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6])+"_s"+num2str(si)+"_c"+string(argv[7])+".csv";
        outfileR[si].open(outfileName.c_str());
    }

    // -----------------------------------------------------------------------80
    int ***sec = new int**[switchNum];
    // create RLearn for recvs
    RLearn*** rLearn = new RLearn**[switchNum];
    for (int i = 0; i < switchNum; ++i)
    {
        rLearn[i] = new RLearn*[actionSize];
    }
    //RLearn* rLearn[switchNum][actionSize];
    for(int si = 0; si < switchNum; si++)
    {
        for(int i = 0; i < actionSize; i++)
        {
            rLearn[si][i] = new RLearn();
            initRLearn(rLearn[si][i]);

        }
    }

    // -----------------------------------------------------------------------80
    floats falsePos;
    floats falsePos0;
    floats haoFalsePos;
    floats haoFalsePos0;
    floats haoFalsePosTotal;
    floats haoFalsePos0Total;
    floats overAggr;

    falsePos.assign(switchNum,0);
    falsePos0.assign(switchNum,0);
    haoFalsePos.assign(switchNum,0);
    haoFalsePos0.assign(switchNum,0);
    haoFalsePosTotal.assign(switchNum,0);
    haoFalsePos0Total.assign(switchNum,0);
    overAggr.assign(switchNum,0);


    vector<string> overBigKeys;
    vector<size_t> overBigKeyNos;
    vector<string> blackKeys;
    vector<size_t> blackKeyNos;
    vector<int> blackActions;
    VBkInfo vBkInfo;

    floatss haoOvers;
    haoOvers.resize(switchNum);
    for(int si = 0; si < switchNum; si++)
    {
        haoOvers[si].assign(actionSize, 0);
    }

    floatss haoOvers0;
    haoOvers0.resize(switchNum);
    for(int si = 0; si < switchNum; si++)
    {
        haoOvers0[si].assign(actionSize, 0);
    }

    floatss haoOversInv;
    haoOversInv.resize(switchNum);
    for(int si = 0; si < switchNum; si++)
    {
        haoOversInv[si].assign(actionSize, 0);
    }

    floatss haoOversInv0;
    haoOversInv0.resize(switchNum);
    for(int si = 0; si < switchNum; si++)
    {
        haoOversInv0[si].assign(actionSize, 0);
    }

    floats haoOversAction;
    haoOversAction.assign(actionSize+1,0);

    floats haoOversActionTtl;
    haoOversActionTtl.assign(actionSize+1,0);

    floats haoOversAction0;
    haoOversAction0.assign(actionSize+1,0);

    floats haoOversActionTtl0;
    haoOversActionTtl0.assign(actionSize+1,0);

    size_ts countNum;
    size_ts countNum0;
    doubles countIP;
    doubles countIP0;
    size_ts countBlack;

    doubles countIPTotal;
    doubles countIP0Total;
    doubles keySum;
    doubles pktSum;
    doubles keySumTotal;
    doubles pktSumTotal;
    doubles countIPTotalOff;
    doubles countIP0TotalOff;
    doubles keySumTotalOff;
    doubles pktSumTotalOff;
    doubles aggrSum;

    countNum.assign(switchNum, 0);
    countNum0.assign(switchNum, 0);
    countIP.assign(switchNum, 0);
    countIP0.assign(switchNum, 0);
    countBlack.assign(switchNum, 0);

    countIPTotal.assign(switchNum, 0);
    countIP0Total.assign(switchNum, 0);
    keySum.assign(switchNum, 0);
    pktSum.assign(switchNum, 0);
    keySumTotal.assign(switchNum, 0);
    pktSumTotal.assign(switchNum, 0);
    countIPTotalOff.assign(switchNum, 0);
    countIP0TotalOff.assign(switchNum, 0);
    keySumTotalOff.assign(switchNum, 0);
    pktSumTotalOff.assign(switchNum, 0);
    aggrSum.assign(switchNum, 0);


    floatss keySums;
    floatss countIPs;
    keySums.resize(switchNum);
    countIPs.resize(switchNum);

    floatss countIPs0;
    countIPs0.resize(switchNum);

    for(int si = 0; si < switchNum; si++)
    {
        keySums[si].assign(actionSize,0);
        countIPs[si].assign(actionSize,0);
        countIPs0[si].assign(actionSize,0);
    }

    floatss keySumsInv;
    floatss countIPsInv;
    keySumsInv.resize(switchNum);
    countIPsInv.resize(switchNum);

    floatss countIPsInv0;
    countIPsInv0.resize(switchNum);

    for(int si = 0; si < switchNum; si++)
    {
        keySumsInv[si].assign(actionSize,0);
        countIPsInv[si].assign(actionSize,0);
        countIPsInv0[si].assign(actionSize,0);
    }

    uint64_t line = 0;
    size_tss slotNums;
    slotNums.resize(switchNum);

    for (int si = 0; si < switchNum; si++)
    {
        slotNums[si].assign(actionSize,0);  // resize each of the contained vectors
    }

    // -----------------------------------------------------------------------80
    // Define a trie
    Trie *trie[3];            // define tree
    for(int si = 0; si < switchNum; si++)
        trie[si] = new Trie();            // define tree

    // -----------------------------------------------------------------------80
    std::ifstream infile;

    strings flows;
    size_ts flowNos;

    flows.resize(20000);
    flowNos.resize(20000);

    // -----------------------------------------------------------------------80
    // load packets from file
    try{
    for (int fi = 0; fi < 12; fi++)
    {
        // --------------------------
        // Open trace file
        string pathname = fp[fi];
        std::string inFileName = argv[3];
        inFileName += pathname;
        infile.open(inFileName.c_str());
        cout<<inFileName.c_str()<<endl;
        if(!infile)
            std::cout << "* TestIP File Error " << std::endl;

        // -------------------------------------------------------------------80
        for(int si = 0; si < switchNum; si++)
        {
            countIPTotalOff[si] += countIP[si];
            countIP0TotalOff[si] += countIP0[si];
            keySumTotalOff[si] += keySum[si];
            pktSumTotalOff[si] += pktSum[si];

            //line = 0;
            countBlack[si] = 0;
            countIP[si] = 0.0f;
            countIP0[si] = 0.0f;
            keySum[si] = 0.0f;
            pktSum[si] = 0.0f;

            for(int ai = 0; ai < actionSize; ai++)
            {
                keySumsInv[si][ai] = 0;
                countIPsInv[si][ai] = 0;
                countIPsInv0[si][ai] = 0;
            }

        }

        size_t  ei;
        bool isEndFlag = 0;
        size_t updateInvDis = 10000; // interval for display
        size_t readNum = strtof(argv[5],NULL);    // interval for feedback

        while(!isEndFlag )
        {
            // init overselection rate each cycle
            for (int si = 0; si < switchNum; si++)
            {
                for(int ai = 0; ai < actionSize; ai++)
                {
                    keySumsInv[si][ai] = 0;
                    countIPsInv[si][ai] = 0;
                    countIPsInv0[si][ai] = 0;
                }
            }

            if(runTime < 10)
            {
                for(int si = 0; si < switchNum; si++)
                {
                    countNum[si] = 0;
                    countNum0[si] = 0;
                    countIP[si] = 0;
                    countIP0[si] = 0;
                    countBlack[si] = 0;

                    countIPTotal[si] = 0;
                    countIP0Total[si] = 0;
                    keySum[si] = 0;
                    pktSum[si] = 0;
                    keySumTotal[si] = 0;
                    pktSumTotal[si] = 0;
                    countIPTotalOff[si] = 0;
                    countIP0TotalOff[si] = 0;
                    keySumTotalOff[si] = 0;
                    pktSumTotalOff[si] = 0;
                    aggrSum[si] = 0;

                    for(int ai = 0; ai < actionSize; ai++)
                    {
                        keySums[si][ai] = 0;
                        countIPs[si][ai] = 0;
                        countIPs0[si][ai] = 0;
                    }

                }

            }

            if(runTime < 10)
            {
                for (int si = 0; si < switchNum; si++)
                {
                    for(int ai = 0; ai < actionSize; ai++)
                        rLearn[si][ai]->clearList();
                }
            }

            // ---------------------------------------------------------------80
            // read file
            size_t flowNum = 0;
            readFile0(infile, flows, flowNos, readNum, isEndFlag, flowNum);
            //sleep(1);
            //flowNum += 1;
            //size_t flowNum = flows.size();

            # pragma omp parallel for \
            shared ( infile,isEndFlag, updateInvDis, line, mask, keySum, pktSum, aggrSum, keySums, countIPs, countIP, countIPs0, countNum, countIP0, countNum0, \
                     countBlack,  trie, vuniquePrefix,vuniqueAggPrefix, cuckooFilter, \
                     cuckooFilterInit0,cuckooBlackKeyTable,cuckooTableKey,cuckooAggrKeyTable,cuckooFilterFlowEst) \
            private ( ei )

            for(ei = 0; ei < flowNum; ei++)
            {
                // read 1 flow
                bool readFlag1 = 0;
                size_t flowNoInt = 1;
                uint32_t flowInt;

                #pragma omp critical
                {
                    flowInt = parseIPV4string(flows[ei].c_str());
                    flowNoInt = flowNos[ei];
                    readFlag1 = 1;

                }

                if(readFlag1)
                {
                    // return position
                    long hBuckhit = -1;
                    int slot_ihit = -1;
                    long hBuck = -1;
                    int slot_i = -1;

                    // -------------------------------------------------------80
                    // lookup  key
                    {
                        int prefix;
                        int flag_look,flag_look0, flag_lookkey[switchNum];
                        bool flag_lookblack[switchNum];

                        bzero(&flag_lookblack,sizeof(flag_lookblack));

                        uint32_t subIP;
                        string flowstr;
                        string flowIpv4 = parsedec2IPV4(flowInt);

                        uint32_t ip = flowInt; //convert IP to int
                        string keyType_cur = "0";
                        string flow_action_str;

                        // ---------------------------------------------------80
                        // look up blackkey
                        int iflowaction;
                        for(int si = 0; si < switchNum; si++)
                        {
                            size_t keySumLocal = 0;
                            size_t aggrSumLocal = 0;
                            size_t countIPLocal = 0;
                            size_t countNumLocal = 0;
                            size_t countIP0Local = 0;
                            size_t countNum0Local = 0;

                            size_ts keySumsLocal;
                            size_ts countIPsLocal;
                            size_ts countIPsLocal0;

                            keySumsLocal.assign(actionSize,0);
                            countIPsLocal.assign(actionSize,0);
                            countIPsLocal0.assign(actionSize,0);

                            if(CUCKOO_BLACK_SIZE > 0)
                            {
                                int mi = 24;
                                //for(int mi = 0; mi <= 32-8; mi++)
                                {

                                    subIP = ip & mask[mi]; // mask
                                    flowstr = parsedec2IPV4(ip);//s.str();
                                    prefix = mi+8;

                                    flag_lookblack[si] = cuckooBlackKeyTable[si].LookUpKeyAction(flowstr,prefix,iflowaction);
                                    if (flag_lookblack[si])
                                    {
                                        countBlack[si] += flowNoInt;
                                        //break;
                                    }
                                }


                            }


                        // ---------------------------------------------------80
                        //for(int si = 0; si < switchNum; si++)
                        //{
                            if (flag_lookblack[si] == 0) // not a blackkey
                            {
                            // look up key
                                flag_look = 0;
                                for(int mi = vuniquePrefix[si].uPres.size()-1; mi >=0; mi--)
                                {

                                    subIP = ip & mask[vuniquePrefix[si].uPres[mi]-8];
                                    flowstr = parsedec2IPV4(subIP);
                                    prefix = vuniquePrefix[si].uPres[mi];


                                    flag_lookkey[si] = cuckooTableKey[si].LookUpKeyActionCount(flowstr,
                                                   prefix,iflowaction,flowNoInt);
                                    if (flag_lookkey[si])
                                    {
                                        keyType_cur = "1";
                                        flow_action_str = num2str(iflowaction);
                                        //#pragma omp atomic
                                        keySumLocal+= (flowNoInt);
                                        for(int ai = 0; ai <actionSize; ai++)
                                        {
                                            if(iflowaction == flowactionunique[si][ai])  // ?
                                            {
                                                //#pragma omp atomic
                                                keySumsLocal[ai] += (flowNoInt);
                                            }

                                        }
                                        break;
                                    }
                                }


                                // -----------------------------------------------80
                                // lookup key from cuckoo filter
                                vector<int> iactions;
                                bool overbigFlag[switchNum];

                                //for(int si = 0; si < switchNum; si++)
                                {
                                    bool isAKey = 0;
                                    for(int mi = 0; mi < vuniqueAggPrefix[si].uPres.size(); mi++)
                                    {
                                        overbigFlag[si] = 0;
                                        subIP = ip & mask[vuniqueAggPrefix[si].uPres[mi]-8];
                                        flowstr = parsedec2IPV4(subIP)+"/"+num2str(vuniqueAggPrefix[si].uPres[mi]);
                                        iactions.clear();
                                        hBuckhit = -1;
                                        slot_ihit = -1;
                                        flag_look = cuckooFilter[si].LookUpKeyActionsCount(flowstr,iactions,flowNoInt, mL0,
                                                    keyCountcs,  keyCountcs0, keyCountDiffs,hBuckhit,slot_ihit);
                                        if(hBuckhit!= -1 && slot_ihit!= -1)
                                        {
                                            hBuck = hBuckhit;
                                            slot_i = slot_ihit;
                                        }
                                        if (flag_look != 0)
                                        {
                                            isAKey = 1;
                                            countNumLocal += flag_look;
                                            countIPLocal+= flag_look*(flowNoInt);

                                            // each action lookups
                                            for(int ai = 0; ai <actionSize; ai++)
                                            {
                                                for(int li = 0; li < iactions.size(); li++)
                                                {
                                                    if(iactions[li] == flowactionunique[si][ai])
                                                    {
                                                        countIPsLocal[ai] += (flowNoInt);
                                                    }

                                                }
                                            }

                                            // overselection count
                                            if(!flag_lookkey[si] && !overbigFlag[si])
                                            {
                                                overbigFlag[si] = 1;
                                                #pragma omp critical
                                                {
                                                    trie[si]->addWordCountNum(DecToBin(flowInt),32, isnonkey, iactions[0], countIPLocal); // ?
                                                }

                                            }
                                        }
                                    }
                                }

                                // ---------------------------------------------------80
                                // look up aggregate key

                                //if(cuckooAggrKeyTable.mm > 1)
                                //for(int si = 0; si < switchNum; si++)
                                {
                                    bool isAggregatekey = 0;
                                    for(int mi = vuniqueAggPrefix[si].uPres.size()-1; mi >=0; mi--)
                                    {
                                        subIP = ip & mask[vuniqueAggPrefix[si].uPres[mi]-8];
                                        flowstr = parsedec2IPV4(subIP);
                                        prefix = vuniqueAggPrefix[si].uPres[mi];
                                        isAggregatekey = cuckooAggrKeyTable[si].LookUpKey(flowstr,prefix);
                                        if (isAggregatekey && !flag_lookkey[si])
                                        {
                                            //#pragma omp atomic
                                            aggrSumLocal += (flowNoInt);
                                            break;

                                        }
                                    }
                                }

                            }
                        //}

                        // -------------------------------------------------------80
                        // lookup key from cuckoo filterInit0

                        //for(int si = 0; si < switchNum; si++)
                        //{
                            vector<int> iactions;
                            for(int mi = 0; mi < vuniquePrefix[si].uPres.size(); mi++)
                            {
                                subIP = ip & mask[vuniquePrefix[si].uPres[mi]-8];
                                flowstr = parsedec2IPV4(subIP)+"/"+num2str(vuniquePrefix[si].uPres[mi]);
                                vector<int>().swap(iactions);

                                flag_look0 =cuckooFilterInit0[si].LookUpKeyActions(flowstr,iactions);
                                if (flag_look0)
                                {
                                    //#pragma omp atomic
                                    countNum0Local += flag_look0;
                                    //#pragma omp atomic
                                    countIP0Local += flag_look0*(flowNoInt);

                                    // each action lookups
                                    for(int ai = 0; ai <actionSize; ai++)
                                    {
                                        for(int li = 0; li < iactions.size(); li++)
                                        {
                                            if(iactions[li] == flowactionunique[si][ai])
                                            {
                                                countIPsLocal0[ai] += (flowNoInt);
                                            }

                                        }
                                    }

                                }

                            }
                       // }
                        // lookup key end

                        // -------------------------------------------------------80


                        #pragma omp critical
                        {
                            //for(int si = 0; si < switchNum; si++)
                           // {
                                //#pragma omp atomic
                                pktSum[si] += flowNoInt;
                                //#pragma omp atomic
                                keySum[si] += keySumLocal;
                                //#pragma omp atomic
                                aggrSum[si] += aggrSumLocal;
                                //#pragma omp atomic
                                countIP[si] += countIPLocal;
                                //#pragma omp atomic
                                countNum[si] += countNumLocal;
                                //#pragma omp atomic
                                countIP0[si] += countIP0Local;
                                //#pragma omp atomic
                                countNum0[si] += countNum0Local;
                                for(int ai = 0; ai <actionSize; ai++)
                                {
                                    //#pragma omp atomic
                                    keySums[si][ai] += keySumsLocal[ai];
                                    //#pragma omp atomic
                                    countIPs[si][ai] += countIPsLocal[ai];

                                    countIPs0[si][ai] += countIPsLocal0[ai];

                                    //#pragma omp atomic
                                    keySumsInv[si][ai] += keySumsLocal[ai];
                                    //#pragma omp atomic
                                    countIPsInv[si][ai] += countIPsLocal[ai];

                                    countIPsInv0[si][ai] += countIPsLocal0[ai];
                                }
                            }
                            size_ts().swap(keySumsLocal);
                            size_ts().swap(countIPsLocal);

                        }

                    }
                }
                else
                {
                    isEndFlag = 1; // end read flag
                }

                // update global variable
                line ++;

                // -----------------------------------------------------------80
                // update rates and write to file
                #pragma omp critical
                {
                    for(int si = 0; si < switchNum; si ++)
                    {
                        //if(size_t(line)%200 == 0)
                        {
                            // update total overselects value
                            countIPTotal[si] = countIPTotalOff[si] + countIP[si];
                            countIP0Total[si] = countIP0TotalOff[si] + countIP0[si];
                            keySumTotal[si] = keySumTotalOff[si] + keySum[si];
                            pktSumTotal[si] = pktSumTotalOff[si] + pktSum[si];

                            // -----------------------------------------------80
                            // overselection rate
                            // ----------------------------------------------
                            if ((pktSum[si] - keySum[si]) != 0)
                            {
                                falsePos[si] = float(countIP[si]-keySum[si])/float(pktSum[si]-keySum[si]);
                                falsePos0[si] = float(countIP0[si]-keySum[si])/float(pktSum[si]-keySum[si]);
                            }


                            if(keySum[si] != 0)
                            {
                                haoFalsePos[si] = float(countIP[si]-keySum[si])/float(keySum[si]);
                                haoFalsePos0[si] = float(countIP0[si]-keySum[si])/float(keySum[si]);
                            }

                            for(int ai = 0; ai < actionSize; ai++)
                            {
                                if(keySums[si][ai] != 0)
                                {
                                    haoOvers[si][ai] = float(countIPs[si][ai]-keySums[si][ai])/float(keySums[si][ai]);
                                }
                                if(keySumsInv[si][ai] != 0)
                                {
                                    haoOversInv[si][ai] = float(countIPsInv[si][ai]-keySumsInv[si][ai])/float(keySumsInv[si][ai]);
                                }

                                if(keySums[si][ai] != 0)
                                {
                                    haoOvers0[si][ai] = float(countIPs0[si][ai]-keySums[si][ai])/float(keySums[si][ai]);
                                }

                                if(keySumsInv[si][ai] != 0)
                                {
                                    haoOversInv0[si][ai] = float(countIPsInv0[si][ai]-keySumsInv[si][ai])/float(keySumsInv[si][ai]);
                                }
                            }

                            if(keySumTotal[si] != 0)
                            {
                                haoFalsePosTotal[si] = (countIPTotal[si]-keySumTotal[si])/(keySumTotal[si]);
                                haoFalsePos0Total[si] = (countIP0Total[si]-keySumTotal[si])/(keySumTotal[si]);
                            }

                            if(keySumTotal[si] != 0)
                                overAggr[si] = aggrSum[si]/keySumTotal[si];

                        }

                    }

                    // -----------------------------------------------------------80
                    // Get the ovs for a receiver
                    for (int ai = 0; ai < flowActions0.size(); ai++)
                    {
                        float countIPsSum = 0;
                        float keySumsSum = 0;

                        float countIPsSum0 = 0;

                        float countIPsSumTtl = 0;
                        float keySumsSumTtl = 0;

                        float countIPsSumTtl0 = 0;


                        for(int si = 0; si < switchNum; si++)
                        {
                            for(int ri = 0; ri < actionSize; ri++)
                            {

                                if (flowactionunique[si][ri] == flowActions0[ai])
                                {
                                    countIPsSum += countIPsInv[si][ri];
                                    keySumsSum += keySumsInv[si][ri];
                                    countIPsSum0 += countIPsInv0[si][ri];

                                    countIPsSumTtl += countIPs[si][ri];
                                    keySumsSumTtl += keySums[si][ri];
                                    countIPsSumTtl0 += countIPs0[si][ri];

                                }
                            }

                        }

                        if(keySumsSum != 0)
                            haoOversAction[ai] = (countIPsSum - keySumsSum)/keySumsSum;

                        if(keySumsSumTtl != 0)
                            haoOversActionTtl[ai] = (countIPsSumTtl - keySumsSumTtl)/keySumsSumTtl;

                        if(keySumsSum != 0)
                            haoOversAction0[ai] = (countIPsSum0 - keySumsSum)/keySumsSum;

                        if(keySumsSumTtl != 0)
                            haoOversActionTtl0[ai] = (countIPsSumTtl0 - keySumsSumTtl)/keySumsSumTtl;
                    }
                }

            }

            //vector<string>().swap(flows);
            //vector<size_t>().swap(flowNos);

            // -----------------------------------------------80
            // Display and write to file
            for(int si = 0; si < switchNum; si++)
            {
                outfile0[si]<<"line,"<<line<<",total,"<<haoFalsePosTotal[si]<<",total0,"<<haoFalsePos0Total[si]<<",false,"<<
                falsePos[si]<<",over,"<<haoFalsePos[si]<<",false0,"<<falsePos0[si]<<",over0,"<<
                haoFalsePos0[si]<<",overaggr,"<<overAggr[si]<<",overcuckoo,"<<(haoFalsePosTotal[si]-overAggr[si])<<",";
                for(int ai = 0; ai < actionSize; ai++)
                    outfile0[si]<<"over_uai,"<<haoOversAction[flowactionunique[si][ai]]<<",";
                for(int ai = 0; ai < actionSize; ai++)
                    outfile0[si]<<"over_utai,"<<haoOversActionTtl[flowactionunique[si][ai]]<<",";
                for(int ai = 0; ai < actionSize; ai++)
                    outfile0[si]<<"over_ai,"<<haoOvers[si][ai]<<",";

                outfile0[si]<<"countIP,"<<countIPTotal[si]<<",countIP0,"<<countIP0Total[si]<<",keysum,"<<keySumTotal[si]<<
                ",pktSum,"<<pktSumTotal[si]<<",aggrSum,"<<aggrSum[si]<<",blackkey_num,"<<countBlack[si]<<",feedback,"<<blackBackSize<<",feedsumportion,"
                <<feedSumPortion[si]<<",finger0,"<<finger0<<",finger,"<<finger<<",time,"<<timeInv<<",runTime,"<<runTime<<",";

                for(int ai = 0; ai < actionSize+1; ai++)
                    outfile0[si]<<"over_uai0,"<<haoOversAction0[ai]<<",";
                for(int ai = 0; ai < actionSize+1; ai++)
                    outfile0[si]<<"over_utai,"<<haoOversActionTtl[ai]<<",";
                for(int ai = 0; ai < actionSize+1; ai++)
                    outfile0[si]<<"over_utai0,"<<haoOversActionTtl0[ai]<<",";
                for(int ai = 0; ai < actionSize; ai++)
                    outfile0[si]<<"over_ai0,"<<haoOvers0[si][ai]<<",";
                for(int ai = 0; ai < actionSize; ai++)
                    outfile0[si]<<"over_inv_ai0,"<<haoOversInv0[si][ai]<<",";
                outfile0[si]<<endl;


                cout<<endl<<"line "<<line<<" total "<<haoFalsePosTotal[si]<<" total0 "<<haoFalsePos0Total[si]<<" false "<<
                falsePos[si]<<" over "<<haoFalsePos[si]<<" false0 "<<falsePos0[si]<<" over0 "<<
                haoFalsePos0[si]<<" overaggr "<<overAggr[si]<<" overcuckoo "<<(haoFalsePosTotal[si]-overAggr[si])<<endl;

                for(int ai = 0; ai < actionSize; ai++)
                    cout<<"over_uai "<<haoOversAction[flowactionunique[si][ai]]<<" ";
                for(int ai = 0; ai < actionSize; ai++)
                    cout<<"over_utai,"<<haoOversActionTtl[flowactionunique[si][ai]]<<",";
                for(int ai = 0; ai < actionSize; ai++)
                    cout<<"over_ai,"<<haoOvers[si][ai]<<",";

                cout<<"countIP "<<countIPTotal[si]<<" countIP0 "<<countIP0Total[si]<<" keysum "<<keySumTotal[si]<<
                " pktSum "<<pktSumTotal[si]<<" aggrSum "<<aggrSum[si]<<" blackkey_num "<<countBlack[si]<<" feedback "<<blackBackSize<<" feedsumportion "
                <<feedSumPortion[si]<<" finger0 "<<finger0<<" finger "<<finger<<" time "<<timeInv<<",runTime,"<<runTime<<endl<<endl;
            }

            // ---------------------------------------------------------------80
            // feed back portionFeedBack% overselections and overselection rates for actions
            //if(size_t(line)%2000000 == 0 )
            if(size_t(runTime)%5 == 0)
            {
                VBkInfo().swap(vBkInfo);
                // -----------------------------------------------------------80
                // the overselects for feedback
                for(int si = 0; si < switchNum; si++)
                {
                    blackKeys.resize(3000);
                    blackKeyNos.resize(3000);
                    blackActions.resize(3000);

                    vector<char> word;
                    size_t i_seq = 0;
                    trie[si]->getLeaf(trie[si]->root,word,blackKeys,blackKeyNos, blackActions, i_seq);
                    trie[si]->deleteChild(trie[si]->root);
                    delete trie[si];
                    trie[si] = new Trie();            // define tree

                    blackKeys.resize(i_seq);
                    blackKeyNos.resize(i_seq);
                    blackActions.resize(i_seq);

                    cout<<"* i_seq: "<<i_seq<<endl;

                    int init = 0;
                    float overTotalSum = accumulate(blackKeyNos.begin(), blackKeyNos.end(), init);

                    // sort balckkeys
                    keySort3(blackKeys,blackKeyNos, blackActions);

                    float portionFeedBack = strtof(argv[4],NULL);
                    blackBackSize = portionFeedBack*blackKeys.size();

                    blackKeys.erase(blackKeys.begin()+blackBackSize, blackKeys.end());
                    blackKeyNos.erase(blackKeyNos.begin()+blackBackSize, blackKeyNos.end());
                    blackActions.erase(blackActions.begin()+blackBackSize, blackActions.end());

                    float feedSum = accumulate(blackKeyNos.begin(), blackKeyNos.end(), init);
                    feedSumPortion[si] = feedSum/overTotalSum;

                    cout<<"* feedback sum: "<<feedSum<<" num: "<<blackBackSize<<endl;

                    // output
                    BkInfo bkInfo;
                    bkInfo.bkActions = blackActions;
                    bkInfo.bkPres.assign(blackKeys.size(),32);
                    bkInfo.bks = blackKeys;
                    bkInfo.bkSizes = blackKeyNos;

                    vBkInfo.push_back(bkInfo);
                }


                // -----------------------------------------------------------80
                // Get instant reward and Q learn
                //cout<<"* Q learn!"<<endl;
                for(int si = 0; si < switchNum; si++)
                {
                    for(int ri = 0; ri < actionSize; ri++)
                    {
                        //cout<<rLearn[ri]<<endl;
                        rLearn[si][ri]->update(slotNums[si][ri],haoOversActionTtl[flowactionunique[si][ri]]);
                        //cout<<"* qleran!"<<endl;

                        // compute ebuse0
                        if(runTime<20)
                        {
                            EPSILON0 = 1.0;
                        }
                        else if(runTime<50)
                        {
                            EPSILON0 = 0.6;
                        }
                        else
                        {
                            EPSILON0 = 0.20;
                        }
                        rLearn[si][ri]->qLearn();
                        //printQList(rLearn[si][ri]);
                    }
                }


                // -----------------------------------------------------------80
                // feedback
                // ------------------------------------

                //time before calling function
                gettimeofday(&gen_start, &tzp);

                // -----------------------------------------------------------80
                // call function
                cout<<"* feedback bks!"<<endl;
                feedbackBlackkeyRL(vBkInfo, rLearn, actionSize,switchNum,slotNums, flowactionunique, line, cuckooBlackKeyTable);
                cout<<"* feedback bks end!"<<endl;

                // ------------------------------
                // time after calling function
                gettimeofday(&gen_end, &tzp);

                // time interval
                timeInv = print_elapsed("Aggr: ", &gen_start,&gen_end, 1);

                strings().swap(blackKeys);
                vector<size_t>().swap(blackKeyNos);
                vector<int>().swap(blackActions);
            }

            //strings().swap(flows);
            //vector<size_t>().swap(flowNos);

        }

        //ifstream().swap(infile);
        infile.clear();
        infile.close();
    } //for fi


    }catch (std::bad_alloc& ba){
        cerr << "bad_alloc caught: " << ba.what() << endl;
    }

    // -----------------------------------------------------------------------80
    /*clear the data structure*/
    mask.clear();

    //cuckooFilter.ClearTable();
    //cuckooBlackKeyTable.ClearTable();
    //cuckooTableKey.ClearTable();
    //cuckooAggrKeyTable.ClearTable();

    // (outfile0.close());
}

bool readFile0(ifstream& infile, vector<string> &flow, vector<size_t> &flow_cnt, size_t readNum, bool& isEndFlag, size_t& i_seq)
{
    //cout<<"* Read File!"<<endl;
    Trie *bTrie = new Trie();
    uint32_t flowInt;
    size_t flowNo;
    size_t in_num = 0;

    double timeInv = 0.0f;
    double nTime = 0.0f;

    while((infile >> flowInt >>flowNo) && timeInv < 1.0f/*&& (in_num < readNum)*/)
    {

        bTrie->addWordCount(DecToBin(flowInt),32, isnonkey, flowNo);
        in_num ++;

        // generate next time
        nTime = nextTime(rateParameter0);
        timeInv += nTime;

        if(nTime< 0)
        {
            cout<<"* readFile0: timeInv < 0. wrong!"<<endl;
            exit(0);
        }

    }

    //cout<<"* Read File end!"<<" timeInv: "<<timeInv<<endl;
    // run time
    runTime += timeInv;
    if(runTime >= INFINITY)
    {
        cout<<"* readFile0:runTime: "<<runTime<<" timeInv: "<<timeInv<<" nTime: "<<nTime<<endl;
        exit(0);
    }
    //cout<<"* timeInv: "<<endl;

    vector<char> word;
    vector<int> flowActions;

    //vector<string>().swap(flow);
    //vector<size_t>().swap(flow_cnt);
    vector<int>().swap(flowActions);
    flowActions.resize(flow.size());

    bTrie->getLeaf(bTrie->root,word,flow,flow_cnt,flowActions, i_seq);

    if(timeInv < 1.0f/*in_num < readNum*/)
    {
        isEndFlag = 1;
        cout<<"* readFile0:timeInv: "<<timeInv<<endl;
    }

    //cout<<"* Flow size: "<<flow.size()<<endl;

    //free(bTrie);
    bTrie->deleteChild(bTrie->root);
    delete bTrie;
    return true;
}


void initRLearn(RLearn* rLearn)
{
    // build a RL for rev1

    // create a qlist and initialize it
    float minState = 20;
    float maxState = 820;   // ovs
    float minAction = 20;
    float maxAction = 820; // #slots

    float state = minState;
    float action = minAction;
    float ovs = 0;
    float ovsTarget = 0.01;

    size_t numS = 8;
    size_t numA = 8;

    //float initState = 500;
    //float initAction = 500;


    //rLearn->rLearnInit();

    // initialize the table
    cout<<"* init qlist..."<<endl;
    rLearn->initQtable(minState, maxState, minAction, maxAction, numS,  numA,  ovsTarget);


    // initialized state
    cout<<"* init rlearn state... "<<endl;
    rLearn->update(state, ovs);

    printQList(rLearn);

    // get suggected action

    // get immedite reward


    // Q leran
    //rLearn.qLearn();
}

void updateBlacklist(vector<string>& overBigKeys, vector<int>& overActions, RLearn* rLearn, int actionSeq
                     , vector<string>& blackkeyPres, vector<int>& blackActionPres, ofstream& blackKeyFileOut,
                     size_t& slotNum, int si, CuckooTable& cuckooBlackKeyTable)
{
    // for each recv
    // get the suggested action
    //slotNum = rLearn->selectActionSuggest();
    cout<<"* suggested action++++++++++++++++++++++++++++++++++: "<<slotNum<<endl;

    // write to file
    outfileR[si]<<slotNum<<",";

    // -----------------------------------------------
    // identify the blackkeys for each action
    vector<string> blackKeysRecv;
    vector<int> blackActRecv;
    blackKeysRecv.clear();
    blackActRecv.clear();

    size_t blackSize = overBigKeys.size();
    for(size_t i = 0; i < blackSize; i++)
    {
        if(i < slotNum)
        {
            if(overActions[i] == actionSeq)
            {
                blackKeysRecv.push_back(overBigKeys[i]);
                blackActRecv.push_back(overActions[i]);
                //overBigKeys.erase();
            }
        }
        else // blacklist is full
        {
            break;
        }

    }

    size_t blackKeysRecvSize = blackKeysRecv.size();
    cout<<"* action: "<<actionSeq<<" #cur blackkeys: "<<blackKeysRecvSize<<endl;

    // --------------------------------------------
    // put current blackkeys to cuckoo table
    for(size_t i = 0; i < blackKeysRecvSize; i++)
    {
        if(i<slotNum)
        {
            bool addFlag = cuckooBlackKeyTable.AddKeyPrefix(blackKeysRecv[i],32, 4);
            if(!addFlag)
            {
                cout<<"1 add fail..."<<endl;
            }
        }
    }

    // -------------------------------------------
    // load the previous blackkeys for each action
    size_t keysNum = blackKeysRecvSize;
    size_t keyPresNum = blackkeyPres.size();


    // record old keys
    vector<string> blackkeyPresOld;
    vector<int> actionPresOld;
    size_t keysPreSeq = 0;

    for(size_t i = 0; i < keyPresNum; i++)
    {
        if(keysNum < slotNum && blackActionPres[i] == actionSeq) // IF with the same action
        {
            bool blackFlag = 0;
            int iflowaction;
            int prefix = 32;
            blackFlag = cuckooBlackKeyTable.LookUpKeyAction(blackkeyPres[i],prefix,iflowaction);

            if(!blackFlag)
            {
                blackKeysRecv.push_back(blackkeyPres[i]);
                blackActRecv.push_back(blackActionPres[i]);
                keysNum ++;
                keysPreSeq ++;

                // put it into the blacklist
                bool addFlag = cuckooBlackKeyTable.AddKeyPrefix(blackkeyPres[i],32, 4);
                if(!addFlag)
                {
                    cout<<"2 add fail..."<<endl;
                }
            }

        }

    }

    size_t pre_action_num = keysPreSeq;

    for(size_t i = keysPreSeq-1; i < keyPresNum; i++)
    {
        if(blackActionPres[i] == actionSeq)
        {
            blackkeyPresOld.push_back(blackkeyPres[i]);
            actionPresOld.push_back(blackActionPres[i]);
            pre_action_num ++;
        }

        if(pre_action_num > 1000)
        break;
    }


    // -----------------------------------------------
    // wirte blackkeys for each action to blacklist
    size_t bkSize = blackKeysRecv.size();

    for(int i = 0; i < bkSize; i++)
        blackKeyFileOut<<blackKeysRecv[i]<<" "<<32<<" "<<blackActRecv[i]<<endl;

    // write other keys
    size_t bkOldSize = blackkeyPresOld.size();

    for(int i = 0; i < bkOldSize; i++)
        blackKeyFileOut<<blackkeyPresOld[i]<<" "<<32<<" "<<actionPresOld[i]<<endl;

    vector<string>().swap(blackKeysRecv);
    vector<int>().swap(blackActRecv);

}
void feedbackBlackkeyRL(VBkInfo& vBkInfo, RLearn** rLearn[],
                        int actionSize, int switchNum, size_tss& slotNums, intss& flowActionUnique,
                         size_t line, CuckooTable* cuckooBlackKeyTable)
{

    selectAction(rLearn, actionSize, switchNum,slotNums);

    // -------------------------------------------
    // load the previous blackkeys for each action
    for(int si = 0; si < switchNum; si++)
    {
        cout<<"* feedbackBlackkeyRL:si: "<<si<<endl;
        string blackFileName = BLACKFILENAME + num2str(si);
        ifstream blackKeyFileIn;
        blackKeyFileIn.open(blackFileName.c_str());

        vector<string> blackkeyPres;
        vector<int> actionPres;
        string blackKeyStr;
        int prefix = 32;
        int actionInt = 0;

        while(blackKeyFileIn>>blackKeyStr>>prefix>>actionInt)
        {
            blackkeyPres.push_back(blackKeyStr);
            actionPres.push_back(actionInt);

        }
        blackKeyFileIn.clear();
        blackKeyFileIn.close();

        // -------------------------------------------
        // Add blackkey to cuckooTable
        //cout<<"* Add balckkey to cuckooTable!"<<endl;
        float loadFactor = 0.9f;
        int slotNo = 4;
        size_t blackKeySize = CUCKOO_BLACK_SIZE;
        size_t bucketSize = int(blackKeySize/(loadFactor*slotNo))+1;
        int fingerprintNew = 12;
        long MaxKickoutNum = 1000;
        //cuckooBlackKeyTable[si].ClearTable();
        cuckooBlackKeyTable[si].CuckooTableInit(bucketSize,fingerprintNew,slotNo, \
                                                MaxKickoutNum);

        // write to history files
        ofstream blackKeyFileOut;
        blackKeyFileOut.open(blackFileName.c_str());


        /*if(line < 80000)
        {
            for(int i = 0; i < actionSize; i++)
            slotNums[i] = CUCKOO_BLACK_SIZE/4;
        }
        else*/


        //printVec(slotNums);

        for(int i = 0; i < actionSize; i++)
        {
            strings overBigKeys = vBkInfo[si].bks;
            ints overActions = vBkInfo[si].bkActions;

            updateBlacklist(overBigKeys, overActions, rLearn[si][i], flowActionUnique[si][i], blackkeyPres, actionPres,
                            blackKeyFileOut, slotNums[si][i], si, cuckooBlackKeyTable[si]);

            // cuckoo black table occupy rate

        }

        cout<<"* cuckooBlackKeyTable[si].occupyRate(): "<<cuckooBlackKeyTable[si].occupyRate()<<endl;

        outfileR[si]<<endl;

        blackKeyFileOut.clear();
        blackKeyFileOut.close();

        vector<string>().swap(blackkeyPres);
        vector<int>().swap(actionPres);
    }


}

void printQList(RLearn* rLearn)
{
    cout<<"Print Qlist!"<<endl;
    for(int i = 0; i < rLearn->_numS; i++)
    {
        cout<<rLearn->qList[i][0].state<<" ";
        for(int j = 0; j < COLNUM; j++)
        {
            // print aValue
            cout<<rLearn->qList[i][j].qValue<<" ";
        }
        cout<<endl;
    }
}

void selectAction(RLearn** rLearn[], int actionSize, int switchNum, size_tss& slotNums)
{
    floatss actionVec;
    actionVec.resize(switchNum*actionSize);
    for(int si = 0; si < actionVec.size(); si ++)
    {
        actionVec[si].assign(COLNUM,0);
    }

    //if()
    vector<vector<vector<QEntry> > > qVecs;
    qVecs.resize(switchNum);

    for(int si = 0; si < switchNum; si++)
    {
        qVecs[si].resize(actionSize);
    }

    //range
    int iMin[switchNum][actionSize], iMax[switchNum][actionSize];
    for(int si = 0; si < switchNum; si++)
    {
        for(int i = 0; i < actionSize; i++)
        {
            cout<<"* si: "<<si<<" ai: "<<i<<endl;
            int stateIndex = rLearn[si][i]->findIndex(rLearn[si][i]->_states,rLearn[si][i]->_state);
            cout<<"* stateIndex: "<<stateIndex<<endl;

            qVecs[si][i] = (rLearn[si][i]->qList[stateIndex]);

            cout<<"* actionVec: "<<endl;
            for(int ci = 0; ci < COLNUM; ci++)
            {
                actionVec[si*actionSize+i][ci]  = rLearn[si][i]->qList[stateIndex][ci].action;
                //cout<<"* actionVec[si*actionSize+i][ci]: "<<actionVec[si*actionSize+i][ci]<<" ";

            }
            //cout<<endl;

            int randNum = (int(rand())%100);
            float randx = float(randNum)/100.0;
            cout<<"* randNum: "<<randNum<<endl;

            // random action
            if(randx < EPSILON0)
            {
                int randa = -1;
                while(rLearn[si][i]->qList[stateIndex][randa].qValue==-100 || randa == -1/* || randa == 0*/)
                {
                    //if(stateQStart > 0 && stateQStart < _numS-1)
                    randa =rand()%COLNUM;
                    cout<<"* randa: "<<randa<<endl;
                    iMin[si][i] = randa;
                    iMax[si][i] = randa + 1;

                }

            }
            else  // max action
            {
                iMin[si][i] = 0;
                iMax[si][i] = COLNUM;
            }
        }
    }

    /*if(maxHaoOver>0.01)
    {
        iMin[indexHaoOver] = COLNUM-1;
        iMax[indexHaoOver] = COLNUM;
    }

    if(minHaoOver < 0.01 && qVecs[indexMinHaoOver][1].qValue != -100)
    {
        iMin[indexMinHaoOver] = 1;
        iMax[indexMinHaoOver] = 2;
    }
    else if (minHaoOver < 0.01 && qVecs[indexMinHaoOver][2].qValue != -100)
    {
        iMin[indexMinHaoOver] = 2;
        iMax[indexMinHaoOver] = 3;
    }*/
    cout<<"* find max actions!"<<endl;

    vector<vector<QSum> > qSums;
    qSums.resize(switchNum);

    for(int si = 0; si < switchNum; si ++)
    {
        for (int i1 = iMin[si][0]; i1 < iMax[si][0]; i1++)
        {

            for(int i2 = iMin[si][1]; i2 < iMax[si][1]; i2++)
            {

                bool isValid = (qVecs[si][0][i1].qValue != -100) && (qVecs[si][1][i2].qValue != -100) ;
                // no qvalue is -100
                if(isValid)
                {

                    QSum qSum;

                    // state
                    qSum.states[0] = qVecs[si][0][i1].state;
                    qSum.states[1] = qVecs[si][1][i2].state;

                    // action
                    qSum.actions[0] = qVecs[si][0][i1].action;
                    qSum.actions[1] = qVecs[si][1][i2].action;

                    // actionsum
                    qSum.actionSum = qVecs[si][0][i1].action + qVecs[si][1][i2].action;

                    //sum
                    qSum.sum = qVecs[si][0][i1].qValue + qVecs[si][1][i2].qValue;

                    // all actions should be constrainted by the blacklist volume
                    if(qSum.actionSum <= CUCKOO_BLACK_SIZE)
                        qSums[si].push_back(qSum);
                }
            }
        }
    }

    // find the peak assignment
    cout<<"* find slotNums!"<<endl;
    for(int si = 0; si < switchNum; si++)
    {
        if(qSums[si].size() == 0) // the old srategy
        {
            for (int i = 0; i < actionSize; i++)
                slotNums[si][i] = qVecs[si][i][0].action;
        }
        else
        {
            findMax(qSums[si], slotNums[si], actionSize, switchNum);
        }
    }

    // find the index of the selected action
    cout<<"* find the index of selected action!"<<endl;

    for(int si = 0; si < switchNum; si ++)
    {
        printVec(slotNums[si]);
        //printVec(actionVec[si]);
        for(int ai = 0; ai < actionSize; ai++)
        {
            rLearn[si][ai]->_actionSuggestIndex = rLearn[si][ai]->findIndex(actionVec[si*actionSize+ai], (float)slotNums[si][ai]);
            cout<<"* rLearn[si][ai]->_actionSuggestIndex: "<<rLearn[si][ai]->_actionSuggestIndex<<endl;
        }
    }

    vector<vector<QSum> >().swap(qSums);
    vector<vector<vector<QEntry> > >().swap(qVecs);
}

void findMax(vector<QSum>& qSums, size_ts& slotNums, int actionSize, int switchNum )
{
    size_t qSumSize = qSums.size();

    if(qSumSize == 0)
    {
        cout<<"No feasible solutions+++++++++++++++++++++++++++++++++"<<endl;
    }

    float qMax = qSums[0].sum;

    size_t index = 0;

    for(size_t i = 0; i < qSumSize; i++)
    {
        if(qSums[i].sum > qMax)
        {
            qMax = qSums[i].sum;
            index = i;
        }
    }

    for(int i = 0; i < actionSize; i++)
        slotNums[i] = qSums[index].actions[i];
}

void printVec(vector<size_t>& vec)
{
    size_t vecSize = vec.size();
    for(size_t i = 0; i < vecSize; i++)
    {
        cout<<vec[i]<<" ";
    }
    cout<<endl;
}

void loadKeys2Filter(string& inFileName, vector<size_t>& mask, VUPrefix& vuniquePrefix,
                     VUPrefix& vuniqueAggPrefix, char mL0[][ACTIONSIZE][20], char * argv[],
                     int& finger, int& finger0,int switchNum, int actionSize, intss& uniqueActs,
                      CuckooTable* cuckooTableKey, CuckooTable* cuckooBlackKeyTable, CuckooTable* cuckooAggrKeyTable,
                      CuckooFilter* cuckooFilter, CuckooFilter* cuckooFilterInit0)
{

    // ------------------------------
    /*load key file*/
    for(int si = 0; si < switchNum; si++)
    {
        cout<<"* switch: "<<si<<"++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;

        string inFileName0 = inFileName + num2str(si+1);
        std::ifstream infile(inFileName0.c_str());
        if(!infile)
            std::cout << "Train File Error " << std::endl;

        int flowPrefixInt;
        string flowStr;
        int actionInt;

        vector<int> keyprefixlengths;
        vector<string> keys;
        vector<int> keyActions;

        while(infile >> flowPrefixInt >> flowStr>>actionInt)
        {
            keyprefixlengths.push_back(flowPrefixInt);
            keys.push_back(flowStr);
            keyActions.push_back(actionInt);
        }
        cout<<"* Key size: "<<keys.size()<<endl;
        infile.clear();
        infile.close();

        // --------------------------------
        // get unique prefix length
        vector<int> uniquePrefix;
        vector<int> uniqueAggPrefix;
        prefixNum(keyprefixlengths, uniquePrefix);

        // --------------------------------
        /*cuckoo table*/
        // load factor
        // m: key number, f: fingerprint width, bc: slots num per bucket,
        // MaxNumKicks: max kickout number in a cuckoo filter
        cout<<"* Init cuckoo table ... ..."<<endl<<endl;
        float a;
        int f,bc;
        a = 0.9;
        bc = 4;
        long m = long(keys.size()/(a*bc))+1;
        f = 12;
        long MaxNumKicks = 1000;
        cuckooTableKey[si].ClearTable();
        cuckooTableKey[si].CuckooTableInit(m,f,bc,MaxNumKicks);

        // --------------------------------
        // Add original key to cuckoo table
        cout<<"* Add key to cuckoo table ... ..."<<endl;
        bool isAddTable;
        for (int i = 0; i < keys.size(); i++)
        {
            isAddTable = cuckooTableKey[si].AddKeyPrefix(keys[i],(keyprefixlengths[i]),keyActions[i]);

            if(isAddTable == 0)
            {
                cout<<"* Flag_add fail"<<endl;
                cout<<"* Order: "<<i<<"  ";
            }

        }

        // -----------------------------------------------
        // init cuckooFilter for flow estimation
        long flowEstSize = FLOW_EST_SIZE;
        m = flowEstSize/(a*bc)+1;
        f = 16;
        //cuckooFilterFlowEst[si].ClearTable();
        //cuckooFilterFlowEst[si].cuckooFilterInit(m,f,bc,MaxNumKicks);

        // init black table
        m = CUCKOO_BLACK_SIZE/(a*bc)+1;
        cuckooBlackKeyTable[si].ClearTable();
        cuckooBlackKeyTable[si].CuckooTableInit(m,f,bc,
                                                MaxNumKicks);
        // -----------------------------------
        // Init blackkey file
        cout<<"* Write blackkey to file!"<<endl;
        ofstream blackKeyFileOut;
        BLACKFILENAME = "../test/blackkeyfile_" + string(argv[2]) + '_' +
                        string(argv[4])+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6]);
        //for(int si =0; si <switchNum; si++)
        {
            string bkFile = BLACKFILENAME + num2str(si);
            blackKeyFileOut.open(bkFile.c_str());
            blackKeyFileOut.clear();
            blackKeyFileOut.close();
        }



        // -----------------------------------
        // Init aggr file
        cout<<"* Write aggr to file!"<<endl;
        ofstream aggrFileOut;
        AGGRFILENAME = "../test/aggrfile" + string(argv[2]) + '_' +
                       string(argv[4])+ "_tstNum_"+ string(argv[5])+"_b"+string(argv[6]);
        aggrFileOut.open(AGGRFILENAME.c_str());
        aggrFileOut.clear();
        aggrFileOut.close();

        // -----------------------------------------------
        // Two counters for estimator
        int2s bigNonCounts;
        long2s timeCounts;
        bigNonCounts= vector<vector<int> > (m, vector<int>(bc, 0));
        timeCounts = vector<vector<long> > (m, vector<long>(bc, 0));

        // ---------------------------------------------
        /* init cuckoo filter without aggregation*/
        cout<<"* Init cuckoo filter ... ..."<<endl;

        float storage = strtof(argv[2],NULL); // storage size
        //int finger = 0;
        long but = 200000/3;
        //cout<<sizeof(mL0)<<endl;
        //return 0;

        initCuckoo(keys,keyprefixlengths,keyActions,storage, finger, mL0, cuckooFilter[si],cuckooFilterInit0[si]);
        finger0 = finger;

        cout<<"* finger: "<<finger<<endl;

        // ------------------------------------------------
        // Init aggregation
        bool isInit = 1;

        initAggregation(keys,keyprefixlengths,keyActions, mask, actionSize, storage,
        isInit, finger,uniqueAggPrefix,mL0,cuckooFilter[si], cuckooAggrKeyTable[si], uniqueActs[si]);

        // output
        UPrefix uPrefix;
        uPrefix.uPres  = uniquePrefix;
        vuniquePrefix.push_back(uPrefix);

        uPrefix.uPres.clear();
        uPrefix.uPres  = uniqueAggPrefix;

        vuniqueAggPrefix.push_back(uPrefix);

        vector<int>().swap(keyprefixlengths);
        vector<string>().swap(keys);
        vector<int>().swap(keyActions);
    }// end si
}

double nextTime(double rateParameter)
{
    double nTime = -log(1.0f - (double) rand() / double(RAND_MAX + 1.0)) / rateParameter;

    if(nTime == INFINITY)
    {
        double nTime = -log(1.0f - (double) rand() / (RAND_MAX + 1.0)) / rateParameter;
        cout<<"* nextTime: INFINITY!"<<endl;
        exit(0);
    }
    return nTime;
}





