/**
cuckooFilterMain.cpp

create by: Yunhong
create time: 05/22/2015
update time: 05/25/2015
update time: 06/09/2015
*/


#include "cuckoo_filter.h"
#include <arpa/inet.h>

using namespace std;

int main(int argc, char * argv[])
{
    //define variables
    std::vector<std::string> key;
    std::vector<string> key_prefix;
    std::vector<std::string> testIP;
    std::vector<std::string> testIP_cnt;

    float a;                        // load factor
    int n,m,f,bc,MaxNumKicks;       // parameters for cuckoo filter
    long lines;
    long numIP,countIP;
    string str;
    int flag_add,flag_look;

    FILE * fp;
    struct timeval gen_start, gen_end; /* gettimeofday stuff */
    struct timezone tzp;

    // get start time
    gettimeofday(&gen_start, &tzp);

    //---------------------------------------------------------------------------------------------
    // load file from csv file
    // load key
    key.clear();
    fp = fopen("/home/yunhong/Documents/research/PCAP/BT_files/BT_Originalkey_file.csv", "r" );
    if( ! fp )
        return 1;
    lines = 300000;
    while(lines > 0)
    {
        readCSVKey( fp, key,key_prefix);
        lines --;

    }

    //-------------------------------------------------------------------------------------
    // load test_ip
    testIP.clear();
    fp = fopen("/home/yunhong/Documents/research/PCAP/BT_files/BT_flow_file.csv", "r" );
    if( ! fp )
        return 1;
    lines = 409280;
    while(lines > 0)
    {
        readCSV( fp, testIP,testIP_cnt );
        lines --;

    }
    numIP = testIP.size();
    //-------------------------------------------------------------------------------------------
    cout<<key[1]<<endl;
    cout<<testIP[1]<<endl;
    if(key[1]==testIP[1])cout<<"equal"<<endl;
    a = 0.90;
    bc = 4;
    n = key.size();
    m = int(n/(a*bc));
    f = 8;
    MaxNumKicks = 1000;

    // init cuckoo filter
    CuckooFilter cuckooFilter(m,f,bc,MaxNumKicks);

    std::stringstream s;
    for(int i = 0; i < n; i++)
    {
        s.str("");
        s<<key[i]<<"/32";
        str = s.str();
        // add keys to cuckoo filter
        flag_add = cuckooFilter.AddKey(str);

        if(flag_add == 0)
        {
            cout << "flag_add"<<"0"<<endl;
            break;
        }

    }

    //----------------------------------------------------------------
    //look up test IPs from cuckoo filter
    int countNum = 0;
    bool isKeyFlag = 1;
    countIP = 0;
    flag_look = 0;
    for(int i = 0; i < numIP; i++)
    {
        s.str("");
        s<<testIP[i]<<"/32";
        str = s.str();

        flag_look = cuckooFilter.LookUpKeyCol(str);
        if (flag_look)
        {
            countNum += flag_look;
            countIP += flag_look*str2num(testIP_cnt[i]);
        }

    }
    //------------------------------------------
    // get the end time
    gettimeofday(&gen_end, &tzp);
    print_elapsed("cuckooFilter: ", &gen_start, &gen_end, 1);
    cout<<endl;

    //---------------------------------------------------------------------------
    cout <<"numIP  "<<testIP.size()<<endl;
    cout<<"keynum  "<<key.size()<<endl;
    cout<<"countIP  "<<countIP<<endl;

    float flase_pos = float(countIP-51434265)/float(70938691-51434265);
    float hao_flase_pos = float(countIP-51434265)/float(51434265);

    cout<<"flase pos:  "<<flase_pos<<"  hao pos:  "<<hao_flase_pos<<"  countNUM:  "<<countNum<<endl;
    cout << endl;


    return 1;
}
