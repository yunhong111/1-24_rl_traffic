/**
otherFun.cpp

create by: Yunhong
create time: 05/22/2015
update time: 05/25/2015
*/

#include "otherFun.h"
#include <arpa/inet.h>

double BIGNONKEYTHLD = 0.00000769230769230769;//0.00000169558403281403;//0.00000070483398121908;//0.00000140966796243816;//0.00117;100 packets
double BIGKEYTHLD = 0.85;
double KEYTHLD = 1;
double GLOBAL_INITKEYTHLD = 0.00012;
int GLOBAL_BIGNONKEYNUM = 3;
double gFPR;
double nonkey_weight_sum = 0;
double key_weight_sum = 0;
double para_k = 0.0000012;
double diff_FPR = 0;
double g_ebsr = 200;
double g_blackkey_weight_sum = 0;
ofstream outfile;
size_t GLOBAL_KEY_SIZE = 0;
size_t FLOW_EST_COUNT = 0;
size_t FLOW_EST_SIZE =6000;
size_t CUCKOO_BLACK_SIZE = 0;
int AGGR_PREFIX = 0;
string BLACKFILENAME;
string AGGRFILENAME;

vector<double> g_vweightThld;
vector<double> g_vnonkey_weight_sum; // nonkey packet ratio
vector<double> g_vnonkey_weight_sum0; // initial nonkey packet ratio
vector<double> g_vkey_weight_sum; // nonkey packet ratio
vector<double> g_vgFPR; // nonkey packet ratio
vector<double> g_vpara_k;
vector<long> g_vcountkey;
vector<long> g_vcountblackkey;
vector<long> g_vcountkey_init;
vector<double> g_vblackkey_weight_sum;

using namespace std;


/** controler for parameters
*/

float ctlPara(double givenFPR, double expFPR, float threshold_old, float threshold_old_old, int action)
{
    float threMax = 1;
    float threMin = 0.1;
    float threshold = threshold_old;
    float ctl_value = 0.05;
    ctl_value = g_vpara_k[action]*(givenFPR - expFPR);
    //para_k = 2;
    if((givenFPR - expFPR)*diff_FPR < 0)
    {
        threshold_old = threshold_old_old;
        g_vpara_k[action] /= 2;
        threshold = threshold_old;
        cout<<"divided by 2"<<" threshold"<<threshold<<endl;
        cout<<endl;
    }
    else{
        while(fabs(ctl_value)>0.1)
        {
            ctl_value = g_vpara_k[action]*(givenFPR - expFPR);
            if (g_vpara_k[action] > 0.000001)
            g_vpara_k[action] -= 0.000001;
            else if (g_vpara_k[action] > 0.0000001)
            g_vpara_k[action] -= 0.0000001;
            else
            break;
        }


        if ((givenFPR - expFPR) > g_ebsr)
        {
            threshold = threshold_old+ctl_value;
        }
        else if ((givenFPR - expFPR)<-g_ebsr)
        {
            threshold = threshold_old+ctl_value;
        }

    }

    if (threshold > threMax)
    {
        threshold = threMax;
    }
    else if (threshold < threMin)
    {
        threshold = threMin;
    }

    return threshold;
}
/**==============================================================
 * print_elapsed (prints timing statistics)
 *==============================================================*/
long print_elapsed(const char* desc, struct timeval* start, struct timeval* end, int niters)
{

    struct timeval elapsed;
    /* calculate elapsed time */
    if(start->tv_usec > end->tv_usec)
    {

        end->tv_usec += 1000000;
        end->tv_sec--;
    }
    elapsed.tv_usec = end->tv_usec - start->tv_usec;
    elapsed.tv_sec  = end->tv_sec  - start->tv_sec;

    long timeInv = (elapsed.tv_sec*1000000 + elapsed.tv_usec);

    return timeInv;

    //printf("\n %s total elapsed time = %ld (sec)",
          // desc, (elapsed.tv_sec*1000000 + elapsed.tv_usec) / (1000000*niters));
    //outfile<<"Time: "<<((elapsed.tv_sec*1000000 + elapsed.tv_usec) / (1000000*niters))<<endl;
}


/**
check whether an integer is a prime
*/
bool IsPrime(int number)
{
    // Given:   num an integer > 1
    // Returns: true if num is prime
    // 			false otherwise.

    int i;

    for (i=2; i<number; i++)
    {
        if (number % i == 0)
        {
            return false;
        }
    }

    return true;
}

/**
convert an integer to a string
*/
std::string num2str(int num)
{
    std::stringstream ss;
    ss << num;
    return ss.str();
}

/**
convert a string to an integer
*/
int str2num(string str)
{
    int numb;
    istringstream ( str ) >> numb;
    return numb;
}

/**

  Read line from a CSV file

  @param[in] fp file pointer to open file
  @param[in] vls reference to vector of strings to hold next line

  */
void readCSV( FILE *fp, std::vector<std::string>& vls, std::vector<string>& vcnt )
{
    //vls.clear();
    if( ! fp )
        return;
    char buf[10000];
    if( ! fgets( buf,999,fp) )
        return;
    std::string s = buf;
    int p,q,t;
    q = -1;
    // loop over columns
    while( 1 )
    {
        p = q;
        q = s.find_first_of(", ",p+1);
        if( q == -1 )
            break;
        vls.push_back( s.substr(p+1,q-p-1) );
        t = s.find_first_of(",\n",p+1);
        vcnt.push_back( s.substr(q+1,t-q-1) );
    }
}

void readCSVKey( FILE *fp, std::vector<std::string>& vls, std::vector<string>& vcnt )
{
    //vls.clear();
    if( ! fp )
        return;
    char buf[10000];
    if( ! fgets( buf,999,fp) )
        return;
    std::string s = buf;
    int p,q,t;
    q = -1;
    // loop over columns
    while( 1 )
    {
        p = q;
        q = s.find_first_of("/",p+1);
        if( q == -1 )
            break;
        vls.push_back( s.substr(p+1,q-p-1) );
        t = s.find_first_of(",\n",p+1);
        vcnt.push_back( s.substr(q+1,t-q-1) );
    }
}

void readCSVKeySlash( FILE *fp, std::vector<std::string>& vls)
{
    //vls.clear();
    if( ! fp )
        return;
    char buf[10000];
    if( ! fgets( buf,999,fp) )
        return;
    std::string s = buf;
    int p,q,t;
    q = -1;
    // loop over columns
    while( 1 )
    {
        p = q;
        q = s.find_first_of(",\n",p+1);
        if( q == -1 )
            break;
        vls.push_back( s.substr(p+1,q-p-1) );
    }
}

void readCSVcnt( FILE *fp, std::vector<std::string>& vls, std::vector<std::string>& vpr, std::vector<string>& vcnt )
{
    //vls.clear();
    if( ! fp )
        return;
    char buf[10000];
    if( ! fgets( buf,999,fp) )
        return;
    std::string s = buf;
    int p,q,r,t;
    q = -1;
    // loop over columns
    //while( 1 )
    {
        //cout<<s<<endl;
        p = q;
        q = s.find_first_of(", ",p+1);
        if( q == -1 )
            return;
        vls.push_back(s.substr(p+1,q-p-1));

        r = s.find_first_of(", ",q+1);
        vpr.push_back(s.substr(q+1,r-q-1));

        t = s.find_first_of(",\n",r+1);
        vcnt.push_back( s.substr(r+1,t-r-1) );
        //cout<<p<<"  "<<q<<"  "<<r<<"  "<<t<<endl;
    }
}
void readCSVaction( FILE *fp, std::vector<std::string>& vls, std::vector<std::string>& vpr, std::vector<string>& action, std::vector<string>& vcnt )
{
    //vls.clear();
    if( ! fp )
        return;
    char buf[10000];
    if( ! fgets( buf,999,fp) )
        return;
    std::string s = buf;
    int p,q,r,t,z;
    q = -1;
    // loop over columns
    //while( 1 )
    {
        //cout<<s<<endl;
        p = q;
        q = s.find_first_of(", ",p+1);
        if( q == -1 )
            return;
        vls.push_back(s.substr(p+1,q-p-1));

        r = s.find_first_of(", ",q+1);
        vpr.push_back(s.substr(q+1,r-q-1));

        t = s.find_first_of(", ",r+1);
        action.push_back( s.substr(r+1,t-r-1) );

        z = s.find_first_of(",\n",t+1);
        vcnt.push_back( s.substr(t+1,z-t-1) );
        //cout<<p<<"  "<<q<<"  "<<r<<"  "<<t<<endl;
    }
}

void readCSVactioncnt( FILE *fp, std::vector<std::string>& vls, std::vector<std::string>& vpr, std::vector<string>& action, std::vector<string>& keytype, std::vector<string>& vcnt )
{
    //vls.clear();
    if( ! fp )
        return;
    char buf[10000];
    if( ! fgets( buf,999,fp) )
        return;
    std::string s = buf;
    int p,q,r,t,z,z1;
    q = -1;
    // loop over columns
    //while( 1 )
    {
        //cout<<s<<endl;
        p = q;
        q = s.find_first_of(", ",p+1);
        if( q == -1 )
            return;
        vls.push_back(s.substr(p+1,q-p-1));

        r = s.find_first_of(", ",q+1);
        vpr.push_back(s.substr(q+1,r-q-1));

        t = s.find_first_of(", ",r+1);
        action.push_back( s.substr(r+1,t-r-1) );

        z = s.find_first_of(", ",t+1);
        keytype.push_back( (s.substr(t+1,z-t-1)) );

        z1 = s.find_first_of(",\n",z+1);
        vcnt.push_back( s.substr(z+1,z1-z-1) );
        //cout<<p<<"  "<<q<<"  "<<r<<"  "<<t<<endl;
    }
}
/**
IPv4 to dec
*/
uint32_t parseIPV4string(const char* ipAddress)
{

    uint32_t ipbytes[4];
    sscanf(ipAddress, "%d.%d.%d.%d", &ipbytes[3], &ipbytes[2], &ipbytes[1], &ipbytes[0]);

    return ipbytes[0] | ipbytes[1] << 8 | ipbytes[2] << 16 | ipbytes[3] << 24;
}

/**
dec to IPv4
*/
string parsedec2IPV4(uint64_t dec)
{
    struct in_addr addr;
    addr.s_addr = dec;
    char* ip_str = inet_ntoa(addr);
    int ipbytes[4];
    sscanf(ip_str, "%d.%d.%d.%d", &ipbytes[3], &ipbytes[2], &ipbytes[1], &ipbytes[0]);
    std::stringstream s;
    char ch='.'; //to temporarily store the '.'
    s << ipbytes[0] << ch << ipbytes[1] << ch << ipbytes[2] << ch << ipbytes[3];
    return s.str();

}

int char_to_ipv4(char *src)
{
    unsigned int ipv4;
    unsigned int octet;
    int count;
    int i;

    ipv4 = 0;
    octet = 0;
    count = 0;
    for ( i = 0; i < 16; i++, src++ )
    {
        if ( *src == '.' || *src == 0 || i == 15 )
        {
            if ( octet > 255 )
            {
                ipv4 = 0;
                break;
            }

            ipv4 <<= 8;
            ipv4 += octet;
            octet = 0;
            count++;
            if ( *src == 0 )
                break;
        }
        else if ( *src < '0' || *src > '9' )
        {
            ipv4 = 0;
            break;
        }
        else
        {
            octet *= 10;

            octet += (*src - '0');

        }

    }

    if ( count != 4 )
        ipv4 = 0;
    return ipv4;
}

/**
dec to binary(32 bits)
*/
string DecToBin(uint64_t dec)
{
    char bin32[]  = "00000000000000000000000000000000";
    for (int pos = 31; pos >= 0; --pos)
    {
        if (dec % 2)
            bin32[pos] = '1';
        dec /= 2;
    }
    return bin32;
}
/**
binary to dec
*/
uint32_t binConv(string binNum)
{
    uint32_t decimal = 0;
    int size_bin = binNum.size();

    for (int counter = 0; counter < size_bin; counter++)
    {
        if (binNum[counter] == '1')
            decimal = (decimal + pow(2.0,counter));
        else
            decimal = (decimal + 0);
    }
    return decimal;
}
/**
IPv4 to binary (32 bits)
*/
string parseIPV42bin(const char* ipAddress)
{
    uint32_t ip = parseIPV4string(ipAddress);
    return DecToBin(ip);
}


/**
binary to IPv4
*/
string parsebin2IPV4(const char* ipAddress)
{
    struct in_addr addr;
    addr.s_addr = binConv(ipAddress);
    char* ip_str = inet_ntoa(addr);
    int ipbytes[4];
    sscanf(ip_str, "%d.%d.%d.%d", &ipbytes[3], &ipbytes[2], &ipbytes[1], &ipbytes[0]);
    std::stringstream s;
    char ch='.'; //to temporarily store the '.'
    s << ipbytes[0] << ch << ipbytes[1] << ch << ipbytes[2] << ch << ipbytes[3];
    return s.str();
}
/**
control real rules to a range of prefix
*/
void ctlPrefix(vector<string> &key, vector<string>& key_cnt,vector<string> &keyprefix,vector<keyType> &keytype,vector<string> &action,vector<string> &prefix,
               vector<string> &key_out,vector<string>& key_cnt_out,vector<string> &keyprefix_out,vector<keyType> &keytype_out,vector<string> &action_out )
{
    string prefixtmp;
    uint32_t keytmp;
    int cnttmp,sumtmp = 0;
    bool equalPrefix =  0;
    for (int i = 0; i < key.size(); i ++)
    {
        equalPrefix =  0;
        // find the shortest prefix length
        for (int pi = 0; pi < prefix.size(); pi++)
        {

            if (keyprefix[i]==prefix[pi]) // don't need to extent
            {
                //cout<<"="<<endl;
                key_out.push_back(key[i]);
                keyprefix_out.push_back(keyprefix[i]);
                key_cnt_out.push_back(key_cnt[i]);
                keytype_out.push_back(keytype[i]);
                action_out.push_back(action[i]);
                prefixtmp = keyprefix[i];
                equalPrefix =  1;
                break;
            }
        }

        if(equalPrefix) continue;

        for (int pi = 0; pi < prefix.size(); pi++)
        {

            if (str2num(keyprefix[i])<str2num(prefix[pi]))
            {
                //cout<<"<"<<endl;
                prefixtmp = prefix[pi];
                //cout<<"tmp: "<<prefixtmp<<endl;
                //cout<<"prefix: "<<keyprefix[i]<<endl;
                //cout<<int(pow(2,str2num(prefixtmp)-str2num(keyprefix[i])) )<<endl;
                break;
            }
        }


        // extent one prefix to the shortest prefix
        sumtmp = 0;
        for(int pi = 0; pi < int(pow(2,str2num(prefixtmp)-str2num(keyprefix[i])) ); pi++)
        {

            keytmp = parseIPV4string(key[i].c_str())+pi*(1<<(32-str2num(prefixtmp)));
            key_out.push_back(parsedec2IPV4(keytmp));
            keyprefix_out.push_back(prefixtmp);
            keytype_out.push_back(keytype[i]);
            action_out.push_back(action[i]);
            cnttmp = (str2num(key_cnt[i])/pow(2,str2num(prefixtmp)-str2num(keyprefix[i])));

            // the last prefix
            if(pi == int(pow(2,str2num(prefixtmp)-str2num(keyprefix[i])))-1)
            {
                key_cnt_out.push_back(num2str((str2num(key_cnt[i])-sumtmp)));

            }
            else // not the lasr prefix
            {
                key_cnt_out.push_back(num2str(cnttmp));
            }
            sumtmp += cnttmp;
        }

    }
}

/**
get random keys from flow
*/
void randKey(long& n, long& ip_num, vector<string>& flow, vector<string>& flow_cnt, vector<string> &key, vector<string> &key_cnt,vector<keyType> keytype)
{
    // randomly generate keys from flows
    // initialize rand
    /* initialize random seed: */
    long lines = n;//100000;
    srand (time(NULL));
    vector<long> index;

    index.clear();
    float factor = 1.05;
    while(index.size()<n)
    {
        index.clear();
        for(int i = 0; i < lines*factor; i++)
        {
            long ith = rand()%ip_num;
            index.push_back(ith);
        }

        std::sort(index.begin(),index.end());

        {
            std::vector<long>::iterator it;
            it = std::unique (index.begin(), index.end());
            index.resize( std::distance(index.begin(),it));
        }
        factor = factor+0.0005;
    }
    cout<<factor<<endl;
    key.clear();
    key_cnt.clear();
    std::ofstream myfile;
    myfile.open ("/home/yunhong/Documents/research/PCAP/BT_key_file_index_prefix.csv"); // write to file: BT_key_file_index_2_3
    for(int i = 0; i < lines; i++)
    {

        key.push_back(flow[index[i]]);
        key_cnt.push_back(flow_cnt[index[i]]);
        keytype[index[i]] = iskey;
        cout << key[i]<<","<<index[i]<<"\n";
        myfile << key[i]<<","<<index[i]<<"\n";

    }
    myfile.close();
}



void clearglobal_value()
{
    g_vweightThld.clear();
    g_vnonkey_weight_sum.clear(); // nonkey packet ratio
    g_vnonkey_weight_sum0.clear(); // initial nonkey packet ratio
    g_vkey_weight_sum.clear(); // nonkey packet ratio
    g_vgFPR.clear(); // nonkey packet ratio
    g_vpara_k.clear();
    g_vcountkey.clear();
    g_vcountblackkey.clear();
    g_vcountkey_init.clear();
    g_vblackkey_weight_sum.clear();

}



void setGlobalValue(int actionSize)
{
    g_vweightThld.assign(actionSize,BIGKEYTHLD);
    g_vnonkey_weight_sum.assign(actionSize,0); // nonkey packet ratio
    g_vkey_weight_sum.assign(actionSize,0); // nonkey packet ratio
    g_vgFPR.assign(actionSize,0);
    g_vpara_k.assign(actionSize,para_k);
    g_vnonkey_weight_sum0.assign(actionSize,0);
    g_vcountkey.assign(actionSize,0);
    g_vcountblackkey.assign(actionSize,0);
    g_vcountkey_init.assign(actionSize,0);
    g_vblackkey_weight_sum.assign(actionSize,0);

}

double gaussrand(double miu, double sigma)
{
	static double V1, V2, S;
	static int phase = 0;
	double X;


	if(phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return (miu+sigma*X);
}
