/**
cuckoo_filter.cpp

create by: Yunhong
create time: 05/22/2015
update time: 05/25/2015
*/


#include "cuckoo_filter.h"
using namespace std;

/**
constructor: init variables
*/
CuckooFilter::CuckooFilter(int m,int f,int bc,int MaxNumKicks)
{
    // Set bucket count
    mm = m;

    // Start with empty hash table
    mL = vector<vector<int> > (m, vector<int>(bc, 0));
    mhigh = vector<vector<bool> > (m, vector<bool>(bc, 0));

    mf = f;
    mbc = bc;

    // maximum kickout number
    mMaxNumKicks = MaxNumKicks;

    // Initialize hash function
    InitHashFunction();

}

CuckooFilter::~CuckooFilter()
{

}

/**
add keys to cuckoo filter
    //-------------------------- Add -------------------------------
    // Syntax:       H.Add(key,value);
    //
    // Inputs:       key is an alphanumeric character array
    //
    //               value is an arbitrary object
    //
    // Description:  Inserts the given key-value pair into H
    //--------------------------------------------------------------

*/
bool CuckooFilter::AddKey(const std::string& key)
{

    //define variables
    long f_key,t_key;
    long h1,h2,ci,ht;
    int I1,I2,It,high_bit_t,high_bit,e_key_valid,e_key_high_bit;
    int ri;
    vector<int> vec1,vec2;

    // compute fingerprint
    f_key = Fingerprint(key);

    // if fingerprint is mm, print it
    if (f_key == 0)
        cout<<"fingerprint is mm"<<endl;

    // Compute the first hash code
    h1 = HashCode(key);

    // compute the second hash code
    h2 = (h1 ^ HashCode(num2str(f_key)));

    high_bit = 0;
    if (h2 >= mm)
    {
        h2 = (h2%mm);
        high_bit = 1;

    }

    // find the first vacant position in table
    vec1 = mL[h1];
    vec2 = mL[h2];
    I1 = find(vec1.begin(),vec1.end(),0)-vec1.begin();//find(~L(:,h1),1); // the first 0 in bucket h1
    I2 = find(vec2.begin(),vec2.end(),0)-vec2.begin();//find(~L(:,h2),1);

    // if find vacant position in h1, store the item and return
    if (I1 != mbc)//~isempty(I1) ----------------------------------------------------
    {
        mL[h1][I1] = f_key;
        mhigh[h1][I1] = 0;
        return true;
    }
    // if find vacant position in h2, store the item and return
    else if (I2 != mbc)//~isempty(I2)
    {
        mL[h2][I2] = f_key;
        mhigh[h2][I2] = high_bit;
        return true;
    }

    // fail to insert keys
    t_key = f_key;
    high_bit_t = high_bit;
    // get a bucket randomly
    if (rand() % 2== 0)
    {
        ci = h1;
        high_bit_t = 0;
    }
    else
    {
        ci = h2;
    }


    for (int ni = 1; ni < mMaxNumKicks; ni++)
    {
        // find a slot randomly
        ri = rand()%mbc;
        //kick out the key
        e_key_valid = mL[ci][ri];
        e_key_high_bit = mhigh[ci][ri];

        // put the key inside
        mL[ci][ri] = t_key;// put f_key to the selected bucket
        mhigh[ci][ri] = high_bit_t;

        // compute the second bucket for e_key
        ht = (ci+e_key_high_bit*(mm)) ^ HashCode(num2str(e_key_valid));
        high_bit_t = 0;
        if (ht >= mm)
        {
            ht = (ht%mm);
            high_bit_t = 1;
        }

        // find a vacant slot for e_key
        vec1 = mL[ht];
        It = find(vec1.begin(),vec1.end(),0)-vec1.begin();//find(~L(:,h1),1); // the first 0 in bucket h1

        if (It != mbc)//~isempty(It)
        {
            mL[ht][It] = e_key_valid;
            mhigh[ht][It] = high_bit_t;
            return true;
        }

        //fail to find, continue to kick out other items
        t_key = e_key_valid;
        //high_bit_t = e_key_high_bit;
        ci = ht;

    }

    return false;
    //             disp('insertion failure');

}
/**
lookup keys from cuckoo filter
    //---------------------- ContainsKey ---------------------------
    // Syntax:       bool = H.ContainsKeys(key);
    //
    // Inputs:       key is an alphanumeric character array
    //
    // Outputs:      bool = {true,false}
    //
    // Description:  Determines if H contains a key-value pair with
    //               the given key
    //--------------------------------------------------------------
    // Algorithm 2: Lookup(x)
    // f = fingerprint(x);
    // i1 = hash(x);
    // i2 = i1 ⊕ hash(f);
    // if bucket[i1] or bucket[i2] has f then
    // return True;
    // return False;
*/
bool CuckooFilter::LookUpKey(const std::string& key)
{

    //define variables
    long f_key,t_key;
    long h1,h2;

    // compute fingerprint
    f_key = Fingerprint(key);

    if (f_key == 0)
        cout<<"fingerprint is mm"<<endl;

    // Compute the first hash code
    h1 = HashCode(key);

    // compute the second hash code
    h2 = (h1 ^ HashCode(num2str(f_key)));
    if (h2 >= mm)
    {
        h2 = (h2%mm);

    }


    // Search for key in list(hk)
    for (int i = 0; i<mbc; i++)
    {
        if ((mL[h1][i]) == f_key | (mL[h2][i]) == f_key)
        {
            return true;
        }

    }

    return false;
}
/**
remove a key from cuckoo filter
*/
bool CuckooFilter::RemoveKey(const std::string& key)
{

    return true;
}
void CuckooFilter::ClearTable()
{

}

/**
universal hash function implementation
    //
    // Compute hash of integer vector
    //
    // Reference: http://en.wikipedia.org/wiki/Universal_hashing
    //            Sections: Hashing integers
    //                      Hashing strings
    //
*/
long CuckooFilter::HashCode(const std::string& key)
{
    //define variables
    int ll;
    long hk;
    vector<int> key_double;

    // Convert character array to integer array
    ll = key.size();

    key_double.reserve(ll);
    for(int i = 0; i < ll; i++)
    {
        key_double[i] = int(key[i])-45;

    }


    hk = (key_double[0]);
    for (int i = 1; i<ll; i++)
    {
        // Could be implemented more efficiently in practice via bit
        // shifts (see reference)
        hk = (mc * hk + key_double[i])%mp;
    }

    hk = ((ma * hk + mb)%mp)%mm;

    return hk;


}
/**
fingerprint calculation
    //
    // Compute hash of integer vector
    //
    // Reference: http://en.wikipedia.org/wiki/Universal_hashing
    //            Sections: Hashing integers
    //                      Hashing strings
    //
*/
long CuckooFilter::Fingerprint(const std::string& key)
{
    //define variables
    int ll;
    long hk;
    vector<int> key_double;
    // Convert character array to integer array
    ll = key.size();

    key_double.reserve(ll);
    for(int i = 0; i < ll; i++)
    {
        key_double[i] = int(key[i])-45;
    }

    hk = (key_double[0]);
    for (int i = 1; i<ll; i++)
    {
        // Could be implemented more efficiently in practice via bit
        // shifts (see reference)
        hk = (mcf * hk + key_double[i])%mpf;

    }

    hk = ((maf * hk + mbf)%mpf)%int(pow(2,mf))+1;
    return hk;
}
/**
init hash funciton
*/
bool CuckooFilter::InitHashFunction()
{
    //define variables
    int ff;
    int pp;
    // Set prime parameter
    ff = 1000; // fudge factor
    pp = ff * max(mm + 1,76);
    pp = pp + ~(pp%2); // make odd
    while (IsPrime(pp) == false)
        pp = pp + 2;

    mp = pp; // sufficiently large prime number

    /* initialize random seed: */
    srand (time(NULL));

    /* generate secret number between 1 and 10: */
    // Randomized parameters
    ma = rand() % pp;
    mb = rand() % pp;
    mc = rand() % pp;
    while(ma == 0)
    {
        ma = rand() % pp; //ma!=0
    }
    while(mc == 0)
    {
        mc = rand() % pp;//randi([1,(pp - 1)]); mc != 0
    }

    // hash for fingerprint
    // Set prime parameter
    pp = ff * max(int(pow(2,mf) + 1),76);
    pp = pp + ~(pp%2); // make odd
    while (IsPrime(pp) == false)
        pp = pp + 2;

    mpf = pp; // sufficiently large prime number

    /* generate secret number between 1 and pp: */
    // Randomized parameters
    maf = rand() % pp;
    mbf = rand() % pp;
    mcf = rand() % pp;
    while(maf == 0)
    {
        maf = rand() % pp; //ma!=0
    }
    while(mcf == 0)
    {
        mcf = rand() % pp;//randi([1,(pp - 1)]); mc != 0
    }
    return true;

}

