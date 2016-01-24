/**
cuckoo_filter.cpp

create by: Yunhong
create time: 05/22/2015
update time: 05/25/2015
*/


#include "cuckoo_table.h"
using namespace std;

/**
constructor: init variables
*/
CuckooTable::CuckooTable()
{


}
void CuckooTable::CuckooTableInit(long m,int f,int bc,long MaxNumKicks)
{
    // Set bucket count
    mm = m;

    // Start with empty hash table


    if(mL.size() == 0)
    {
        mL = vector<vector<uint32_t> > (m, vector<uint32_t>(bc, 0));
        maction = vector<vector<int> > (m, vector<int>(bc, 0));
        mprefix= vector<vector<int> > (m, vector<int>(bc, 0));
        mhigh = vector<vector<bool> > (m, vector<bool>(bc, 0));
        mcount = vector<vector<size_t> > (m, vector<size_t>(bc, 0));
    }
    else
    {
        for(size_t i = 0; i < m; i++)
        {
            for(int j = 0; j < bc; j++)
            {
                mL[i][j] = 0;
                maction[i][j] = 0;
                mhigh[i][j] = 0;
                mprefix[i][j] = 0;
                mcount[i][j] = 0;
            }
        }

    }

    mf = f;
    mbc = bc;

    // maximum kickout number
    mMaxNumKicks = MaxNumKicks;

    // Initialize hash function
    InitHashFunction();

}

CuckooTable::~CuckooTable()
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
bool CuckooTable::AddKey(const string& key,int action)
{

    //define variables
    uint32_t f_key,t_key;
    long h1,h2,ci,ht;
    int I1,I2,It,high_bit_t,high_bit,e_key_valid,e_key_high_bit,e_action,t_action;
    int ri;
    vector<uint32_t> vec1,vec2;

    // compute fingerprint
    //f_key = Fingerprint(key);
    f_key = parseIPV4string(key.c_str());

    // if fingerprint is mm, print it
    //if (f_key == 0)
       // cout<<"fingerprint is mm"<<endl;

    // Compute the first hash code
    h1 = HashCode(key);

    // compute the second hash code
    //h2 = (h1 ^ HashCode(num2str(f_key)));
    h2 = HashCode((key+num2str(mm*mbc)));

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
        maction[h1][I1] = action;
        return true;
    }
    // if find vacant position in h2, store the item and return
    else if (I2 != mbc)//~isempty(I2)
    {
        mL[h2][I2] = f_key;
        mhigh[h2][I2] = high_bit;
        maction[h2][I2] = action;
        return true;
    }

    // fail to insert keys
    t_key = f_key;
    high_bit_t = high_bit;
    t_action = action;
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


    for (uint16_t ni = 1; ni < mMaxNumKicks; ni++)
    {
        // find a slot randomly
        ri = rand()%mbc;
        //kick out the key
        e_key_valid = mL[ci][ri];
        e_key_high_bit = mhigh[ci][ri];
        e_action = maction[ci][ri];

        // put the key inside
        mL[ci][ri] = t_key;// put f_key to the selected bucket
        mhigh[ci][ri] = high_bit_t;
        maction[ci][ri] = t_action;

        // compute the second bucket for e_key
        //ht = (ci+e_key_high_bit*(mm)) ^ HashCode(num2str(e_key_valid));
        h1 = HashCode((parsedec2IPV4(e_key_valid)));
        h2 = HashCode((parsedec2IPV4(e_key_valid)+num2str(mm*mbc)));
        if(ci == h1)
            ht = h2;
        else
            ht = h1;

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
            maction[ht][It] = e_action;
            return true;
        }

        //fail to find, continue to kick out other items
        t_key = e_key_valid;
        t_action = e_action;
        //high_bit_t = e_key_high_bit;
        ci = ht;

    }

    return false;
    //             disp('insertion failure');

}

bool CuckooTable::AddKeyPrefix(const string& key,int prefix, int action)
{

    //define variables
    uint32_t f_key,t_key;
    long h1,h2,ci,ht;
    int I1,I2,It,high_bit_t,high_bit,e_key_valid,e_key_high_bit,e_action,t_action,t_prefix,e_prefix;
    int ri;
    vector<uint32_t> vec1,vec2;

    // compute fingerprint
    //f_key = Fingerprint(key);
    f_key = parseIPV4string(key.c_str());

    // if fingerprint is mm, print it
    //if (f_key == 0)
       // cout<<"fingerprint is mm"<<endl;

    // Compute the first hash code
    h1 = HashCode(key);

    // compute the second hash code
    //h2 = (h1 ^ HashCode(num2str(f_key)));
    h2 = HashCode((key+num2str(mm*mbc)));

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
        maction[h1][I1] = action;
        mprefix[h1][I1] = prefix;
        return true;
    }
    // if find vacant position in h2, store the item and return
    else if (I2 != mbc)//~isempty(I2)
    {
        mL[h2][I2] = f_key;
        mhigh[h2][I2] = high_bit;
        maction[h2][I2] = action;
        mprefix[h2][I2] = prefix;
        return true;
    }

    // fail to insert keys
    t_key = f_key;
    high_bit_t = high_bit;
    t_action = action;
    t_prefix = prefix;
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


    for (uint16_t ni = 1; ni < mMaxNumKicks; ni++)
    {
        // find a slot randomly
        ri = rand()%mbc;
        //kick out the key
        e_key_valid = mL[ci][ri];
        e_key_high_bit = mhigh[ci][ri];
        e_action = maction[ci][ri];
        e_prefix = mprefix[ci][ri];

        // put the key inside
        mL[ci][ri] = t_key;// put f_key to the selected bucket
        mhigh[ci][ri] = high_bit_t;
        maction[ci][ri] = t_action;
        mprefix[ci][ri] = t_prefix;

        // compute the second bucket for e_key
        //ht = (ci+e_key_high_bit*(mm)) ^ HashCode(num2str(e_key_valid));
        h1 = HashCode((parsedec2IPV4(e_key_valid)));
        h2 = HashCode((parsedec2IPV4(e_key_valid)+num2str(mm*mbc)));
        if(ci == h1)
            ht = h2;
        else
            ht = h1;

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
            maction[ht][It] = e_action;
            mprefix[ht][It] = e_prefix;
            return true;
        }

        //fail to find, continue to kick out other items
        t_key = e_key_valid;
        t_action = e_action;
        t_prefix = e_prefix;
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
bool CuckooTable::LookUpKey(const string& key, int &prefix)
{

    //define variables
    uint32_t f_key,t_key;
    long h1,h2;

    // compute fingerprint
    //f_key = Fingerprint(key);
    f_key = parseIPV4string(key.c_str());

    //if (f_key == 0)
        //cout<<"fingerprint is mm"<<endl;

    // Compute the first hash code
    h1 = HashCode(key);

    // compute the second hash code
    // h2 = (h1 ^ HashCode(num2str(f_key)));
    h2 = HashCode(key+num2str(mm*mbc));
    if (h2 >= mm)
    {
        h2 = (h2%mm);

    }


    // Search for key in list(hk)
    for (int i = 0; i<mbc; i++)
    {
        if ((mL[h1][i] == f_key && maction[h1][i] == prefix) | (mL[h2][i] == f_key&& maction[h2][i] == prefix))
        {
            return true;
        }

    }

    return false;
}

bool CuckooTable::LookUpKeyAction(const string& key, int &prefix, int& action)
{

    //define variables
    uint32_t f_key,t_key;
    long h1,h2;

    // compute fingerprint
    //f_key = Fingerprint(key);
    f_key = parseIPV4string(key.c_str());

    //if (f_key == 0)
        //cout<<"fingerprint is mm"<<endl;

    // Compute the first hash code
    h1 = HashCode(key);

    // compute the second hash code
    // h2 = (h1 ^ HashCode(num2str(f_key)));
    h2 = HashCode(key+num2str(mm*mbc));

    if (h2 >= mm)
    {
        h2 = (h2%mm);

    }



    // Search for key in list(hk)
    for (int i = 0; i<mbc; i++)
    {

        if ((mL[h1][i] == f_key && mprefix[h1][i] == prefix))
        {

            action = maction[h1][i];

            return true;
        }
        else if((mL[h2][i] == f_key&& mprefix[h2][i] == prefix))
        {

            action = maction[h2][i];

            return true;
        }
    }

    return false;
}

bool CuckooTable::LookUpKeyActionCount(const string& key, int &prefix, int& action, size_t& flowNo)
{

    //define variables
    uint32_t f_key,t_key;
    long h1,h2;

    // compute fingerprint
    //f_key = Fingerprint(key);
    f_key = parseIPV4string(key.c_str());

    //if (f_key == 0)
        //cout<<"fingerprint is mm"<<endl;

    // Compute the first hash code
    h1 = HashCode(key);

    // compute the second hash code
    // h2 = (h1 ^ HashCode(num2str(f_key)));
    h2 = HashCode(key+num2str(mm*mbc));
    if (h2 >= mm)
    {
        h2 = (h2%mm);

    }


    // Search for key in list(hk)
    for (int i = 0; i<mbc; i++)
    {
        if ((mL[h1][i] == f_key && mprefix[h1][i] == prefix))
        {
            action = maction[h1][i];
            mcount[h1][i] += flowNo;
            return true;
        }
        else if((mL[h2][i] == f_key&& mprefix[h2][i] == prefix))
        {
            action = maction[h2][i];
            mcount[h2][i] += flowNo;
            return true;
        }
    }

    return false;
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
int CuckooTable::LookUpKeyCol(const string& key,int & iaction)
{

    //define variables
    long f_key,t_key;
    long h1,h2;
    int returnValue = 0;
    bool flag;
    vector<int> action;

    // compute fingerprint
    //f_key = Fingerprint(key);
    f_key = parseIPV4string(key.c_str());

    //if (f_key == 0)
        //cout<<"fingerprint is mm"<<endl;

    // Compute the first hash code
    h1 = HashCode(key);

    // compute the second hash code
    //h2 = (h1 ^ HashCode(num2str(f_key)));
    h2 = HashCode(key+num2str(mm*mbc));

    if (h2 >= mm)
    {
        h2 = (h2%mm);

    }
    // Search for key in list(hk)
    for (int i = 0; i<mbc; i++)
    {
        if((mL[h1][i]) == f_key)
        {
            action.push_back(maction[h1][i]);
            break;
        }
        else if ((mL[h2][i]) == f_key)
        {

            action.push_back(maction[h2][i]);
            break;
        }

    }
    returnValue = action.size();

    if (action.size()>0)
    {
        iaction = action[0];
    return 1;
    }
    else
        return 0;
}
/**
remove a key from cuckoo filter
*/
bool CuckooTable::RemoveKey(const string& key)
{

    return true;
}
void CuckooTable::ClearTable()
{

    mL.clear();
    maction.clear();
    mhigh.clear();

    vector<vector<uint32_t> >().swap(mL);
    vector<vector<int> >().swap(maction);
    vector<vector<int> >().swap(mprefix);
    vector<vector<bool> >().swap(mhigh);
    vector<vector<size_t> >().swap(mcount);

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
long CuckooTable::HashCode(const string& key)
{
    //define variables
    int ll = 0;
    long hk = 0;
    vector<int> key_double;

    // Convert character array to integer array
    ll = key.size();

    key_double.reserve(ll);
    for(int i = 0; i < ll; i++)
    {
        key_double[i] = int(key[i])-45;

    }

    if(key.size()>0)
    hk = (key_double[0]);

    for (int i = 1; i<ll; i++)
    {
        // Could be implemented more efficiently in practice via bit
        // shifts (see reference)
        hk = (mc * hk + key_double[i])%mp;
    }
    if(mp!=0 && mm!= 0)
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
long CuckooTable::Fingerprint(const string& key)
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

    hk = ((maf * hk + mbf)%mpf)%long(pow(2,mf))+1;
    return hk;
}
/**
init hash funciton
*/
bool CuckooTable::InitHashFunction()
{
    //define variables
    long ff;
    long pp;
    // Set prime parameter
    ff = 1000; // fudge factor
    pp = ff * max(long(mm + 1),long(76));
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
    pp = ff * max(long(pow(2,mf) + 1),long(76));
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

void CuckooTable::clearCount()
{
    mcount = vector<vector<size_t> > (mm, vector<size_t>(mbc, 0));
}

void CuckooTable::returnKey(vector<string>& keys, vector<size_t>& keyNos,
                            vector<int>& keyPrefixes, vector<int>& keyActions)
{
    // Search for key in list(hk)
    for(size_t h = 0; h < mm; h++)
    {
        for (int i = 0; i<mbc; i++)
        {
            if((mL[h][i]) != 0)
            {
                keys.push_back(parsedec2IPV4(mL[h][i]));
                keyNos.push_back(mcount[h][i]);
                keyPrefixes.push_back(mprefix[h][i]);
                keyActions.push_back(maction[h][i]);
            }

        }
    }
}

float CuckooTable::occupyRate()
{
    // Search for key in list(hk)
    float num = 0.0f;
    for(size_t h = 0; h < mm; h++)
    {
        for (int i = 0; i<mbc; i++)
        {
            if((mL[h][i]) != 0)
            {
                num += 1.0f;
            }

        }
    }

    float ocyRate = num/float(mm*mbc);
    return ocyRate;
}

