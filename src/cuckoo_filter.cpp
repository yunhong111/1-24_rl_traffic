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
CuckooFilter::CuckooFilter()
{
}
void CuckooFilter::cuckooFilterInit(long m,int f,int bc,long MaxNumKicks)
{
    // Set bucket count
    mm = m;

    // Start with empty hash table
    //mL0 = vector<vector<string> > (m, vector<string>(bc, '\0'));
    if(mL.size() == 0)
    {
        mL = vector<vector<int> > (m, vector<int>(bc, 0));
        maction = vector<vector<int> > (m, vector<int>(bc, 0));
        mhigh = vector<vector<bool> > (m, vector<bool>(bc, 0));
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
            }
        }

    }

    //mcount = vector<vector<size_t> > (m, vector<size_t>(bc, 0));
   // mcount0 = vector<vector<size_t> > (m, vector<size_t>(bc, 0));
   // mcountdiff = vector<vector<size_t> > (m, vector<size_t>(bc, 0));

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
bool CuckooFilter::AddKey(const string& key,int action)
{

    //define variables
    long f_key,t_key;
    long h1,h2,ci,ht;
    int I1,I2,It,high_bit_t,high_bit,e_key_valid,e_key_high_bit,e_action,t_action;
    string e_key_valid0,t_key0;
    int ri;
    vector<int> vec1,vec2;

    // compute fingerprint
    f_key = Fingerprint(key);
    //cout<<key<<" "<<"fkey: "<<f_key<<" ";

    // if fingerprint is mm, print it
    //if (f_key == 0)
        //cout<<"fingerprint is mm"<<endl;

    // Compute the first hash code
    h1 = HashCode(key);
    //cout<<"h1: "<<h1<<" ";

    // compute the second hash code
    h2 = (h1 ^ HashCodeXor(num2str(f_key)));

    high_bit = 0;
    if (h2 >= mm)
    {
        h2 = (h2%mm);
        high_bit = 1;

    }
    //cout<<"h2: "<<h2<<" ";

    // find the first vacant position in table
    vec1 = mL[h1];
    vec2 = mL[h2];
    I1 = find(vec1.begin(),vec1.end(),0)-vec1.begin();//find(~L(:,h1),1); // the first 0 in bucket h1
    I2 = find(vec2.begin(),vec2.end(),0)-vec2.begin();//find(~L(:,h2),1);

    // if find vacant position in h1, store the item and return
    if (I1 != mbc)//~isempty(I1) ----------------------------------------------------
    {
        mL[h1][I1] = f_key;
        //mL0[h1][I1] = key.c_str();
        mhigh[h1][I1] = 0;
        maction[h1][I1] = action;
        return true;
    }
    // if find vacant position in h2, store the item and return
    else if (I2 != mbc)//~isempty(I2)
    {
        mL[h2][I2] = f_key;
        //mL0[h2][I2] = key.c_str();
        mhigh[h2][I2] = high_bit;
        maction[h2][I2] = action;
        return true;
    }

    // fail to insert keys
    t_key = f_key;
    t_key0 = key.c_str();
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
        //e_key_valid0 = mL0[ci][ri];
        e_key_high_bit = mhigh[ci][ri];
        e_action = maction[ci][ri];

        // put the key inside
        mL[ci][ri] = t_key;// put f_key to the selected bucket
        //mL0[ci][ri] = t_key0.c_str();
        mhigh[ci][ri] = high_bit_t;
        maction[ci][ri] = t_action;

        // compute the second bucket for e_key
        ht = (ci+e_key_high_bit*(mm)) ^ HashCodeXor(num2str(e_key_valid));
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
           // mL0[ht][It] = e_key_valid0.c_str();
            mhigh[ht][It] = high_bit_t;
            maction[ht][It] = e_action;
            return true;
        }

        //fail to find, continue to kick out other items
        t_key = e_key_valid;
        t_key0 = e_key_valid0;
        t_action = e_action;
        //high_bit_t = e_key_high_bit;
        ci = ht;

    }
    //cout<<"key: "<<key<<" f_key: "<<f_key<<" h1: "<<h1<<" h2: "<<h2<<" ht: "<<ht<<endl;
    for(int i = 0; i < mbc; i++)
    {
        cout<<"h1: "<<mL[h1][i]<<" ";
         cout<<"h2: "<<mL[h2][i]<<" ";
    }
    return false;
    //             disp('insertion failure');

}

bool CuckooFilter::AddKeyCount(const string& key,int action, char mL0[][ACTIONSIZE][20])
{

    //define variables
    long f_key,t_key;
    long h1,h2,ci,ht;
    int I1,I2,It,high_bit_t,high_bit,e_key_valid,e_key_high_bit,e_action,t_action;
    char e_key_valid0[20],t_key0[20];
    int ri;
    vector<int> vec1,vec2;

    // compute fingerprint
    f_key = Fingerprint(key);
    //cout<<key<<" "<<"fkey: "<<f_key<<" ";

    // if fingerprint is mm, print it
    //if (f_key == 0)
        //cout<<"fingerprint is mm"<<endl;

    // Compute the first hash code
    h1 = HashCode(key);
    //cout<<"h1: "<<h1<<" ";

    // compute the second hash code
    h2 = (h1 ^ HashCodeXor(num2str(f_key)));

    high_bit = 0;
    if (h2 >= mm)
    {
        h2 = (h2%mm);
        high_bit = 1;

    }
    //cout<<"h2: "<<h2<<" ";

    // find the first vacant position in table
    vec1 = mL[h1];
    vec2 = mL[h2];
    I1 = find(vec1.begin(),vec1.end(),0)-vec1.begin();//find(~L(:,h1),1); // the first 0 in bucket h1
    I2 = find(vec2.begin(),vec2.end(),0)-vec2.begin();//find(~L(:,h2),1);

    // if find vacant position in h1, store the item and return
    if (I1 != mbc)//~isempty(I1) ----------------------------------------------------
    {
        mL[h1][I1] = f_key;
        bzero(&mL0[h1][I1], sizeof(mL0[h1][I1]));
        strcpy(mL0[h1][I1], key.c_str());
        mhigh[h1][I1] = 0;
        maction[h1][I1] = action;
        //cout<<"1: "<<key<<" "<<mL0[h1][I1]<<"        ";
        return true;
    }
    // if find vacant position in h2, store the item and return
    else if (I2 != mbc)//~isempty(I2)
    {
        mL[h2][I2] = f_key;
        bzero(&mL0[h2][I2], sizeof(mL0[h2][I2]));
        strcpy(mL0[h2][I2], key.c_str());
        //mL0[h2][I2] = key.c_str();
        mhigh[h2][I2] = high_bit;
        maction[h2][I2] = action;
        //if(strcmp(key.c_str(),mL0[h2][I2] ) == 0)
        //cout<<"2: "<<key<<" "<<mL0[h2][I2] <<"        ";
        return true;
    }

    // fail to insert keys
    t_key = f_key;
    bzero(&t_key0,sizeof(t_key0));
    strcpy(t_key0, key.c_str());
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
        bzero(&e_key_valid0,sizeof(e_key_valid0));
        strcpy(e_key_valid0, mL0[ci][ri]);
        e_key_high_bit = mhigh[ci][ri];
        e_action = maction[ci][ri];
        //cout<<"2: "<<e_key_valid0<<" "<<mL0[ci][ri]<<"        ";

        // put the key inside
        mL[ci][ri] = t_key;// put f_key to the selected bucket
        bzero(&mL0[ci][ri], sizeof(mL0[ci][ri]));
        strcpy(mL0[ci][ri], t_key0);
        //mL0[ci][ri] = t_key0.c_str();
        mhigh[ci][ri] = high_bit_t;
        maction[ci][ri] = t_action;

        //cout<<"3: "<<t_key0<<" "<<mL0[ci][ri]<<"        ";

        // compute the second bucket for e_key
        ht = (ci+e_key_high_bit*(mm)) ^ HashCodeXor(num2str(e_key_valid));
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
            bzero(&mL0[ht][It], sizeof(mL0[ht][It]));
            strcpy(mL0[ht][It], e_key_valid0);
            //mL0[ht][It] = e_key_valid0.c_str();
            mhigh[ht][It] = high_bit_t;
            maction[ht][It] = e_action;
            //cout<<"ke: "<<e_key_valid0<<" "<<mL0[ht][It] <<"        ";
            //cout<<"k: "<<key<<" "<<mL0[ht][It] <<"        ";
            return true;
        }

        //fail to find, continue to kick out other items
        t_key = e_key_valid;
        bzero(&t_key0,sizeof(t_key0));
        strcpy(t_key0, e_key_valid0);
        //t_key0 = e_key_valid0;
        t_action = e_action;
        //high_bit_t = e_key_high_bit;
        ci = ht;
        //cout<<"4: "<<t_key0<<" "<<e_key_valid0<<"        ";

    }
    //cout<<"key: "<<key<<" f_key: "<<f_key<<" h1: "<<h1<<" h2: "<<h2<<" ht: "<<ht<<endl;
    for(int i = 0; i < mbc; i++)
    {
        cout<<"h1: "<<mL[h1][i]<<" ";
         cout<<"h2: "<<mL[h2][i]<<" ";
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
bool CuckooFilter::LookUpKey(const string& key)
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
    h2 = (h1 ^ HashCodeXor(num2str(f_key)));
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
    /*cout<<"key: "<<key<<" size: "<<key.length()<<" f_key: "<<f_key<<" h1: "<<h1<<" h2: "<<h2<<endl;
    for(int i = 0; i < mbc; i++)
    {
        cout<<"h1: "<<mL[h1][i]<<" ";
         cout<<"h2: "<<mL[h2][i]<<" ";
    }*/
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
int CuckooFilter::LookUpKeyCol(const string& key,int & action_out)
{

    //define variables
    long f_key,t_key;
    long h1,h2;
    int returnValue = 0;
    vector<int> action;
    bool flag = 0;

    // compute fingerprint
    f_key = Fingerprint(key);

    if (f_key == 0)
        cout<<"fingerprint is mm"<<endl;

    // Compute the first hash code
    h1 = HashCode(key);

    // compute the second hash code
    h2 = (h1 ^ HashCodeXor(num2str(f_key)));
    if (h2 >= mm)
    {
        h2 = (h2%mm);

    }
    // Search for key in list(hk)
    for (int i = 0; i<mbc; i++)
    {
        flag = 0;

        if((mL[h1][i]) == f_key)
        {
            for (int ai = 0; ai < action.size(); ai++)
            {
                if (maction[h1][i] == action[ai]) {flag = 1; break;}
            }
            if (flag == 0)
            action.push_back(maction[h1][i]);
        }
        flag = 0;
        if ((mL[h2][i]) == f_key)
        {
            for (int ai = 0; ai < action.size(); ai++)
            {
                if (maction[h2][i] == action[ai]) {flag = 1; break;}
            }
            if (flag == 0)
            action.push_back(maction[h2][i]);
        }

    }
    returnValue = action.size();
    if (returnValue > 0)
    action_out = action[0];
    return returnValue;
}

int CuckooFilter::LookUpKeyActions(const string& key,vector<int>& action)
{

    //define variables
    long f_key,t_key;
    long h1,h2;
    int returnValue = 0;
    bool flag = 0;

    // compute fingerprint
    f_key = Fingerprint(key);

    if (f_key == 0)
        cout<<"fingerprint is mm"<<endl;

    // Compute the first hash code
    h1 = HashCode(key);

    // compute the second hash code
    h2 = (h1 ^ HashCodeXor(num2str(f_key)));

    if (h2 >= mm)
    {
        h2 = (h2%mm);

    }
    // Search for key in list(hk)
    for (int i = 0; i<mbc; i++)
    {
        flag = 0;

        if((mL[h1][i]) == f_key)
        {
            for (int ai = 0; ai < action.size(); ai++)
            {
                if (maction[h1][i] == action[ai]) {flag = 1; break;}
            }
            if (flag == 0)
            action.push_back(maction[h1][i]);

        }
        flag = 0;
        if ((mL[h2][i]) == f_key)
        {
            for (int ai = 0; ai < action.size(); ai++)
            {
                if (maction[h2][i] == action[ai]) {flag = 1; break;}
            }
            if (flag == 0)
            action.push_back(maction[h2][i]);

        }

    }
    returnValue = action.size();
    return returnValue;
}

int CuckooFilter::LookUpKeyActionsCount(const string& key,vector<int>& action, size_t& flowNo
                              ,char mL0[][ACTIONSIZE][20], vector<vector<size_t> >& mcount, vector<vector<size_t> >& mcount0,
                            vector<vector<size_t> >& mcountdiff, long& hBuck, int& slot_i)
{

    //define variables
    long f_key,t_key;
    long h1,h2;
    int returnValue = 0;
    bool flag = 0;

    // compute fingerprint
    f_key = Fingerprint(key);

    if (f_key == 0)
        cout<<"fingerprint is mm"<<endl;

    // Compute the first hash code
    h1 = HashCode(key);

    // compute the second hash code
    h2 = (h1 ^ HashCodeXor(num2str(f_key)));

    if (h2 >= mm)
    {
        h2 = (h2%mm);

    }
    // Search for key in list(hk)
    for (int i = 0; i<mbc; i++)
    {
        flag = 0;

        if((mL[h1][i]) == f_key )
        {
            for (int ai = 0; ai < action.size(); ai++)
            {
                if (maction[h1][i] == action[ai]) {flag = 1; break;}
            }
            if (flag == 0)
            action.push_back(maction[h1][i]);

            //mcount[h1][i] += flowNo;
            //cout<<mL0[h1][i]<<" "<<key.c_str()<<endl;

        }
        flag = 0;
        if ((mL[h2][i]) == f_key )
        {
            for (int ai = 0; ai < action.size(); ai++)
            {
                if (maction[h2][i] == action[ai]) {flag = 1; break;}
            }
            if (flag == 0)
            action.push_back(maction[h2][i]);

            //mcount[h2][i] += flowNo;

        }

    }

    // Search for key in list(hk)
    for (int i = 0; i<mbc; i++)
    {

        if((strcmp(mL0[h1][i] ,key.c_str())==0))
        {

            //mcount0[h1][i] += flowNo;
            //mcountdiff[h1][i] = mcount[h1][i] - mcount0[h1][i];
            //cout<<"match key"<<endl;
            //cout<<mL0[h1][i]<<" "<<key.c_str()<<endl;
            hBuck = h1;
            slot_i = i;
            break;

        }

        if ((strcmp(mL0[h2][i],key.c_str())==0))
        {
            //mcount0[h2][i] += flowNo;
            //mcountdiff[h1][i] = mcount[h2][i] - mcount0[h2][i];
            //cout<<"match key"<<endl;
            //cout<<mL0[h2][i]<<" "<<key.c_str()<<endl;
            hBuck = h2;
            slot_i = i;
            break;

        }

    }

    returnValue = action.size();
    /*if(!returnValue && hBuck>=0)
    {
        cout<<"wrong!"<<endl;
    }*/

    return returnValue;
}

int CuckooFilter::LookUpKeyBack(const string& key, long& posx, int posy)
{

    //define variables
    long f_key,t_key;
    long h1,h2;
    int returnValue = 0;

    // compute fingerprint
    f_key = Fingerprint(key);

    if (f_key == 0)
        cout<<"fingerprint is mm"<<endl;

    // Compute the first hash code
    h1 = HashCode(key);

    // compute the second hash code
    h2 = (h1 ^ HashCodeXor(num2str(f_key)));
    if (h2 >= mm)
    {
        h2 = (h2%mm);

    }
    // Search for key in list(hk)
    for (int i = 0; i<mbc; i++)
    {

        if((mL[h1][i]) == f_key)
        {
            posx = h1;
            posy = i;
            return 1;
        }

        if ((mL[h2][i]) == f_key)
        {
            posx = h2;
            posy = i;
            return 1;
        }

    }
    return 0;
}
/**
remove a key from cuckoo filter
*/
bool CuckooFilter::RemoveKey(const string& key)
{
    //define variables
    long f_key,t_key;
    long h1,h2;
    int returnValue = 0;

    // compute fingerprint
    f_key = Fingerprint(key);

    if (f_key == 0)
        cout<<"fingerprint is mm"<<endl;

    // Compute the first hash code
    h1 = HashCode(key);

    // compute the second hash code
    h2 = (h1 ^ HashCodeXor(num2str(f_key)));
    if (h2 >= mm)
    {
        h2 = (h2%mm);

    }
    // Search for key in list(hk)
    for (int i = 0; i<mbc; i++)
    {

        if((mL[h1][i]) == f_key)
        {
            mL[h1][i] = 0;
            maction[h1][i] = 0;
            return true;
        }

        if ((mL[h2][i]) == f_key)
        {
            mL[h2][i] = 0;
            maction[h2][i] = 0;
            return true;
        }

    }
    return false;
}

/**
remove a key from cuckoo filter
*/
bool CuckooFilter::RemovePos(long& xPos, int& yPos)
{
    mL[xPos][yPos] = 0;
    maction[xPos][yPos] = 0;
    return true;

}
void CuckooFilter::ClearTable()
{
    //mL.clear();
    //maction.clear();
    //mhigh.clear();

    //vector<vector<int> >().swap(mL);
    //ector<vector<int> >().swap(maction);
    //vector<vector<bool> >().swap(mhigh);

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
/*long CuckooFilter::HashCode(const string& key)
{
    //define variables
    int ll = 0;
    long hk= 0;
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

    //if(mp!=0 && mm!= 0)
    hk = ((ma * hk + mb)%mp)%mm;

    return hk;


}*/

long CuckooFilter::HashCode(const string& key)
{
    //define variables
    long hk;
    hk = SHA1((const unsigned char *)key.c_str(), key.length());
    hk = hk % mm;
    //cout<<"hk: "<<hk<<" "<<endl;

    return hk;
}

long CuckooFilter::HashCodeXor(const string& key)
{
    //define variables
    long hk;
    hk = SHA1((const unsigned char *)key.c_str(), key.length());
    hk = hk % mm;
    //cout<<"hk: "<<hk<<" "<<endl;

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
/*long CuckooFilter::Fingerprint(const string& key)
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
}*/

long CuckooFilter::Fingerprint(const string& key)
{
    //define variables
    long hk;
    hk = BOB((const unsigned char *)key.c_str(), key.length());
    long mask = long(pow(2,mf));
    hk = hk % mask +1;
    //cout<<key<<" "<<hk<<endl;
    //cout<<"mf: "<<mf<<" "<<endl;

    return hk;
}
/**
init hash funciton
*/
bool CuckooFilter::InitHashFunction()
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

void CuckooFilter::returnKey(vector<string>& keys,vector<int>& keyActions,
                              char mL0[][ACTIONSIZE][20], vector<vector<size_t> >& mcount, vector<vector<size_t> >& mcount0,
                            vector<vector<size_t> >& mcountdiff,vector<size_t>& keyCounts, vector<size_t>& keyCounts0,
                            vector<size_t>& keyCountDiffs)
{
    // Search for key in list(hk)
    for(size_t h = 0; h < mm; h++)
    {
        for (int i = 0; i<mbc; i++)
        {
            if((mL[h][i]) != 0)
            {
                keys.push_back((mL0[h][i]));
                keyActions.push_back(maction[h][i]);
                keyCounts.push_back(mcount[h][i]);
                keyCounts0.push_back(mcount0[h][i]);
                keyCountDiffs.push_back(mcountdiff[h][i]);
            }

        }
    }
}

