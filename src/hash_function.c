#include "hash_function.h"
#include "sha1.h"


#define mix(a,b,c) \
{ \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8); \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12);  \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}




unsigned int 
BOB(const unsigned char * str, unsigned int len)
{
	//register ub4 a,b,c,len;
	unsigned int a,b,c;
	unsigned int initval = 0;
	/* Set up the internal state */
	//len = length;
	a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
	c = initval;         /* the previous hash value */

	/*---------------------------------------- handle most of the key */
	while (len >= 12)
	{
		a += (str[0] +((unsigned int)str[1]<<8) +((unsigned int)str[2]<<16) +((unsigned int)str[3]<<24));
		b += (str[4] +((unsigned int)str[5]<<8) +((unsigned int)str[6]<<16) +((unsigned int)str[7]<<24));
		c += (str[8] +((unsigned int)str[9]<<8) +((unsigned int)str[10]<<16)+((unsigned int)str[11]<<24));
		mix(a,b,c);
		str += 12; len -= 12;
	}

	/*------------------------------------- handle the last 11 bytes */
	c += len;
	switch(len)              /* all the case statements fall through */
	{
		case 11: c+=((unsigned int)str[10]<<24);
		case 10: c+=((unsigned int)str[9]<<16);
		case 9 : c+=((unsigned int)str[8]<<8);
		/* the first byte of c is reserved for the length */
		case 8 : b+=((unsigned int)str[7]<<24);
		case 7 : b+=((unsigned int)str[6]<<16);
		case 6 : b+=((unsigned int)str[5]<<8);
		case 5 : b+=str[4];
		case 4 : a+=((unsigned int)str[3]<<24);
		case 3 : a+=((unsigned int)str[2]<<16);
		case 2 : a+=((unsigned int)str[1]<<8);
		case 1 : a+=str[0];
		/* case 0: nothing left to add */
	}
	mix(a,b,c);
	/*-------------------------------------------- report the result */
	return c;
}


/*
unsigned int
PJWHash (const unsigned char *str, unsigned int len)
{
    unsigned int BitsInUnignedInt = (unsigned int) (sizeof (unsigned
							    int) * 8);
    unsigned int ThreeQuarters = (unsigned int) ((BitsInUnignedInt * 3) / 4);
    unsigned int OneEighth = (unsigned int) (BitsInUnignedInt / 8);
    unsigned int HighBits = (unsigned int) (0xFFFFFFFF) <<
	(BitsInUnignedInt - OneEighth);
    unsigned int hash = 0;
    unsigned int test = 0;
   unsigned int i = 0;
	 for( i = 0; i < len; i++)
    {
		hash = (hash << OneEighth) + (*str++);
		if ((test = hash & HighBits) != 0)
		{
			hash = ((hash ^ (test >> ThreeQuarters)) & (~HighBits));
		}
	}
    return hash ;
}
*/



unsigned int 
BOB1(const unsigned char * str, unsigned int len)
{
	//register ub4 a,b,c,len;
	unsigned int a,b,c;
	unsigned int initval = 2;
	/* Set up the internal state */
	//len = length;
	a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
	c = initval;         /* the previous hash value */

	/*---------------------------------------- handle most of the key */
	while (len >= 12)
	{
		a += (str[0] +((unsigned int)str[1]<<8) +((unsigned int)str[2]<<16) +((unsigned int)str[3]<<24));
		b += (str[4] +((unsigned int)str[5]<<8) +((unsigned int)str[6]<<16) +((unsigned int)str[7]<<24));
		c += (str[8] +((unsigned int)str[9]<<8) +((unsigned int)str[10]<<16)+((unsigned int)str[11]<<24));
		mix(a,b,c);
		str += 12; len -= 12;
	}

	/*------------------------------------- handle the last 11 bytes */
	c += len;
	switch(len)              /* all the case statements fall through */
	{
		case 11: c+=((unsigned int)str[10]<<24);
		case 10: c+=((unsigned int)str[9]<<16);
		case 9 : c+=((unsigned int)str[8]<<8);
		/* the first byte of c is reserved for the length */
		case 8 : b+=((unsigned int)str[7]<<24);
		case 7 : b+=((unsigned int)str[6]<<16);
		case 6 : b+=((unsigned int)str[5]<<8);
		case 5 : b+=str[4];
		case 4 : a+=((unsigned int)str[3]<<24);
		case 3 : a+=((unsigned int)str[2]<<16);
		case 2 : a+=((unsigned int)str[1]<<8);
		case 1 : a+=str[0];
		/* case 0: nothing left to add */
	}
	mix(a,b,c);
	/*-------------------------------------------- report the result */
	return c;
}

unsigned int
SHA1(const unsigned char *str, unsigned int len)
{
	SHA1Context sha1Text;
	unsigned int *p;
	int i ;
	unsigned int hash;
	unsigned char messageDegest[20];
	SHA1Reset(&sha1Text);
	SHA1Input(&sha1Text, str, len);
	SHA1Result(&sha1Text, messageDegest);
	p = (unsigned int *)messageDegest;
	hash = p[0];
	for(i = 1; i < 5; i++){
		hash ^= p[i];
	}
	
	return hash;
}

