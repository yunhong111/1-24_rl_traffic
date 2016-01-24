/*a total of 18 hash functions*/
#ifndef HASH_FUNCTION_H
#define HASH_FUNCTION_H

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

using namespace std;

unsigned int 
BOB(const unsigned char * str, unsigned int len);
/*
unsigned int
PJWHash (const unsigned char *str, unsigned int len);
*/
unsigned int
SHA1(const unsigned char *str, unsigned int len);

unsigned int 
BOB1(const unsigned char * str, unsigned int len);


#endif
