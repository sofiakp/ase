#ifndef SWAK_H
#define SWAK_H

#ifndef Assert
#define Assert(cond) do {if (!(cond)) throw runtime_error("Assertion failed: " #cond); } while(false)
#endif

#ifndef AssertMsg
#define AssertMsg(cond, msg) do {if (!(cond)) throw runtime_error(string("Assertion failed: " #cond "\n") + msg); } while(false)
#endif

namespace std
{
};
using namespace std;

namespace __gnu_cxx
{
};
using namespace __gnu_cxx;

#define _FILE_OFFSET_BITS 64

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <vector>
#include <list>
#include <ext/slist>
#include <string>
#include <sstream>
#include <iostream>
#include <streambuf>
#include <fstream>
#include <iomanip>
#include <deque>
#include <algorithm>
#include <limits>
#include <set>
#include <map>
#include <queue>
#include <ext/functional>
#include <ext/hash_set>
#include <ext/hash_map>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include <dirent.h>
#include <pthread.h>
#include <regex.h>

typedef unsigned char uint8;
typedef signed char int8;
typedef unsigned short uint16;
typedef signed short int16;
typedef unsigned int uint32;
typedef signed int int32;
typedef unsigned long long uint64;
typedef signed long long int64;
typedef unsigned int uword;
typedef signed int word;


typedef map<string, string> MapSS;
typedef map<string, int64> MapSI;
typedef vector<string> VecS;
typedef vector<int32> VecI;
typedef vector<int64> VecL;
typedef vector<float> VecF;
typedef vector<double> VecD;

#define UWORD_BITS ((sizeof(uword) / sizeof(uint8)) * 8)
#define WORD_BITS ((sizeof(word) / sizeof(uint8)) * 8)
#define INT_BITS ((sizeof(int) / sizeof(uint8)) * 8)
#define CHAR_BITS ((sizeof(char) / sizeof(uint8)) * 8)

#define UINT8_MIN 0
#define UINT8_MAX 255
#define INT8_MIN -128
#define INT8_MAX 127
#define UINT16_MIN 0
#define UINT16_MAX 65535
#define INT16_MIN -32768
#define INT16_MAX 32767
#define UINT32_MIN uint32(0U)
#define UINT32_MAX uint32(4294967295U)
#define INT32_MIN int32(-2147483647 - 1)
#define INT32_MAX int32(2147483647)
#define UINT_MIN 0

#ifndef UINT_MAX
#define UINT_MAX 4294967295
#endif

#ifndef INT_MIN
#define INT_MIN int(-2147483647 - 1)
#endif

#ifndef INT_MAX
#define INT_MAX int(2147483647)
#endif

#define SIZE_T_MIN 0

#ifndef SIZE_T_MAX
#define SIZE_T_MAX 4294967295
#endif

#define INT64_CONST(x) (x##LL)

////////////////

// For some var, you want to cout << "var: " << var << " ";
#ifndef pprint
#define pprint(x)  #x << ": " << x << " "
#endif

#ifndef ppout
#define ppout(x) do {cout << #x << ": " << x << " " << endl;} while(false)
#endif

#ifndef pperr
#define pperr(x) do {cerr << #x << ": " << x << " " << endl;} while(false)
#endif

// Same, but with strings instead of <<
#ifndef sprint
#define sprint(x) (string(" " #x "=") + ToStr(x) + string(" "))
#endif

// STL loop helpers
#define StlFor(i, container) for (int64 i = 0; i < container.size(); ++i)
#define StlForMap(type1, type2, iter, container) for (map< type1, type2 >::iterator iter = container.begin(); iter != container.end(); ++iter)
#define StlForMapConst(type1, type2, iter, container) for (map< type1, type2 >::const_iterator iter = container.begin(); iter != container.end(); ++iter)
#define StlForIter(iter, cont_type, container) for (cont_type::iterator iter = container.begin(); iter != container.end(); ++iter)
#define StlForIterConst(iter, cont_type, container) for (cont_type::const_iterator iter = container.begin(); iter != container.end(); ++iter)

#include "StringUtil.h"

#endif
