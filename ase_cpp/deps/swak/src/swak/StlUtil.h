#ifndef SWAK_STL_UTIL_H_
#define SWAK_STL_UTIL_H_

#include "System.h"

namespace Vector
{
  template <typename T>
  string ToStr(const vector<T> &vec, string separator = ", ", string delim_l = "[", string delim_r = "]")
  {
    ostringstream ss;
    ss << delim_l;
    if (!vec.empty())
    {
      ss << vec[0];
      for (int i = 1; i < vec.size(); ++i)
        ss << separator << vec[i];
    }
    ss << delim_r;
    return ss.str();
  }
};

namespace Map
{
  template <typename K, typename V> bool Contains(const map<K, V> &the_map, const K &key)
  {
    return the_map.find(key) != the_map.end();
  }
};

namespace __gnu_cxx
{
  template <typename T>
  ostream& operator<< (ostream &os, const vector<T> &vec)
  {
    os << Vector::ToStr(vec);
    return os;
  }
};


#endif
