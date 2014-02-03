#ifndef SWAK_SAFEVEC
#define SWAK_SAFEVEC

//#include "Swak.h"
#include <vector>
#include <stdexcept>
#include <sstream>

using namespace std;

template<typename T, class Allocator = allocator<T> >
class SafeVec : public vector<T>
{
public:
  explicit SafeVec() : vector<T, Allocator>() {};
  explicit SafeVec(size_t n, const T& value=T()) : vector<T, Allocator>(n, value) {} ;

  /*
  template <class InputIterator>
  SafeVector ( InputIterator first, InputIterator last, const Allocator& = Allocator() );
  SafeVector ( const vector<T,Allocator>& x );
  */
  
  T & operator[] ( size_t n )
  {
    if (n >= this->size()) 
    {
      stringstream msg;
      msg << "SafeVec index out of bounds in T & operator[] i=" << n << " size()=" << this->size();
      throw runtime_error(msg.str());
    }
    return vector<T>::operator[](n);
  }

  const T & operator[] ( size_t n ) const
  {
    if (n >= this->size()) 
    {
      stringstream msg;
      msg << "SafeVec index out of bounds in const T & operator[] i=" << n << " size()=" << this->size();
      throw runtime_error(msg.str());
    }
    return vector<T>::operator[](n);
  }
};


/*
typedef size_t size_type;
#define const_reference const T &
#define reference T &
#define iterator typename T::iterator
#define const_iterator typename T::const_iterator

template<typename T>
class SafeVec
{
  public:
    vector<T> vec;

    //class iterator : vector<T>::iterator { };
    //class const_iterator : vector<T>::const_iterator { };

    template <class InputIterator>
    void assign ( InputIterator first, InputIterator last ) {vec.assign(first, last); }
    void assign ( size_type n, const T& u ) { vec.assign(n, u); }

    const_reference at ( size_type n ) const { return vec.at(n); }
          reference at ( size_type n ) { return vec.at(n); }
  
    reference back () { return vec.back(); }
    const_reference back () const { return vec.back(); }

    iterator begin () { return vec.begin(); }
    const_iterator begin () const { return vec.begin(); }

    size_type capacity () const { return vec.capacity(); }

    void clear ( ) { return vec.clear(); }

    bool empty () const { return vec.empty(); } 

    iterator end () { return vec.end(); }
    const_iterator end () const {return vec.end(); }

    iterator erase ( iterator position ) { return vec.erase(position); }
    iterator erase ( iterator first, iterator last ) { return vec.erase(first, last); }
};

*/

#endif
