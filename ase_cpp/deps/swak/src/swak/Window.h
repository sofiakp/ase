#ifndef SWAK_WINDOW
#define SWAK_WINDOW

#include "Swak.h"

namespace Swak
{
  template <typename T, typename Container = vector<T>, bool fast=false >
    class Window
    {
      private:
        Container * items;
        int start, end;

      public:
        Window()
        {
        }

        Window(Container *items_, int start_, int end_) 
        {
          AssertMsg(items_ != NULL, "In Window::Window()");
          Init(items_, start_, end_);
        }

        void Init(Container *items_, int start_, int end_) 
        {
          AssertMsg(items_ != NULL, "In Window::Init()");
          items = (items_); start = (start_); end = (end_);
          AssertMsg(start >= 0 && end <= items->size(), "In Window::Init()");
        }

        int size() const
        {
          return end - start;
        }

        int Start() const
        {
          return start;
        }

        int End() const
        {
          return end;
        }

        T & operator[] ( size_t n )
        {
          if (!fast && start + n >= items->size())
          {
            stringstream msg;
            msg << "Window index out of bounds in T & operator[] start=" << start << " i=" << n << " size()=" << this->size();
            throw runtime_error(msg.str());
          }

          return (*items)[start + n];
        }

        const T & operator[] ( size_t n ) const
        {
          if (!fast && start + n >= items->size())
          {
            stringstream msg;
            msg << "Window index out of bounds in const T & operator[] start=" << start << " i=" << n << " size()=" << this->size();
            throw runtime_error(msg.str());
          }
          return (*items)[start + n];
        }

        T & Last()
        {
          return operator[](size() - 1);
        }

        const T & Last() const
        {
          return operator[](size() - 1);
        }

    };

};

#endif
