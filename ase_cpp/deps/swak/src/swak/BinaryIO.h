#ifndef SWAK_BIN_IO_H
#define SWAK_BIN_IO_H

#include "swak/Swak.h"

namespace Swak
{
  class BinWriter
  {
    private:
      ofstream fout;

    public:
      BinWriter()
      {
      }

      // For ints, floats, chars, etc
      template <typename T>
      void Write(const T &data)
      {
        fout.write(reinterpret_cast<const char*>(&data), sizeof(T));
        AssertMsg(fout.good(), "Error writing to binary file");
      }

      void Write(const string &str)
      {
        size_t size = str.size();
        fout.write(reinterpret_cast<const char*>(&size), sizeof(size));
        fout.write(str.c_str(), size);
        AssertMsg(fout.good(), "Error writing to binary file");
      }

      void Init(const string &file, const string &magic)
      {
        fout.open(file.c_str(), ios::out | ios::binary);
        AssertMsg(fout.good(), "Couldn't open " + file + " for binary writing");
        Write(magic);
      }
      
      void Close()
      {
        fout.close();
      }

      ~BinWriter()
      {
        Close();
      }
  };
  
  class BinReader
  {
    private:
      ifstream fin;
      char * buf;
      static const size_t default_buffer_size = 1024;

    public:
      BinReader()
      {
        buf = (char*)malloc(default_buffer_size);
      }

      // For ints, floats, chars, etc
      template <typename T>
      void Read(T &data)
      {
        fin.read(reinterpret_cast<char*>(&data), sizeof(T));
        AssertMsg(fin.good(), "Error reading from binary file");
      }

      void Read(string &str)
      {
        size_t size;
        fin.read(reinterpret_cast<char*>(&size), sizeof(size));

        fin.read(buf, size);
        str.assign(buf, size);

        AssertMsg(fin.good(), "Error reading from binary file");
      }

      void Init(const string &file, const string &expected_magic, size_t max_string_size=1024)
      {
        if (max_string_size > default_buffer_size)
          buf = (char*)realloc(buf, max_string_size);

        AssertMsg(buf != NULL, "Could not allocate buffer for binary reader");

        fin.open(file.c_str(), ios::in | ios::binary);
        AssertMsg(fin.good(), "Couldn't open " + file + " for binary reading");

        string magic;
        Read(magic);
        AssertMsg(magic == expected_magic, "When reading binary file " + file + " expected magic " + expected_magic + " in file, but found " + magic);
      }

      void Close()
      {
        fin.close();
      }

      ~BinReader()
      {
        Close();
        free(buf);
      }
  };

};

#endif
