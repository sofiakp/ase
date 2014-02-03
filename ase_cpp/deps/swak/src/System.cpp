#include "swak/System.h"

void SplitPath(const string &pathname, string &dir, string &file)
{
  size_t pos = pathname.find_last_of("/\\");
  if (pos == string::npos)
  {
    dir.clear();
    file = pathname;
  }
  else
  {
    file = pathname.substr(pos + 1, pathname.length() - (pos + 1));
    size_t dend = pathname.find_last_not_of("/\\", pos);
    if (dend == string::npos)
      dir = pathname.substr(0, pos + 1);
    else
      dir = pathname.substr(0, dend + 1);
  }
}

string PathDir(const string &pathname)
{
  size_t pos = pathname.find_last_of("/\\");
  if (pos == string::npos)
    return string();
  else
  {
    size_t dend = pathname.find_last_not_of("/\\", pos);
    if (dend == string::npos)
      return pathname.substr(0, pos + 1);
    else
      return pathname.substr(0, dend + 1);
  }
}

string PathFile(const string &pathname)
{
  size_t pos = pathname.find_last_of("/\\");
  if (pos == string::npos)
    return pathname;
  else
    return pathname.substr(pos + 1, pathname.length() - (pos + 1));
}

string JoinPath(const string &dir, const string &file)
{
  if (!file.empty() && ((file[0] == '/') || (file[0] == '\\')))
    return file;

  if (!dir.empty() && 
      (dir[dir.length() - 1] != '/') && (dir[dir.length() - 1] != '\\'))
    return dir + "/" + file;
  else
    return dir + file;
}

string Basename(const string &pathname)
{
  string dir, file;
  SplitPath(pathname, dir, file);
  return file;
}

bool GetDirList(const string &dir, vector<string> &files, vector<string> &subdirs)
{
  DIR *d = opendir(dir.c_str());
  if (!d)
    return false;
  
  string name = dir + '/';
  int base_dir_len = name.length();
  
  struct dirent *de;
  struct stat stat_buf;
  while ((de = readdir(d)))
  {
    if ((strcmp(de->d_name, ".") == 0) || (strcmp(de->d_name, "..") == 0))
      continue;
    name.resize(base_dir_len);
    name += de->d_name;
    if (!stat(name.c_str(), &stat_buf))
    {
      if (S_ISREG(stat_buf.st_mode))
        files.push_back(de->d_name);
      else if (S_ISDIR(stat_buf.st_mode))
        subdirs.push_back(de->d_name);
    }
  }
  
  closedir(d);
  return true;
}

bool MakeDir(const string &dir)
{
  return (mkdir(dir.c_str(), 0777) == 0);
}

bool RemoveDir(const string &dir)
{
  return (rmdir(dir.c_str()) == 0);
}

bool RemoveRecursiveDir(const string &dir)
{
  vector<string> files, subdirs;
  if (!GetDirList(dir, files, subdirs))
    return false;

  for (int i = 0; i < subdirs.size(); ++i)
    if (!RemoveRecursiveDir(JoinPath(dir, subdirs[i])))
      return false;

  for (int i = 0; i < files.size(); ++i)
    if (!DeleteFile(JoinPath(dir, files[i])))
      return false;

  return RemoveDir(dir);
}

bool FileReadable(const string &pathname)
{
  ifstream f(pathname.c_str(), ifstream::in);
  return f.is_open();
}

bool CreateFile(const string &pathname)
{
  ofstream f(pathname.c_str(), ofstream::out);
  bool ok = f.is_open();
  f.close();

  return ok;
}

bool DeleteFile(const string &pathname)
{
  return (remove(pathname.c_str()) == 0);
}

bool RenameFile(const string &from, const string &to)
{
  return (rename(from.c_str(), to.c_str()) == 0);
}

bool LinkFile(const string &from, const string &to)
{
  return (link(from.c_str(), to.c_str()) == 0);
}

bool ResizeFile(const string &pathname, uint64 size)
{
  if (!FileReadable(pathname))
    CreateFile(pathname);

  return (truncate(pathname.c_str(), size) == 0);
}

bool CopyFile(const string &pathname, const string &toname)
{
  ofstream os(toname.c_str(), ios::binary);
  if (!os.good())
    return false;

  return CopyFile(pathname, os);
}

bool CopyFile(const string &pathname, ostream &os)
{
  if (!os.good())
    return false;

  ifstream is(pathname.c_str(), ios::binary);
  if (!is.good())
    return false;

  is.seekg (0, ios::end);
  int64 length = is.tellg();
  is.seekg (0, ios::beg);

  const int64 megabyte = 1048576;
  char * buffer = new char[megabyte]; // 1mb
  while (length > 0)
  {
    int num_bytes = min(length, megabyte);
    is.read(buffer, num_bytes);
    os.write(buffer, num_bytes);
    length -= num_bytes;
  }

  delete [] buffer;

  return true;
}

struct PathNameParts
{
  int d_start, d_end; // Directory start/end
  int f_start, f_end; // File start/end
  int e_start, e_end; // Extension start/end
  int c_start, c_end; // Compressed extension start/end

  PathNameParts(const string &pathname)
  {
    d_start = d_end = f_start = 0;
    f_end = e_start = e_end = c_start = c_end = pathname.length();

    size_t spos = pathname.find_last_of("/\\");
    if (spos != string::npos)
    {
      f_start = spos + 1;
      size_t dend = pathname.find_last_not_of("/\\", spos);
      if (dend == string::npos)
        d_end = spos + 1;
      else
        d_end = dend + 1;
    }

    if (f_start < pathname.length())
    {
      size_t epos = pathname.find_last_of('.');
      if ((epos != string::npos) && (epos > f_start))
      {
        string temp = pathname.substr(epos + 1, pathname.length() - (epos + 1));
        for (int i = 0; i < temp.length(); ++i)
          temp[i] = tolower(temp[i]);
        if ((temp == "gz") || (temp == "bz") | (temp == "bz2"))
        {
          c_start = epos + 1;
          e_start = e_end = epos;

          size_t eepos = pathname.find_last_of('.', epos - 1);
          if ((eepos == string::npos) || (epos <= f_start))
            f_end = epos;
          else
          {
            f_end = eepos;
            e_start = eepos + 1;
          }
        }
        else
        {
          f_end = epos;
          e_start = epos + 1;
        }
      }
    }
  }
};

bool FileHasCompressedExt(const string &pathname)
{
  PathNameParts parts(pathname);
  return (parts.c_start < parts.c_end);
}

bool FileHasExt(const string &pathname)
{
  PathNameParts parts(pathname);
  return (parts.e_start < parts.e_end) || (parts.c_start < parts.c_end);
}

string FileExt(const string &pathname)
{
  PathNameParts parts(pathname);
  if (parts.c_start < parts.c_end)
    return pathname.substr(parts.c_start, parts.c_end - parts.c_start);
  else
    return pathname.substr(parts.e_start, parts.e_end - parts.e_start);
}

string FileBase(const string &pathname)
{
  PathNameParts parts(pathname);
  if (parts.c_start < parts.c_end)
    return pathname.substr(0, parts.e_end);
  else
    return pathname.substr(0, parts.f_end);
}

bool FileHasUncompressedExt(const string &pathname)
{
  PathNameParts parts(pathname);
  return (parts.e_start < parts.e_end);
}

string UncompressedFileExt(const string &pathname)
{
  PathNameParts parts(pathname);
  return pathname.substr(parts.e_start, parts.e_end - parts.e_start);
}

string UncompressedFileBase(const string &pathname)
{
  PathNameParts parts(pathname);
  return pathname.substr(0, parts.f_end);
}

istream *InFileStream(const string &pathname)
{
  ifstream *stream = new ifstream(pathname.c_str(), ifstream::binary);
  if (!stream->is_open())
  {
    delete stream;
    return NULL;
  }
  else
    return stream;
}

ostream *OutFileStream(const string &pathname)
{
  ofstream *stream = new ofstream(pathname.c_str(), ofstream::binary);
  if (!stream->is_open())
  {
    delete stream;
    return NULL;
  }
  else
    return stream;
}

string TempPathname()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  uint32 temp = (uint32(tv.tv_sec) * 1000 + int32(tv.tv_usec / 1000) + (uint32(getpid()) << 16)) ^ ((uint32(rand()) << 15) + uint32(rand()));

  string fname = "tmp_";
  for (int i = 7; i >= 0; --i)
  {
    int c = (temp >> (4 * i)) & 0xf;
    fname += (c < 10) ? ('0' + c) : ('a' + (c - 10));
  }

  const char *tmp_dir = getenv("TMP_DIR");
  if (!tmp_dir)
    tmp_dir = "/tmp";

  return JoinPath(tmp_dir, fname);
}

double TimeStamp()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return double(int64(tv.tv_sec) * 1000 + int64(tv.tv_usec / 1000)) * 0.001;
}

int64 TimeStampInt()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return int64(tv.tv_sec) * 1000 + int64(tv.tv_usec / 1000);
}

istream * OpenFileOrStdin(const string &filename, ifstream &file_stream)
{
  if (filename == "-")
    return &cin;
  else
  {
    file_stream.open(filename.c_str());
    AssertMsg(file_stream.good(), "Couldn't open file for reading: " + filename);
    return &file_stream;
  }
}

ostream * OpenFileOrStdout(const string &filename, ofstream &file_stream)
{
  if (filename == "-")
    return &cout;
  else
  {
    file_stream.open(filename.c_str());
    AssertMsg(file_stream.good(), "Couldn't open file for writing: " + filename);
    return &file_stream;
  }
}

istream * OpenFileOrStdin(const string &filename)
{
  if (filename == "-")
    return &cin;
  else
  {
    ifstream * stream_ptr = new ifstream();
    return OpenFileOrStdin(filename, *stream_ptr);
  }
}

ostream * OpenFileOrStdout(const string &filename)
{
  if (filename == "-")
    return &cout;
  else
  {
    ofstream * stream_ptr = new ofstream();
    return OpenFileOrStdout(filename, *stream_ptr);
  }
}

/////

istream *InFileStream(const string &pathname, ifstream &file_stream)
{
  file_stream.open(pathname.c_str());
  AssertMsg(file_stream.good(), "Couldn't open file for reading: " + pathname);
  return &file_stream;
}

ostream *OutFileStream(const string &pathname, ofstream &file_stream)
{
  file_stream.open(pathname.c_str());
  AssertMsg(file_stream.good(), "Couldn't open file for writing: " + pathname);
  return &file_stream;
}

/////

istream * OpenFile(const string &filename, ifstream &file_stream)
{
  return InFileStream(filename, file_stream);
}

ostream * OpenFile(const string &filename, ofstream &file_stream)
{
  return OutFileStream(filename, file_stream);
}
