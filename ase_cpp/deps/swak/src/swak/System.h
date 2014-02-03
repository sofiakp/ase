#ifndef SWAK_SYSTEM_H
#define SWAK_SYSTEM_H

#include "Swak.h"

void SplitPath(const string &pathname, string &dir, string &file);
string PathDir(const string &pathname);
string PathFile(const string &pathname);
string JoinPath(const string &dir, const string &file);
string Basename(const string &pathname);

bool GetDirList(const string &dir, vector<string> &files, vector<string> &subdirs);
bool MakeDir(const string &dir);
bool RemoveDir(const string &dir);
bool RemoveRecursiveDir(const string &dir);
 
bool FileReadable(const string &pathname);
bool CreateFile(const string &pathname);
bool DeleteFile(const string &pathname);
bool RenameFile(const string &from, const string &to);
bool LinkFile(const string &from, const string &to);
bool ResizeFile(const string &pathname, uint64 size);
bool CopyFile(const string &pathname, const string &toname);
bool CopyFile(const string &pathname, ostream &os);

bool FileHasCompressedExt(const string &pathname); // Does it end in .gz .bz2 .zip etc?
bool FileHasExt(const string &pathname);  // Does it have any extension at all? (of any kind, compressed or not)
string FileExt(const string &pathname);   // Gets the last extension on the file of any kind
string FileBase(const string &pathname);  // Gets everything up to the last extension of any kind
bool FileHasUncompressedExt(const string &pathname);     // Does it end in a non-zipped extension?
string UncompressedFileExt(const string &pathname);      // Gets txt out of blah.txt and blah.txt.gz
string UncompressedFileBase(const string &pathname);     // Gets blah out of blah.txt and blah.txt.gz

istream *InFileStream(const string &pathname);
ostream *OutFileStream(const string &pathname);
istream *InFileStream(const string &pathname, ifstream &file_stream);
ostream *OutFileStream(const string &pathname, ofstream &file_stream);

istream * OpenFileOrStdin(const string &filename, ifstream &file_stream);
ostream * OpenFileOrStdout(const string &filename, ofstream &file_stream);
istream * OpenFileOrStdin(const string &filename);
ostream * OpenFileOrStdout(const string &filename);

istream * OpenFile(const string &pathname, ifstream &file_stream);
ostream * OpenFile(const string &pathname, ofstream &file_stream);

string TempPathname();

double TimeStamp();
int64 TimeStampInt();

#endif
