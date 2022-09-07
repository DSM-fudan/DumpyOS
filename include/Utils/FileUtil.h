//
// Created by wzy on 2021/8/7.
//

#ifndef MULGIFT_FILEUTIL_H
#define MULGIFT_FILEUTIL_H
#include <string>
#include <vector>
using namespace std;

class FileUtil {

public:
    static bool checkFileExists(const char *name);

    static void Getfilepath(const char *path, const char *filename, char *filepath);

    static bool checkDirClean(const char *path);

    static long getFileSize(const char *fname);

    static void getFiles(const string& path, vector<string>& files );

    static float *readSeries(FILE *f);

    static void writeSeries(FILE *f, float *ts);

    static int FileRemove(const char *fname);

    static void mergeFiles(const string sources[], const string &dest, int num);

    static void deleteFiles(const string *fs, int num);

    static vector<float> *readSeriesVector(FILE *f);

    static void generateQueryFile(const string& data_file, int num);

    static vector<float> *readSeriesVector(FILE *f, int offsetInc);

    static float *readSeries(FILE *f, int num);

    static float **readSeries(FILE *f, const vector<int> &offsets, const vector<int> &cnts);

    static bool createDir(const char *path);

    static float *readSeries(FILE *f, const vector<int> &offsets);

    static float *readSeriesOffset(FILE *f, int offset);

    static void renameFile(const string &old_file, const string &new_file);
};


#endif //MULGIFT_FILEUTIL_H
