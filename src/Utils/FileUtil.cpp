//
// Created by wzy on 2021/8/7.
//

#include "../../include/Utils/FileUtil.h"
#include "../../include/DataStructures/TimeSeries.h"
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <unordered_set>
#include <cstdio>
#include <filesystem>
#include <dirent.h>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <sys/stat.h>

bool FileUtil::checkFileExists(const char *name) {
    ifstream f(name);
    return f.good();
}

int FileUtil::FileRemove(const char* fname)
{
    return remove(fname);
}

void FileUtil::deleteFiles(const string fs[], int num){
    for(int i=0;i<num;++i)
        FileRemove(fs[i].c_str());
}

long FileUtil::getFileSize(const char* fname)
{
    struct stat statbuf{};
    if(stat(fname,&statbuf)==0)
        return statbuf.st_size;
    return -1;
}

void FileUtil::Getfilepath(const char *path, const char *filename,  char *filepath)
{
    strcpy(filepath, path);
    if(filepath[strlen(path) - 1] != '/')
        strcat(filepath, "/");
    strcat(filepath, filename);
}

bool FileUtil::createDir(const char * path){
    return !mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);///home/newdir
}

bool FileUtil::checkDirClean(const char* path)
{
    if(!checkFileExists(path))
        createDir(path);
    DIR *dir;
    struct dirent *dirinfo;
    struct stat statbuf{};
    char filepath[256] = {0};
    lstat(path, &statbuf);

    if (S_ISREG(statbuf.st_mode))//判断是否是常规文件
        {
        remove(path);
        }
    else if (S_ISDIR(statbuf.st_mode))//判断是否是目录
        {
        if ((dir = opendir(path)) == nullptr)
            return true;
        while ((dirinfo = readdir(dir)) != nullptr)
        {
            Getfilepath(path, dirinfo->d_name, filepath);
            if (strcmp(dirinfo->d_name, ".") == 0 || strcmp(dirinfo->d_name, "..") == 0)//判断是否是特殊目录
                continue;
            checkDirClean(filepath);
            rmdir(filepath);
        }
        closedir(dir);
        }
    return false;
}

void FileUtil::mergeFiles(const string sources[], const string& dest, int num ){
    FILE * out = fopen(dest.c_str(), "wb");
    for(int i=0;i<num;++i){
        long size = getFileSize(sources[i].c_str());
        char * tmp = new char[size];
        FILE *in = fopen(sources[i].c_str(), "rb");
        fread(tmp, 1, size, in);
        fclose(in);
        fwrite(tmp, 1, size, out);
        free(tmp);
    }
    fclose(out);
}

void FileUtil::getFiles(const string& path, vector<string>& files )
{
    files.clear();
    for(auto &entry:filesystem::directory_iterator(path))
        if(entry.is_directory())    continue;
        else files.push_back(entry.path());
}

float * FileUtil::readSeries(FILE* f){
    auto *ts = new float[Const::tsLength];
    fread(ts, sizeof(float), Const::tsLength, f);
    return ts;
}

float * FileUtil::readSeries(FILE* f, int num){
    long size = Const::tsLength * (long)num;
    auto *ts = new float[size];
    fread(ts, sizeof(float), size, f);
    return ts;
}

float * FileUtil::readSeriesOffset(FILE* f, int offset){
    long off = Const::tsLengthBytes * (long)offset;
    auto *ts = new float[Const::tsLength];
    fseek(f, off, SEEK_SET);
    fread(ts, sizeof(float), Const::tsLength, f);
    return ts;
}

float** FileUtil::readSeries(FILE* f, const vector<int>& offsets, const vector<int>& cnts){
    int num = offsets.size();
    if(num <= 0)    return nullptr;
    auto**seriesSet = new float*[num];
    float *ts;  int offset; int cnt;
    fseek(f, (long)offsets[0] * Const::tsLengthBytes, SEEK_SET);
    for(int i=0;i<num;++i){
        offset = offsets[i];
        cnt = cnts[i];
        ts = new float[Const::tsLength * cnt];
        fread(ts, sizeof(float), Const::tsLength * cnt, f);
        seriesSet[i] = ts;
        if(i < num - 1)
            fseek(f, (long)(offsets[i + 1] - offset - cnt) * Const::tsLengthBytes, SEEK_CUR);
    }
    return seriesSet;
}

float* FileUtil::readSeries(FILE* f, const vector<int>& offsets){
    int num = offsets.size();
    if(num <= 0)    return nullptr;
    auto*seriesSet = new float[num * Const::tsLength];
    int before  = -1;
    fseek(f, (long)offsets[0] * Const::tsLengthBytes, SEEK_SET);
    for(int i=0;i<num;++i){
        fseek(f, (long )(offsets[i] - before - 1) * Const::tsLengthBytes, SEEK_CUR);
        before = offsets[i];
        fread(seriesSet + i * Const::tsLength, sizeof(float), Const::tsLength, f);
    }
    return seriesSet;
}

vector<float>* FileUtil::readSeriesVector(FILE* f){
    auto *ts = new vector<float>(Const::tsLength);
    fread(&((*ts)[0]), sizeof(float), Const::tsLength, f);
    return ts;
}

vector<float>* FileUtil::readSeriesVector(FILE* f, int offset){
    auto *ts = new vector<float>(Const::tsLength);
    fseek(f, offset, SEEK_SET);
    fread(&((*ts)[0]), sizeof(float), Const::tsLength, f);
    return ts;
}

void FileUtil::writeSeries(FILE* f, float *ts){
    fwrite(ts, sizeof(float), Const::tsLength, f);
}
#include <cmath>
void FileUtil::generateQueryFile(const string& data_file, int num){
    FILE * inf = fopen(data_file.c_str(), "rb");
    int lastIndex= data_file.rfind(".bin");
    string fn = data_file.substr(0, lastIndex) + "_query_in.bin";
    FILE * outf = fopen(fn.c_str(), "wb");
    unordered_set<int>mask;
    int n;
    srand(time(nullptr));
    for(int i=0;i<num;++i) {
        while (true){
           n = rand() % 100000000;
           if(mask.find(n) == mask.end()){
               mask.insert(n);
               break;
           }
        }
        cout << n << ",";
        float *ts;
        while (true){
            ts = readSeriesOffset(inf, n);
            if(!isnan(ts[0])) break;
        }
        cout << ts[0] << endl;
        fwrite(ts, sizeof(float), Const::tsLength, outf);
    }
    fclose(inf);
    fclose(outf);
}