//
// Created by caucher on 2021/12/3.
//

#include "../../include/Expr/DNATranslator.h"
#include "../../include/Utils/TimeSeriesUtil.h"
#include <iostream>
#include <cstdio>
#include <fstream>
#include <cmath>

void DNATranslator::translate(const string &output, int len){
    string fn = "/mnt/g/Series4Similarity_Search/GCA_018340775.1_ASM1834077v1_genomic.fna.gz";
    ifstream fin(fn, ios::in);
    FILE *of = fopen(output.c_str(), "wb");
    string s;
    getline(fin, s);
    int step = 30;
    float cur = 0; char c;
    long num = 1;
    float data[len];
    for(float &i:data)    i=0;

    // read the 1st series
    for(int i=0;i<len && fin;){
        fin >> c;
        switch (tolower(c)) {
            case 'a':   cur += 1; break;
            case 'g':   cur += 2; break;
            case 'c':   cur -= 1; break;
            case 't':   cur -= 2; break;
            default:    continue;
        }
        data[i] = cur;
        ++i;
    }
    if(fin) z_normalize_and_save(data, len, of);

    // read other series
    while(fin){
        float tmp = data[step];
        for(int i = step; i<len;++i)
            data[i-step] = data[i] - tmp;
        cur = data[len - step - 1];
        for(int i = len - step; i < len && fin;){
            fin >> c;
            switch (tolower(c)) {
                case 'a':   cur += 1; break;
                case 'g':   cur += 2; break;
                case 'c':   cur -= 1; break;
                case 't':   cur -= 2; break;
                default:    continue;
            }
            data[i] = cur;
            ++i;
        }
        if(fin) z_normalize_and_save(data, len, of);
        if(++num % 200000 == 0) cout << num <<endl;
    }

    fin.close();
    fclose(of);
}

void DNATranslator::translateNonOverlap(const string &output, int len) {
    string fn = "/mnt/c/Series4Similarity_Search/dna/GRCh38_latest_genomic.fna";
    ifstream fin(fn, ios::in);
    FILE *of = fopen(output.c_str(), "ab");
    string s;
    getline(fin, s);
    float cur = 0; char c;
    long num = 1;
    float data[len];
    for(float &i:data)    i=0;

    while(fin){
        cur  = 0;
        for(int i=0;i<len && fin;){
            fin >> c;
            switch (tolower(c)) {
                case 'a':   cur += 1; break;
                case 'g':   cur += 2; break;
                case 'c':   cur -= 1; break;
                case 't':   cur -= 2; break;
                default:    continue;
            }
            data[i] = cur;
            ++i;
        }
        if(fin) z_normalize_and_save(data, len, of);
        if(++num % 200000 == 0) cout << num <<endl;
    }
    cout << "Total Series Number is " << num << endl;

    fin.close();
    fclose(of);

}

void DNATranslator::z_normalize_and_save(const float *ts, int size, FILE *f) {
    int i;
    double mean = 0;//gsl_stats_mean(ts, 1, size);
    double std = 0;//gsl_stats_sd(ts, 1, size);
    for (i = 0; i < size; i++) {
        mean += ts[i];
    }
    mean /= size;

    for (i = 0; i < size; i++) {
        std += (ts[i] - mean) * (ts[i] - mean);
    }
    std /= size;
    std = sqrt(std);
    float tmp[size];
    for (i = 0; i < size; i++) {
        tmp[i] = (ts[i] - mean) / std;
    }

    fwrite(tmp, sizeof(float), size, f);
}

void DNATranslator::z_normalize(float *ts, int size) {
    int i;
    double mean = 0;//gsl_stats_mean(ts, 1, size);
    double std = 0;//gsl_stats_sd(ts, 1, size);
    for (i = 0; i < size; i++) {
        mean += ts[i];
    }
    mean /= size;

    for (i = 0; i < size; i++) {
        std += (ts[i] - mean) * (ts[i] - mean);
    }
    std /= size;
    std = sqrt(std);
    for (i = 0; i < size; i++) {
        ts[i] = (ts[i] - mean) / std;
    }

}

void DNATranslator::readSIFT(){
    FILE *f = fopen("/mnt/g/Series4Similarity_Search/siftsmall_base.fvecs", "rb");
    float data[128]; int d;
    for(int i=0;i<10;++i){
        fread(&d, sizeof(int ),1 ,f);
        assert(d == 128);
        fread(data, sizeof(float), 128, f);
//        z_normalize(data, 128);
        cout << TimeSeriesUtil::timeSeries2Line(data);
    }
    fclose(f);
}