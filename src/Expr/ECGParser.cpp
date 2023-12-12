//
// Created by pengwang5 on 2022/2/13.
//

#include "../../include/Expr/ECGParser.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Expr/DNATranslator.h"
#include <cmath>

void ECGParser::generateECG(const string &data_dir, const string &output, int ts_length) {
    vector<string>files;
    FileUtil::getFiles(data_dir, files);
    auto *series = new float [ts_length];
    int cur = 0;
    float t;
    FILE *outf = fopen(output.c_str(), "wb");
    for(string& file:files){
        cur = 0;
        FILE *f = fopen(file.c_str(), "rb");
        long size= FileUtil::getFileSize(file.c_str());
        int num = size / sizeof(float);
        for(int i=0;i<num;++i){
            fread(&t, sizeof(float), 1, f);
            if(isnan(t))    continue;
            series[cur++] = t;
            if(cur >= ts_length){
                DNATranslator::z_normalize(series, ts_length);
                fwrite(series, sizeof(float ), ts_length, outf);
                cur = 0;
            }
        }
        fclose(f);
    }
    fclose(outf);
}