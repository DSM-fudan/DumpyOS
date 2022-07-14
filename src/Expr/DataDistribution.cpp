//
// Created by wzy on 2021/9/21.
//

#include "../../include/Expr/DataDistribution.h"
#include "../../include/Expr/Recall.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Utils/SaxUtil.h"
#include "../../include/DataStructures/TimeSeries.h"
#include "../../include/Const.h"
#include <string>
#include <algorithm>
#include <iostream>
using namespace std;
void DataDistribution::getInvSaxHeadDistribution(const string &fn, const string &outfn){
    long size = FileUtil::getFileSize(fn.c_str());
    long sax_char_nums = size / sizeof(unsigned short);
    auto *saxes = new unsigned short [sax_char_nums];
    FILE *f = fopen(fn.c_str(), "rb");
    fread(saxes, sizeof(unsigned short ), sax_char_nums, f);
    fclose(f);
    long series_num = sax_char_nums / Const::segmentNum;
    int results[Const::vertexNum];
    for(int & result : results)
        result = 0;
    for(int i=0;i<series_num;++i)
        results[SaxUtil::invSaxHeadFromSax(saxes + i * Const::segmentNum, Const::bitsCardinality, Const::segmentNum)]++;
    sort(results, results + Const::vertexNum);

    f = fopen(outfn.c_str(), "wb");
    fwrite(results, sizeof(int ), Const::vertexNum, f);
    fclose(f);

}

void DataDistribution::getPrecisionLostDistribution(const string & query_file, const string & result_file){
    int k = 500;
    double max = 2000 * k;
    int reserve[]{0,0,0,0,0,0};
    float query[Const::tsLength];
    FILE *qf = fopen(query_file.c_str(), "rb");
    unsigned short query_sax[Const::segmentNum];
    vector<unsigned short *>saxes(k);
    for(int i=0;i<k;++i)
        saxes[i] = new unsigned short[Const::segmentNum];

    for(int i=0;i<2000;++i){
        fread(query, sizeof(float ), Const::tsLength, qf);
        SaxUtil::saxFromTs(query, query_sax, Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
        vector<float*>* results = Recall::getResult(result_file, i, k);
        for(int j=0;j<k;++j)
            SaxUtil::saxFromTs((*results)[j], saxes[j], Const::tsLengthPerSegment, Const::segmentNum,
                               Const::cardinality);

        vector<vector<unsigned short *> >temp_reserve(6,vector<unsigned short *>());
        for(int j=1;j<=6;++j){
            // first determine where to fetch saxes of results
            vector<unsigned short*>* temp;
            if(j == 1)  temp = &saxes;
            else temp = &(temp_reserve[j - 2]);

            int query_head = SaxUtil::invSaxHeadkFromSax(query_sax, Const::bitsCardinality, Const::segmentNum, j);

            for(auto * n : *temp){
                int head = SaxUtil::invSaxHeadkFromSax(n, Const::bitsCardinality, Const::segmentNum, j);
                if(head == query_head)
                    temp_reserve[j-1].push_back(n);
            }
            reserve[j - 1] += temp_reserve[j-1].size();
        }
    }
    for(int j : reserve)
        cout << j / max <<endl;

}