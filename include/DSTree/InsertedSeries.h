//
// Created by wzy on 2021/8/7.
//

#ifndef MULGIFT_INSERTEDSERIES_H
#define MULGIFT_INSERTEDSERIES_H
#include <unordered_map>
#include "../DataStructures/TimeSeries.h"
#include "../Utils/MathUtil.h"
#include "../Const.h"

using namespace std;


class InsertedSeries {
public:
    float* ts;
    unordered_map<int, unordered_map<int, double>> segmentStdevs;
    double *prefixSum;

    explicit InsertedSeries(float* _ts){
        ts = _ts;
        prefixSum = new double [Const::tsLength];
        prefixSum[0] = 0;
        for(int i=1; i <= Const::tsLength; ++i)
            prefixSum[i] = prefixSum[i-1] + ts[i - 1];
    }

    double getMean(int start, int end) const{  //[,)
        return (prefixSum[end] - prefixSum[start]) / (double)(end - start);
    }

    double getStdEv(int start, int end){ //[,)  TODO:bottleNeck
        if(segmentStdevs.find(start) != segmentStdevs.end() && segmentStdevs[start].find(end) != segmentStdevs[start].end())
            return segmentStdevs[start][end];
        auto tmp = segmentStdevs[start];
        double res = MathUtil::deviation(ts, start, end, getMean(start, end));
        tmp[end] = res;
        return res;
    }

    ~InsertedSeries(){
        delete[] prefixSum;
        auto iter = segmentStdevs.begin();
        while (iter!=segmentStdevs.end()){
            unordered_map<int,double>().swap(iter->second);
            iter++;
        }
        unordered_map<int, unordered_map<int, double>>().swap(segmentStdevs);
    }
};


#endif //MULGIFT_INSERTEDSERIES_H
