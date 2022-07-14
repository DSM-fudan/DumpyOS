//
// Created by wzy on 2021/8/7.
//

#ifndef MULGIFT_PQITEMSERIES_H
#define MULGIFT_PQITEMSERIES_H
#include <vector>
using namespace std;

class PqItemSeries {
public:
    float* ts{};
    double dist{};
    bool needFree= false, needDeepCopy = false;

    PqItemSeries()= default;

    PqItemSeries(float*t, float*q, bool free=false, bool ndc = false);

    PqItemSeries(float*t, double d, bool free= false, bool ndc = false){
        ts = t;
        dist = d;
        needDeepCopy = ndc;
        needFree = free;
    }

    PqItemSeries(const vector<float>&t, double d, bool free= false, bool ndc = false){
        ts = new float[t.size()];
        for(int i=0;i<t.size();++i)
            ts[i] = t[i];
        dist = d;
        needDeepCopy = ndc;
        needFree = free;
    }

    ~PqItemSeries(){
        if(needFree)
            delete[] ts;
    }

    void copyData();

};

struct PqItemSeriesMaxHeap{
    bool operator()(const PqItemSeries* x, const PqItemSeries* y){
        return x->dist < y->dist;
    }
};

struct PqItemSeriesMinHeap{
    bool operator()(const PqItemSeries* x, const PqItemSeries* y){
        return x->dist > y->dist;
    }
};

struct PqItemSeriesMaxHeap2{
    bool operator()(const PqItemSeries& x, const PqItemSeries& y){
        return x.dist < y.dist;
    }
};

#endif //MULGIFT_PQITEMSERIES_H
