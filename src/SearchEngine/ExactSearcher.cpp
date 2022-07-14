//
// Created by wzy on 2021/8/9.
//
#include <iostream>
#include <cmath>
#include "../../include/Searchers/ExactSearcher.h"
#include "../../include/Searchers/DSTreeExactSearcher.h"
#include "../../include/Utils/SaxUtil.h"
#include "../../include/Const.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Utils/TimeSeriesUtil.h"
#include "../../include/Const.h"
//void ExactSearcher::processSingleVertex(Vertex* v, Graph* g , TimeSeries* queryTs, int k,
//                                        vector<PqItemSeries *> &heap, bool* hasProcessed){
//    if(v->level == huge)
//        exactKnnSearch2(v->g, queryTs, k, &heap);
//    else if (v->level == big) {
//        DSTreeExactSearcher::exactSearchKnnWithBsf(v->root, queryTs->ts, k , heap);
//    } else if (v->level == small) {
//        v->smallVertexKnnSearch(*queryTs, heap, k, -1);
//    } else {
//        v->middleVertexApproxKnnSearchWithoutBufferNew(*queryTs, k, heap, -1, g->rowDataFileName);
//    }
//    hasProcessed[v->id] = true;
//}
//
//vector<PqItemSeries *> * ExactSearcher::exactKnnSearch(Graph *g, float *query, int k, vector<PqItemSeries *> *heap) {
//    if(heap == nullptr) heap = new vector<PqItemSeries*>;
//
//    auto* queryTs = new TimeSeries(query);
//    int targetId = SaxUtil::invSaxHeadFromSax(queryTs->sax, TimeSeries::bitsCardinality, TimeSeries::segmentNum);
//
//    bool* hasProcessed = new bool[Const::vertexNum];
//    for(int i=0;i<Const::vertexNum;++i)
//        hasProcessed[i] = false;
//
//    //        System.out.println("Start process the destination vertex.");
//    // 0
//    Vertex* targetVertex = g->VERTEXES[targetId];
//    if(targetVertex->cnt > 0)
//        processSingleVertex(targetVertex, g, queryTs, k, *heap, hasProcessed);
//
//    //        System.out.println("Start process the neighborhood.");
//    // nnList
//    vector<Vertex *> &nnList = g->NeighborsList[targetId];
//    for (auto v : nnList) {
//        if (v->cnt > 0)
//            processSingleVertex(v, g, queryTs, k, *heap, hasProcessed);
//    }
//
//    //        System.out.println("Start process other vertices.");
//    // others
//    for(int i=0;i<Const::vertexNum; ++i){
//        Vertex* v = g->VERTEXES[i];
//        if(v->cnt > 0 && !hasProcessed[i])
//            processSingleVertex(v, g ,queryTs, k, *heap, hasProcessed);
//    }
//
//    sort_heap(heap->begin(),  heap->end(), PqItemSeriesMaxHeap());    // ascending order
//    delete[] hasProcessed;
//    delete queryTs;
//
//    return heap;
//}
//
//vector<PqItemSeries *> * ExactSearcher::exactKnnSearch2(Graph *g, TimeSeries *queryTs, int k, vector<PqItemSeries *> *heap) {
//    if(heap == nullptr) heap = new vector<PqItemSeries*>;
//
//    int targetId = SaxUtil::invSaxHead2FromSax(queryTs->sax, TimeSeries::bitsCardinality, TimeSeries::segmentNum);
//
//    bool* hasProcessed = new bool[Const::vertexNum];
//    for(int i=0;i<Const::vertexNum;++i)
//        hasProcessed[i] = false;
//
//    //        System.out.println("Start process the destination vertex.");
//    // 0
//    Vertex* targetVertex = g->VERTEXES[targetId];
//    if(targetVertex->cnt > 0)
//        processSingleVertex(targetVertex, g, queryTs, k, *heap, hasProcessed);
//
//    //        System.out.println("Start process the neighborhood.");
//    // nnList
//    vector<Vertex *> &nnList = g->NeighborsList[targetId];
//    for (auto v : nnList) {
//        if (v->cnt > 0)
//            processSingleVertex(v, g, queryTs, k, *heap, hasProcessed);
//    }
//
//    //        System.out.println("Start process other vertices.");
//    // others
//    for(int i=0;i<Const::vertexNum; ++i){
//        Vertex* v = g->VERTEXES[i];
//        if(v->cnt > 0 && !hasProcessed[i])
//            processSingleVertex(v, g ,queryTs, k, *heap, hasProcessed);
//    }
//
//    delete[] hasProcessed;
//
//    return heap;
//}
//vector<PqItemSeries *> & ExactSearcher::groundTruthKnnSearch(const string &fn, const string &queryfn, int k,
//                                                             int query_num, int series_num) {
//    FILE  *qf = fopen(queryfn.c_str(), "rb");
//    auto *queries = new float [query_num * Const::tsLength];
//    fread(queries, sizeof(float ), (long)query_num * Const::tsLength, qf);
//    fclose(qf);
//
//    FILE *f = fopen(fn.c_str(), "rb");
//    int seriesNum = FileUtil::getFileSize(fn.c_str()) / Const::tsLengthBytes;
//    auto *heap = new vector<PqItemSeries*>;
//    double bsfMin = numeric_limits<double>::max();
//    for(int i=0;i<seriesNum;++i){
//        auto ts = FileUtil::readSeries(f);
//        double dist = TimeSeriesUtil::euclideanDist(ts, queries, Const::tsLength, bsfMin);
//        if(dist >= bsfMin)  {
//            delete ts;
//            continue;
//        }
//        if(heap->size() < k){
//            heap->push_back(new PqItemSeries(ts, dist, true));
//            push_heap(heap->begin(),heap->end(),PqItemSeriesMaxHeap());
//        } else{
//            pop_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
//            delete heap->back();
//            heap->pop_back();
//            heap->push_back(new PqItemSeries(ts, dist, true));
//            push_heap(heap->begin(),heap->end(),PqItemSeriesMaxHeap());
//        }
//        if(heap->size() >= k)
//            bsfMin = (*heap)[0]->dist;
//    }
//    fclose(f);
//    sort_heap(heap->begin(),  heap->end(), PqItemSeriesMaxHeap());    //ascending order
//    return *heap;
//}

class TempSeries {
public:
    float * ts;
    double dist;


    TempSeries(float *t, double d){
        ts = t;
        dist = d;
    }

    ~TempSeries(){
        delete[] ts;
    }

};

struct TempSeriesMaxHeap{
    bool operator()(const TempSeries* x, const TempSeries* y){
        return x->dist < y->dist;
    }
};

struct TempSeriesMinHeap{
    bool operator()(const TempSeries* x, const TempSeries* y){
        return x->dist > y->dist;
    }
};

void ExactSearcher::groundTruthKnnSearch(const string &fn, const string &queryfn, int k, int query_num, const string &output,
                                         int series_num) {
    FILE  *qf = fopen(queryfn.c_str(), "rb");
    auto *queries = new float [query_num * Const::tsLength];
    fread(queries, sizeof(float ), (long)query_num * Const::tsLength, qf);
    fclose(qf);

    cout << "k:" << k <<endl;
    FILE *f = fopen(fn.c_str(), "rb");
    if(series_num == -1)
        series_num = FileUtil::getFileSize(fn.c_str()) / Const::tsLengthBytes;
    vector<vector<TempSeries*>>heaps(query_num, vector<TempSeries*>());
    vector<double>bsfMins(query_num, 500000);

    long bulk = 1000000;
    auto* ts = new float[Const::tsLength * bulk];
    long cur = 0;


    while (cur < series_num){
        cout << cur <<endl;
        fread(ts, sizeof(float ), Const::tsLength * bulk, f);

        for(long i=0;i<bulk;++i) {
            if(isnan(ts[Const::tsLength * i])) continue;
            for (int j = 0; j < query_num; ++j) {
                double dist = TimeSeriesUtil::euclideanDist(ts + Const::tsLength * i, queries + j * Const::tsLength, Const::tsLength, bsfMins[j]);
                if (dist >= bsfMins[j])
                    continue;
                auto *tmp = new float[Const::tsLength];
                copy(ts + Const::tsLength * i, ts + Const::tsLength * (i + 1), tmp);
                auto *t = new TempSeries(tmp, dist);
                if (heaps[j].size() < k) {
                    heaps[j].push_back(t);
                    push_heap(heaps[j].begin(), heaps[j].end(), TempSeriesMaxHeap());
                } else {
                    pop_heap(heaps[j].begin(), heaps[j].end(), TempSeriesMaxHeap());
                    delete heaps[j].back();
                    heaps[j].pop_back();
                    heaps[j].push_back(t);
                    push_heap(heaps[j].begin(), heaps[j].end(), TempSeriesMaxHeap());
                }
                if (heaps[j].size() >= k)
                    bsfMins[j] = heaps[j][0]->dist;
            }
        }

        cur += bulk;
    }
    fclose(f);

    delete[] ts;
    for(int j=0; j < query_num; ++j)
        sort(heaps[j].begin(),  heaps[j].end(), TempSeriesMaxHeap());
    f = fopen(output.c_str(), "wb");
    for(int j=0; j < query_num; ++j){
        cout << j << ": " << heaps[j].size() << endl;
        for(auto *tmps:heaps[j]){
            fwrite(tmps->ts, sizeof(float ), Const::tsLength, f);
        }
    }
    fclose(f);
}


void ExactSearcher::groundTruthKnnSearchDTW(const string &fn, const string &queryfn, int k, int query_num, const string &output,
                                         int series_num) {
    FILE  *qf = fopen(queryfn.c_str(), "rb");
    auto *queries = new float [query_num * Const::tsLength];
    fread(queries, sizeof(float ), (long)query_num * Const::tsLength, qf);
    fclose(qf);

    cout << "k:" << k <<endl;
    FILE *f = fopen(fn.c_str(), "rb");
    if(series_num == -1)
        series_num = FileUtil::getFileSize(fn.c_str()) / Const::tsLengthBytes;
    vector<vector<TempSeries*>>heaps(query_num, vector<TempSeries*>());
    vector<double>bsfMins(query_num, 500000);

    long bulk = 1000000;
    auto* ts = new float[Const::tsLength * bulk];
    long cur = 0;


    while (cur < series_num){
        cout << cur <<endl;
        fread(ts, sizeof(float ), Const::tsLength * bulk, f);

        for(long i=0;i<bulk;++i) {
            if(isnan(ts[Const::tsLength * i])) continue;
            for (int j = 0; j < query_num; ++j) {
                double dist = TimeSeriesUtil::dtw(ts + Const::tsLength * i, queries + j * Const::tsLength, Const::tsLength,Const::dtw_window_size, bsfMins[j]);
                if (dist >= bsfMins[j])
                    continue;
                auto *tmp = new float[Const::tsLength];
                copy(ts + Const::tsLength * i, ts + Const::tsLength * (i + 1), tmp);
                auto *t = new TempSeries(tmp, dist);
                if (heaps[j].size() < k) {
                    heaps[j].push_back(t);
                    push_heap(heaps[j].begin(), heaps[j].end(), TempSeriesMaxHeap());
                } else {
                    pop_heap(heaps[j].begin(), heaps[j].end(), TempSeriesMaxHeap());
                    delete heaps[j].back();
                    heaps[j].pop_back();
                    heaps[j].push_back(t);
                    push_heap(heaps[j].begin(), heaps[j].end(), TempSeriesMaxHeap());
                }
                if (heaps[j].size() >= k)
                    bsfMins[j] = heaps[j][0]->dist;
            }
        }

        cur += bulk;
    }
    fclose(f);

    delete[] ts;
    for(int j=0; j < query_num; ++j)
        sort(heaps[j].begin(),  heaps[j].end(), TempSeriesMaxHeap());
    f = fopen(output.c_str(), "wb");
    for(int j=0; j < query_num; ++j){
        cout << j << ": " << heaps[j].size() << endl;
        for(auto *tmps:heaps[j]){
            fwrite(tmps->ts, sizeof(float ), Const::tsLength, f);
        }
    }
    fclose(f);
}
