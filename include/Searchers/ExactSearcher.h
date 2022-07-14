//
// Created by wzy on 2021/8/9.
//

#ifndef MULGIFT_EXACTSEARCHER_H
#define MULGIFT_EXACTSEARCHER_H

#include <vector>
#include <string>
#include "../DataStructures/PqItemSeries.h"
using namespace std;

class ExactSearcher {

//    static void processSingleVertex(Vertex *v, Graph *g, TimeSeries *queryTs, int k, vector<PqItemSeries *> &heap,
//                                    bool *hasProcessed);

public:
//    static vector<PqItemSeries *> *exactKnnSearch(Graph *g, float *query, int k, vector<PqItemSeries *> *heap);

    static vector<PqItemSeries *> &
    groundTruthKnnSearch(const string &fn, const string &queryfn, int k, int query_num, int series_num);

    static void
    groundTruthKnnSearch(const string &fn, const string &queryfn, int k, int query_num, const string &output,
                         int series_num);

//    static vector<PqItemSeries *> *exactKnnSearch2(Graph *g, TimeSeries *queryTs, int k, vector<PqItemSeries *> *heap);
    static void groundTruthKnnSearchDTW(const string &fn, const string &queryfn, int k, int query_num, const string &output,
                                 int series_num);
};


#endif //MULGIFT_EXACTSEARCHER_H
