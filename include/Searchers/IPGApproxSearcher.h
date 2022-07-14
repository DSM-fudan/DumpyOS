//
// Created by wzy on 2021/11/9.
//

#ifndef MULGIFT_IPGAPPROXSEARCHER_H
#define MULGIFT_IPGAPPROXSEARCHER_H
#include "../DataStructures/PqItemSeries.h"
#include "../DataStructures/IPGNode.h"
#include <vector>

using namespace std;
class IPGApproxSearcher {

public:
    static vector<PqItemSeries *> *approxKnnSearchModel(IPGNode *root, float *query, int k, int threshold);

    static void
    approxKnnSearchModelsub(IPGNode *root, TimeSeries *queryTs, unsigned short *sax, int k, vector<PqItemSeries *> *heap);

    static vector<PqItemSeries *> *approxKnnSearchDynamicModel(IPGNode *root, float *query, int k, int threshold);

    static void approxKnnSearchDynamicModelsub(IPGNode *root, TimeSeries *queryTs, int k, vector<PqItemSeries *> *heap);

    static vector<PqItemSeries *> *approxKnnSearchDynamicModelRange(IPGNode *root, float *query, int k, int threshold);

    static void
    approxKnnSearchDynamicModelsubRange(IPGNode *root, TimeSeries *queryTs, int k, vector<PqItemSeries *> *heap);

    static vector<PqItemSeries *> *approxKnnSearchModelPaaMuAsSplit(IPGNode *root, float *query, int k);

    static void
    approxKnnSearchModelsubPaaMuAsSplit(IPGNode *root, TimeSeries *queryTs, int k,vector<PqItemSeries *> *heap);
};
#endif //MULGIFT_IPGAPPROXSEARCHER_H
