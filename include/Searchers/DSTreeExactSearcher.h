//
// Created by wzy on 2021/8/9.
//

#ifndef MULGIFT_DSTREEEXACTSEARCHER_H
#define MULGIFT_DSTREEEXACTSEARCHER_H
#include "../../include/DataStructures/PqItemSeries.h"
#include "../DSTree/PqItemNode.h"


class DSTreeExactSearcher {

//    static void processPqInKnnSearch(InsertedSeries *q, int k, vector<PqItemSeriesVector *> &heap, DSTreeNode *bsfDSTreeNode,
//                                     vector<PqItemNode *> &pq);
//
//    static void
//    processPqInKnnSearchWithThreshold(InsertedSeries *q, int k, vector<PqItemSeriesVector *> &heap, DSTreeNode *bsfDSTreeNode,
//                                      vector<PqItemNode *> &pq, int threshold);

public:
    static void exactSearchKnnWithBsf(DSTreeNode *root, float *queryTs, int k, vector<PqItemSeries *> &heap);

//    static void
//    exactSearchKnnWithBsfAndThreshold(DSTreeNode *root, float *queryTs, int k, vector<PqItemSeriesVector *> &heap,
//                                      int threshold);

    static vector<PqItemSeries *> & exactSearchKnn(DSTreeNode *root, float *queryTs, int k);

//    static void exactSearchKnnWithBsfAndThreshold2(DSTreeNode* root, float* queryTs, int k, vector<PqItemSeriesVector *> &heap, int threshold);
//
//    static void processPqInKnnSearchWithThreshold2(InsertedSeries *q, int k, vector<PqItemSeriesVector *> &heap,
//                                                   DSTreeNode *bsfDSTreeNode, vector<PqItemNode *> &pq, int threshold);

    static void processPqInKnnSearch(InsertedSeries *q, int k, vector<PqItemSeries *> &heap, DSTreeNode *bsfDSTreeNode,
                              vector<PqItemNode *> &pq);

    static void
    processPqInKnnSearchWithThreshold(InsertedSeries *q, int k, vector<PqItemSeries *> &heap, DSTreeNode *bsfDSTreeNode,
                                      vector<PqItemNode *> &pq, int threshold);

    static void
    exactSearchKnnWithBsfAndThreshold(DSTreeNode *root, float *queryTs, int k, vector<PqItemSeries *> &heap,
                                      int threshold);
};


#endif //MULGIFT_DSTREEEXACTSEARCHER_H
