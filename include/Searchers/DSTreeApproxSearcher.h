//
// Created by wzy on 2021/8/9.
//

#ifndef MULGIFT_DSTREEAPPROXSEARCHER_H
#define MULGIFT_DSTREEAPPROXSEARCHER_H
#include "../DSTree/DSTreeNode.h"
#include "../DataStructures/PqItemSeries.h"

class DSTreeApproxSearcher {

public:
    static void approxKnnSearch(DSTreeNode *root, float *queryTs, int k, vector<PqItemSeries *> &heap);

    static vector<PqItemSeries *> *
    approxKnnSearchWithThreshold(DSTreeNode *root, float *queryTs, int k, int threshold);

    static void approxKnnSearchWithThresholdSub(DSTreeNode *root, int k, int threshold, vector<PqItemSeries *> &heap,
                                         InsertedSeries *q);
};


#endif //MULGIFT_DSTREEAPPROXSEARCHER_H
