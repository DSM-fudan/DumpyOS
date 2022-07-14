//
// Created by wzy on 2021/11/20.
//

#ifndef MULGIFT_ISAXSEARCHER_H
#define MULGIFT_ISAXSEARCHER_H


#include "../DataStructures/PqItemSeries.h"
#include "../DataStructures/iSAXNode.h"

class iSAXSearcher {
public:
    static vector<PqItemSeries *> * approxKnnSearch(iSAXRoot* root, float* query, int k);


    static vector<PqItemSeries *> *iSAXApproxKnnSearch(iSAXRoot *root, float *query, int k, int threshold);

    static void iSAXApproxKnnSearchSub(iSAXNode *root, unsigned short *sax, TimeSeries *queryTs, int k, int &threshold,
                                vector<PqItemSeries *> *heap);
};


#endif //MULGIFT_ISAXSEARCHER_H
