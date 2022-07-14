//
// Created by wzy on 2022/7/2.
//

#ifndef FADAS_TARSEARCHER_H
#define FADAS_TARSEARCHER_H


#include <vector>
#include "TARGNode.h"
#include "../../include/DataStructures/PqItemSeries.h"


class TARSearcher {
public:
    static vector<PqItemSeries *> * approxSearch(TARGNode *root, float *query, int k, const string &index_dir);

    static vector<PqItemSeries *> *
    approxSearchDTW(TARGNode *root, float *query, int k, const string &index_dir, int window_size);

    static void exactSearchLocal(TARLNode *l_root, vector<PqItemSeries *> *heap, int k, float *query);

    static vector<PqItemSeries *> *
    approxIncSearch(TARGNode *root, float *query, int k, const string &index_dir, int *node_num);

    static void incSearchLocal(TARLNode *l_root, vector<PqItemSeries *> *heap, int k, float *query, int *node_num,
                               const string &query_invsax_str);

    static vector<PqItemSeries *> *exactSearch(TARGNode *root, float *query, int k, const string &index_dir);

    static void
    exactSearchLocalDTW(TARLNode *l_root, vector<PqItemSeries *> *heap, int k, float *query, const string &index_dir);

    static vector<PqItemSeries *> *exactSearchDTW(TARGNode *root, float *query, int k, const string &index_dir);

    static void approxIncSearchSub(TARGNode *root, float *query, int k, const string &index_dir, int *node_num,
                            vector<PqItemSeries *> *heap, string &query_invsax_str);
};


#endif //FADAS_TARSEARCHER_H
