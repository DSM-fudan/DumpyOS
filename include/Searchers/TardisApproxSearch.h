//
// Created by wzy on 2021/9/24.
//

#ifndef MULGIFT_TARDISAPPROXSEARCH_H
#define MULGIFT_TARDISAPPROXSEARCH_H


#include <vector>
#include "../Tardis/TardisTreeNode.h"
#include "../DataStructures/PqItemSeries.h"

class TardisApproxSearch {
public:
    static std::vector<PqItemSeries *> *tardisApproxKnnSearch(TardisTreeNode *root, float *query, int k, int threshold);

    static vector<PqItemSeries *> *tardisApproxKnnSearchModel(TardisTreeNode *root, float *query, int k);

    static void tardisApproxKnnSearchSub(TardisTreeNode *root, unsigned short *sax, TimeSeries *queryTs, int k, int threshold,
                                         vector<PqItemSeries *> *heap);

    static void
    tardisApproxKnnSearchModelsub(TardisTreeNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                  vector<PqItemSeries *> *heap);

    static void
    tardisApproxKnnIncSearchModelSub(TardisTreeNode *root, unsigned short *sax, TimeSeries *queryTs, int k,
                                     int &node_num,
                                     vector<PqItemSeries *> *heap);

    static vector<PqItemSeries *> *tardisApproxKnnIncSearchModel(TardisTreeNode *root, float *query, int k, int node_num);

    static void
    tardisApproxKnnSearchSubGraph(TardisTreeNode *root, unsigned short *sax, TimeSeries *queryTs, int k, int &threshold,
                                  vector<PqItemSeries *> *heap, vector<vector<int>> *g);

    static vector<PqItemSeries *> *
    tardisApproxKnnSearchGraph(TardisTreeNode *root, float *query, int k, int threshold, vector<vector<int>> *g);
};


#endif //MULGIFT_TARDISAPPROXSEARCH_H
