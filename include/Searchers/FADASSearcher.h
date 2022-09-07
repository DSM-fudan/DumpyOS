//
// Created by pengwang5 on 2022/1/16.
//

#ifndef MULGIFT_FADASSEARCHER_H
#define MULGIFT_FADASSEARCHER_H
#include <vector>
#include "../DataStructures/FADASNode.h"

class FADASSearcher {
public:
    static void approxSearchInterNode(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                      vector<PqItemSeries *> *heap, const string &index_dir);

    static vector<PqItemSeries *> *
    approxSearch(FADASNode *root, float *query, int k, vector<vector<int>> *g, const string &index_dir);

    static void approxSearchInterNodeDTW(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                                        vector<PqItemSeries *> *heap, const string &index_dir);

    static vector<PqItemSeries *> * approxSearchDTW(FADASNode *root, float *query, int k, vector<vector<int>> *g,
                                                            const string &index_dir);

    static vector<PqItemSeries *> *exactSearch(FADASNode *root, float *query, int k, vector<vector<int>> *g);

    static vector<PqItemSeries*>* exactSearchDTW(FADASNode* root, float *query, int k, vector<vector<int>> *g);

    static vector<PqItemSeries *> *exactSearchPos(FADASNode *root, float *query, int k, vector<vector<int>> *g);

    static vector<PqItemSeries *> * approxIncSearchDTW(FADASNode *root, float *query, int k, const string &index_dir,
                                                                      int node_num);

    static void approxIncSearchInterNodeDTW(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                                    vector<PqItemSeries *> *heap, const string &index_dir,int &node_num);

    static void
    approxSearchInterNodePos(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                             vector<PqItemSeries *> *heap,
                             const string &index_dir);

    static vector<PqItemSeries *> *
    approxSearchPos(FADASNode *root, float *query, int k, vector<vector<int>> *g, const string &index_dir);

    static vector<PqItemSeries *> *
    approxIncSearch(FADASNode *root, float *query, int k, const string &index_dir, int node_num, float * query_reordered, int*ordering);

    static
    vector<PqItemSeries *> *ngSearch(FADASNode *root, float *query, float *query_reordered, int *ordering, int k,
                                     vector<vector<int>> *g, int nprobes);

    static void
    approxIncSearchInterNode(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                             vector<PqItemSeries *> *heap, const string &index_dir,
                             int &node_num, float *query_reordered, int *ordering);

    static void approxIncSearchInterNodeFuzzy(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                       vector<PqItemSeries *> *heap, const string &index_dir, int &node_num,
                                       unordered_set<float *, createhash, isEqual> *hash_set);

    static vector<PqItemSeries *> *
    approxIncSearchFuzzy(FADASNode *root, float *query, int k, const string &index_dir, int node_num);

    static vector<PqItemSeries *> *
    approxSearchLessPack(FADASNode *root, float *query, int k, vector<vector<int>> *g, const string &index_dir);

    static void approxSearchInterNodeLessPack(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                                      vector<PqItemSeries *> *heap, const string &index_dir);

    static vector<PqItemSeries *> *
    ngSearchNaive(FADASNode *root, float *query, float *query_reordered, int *ordering, int k, vector<vector<int>> *g,
                  int nprobes);

    static vector<PqItemSeries *> *
    ngSearchIdLevelNaive(FADASNode *root, float *query, float *query_reordered, int *ordering, int k,
                         vector<vector<int>> *g, int nprobes);
};


#endif //MULGIFT_FADASSEARCHER_H
