//
// Created by wzy on 2021/11/9.
//

#ifndef MULGIFT_IPGPARTITION_H
#define MULGIFT_IPGPARTITION_H


#include <unordered_map>
#include <vector>
#include "IPGNode.h"

class IPGPartition {

    static void savePartition(unordered_map<int, IPGNode *> *nodes_map, const string &output);

    static void outputPartition(unordered_map<int, IPGNode *> *nodes_map, const string &output, int k);

public:
//    static void bubble(unordered_map<int, IPGNode *> *nodes_map);

    static int growingGraph(unordered_map<int, IPGNode *> *nodes_map, vector<vector<int>> *g, double filling_factor,
                            int seg_num);

    static void zOrderCurve(unordered_map<int, IPGNode *> *nodes_map);

    static int findFirstLE(vector<IPGNode *> &nodes, int start, int end, int target);

    static int balanceMatch(unordered_map<int, IPGNode *> *nodes_map);
};


#endif //MULGIFT_IPGPARTITION_H
