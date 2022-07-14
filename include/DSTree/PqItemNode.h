//
// Created by wzy on 2021/8/7.
//

#ifndef MULGIFT_PQITEMNODE_H
#define MULGIFT_PQITEMNODE_H
#include "DSTreeNode.h"


class PqItemNode {
public:
    PqItemNode() {;}

    PqItemNode(DSTreeNode*x, double d) {
        node = x;
        dist = d;
    }

    DSTreeNode* node{};
    double dist{};
};

struct PqItemNodeComparator{
    bool operator()(PqItemNode*& x, PqItemNode*& y){
        return x->dist > y->dist;
    }
};


class PqItemNodeSeries{
public:
    DSTreeNode* node;
    int offset;
    double dist;

    PqItemNodeSeries(DSTreeNode* _node, int i, double minDist) {
        node = _node;
        offset = i;
        dist = minDist;
    }

};

struct PqItemNodeSeriesMaxHeap{
    bool operator()(PqItemNodeSeries*& x, PqItemNodeSeries*& y){
        return x->dist < y->dist;
    }
};

#endif //MULGIFT_PQITEMNODE_H
