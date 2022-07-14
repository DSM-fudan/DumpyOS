//
// Created by wzy on 2021/11/20.
//

#include <algorithm>
#include "../../include/Searchers/iSAXSearcher.h"
#include "../../include/DataStructures/TimeSeries.h"
#include "../../include/Utils/MathUtil.h"

extern int _search_num;
vector<PqItemSeries *> * iSAXSearcher::approxKnnSearch(iSAXRoot* root, float* query, int k){
    auto* queryTs = new TimeSeries(query);

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    iSAXNode *cur = root->route2Target(sax, false);
    auto *heap = new  vector<PqItemSeries *>();
    cur->exactSearchKnn(k, queryTs, *heap);
    _search_num = cur->size;

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

static float *t_paa;
extern int layer;

struct cand_isax{
    iSAXNode* node;
    double lb;

    cand_isax(iSAXNode*n, double _){
        node = n;
        lb = _;
    }
    cand_isax(){;}
};


bool comp_iSAX(const iSAXNode* x, const iSAXNode* y){
    if(x == nullptr)    return false;
    if(y == nullptr)    return true;
    return SaxUtil::LowerBound_Paa_iSax(t_paa, x->sax, x->layer) < SaxUtil::LowerBound_Paa_iSax(t_paa, y->sax, y->layer);
}

bool comp_iSAX_cand(const cand_isax& x, const cand_isax& y){
    return x.lb < y.lb;
}

vector<PqItemSeries *> * iSAXSearcher::iSAXApproxKnnSearch(iSAXRoot* root, float* query, int k, int threshold){
    auto* queryTs = new TimeSeries(query);
    auto*heap = new vector<PqItemSeries*>();
    t_paa = queryTs->paa;

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    int id1 = SaxUtil::invSaxHeadFromSax(sax, Const::bitsCardinality, Const::segmentNum);
    iSAXNode*node = root->children[id1];
    int rest = threshold;
    if(node == nullptr || node->size < threshold){
        layer = 1;
        vector<cand_isax>candidates;
        candidates.reserve(1 << Const::segmentNum);
        for(iSAXNode*child:root->children) {
            if (child != nullptr ) {
                int child_id = 0;
                for(unsigned short i : child->sax)
                    child_id = (child_id << 1) + i;
                candidates.emplace_back(child, MathUtil::bitDiffNum(id1,child_id, Const::segmentNum));
            }
        }
        sort(candidates.begin(), candidates.end(), comp_iSAX_cand);
        for(auto & candidate : candidates){
            if(candidate.node->size <= rest){
                candidate.node->exactSearchKnnInMemory(k, queryTs, *heap);
                rest -= candidate.node->size;
            }else{
                candidate.node->exactSearchKnnInMemoryEntry(k, queryTs, *heap, rest);
                break;
            }
        }
    }else   iSAXApproxKnnSearchSub(node, sax, queryTs, k, threshold, heap);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}


void iSAXSearcher::iSAXApproxKnnSearchSub( iSAXNode*root,  unsigned short *sax, TimeSeries* queryTs, int k, int& threshold, vector<PqItemSeries*>*heap){
    if(root==nullptr || threshold <=0)  return;
    if(root->isLeafNode()) { root->exactSearchKnnInMemory(k, queryTs, *heap, threshold);    return; }
    iSAXNode *node = root->route1Step(sax), *parent = root;
    while (node!= nullptr && !node->isLeafNode() && node->size > threshold) {
        parent = node;
        node = node->route1Step(sax);
    }

    layer = node->layer;
    node->exactSearchKnnInMemoryEntry(k, queryTs, *heap, threshold);
    if(threshold <= 0)  return;
    if(parent->left == node)
        iSAXApproxKnnSearchSub(parent->right, sax, queryTs, k, threshold, heap);
    else
        iSAXApproxKnnSearchSub(parent->left, sax, queryTs, k, threshold, heap);

}


