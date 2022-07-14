//
// Created by wzy on 2021/11/9.
//

#include "../../include/Searchers/IPGApproxSearcher.h"
#include "../../include/Utils/MathUtil.h"
#include "../../include/Const.h"
#include <iostream>
#include <iomanip>
#include <unordered_set>

static double *t_paa;
extern int layer;
extern vector<IPGNode*>c_nodes;

vector<PqItemSeries *> * IPGApproxSearcher::approxKnnSearchModel(IPGNode* root, float* query, int k, int threshold){
    auto* queryTs = new TimeSeries(query);
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    int pid = -1;
    IPGNode *cur = (*root->children)[SaxUtil::invSaxHeadFromSax(sax, Const::bitsCardinality, Const::segmentNum)];
    if(cur!= nullptr && !cur->isLeafNode())
        approxKnnSearchModelsub(cur, queryTs, sax, k, heap);
    else {
//        layer = 1; delete queryTs; return heap;
        IPGNode *node;
        if (cur != nullptr) { pid = cur->partitionId; node = cur; }
        else {   //cur is nullptr
            double min_dist = numeric_limits<double>::max(), max_size = 0;
            for(auto &iter:*root->children){
                if(iter.second == nullptr)  continue;
                double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, iter.second->sax, 1);
                if(dist < min_dist){
                    min_dist = dist;
                    pid = iter.second->partitionId;
                    max_size = iter.second->size;
                    node = iter.second;
                }else if(dist == min_dist && iter.second->size > max_size){
                    max_size = iter.second->size;
                    pid = iter.second->partitionId;
                    node = iter.second;
                }
            }
            // we only concern whether the nearest node is a leaf or an internal node
            if(!node->isLeafNode()){
                approxKnnSearchModelsub(node, queryTs, sax, k, heap);
                goto FINISH;
            }
        }
        layer = 1;
//        if(pid == -1)   node->exactSearchKnn(k, queryTs, *heap);
//        else{
            assert(node->partitionId != -1);
            for(auto &iter:*root->children){
                if(iter.second!= nullptr && iter.second->partitionId == pid){
                    iter.second->exactSearchKnn(k, queryTs, *heap);
                }
            }
//        }
    }

FINISH:
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void IPGApproxSearcher::approxKnnSearchModelsub(IPGNode* root, TimeSeries* queryTs, unsigned short *sax, int k, vector<PqItemSeries*>*heap){
    IPGNode *cur = root->route(sax), *parent = root;
    while (cur!= nullptr && !cur->isLeafNode()) {
        parent = cur;
        cur = cur->route(sax);
    }
    int pid = -1;
    layer = parent->layer+1;
    if(cur != nullptr)
        pid = cur->partitionId;
    else{
        double min_dist = numeric_limits<double>::max(), max_size = 0;
        IPGNode *node;
        for(auto &iter:*parent->children){
            if(iter.second == nullptr)  continue;
            double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, iter.second->sax, iter.second->bits_cardinality);
            if(dist < min_dist){
                min_dist = dist;
                pid = iter.second->partitionId;
                max_size = iter.second->size;
                node = iter.second;
            }else if(dist == min_dist && iter.second->size > max_size){
                max_size = iter.second->size;
                pid = iter.second->partitionId;
                node = iter.second;
            }
        }
        // we only concern whether the nearest node is a leaf or an internal node
        if(!node->isLeafNode()){
            approxKnnSearchModelsub(node, queryTs, sax, k, heap);
            return;
        }
    }

//    if(pid == -1)   cur->exactSearchKnn(k, queryTs, *heap);
//    else{
        assert(pid != -1);
        for(auto &iter:*parent->children){
            if(iter.second!= nullptr && iter.second->partitionId == pid){
                iter.second->exactSearchKnn(k, queryTs, *heap);
            }
        }
//    }
}

vector<PqItemSeries *> * IPGApproxSearcher::approxKnnSearchDynamicModel(IPGNode* root, float* query, int k, int threshold){
    auto* queryTs = new TimeSeries(query);
    auto*heap = new vector<PqItemSeries*>();

    int pid = -1;
    IPGNode *cur = (*root->children)[SaxUtil::invSaxHeadFromPaa(queryTs->paa, Const::bitsCardinality, Const::segmentNum)];
    if(cur!= nullptr && !cur->isLeafNode())
        approxKnnSearchDynamicModelsub(cur, queryTs, k, heap);
    else {
//        layer = 1; delete queryTs; return heap;
        IPGNode *node;
        if (cur != nullptr) { pid = cur->partitionId; node = cur; }
        else {   //cur is nullptr
            double min_dist = numeric_limits<double>::max(), max_size = 0;
            for(auto &iter:*root->children){
                if(iter.second == nullptr)  continue;
                double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, iter.second->sax, 1);
                if(dist < min_dist){
                    min_dist = dist;
                    pid = iter.second->partitionId;
                    max_size = iter.second->size;
                    node = iter.second;
                }else if(dist == min_dist && iter.second->size > max_size){
                    max_size = iter.second->size;
                    pid = iter.second->partitionId;
                    node = iter.second;
                }
            }
            // we only concern whether the nearest node is a leaf or an internal node
            if(!node->isLeafNode()){
                approxKnnSearchDynamicModelsub(node, queryTs, k, heap);
                goto FINISH;
            }
        }
        layer = 1;
//        if(pid == -1)   node->exactSearchKnn(k, queryTs, *heap);
//        else{
        assert(node->partitionId != -1);
        for(auto &iter:*root->children){
            if(iter.second->partitionId == pid){
                iter.second->exactSearchKnn(k, queryTs, *heap);
            }
        }
//        }
    }

    FINISH:
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void IPGApproxSearcher::approxKnnSearchDynamicModelsub(IPGNode* root, TimeSeries* queryTs, int k,vector<PqItemSeries*>*heap){
    IPGNode *cur = root->route(queryTs->paa), *parent = root;
    while (cur!= nullptr && !cur->isLeafNode()) {
        parent = cur;
        cur = cur->route(queryTs->paa);
    }
    int pid = -1;
    layer = parent->layer+1;
    if(cur != nullptr)
        pid = cur->partitionId;
    else{
        double min_dist = numeric_limits<double>::max(), max_size = 0;
        IPGNode *node;
        for(auto &iter:*parent->children){
            if(iter.second == nullptr)  continue;
            double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, iter.second->sax, iter.second->bits_cardinality);
            if(dist < min_dist){
                min_dist = dist;
                pid = iter.second->partitionId;
                max_size = iter.second->size;
                node = iter.second;
            }else if(dist == min_dist && iter.second->size > max_size){
                max_size = iter.second->size;
                pid = iter.second->partitionId;
                node = iter.second;
            }
        }
        // we only concern whether the nearest node is a leaf or an internal node
        if(!node->isLeafNode()){
            approxKnnSearchDynamicModelsub(node, queryTs, k, heap);
            return;
        }
    }

//    if(pid == -1)   cur->exactSearchKnn(k, queryTs, *heap);
//    else{
    assert(pid != -1);
    for(auto &iter:*parent->children){
        if(iter.second->partitionId == pid){
            iter.second->exactSearchKnn(k, queryTs, *heap);
        }
    }
//    }
}

vector<PqItemSeries *> * IPGApproxSearcher::approxKnnSearchDynamicModelRange(IPGNode* root, float* query, int k, int threshold){
    auto* queryTs = new TimeSeries(query);
    auto*heap = new vector<PqItemSeries*>();

    int pid = -1;
    IPGNode *cur = root->routeIn1stLayer(queryTs->paa);
    if(cur!= nullptr && !cur->isLeafNode())
        approxKnnSearchDynamicModelsubRange(cur, queryTs, k, heap);
    else {
//        layer = 1; delete queryTs; return heap;
        IPGNode *node;
        if (cur != nullptr) { pid = cur->partitionId; node = cur; }
        else {   //cur is nullptr
            double min_dist = numeric_limits<double>::max(), max_size = 0;
            for(auto &iter:*root->children){
                if(iter.second == nullptr)  continue;
                double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, iter.second->sax, 1);
                if(dist < min_dist){
                    min_dist = dist;
                    pid = iter.second->partitionId;
                    max_size = iter.second->size;
                    node = iter.second;
                }else if(dist == min_dist && iter.second->size > max_size){
                    max_size = iter.second->size;
                    pid = iter.second->partitionId;
                    node = iter.second;
                }
            }
            // we only concern whether the nearest node is a leaf or an internal node
            if(!node->isLeafNode()){
                approxKnnSearchDynamicModelsubRange(node, queryTs, k, heap);
                goto FINISH;
            }
        }
        layer = 1;
//        if(pid == -1)   node->exactSearchKnn(k, queryTs, *heap);
//        else{
        assert(node->partitionId != -1);
        for(auto &iter:*root->children){
            if(iter.second->partitionId == pid){
                iter.second->exactSearchKnn(k, queryTs, *heap);
            }
        }
//        }
    }

    FINISH:
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void IPGApproxSearcher::approxKnnSearchDynamicModelsubRange(IPGNode* root, TimeSeries* queryTs, int k,vector<PqItemSeries*>*heap){
    IPGNode *cur = root->routeBelow1stLayer(queryTs->paa), *parent = root;
    while (cur!= nullptr && !cur->isLeafNode()) {
        parent = cur;
        cur = cur->routeBelow1stLayer(queryTs->paa);
    }
    int pid = -1;
    layer = parent->layer+1;
    if(cur != nullptr)
        pid = cur->partitionId;
    else{
        double min_dist = numeric_limits<double>::max(), max_size = 0;
        IPGNode *node;
        for(auto &iter:*parent->children){
            if(iter.second == nullptr)  continue;
            double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, iter.second->sax, iter.second->bits_cardinality);
            if(dist < min_dist){
                min_dist = dist;
                pid = iter.second->partitionId;
                max_size = iter.second->size;
                node = iter.second;
            }else if(dist == min_dist && iter.second->size > max_size){
                max_size = iter.second->size;
                pid = iter.second->partitionId;
                node = iter.second;
            }
        }
        // we only concern whether the nearest node is a leaf or an internal node
        if(!node->isLeafNode()){
            approxKnnSearchDynamicModelsub(node, queryTs, k, heap);
            return;
        }
    }

//    if(pid == -1)   cur->exactSearchKnn(k, queryTs, *heap);
//    else{
    assert(pid != -1);
    for(auto &iter:*parent->children){
        if(iter.second->partitionId == pid){
            iter.second->exactSearchKnn(k, queryTs, *heap);
        }
    }
//    }
}

vector<PqItemSeries *> * IPGApproxSearcher::approxKnnSearchModelPaaMuAsSplit(IPGNode* root, float* query, int k){
    auto* queryTs = new TimeSeries(query);
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    int pid = -1;
    IPGNode *cur = (*root->children)[SaxUtil::invSaxHeadFromSax(sax, Const::bitsCardinality, Const::segmentNum)];
    if(cur!= nullptr && !cur->isLeafNode())
        approxKnnSearchModelsubPaaMuAsSplit(cur, queryTs, k, heap);
    else {
//        layer = 1; delete queryTs; return heap;
        IPGNode *node;
        if (cur != nullptr) { pid = cur->partitionId; node = cur; }
        else {   //cur is nullptr
            double min_dist = numeric_limits<double>::max(), max_size = 0;
            for(auto &iter:*root->children){
                if(iter.second == nullptr)  continue;
                double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, iter.second->sax, 1);
                if(dist < min_dist){
                    min_dist = dist;
                    pid = iter.second->partitionId;
                    max_size = iter.second->size;
                    node = iter.second;
                }else if(dist == min_dist && iter.second->size > max_size){
                    max_size = iter.second->size;
                    pid = iter.second->partitionId;
                    node = iter.second;
                }
            }
            // we only concern whether the nearest node is a leaf or an internal node
            if(!node->isLeafNode()){
                approxKnnSearchModelsub(node, queryTs, sax, k, heap);
                goto FINISH;
            }
        }
        layer = 1;
//        if(pid == -1)   node->exactSearchKnn(k, queryTs, *heap);
//        else{
        assert(node->partitionId != -1);
        for(auto &iter:*root->children){
            if(iter.second->partitionId == pid){
                iter.second->exactSearchKnn(k, queryTs, *heap);
            }
        }
//        }
    }

    FINISH:
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void IPGApproxSearcher::approxKnnSearchModelsubPaaMuAsSplit(IPGNode* root, TimeSeries* queryTs,int k,vector<PqItemSeries*>*heap){
    IPGNode *cur = root->routePaaMu(queryTs->paa), *parent = root;
    while (cur!= nullptr && !cur->isLeafNode()) {
        parent = cur;
        cur = cur->routePaaMu(queryTs->paa);
    }
    int pid = -1;
    layer = parent->layer+1;
    if(cur != nullptr)
        pid = cur->partitionId;
    else{
        double min_dist = numeric_limits<double>::max(), max_size = 0;
        IPGNode *node;
        for(auto &iter:*parent->children){
            if(iter.second == nullptr)  continue;
            double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, iter.second->sax, iter.second->bits_cardinality);
            if(dist < min_dist){
                min_dist = dist;
                pid = iter.second->partitionId;
                max_size = iter.second->size;
                node = iter.second;
            }else if(dist == min_dist && iter.second->size > max_size){
                max_size = iter.second->size;
                pid = iter.second->partitionId;
                node = iter.second;
            }
        }
        // we only concern whether the nearest node is a leaf or an internal node
        if(!node->isLeafNode()){
            approxKnnSearchModelsubPaaMuAsSplit(node, queryTs, k, heap);
            return;
        }
    }

//    if(pid == -1)   cur->exactSearchKnn(k, queryTs, *heap);
//    else{
    assert(pid != -1);
    for(auto &iter:*parent->children){
        if(iter.second!= nullptr && iter.second->partitionId == pid){
            iter.second->exactSearchKnn(k, queryTs, *heap);
        }
    }
//    }
}
