//
// Created by wzy on 2021/9/24.
//

#include "../../include/Searchers/TardisApproxSearch.h"
#include "../../include/Utils/MathUtil.h"
#include <algorithm>
#include <unordered_set>
#include <random>
#include <chrono>

static float *t_paa;
extern int layer;
extern int actual_node_number;

struct cand_tardis{
    TardisTreeNode* node;
    double lb;

    cand_tardis(){;}
    cand_tardis(TardisTreeNode*n, double l){
        node = n;
        lb = l;
    }
};

bool comp_tardis(const TardisTreeNode* x, const TardisTreeNode* y){
    if(x == nullptr)    return false;
    if(y == nullptr)    return true;
    return SaxUtil::LowerBound_Paa_iSax(t_paa, x->sax, x->layer) < SaxUtil::LowerBound_Paa_iSax(t_paa, y->sax, y->layer);
}

bool comp_tardis_cand(const cand_tardis& x, const cand_tardis& y){
    return x.lb < y.lb;
}

bool comp_tardis_hamming(const cand_tardis& x, const cand_tardis& y){
    return x.lb < y.lb;
}

bool comp_tardis_min_heap(const cand_tardis& x, const cand_tardis& y){
    return x.lb > y.lb;
}


vector<PqItemSeries *> * TardisApproxSearch::tardisApproxKnnSearch(TardisTreeNode* root, float* query, int k, int threshold){
    auto* queryTs = new TimeSeries(query);
    t_paa = queryTs->paa;
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    tardisApproxKnnSearchSub(root, sax, queryTs, k, threshold, heap);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void TardisApproxSearch::tardisApproxKnnSearchSub(TardisTreeNode* root, unsigned short *sax, TimeSeries* queryTs, int k, int threshold, vector<PqItemSeries*>*heap){
    if(root->isLeafNode())  return;
    TardisTreeNode *cur = root->route(sax), *parent = root;
    while (cur!= nullptr && !cur->isLeafNode() && cur->size > threshold) {
        parent = cur;
        cur = cur->route(sax);
    }
    int rest = threshold;

    if(cur != nullptr) {
        cur->exactSearchKnnInMemory(k, queryTs, *heap);
        rest -= cur->size;
    }

    int query_id = SaxUtil::invSaxHeadkFromSax(sax, Const::bitsCardinality, Const::segmentNum, parent->layer + 1);

    vector<cand_tardis>candidates;
    candidates.reserve((*parent->children).size());
    for(auto &iter:(*parent->children))
        if(iter.second != nullptr && iter.second != cur)
//            candidates.emplace_back(iter.second, MathUtil::bitDiffNum(iter.first, query_id, Const::segmentNum));
            candidates.emplace_back(iter.second, SaxUtil::LowerBound_Paa_iSax(t_paa, iter.second->sax, iter.second->layer));

//    make_heap(candidates.begin(), candidates.end(), comp_tardis_min_heap);
    sort(candidates.begin(), candidates.end(), comp_tardis_cand);

    for(int i=0;i<candidates.size() && rest > 0;++i){
//        cur = candidates.front().node;
        if(rest >= candidates[i].node->size){
            candidates[i].node->exactSearchKnnInMemory(k, queryTs, *heap);
            rest -= candidates[i].node->size;
        }
        else {
            if(candidates[i].node->isLeafNode())
                candidates[i].node->exactSearchKnnInMemory(k, queryTs, *heap, rest);
            else
                tardisApproxKnnSearchSub(candidates[i].node, sax, queryTs, k, rest, heap);
            return;
        }
//        pop_heap(candidates.begin(), candidates.end(), comp_tardis_min_heap);
//        candidates.pop_back();
    }
}

vector<PqItemSeries *> * TardisApproxSearch::tardisApproxKnnSearchGraph(TardisTreeNode* root, float* query, int k,
                                                                        int threshold, vector<vector<int>>*g){
    auto* queryTs = new TimeSeries(query);
    t_paa = queryTs->paa;
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    tardisApproxKnnSearchSubGraph(root, sax, queryTs, k, threshold, heap, g);
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}


void TardisApproxSearch::tardisApproxKnnSearchSubGraph(TardisTreeNode* root, unsigned short *sax, TimeSeries* queryTs,
                                                       int k, int &threshold, vector<PqItemSeries*>*heap, vector<vector<int>>*g){
    if(root->isLeafNode())  return;
    TardisTreeNode *cur = root->route(sax), *parent = root;
    while (cur!= nullptr && !cur->isLeafNode() && cur->size > threshold) {
        parent = cur;
        cur = cur->route(sax);
    }

//    bool flag[1 << Const::segmentNum];

    int rest = threshold;
    int query_id = SaxUtil::invSaxHeadkFromSax(sax, Const::bitsCardinality, Const::segmentNum, parent->layer + 1);
    if(cur != nullptr) {
        cur->exactSearchKnnInMemory(k, queryTs, *heap);
        rest -= cur->size;
    }
    TardisTreeNode * tar = cur;
//    flag[query_id] = true;

    int bound = Const::neighborNum;
    int i;
    for(i=0;i<bound && rest > 0;++i) {
        int nei =(*g)[query_id][i];
        cur = (*parent->children)[nei];
//        flag[(*g)[query_id][i]] = true;
        if(cur == nullptr || cur == tar)  continue;
        if (rest >= cur->size) {
            cur->exactSearchKnnInMemory(k, queryTs, *heap);
            rest -= cur->size;
        } else {
            if (cur->isLeafNode())
                cur->exactSearchKnnInMemory(k, queryTs, *heap, rest);
            else
                tardisApproxKnnSearchSubGraph(cur, sax, queryTs, k, rest, heap, g);
            break;
        }
    }
    actual_node_number = rest;

//    for(int i=0;i<(1<<Const::segmentNum) && rest > 0;++i){
//        if(!flag[i] &&  (*parent->children)[i] != nullptr){
//            cur = (*parent->children)[i];
//            if (rest >= cur->size) {
//                cur->exactSearchKnnInMemory(k, queryTs, *heap);
//                rest -= cur->size;
//            } else {
//                if (cur->isLeafNode())
//                    cur->exactSearchKnnInMemory(k, queryTs, *heap, rest);
//                else
//                    tardisApproxKnnSearchSubGraph(cur, sax, queryTs, k, rest, heap, g);
//                return;
//            }
//        }
//    }




}


vector<PqItemSeries *> * TardisApproxSearch::tardisApproxKnnSearchModel(TardisTreeNode *root, float *query, int k) {
    auto* queryTs = new TimeSeries(query);
    auto*heap = new vector<PqItemSeries*>();
    t_paa = queryTs->paa;

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    tardisApproxKnnSearchModelsub(root, queryTs, sax, k, heap);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

// deprecated: the fine-grained sorting is not viable
void TardisApproxSearch::tardisApproxKnnSearchModelsub(TardisTreeNode* root, TimeSeries* queryTs, unsigned short *sax, int k, vector<PqItemSeries*>*heap){
    TardisTreeNode *cur = root->route(sax), *parent = root;
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
        TardisTreeNode *node;
        for(auto &iter:*parent->children){
            if(iter.second == nullptr)  continue;
            double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, iter.second->sax, layer);
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
            tardisApproxKnnSearchModelsub(node, queryTs, sax, k, heap);
            return;
        }
    }

    for(auto &iter:*parent->children){
        if(iter.second!= nullptr && iter.second->partitionId == pid && iter.second->isLeafNode()){
            iter.second->exactSearchKnnLeaf(k, queryTs, *heap);
        }
    }
}

vector<PqItemSeries *> * TardisApproxSearch::tardisApproxKnnIncSearchModel(TardisTreeNode* root, float* query, int k, int node_num){
    auto* queryTs = new TimeSeries(query);
    t_paa = queryTs->paa;
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    tardisApproxKnnIncSearchModelSub(root, sax, queryTs, k, node_num, heap);

    actual_node_number = node_num;
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void TardisApproxSearch::tardisApproxKnnIncSearchModelSub(TardisTreeNode* root, unsigned short *sax, TimeSeries* queryTs, int k, int& node_num, vector<PqItemSeries*>*heap){
    if(root->isLeafNode() || node_num <= 0)  return;
    TardisTreeNode *cur = root->route(sax), *parent = root;
    while (cur!= nullptr && !cur->isLeafNode() && cur->getLeafNodeNumber(cur) > node_num) {
        parent = cur;
        cur = cur->route(sax);
    }

    if(cur != nullptr) {
        if(cur->isLeafNode()){
            for(auto &iter:(*parent->children)){
                if(iter.second!= nullptr && iter.second->partitionId == cur->partitionId){
                    iter.second->exactSearchKnnLeaf(k, queryTs, *heap);
                }
            }
            node_num--;
        }else{
            cur->exactSearchKnn(k, queryTs, *heap, node_num);
        }
    }

    if(node_num <= 0)   return;

    int max_pid = 0;
    for(auto &iter:(*parent->children)){
        if(iter.second!= nullptr)
            max_pid = max(max_pid, iter.second->partitionId);
    }

//    vector<int>candidates;
//    for(int i=0;i<=max_pid;++i)
//        if(cur!= nullptr && i == cur->partitionId)  continue;
//        else    candidates.push_back(i);
//    for(auto &iter:(*parent->children)){
//        if(iter.second!= nullptr && iter.second->partitionId == -1){
//            candidates.push_back(-iter.second->id);
//        }
//    }
//
//    auto seed = chrono::system_clock::now().time_since_epoch().count();
//    shuffle(candidates.begin(), candidates.end(), default_random_engine(seed));
//
//    for(int i=0;i<candidates.size() && node_num > 0;++i){
//        if(candidates[i] < 0){
//            assert((*parent->children)[-candidates[i]]->id == -candidates[i]);
//            tardisApproxKnnIncSearchModelSub((*parent->children)[-candidates[i]], sax, queryTs, k, node_num, heap);
//        }else{
//            for(auto &iter:(*parent->children)){
//                if(iter.second!= nullptr && iter.second->partitionId == candidates[i]){
//                    iter.second->exactSearchKnnLeaf(k, queryTs, *heap);
//                }
//            }
//            node_num--;
//        }
//    }



    for(int i=0;i<max_pid && node_num > 0;++i){
        if(cur!= nullptr && i == cur->partitionId)   continue;
        for(auto &iter:(*parent->children)){
            if(iter.second!= nullptr && iter.second->partitionId == i){
                iter.second->exactSearchKnnLeaf(k, queryTs, *heap);
            }
        }
        --node_num;
    }

    for(auto &iter:(*parent->children)){
        if(iter.second!= nullptr && iter.second->partitionId == -1){
            tardisApproxKnnIncSearchModelSub(iter.second, sax, queryTs, k, node_num, heap);
        }
        if(node_num <=0)    return;
    }
}

