//
// Created by wzy on 2022/7/2.
//

#include <unordered_set>
#include <set>
#include "../../include/TAR/TARSearcher.h"
#include "../../include/Utils/SaxUtil.h"
#include "../../include/DataStructures/TimeSeries.h"

static int target_pack;

vector<PqItemSeries *> * TARSearcher::approxSearch(TARGNode *root, float *query, int k, const string &index_dir) {
    auto sax = SaxUtil::saxFromTs(query, Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
    unsigned query_invsax[Const::bitsCardinality];
    SaxUtil::invSaxFromSax(sax, query_invsax, Const::bitsCardinality, Const::segmentNum);
    string query_invsax_str = SaxUtil::invSax2String(query_invsax);

    auto*heap = new vector<PqItemSeries*>();

    int pid;
    TARGNode* target = root->route(query_invsax_str, &pid);
    if(pid == -1){
        pid = target->children.begin()->second->pid;
    }
    {
        string file_name = Const::tardisfn + to_string(pid);
        TARLNode* l_root = TARLNode::loadFromDisk(file_name);
        TARLNode* l_target = l_root->route2Leaf(query_invsax_str);
        target_pack = pid;
        l_target->search(k, query, *heap);
        l_root->deleteDescendants();
        delete l_root;
    }

    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

vector<PqItemSeries *> * TARSearcher::approxSearchDTW(TARGNode *root, float *query, int k, const string &index_dir, int window_size) {
    auto sax = SaxUtil::saxFromTs(query, Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
    unsigned query_invsax[Const::bitsCardinality];
    SaxUtil::invSaxFromSax(sax, query_invsax, Const::bitsCardinality, Const::segmentNum);
    string query_invsax_str = SaxUtil::invSax2String(query_invsax);

    auto*heap = new vector<PqItemSeries*>();

    int pid;
    TARGNode* target = root->route(query_invsax_str, &pid);
    if(pid == -1){
        pid = target->children.begin()->second->pid;
    }
    {
        string file_name = Const::tardisfn + to_string(pid);
        TARLNode* l_root = TARLNode::loadFromDisk(file_name);
        TARLNode* l_target = l_root->route2Leaf(query_invsax_str);
        l_target->search_dtw(k, query, *heap, window_size);
        l_root->deleteDescendants();
        delete l_root;
    }

    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

vector<PqItemSeries *> * TARSearcher::approxIncSearch(TARGNode *root, float *query, int k, const string &index_dir,
                                                      int *node_num) {
    auto sax = SaxUtil::saxFromTs(query, Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
    unsigned query_invsax[Const::bitsCardinality];
    SaxUtil::invSaxFromSax(sax, query_invsax, Const::bitsCardinality, Const::segmentNum);
    string query_invsax_str = SaxUtil::invSax2String(query_invsax);

    auto*heap = new vector<PqItemSeries*>();

    int pid;
    TARGNode* target = root->route(query_invsax_str, &pid);
    if(pid == -1){
        for(auto &iter: target->children){
            if(iter.second->pid == -1){
                approxIncSearchSub(iter.second, query, k, index_dir, node_num, heap, query_invsax_str);
                assert(*node_num < 1);
                sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
                return heap;
            }else
                pid = iter.second->pid;
            break;
        }
    }
    {
        string file_name = Const::tardisfn + to_string(pid);
        TARLNode* l_root = TARLNode::loadFromDisk(file_name);
        incSearchLocal(l_root, heap, k, query, node_num, query_invsax_str);
        l_root->deleteDescendants();
        delete l_root;
        if(*node_num > 0){
            TARGNode* parent = target->ancestor;
            for(auto &iter: parent->children){
                if(iter.second != target){
                    if(iter.second->pid == -1){
                        approxIncSearchSub(iter.second, query, k, index_dir, node_num, heap, query_invsax_str);
                    }else {
                        file_name = Const::tardisfn + to_string(iter.second->pid);
                        l_root = TARLNode::loadFromDisk(file_name);
                        incSearchLocal(l_root, heap, k, query, node_num, query_invsax_str);
                        l_root->deleteDescendants();
                        delete l_root;
                    }
                    if(*node_num < 1)    break;
                }
            }
        }
    }

    assert(*node_num < 1);
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void TARSearcher::approxIncSearchSub(TARGNode *root, float *query, int k, const string &index_dir,
                                  int *node_num,vector<PqItemSeries *> *heap, string& query_invsax_str){
    if(*node_num < 1)   return;
    int pid;
    TARGNode* target = root->route(query_invsax_str, &pid);
    if(pid == -1){
        for(auto &iter: target->children){
            if(iter.second->pid == -1){
                approxIncSearch(iter.second, query, k, index_dir, node_num);
            }else
                pid = iter.second->pid;
            break;
        }
    }
    {
        string file_name = Const::tardisfn + to_string(pid);
        TARLNode* l_root = TARLNode::loadFromDisk(file_name);
        incSearchLocal(l_root, heap, k, query, node_num, query_invsax_str);
        l_root->deleteDescendants();
        delete l_root;
        if(*node_num > 0){
            TARGNode* parent = target->ancestor;
            for(auto &iter: parent->children){
                if(iter.second != target){
                    if(iter.second->pid == -1){
                        approxIncSearch(iter.second, query, k, index_dir, node_num);
                    }else {
                        file_name = Const::tardisfn + to_string(iter.second->pid);
                        l_root = TARLNode::loadFromDisk(file_name);
                        incSearchLocal(l_root, heap, k, query, node_num, query_invsax_str);
                        l_root->deleteDescendants();
                        delete l_root;
                    }
                    if(*node_num < 1)    break;
                }
            }
        }
    }

}

struct PqItemTAR{
    TARGNode* node{};
    double dist{};

    PqItemTAR(TARGNode* n, double d){node = n; dist = d;}
    PqItemTAR(){node = nullptr; dist = 0;}

    bool operator <(const PqItemTAR & pt) const{
        if(node == pt.node) return false;
        else {
            if(dist != pt.dist) return dist < pt.dist;
            else return node < pt.node;
        }
    }
};

extern int LOADED_NODE_CNT;
vector<PqItemSeries *> * TARSearcher::exactSearch(TARGNode *root, float *query, int k, const string &index_dir){
    auto heap = approxSearch(root, query, k, index_dir);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);

    unordered_set<int>visited;
    set<PqItemTAR>pq;
    pq.insert(PqItemTAR(root, 0));

    PqItemTAR cur;
    while(!pq.empty()){
        cur = *pq.begin();
        if(cur.dist > bsf)  break;
        pq.erase(pq.begin());
        if(!cur.node->children.empty()){
            for(auto &iter:cur.node->children)
                // leaf node
                if(iter.second->children.empty()){
                    if(visited.contains(iter.second->pid) || iter.second->pid == target_pack)   continue;
                    pq.insert(PqItemTAR(iter.second, cur.dist));
                    visited.insert(iter.second->pid);
                } else{ // internal node
                    string invsax_str = iter.second->invSAX;
                    auto invsax = SaxUtil::str2Invsax(invsax_str);
                    auto sax = SaxUtil::invSax2Sax(invsax, iter.second->layer);
                    double lb_dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, sax, iter.second->layer);
                    if(lb_dist < bsf){
                        pq.insert(PqItemTAR(iter.second, lb_dist));
                    }
                    delete []sax;
                    delete[] invsax;
                }
        }else{
            string file_name = Const::tardisfn + to_string(cur.node->pid);
            TARLNode* l_root = TARLNode::loadFromDisk(file_name);
            exactSearchLocal(l_root, heap, k ,query);
            LOADED_NODE_CNT += l_root->getLeafNodeNbr();
            bsf = (*heap)[0]->dist;
            l_root->deleteDescendants();
            delete l_root;
        }
    }
    pq.clear();
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

vector<PqItemSeries *> * TARSearcher::exactSearchDTW(TARGNode *root, float *query, int k, const string &index_dir){
    auto heap = approxSearch(root, query, k, index_dir);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);
    auto* lowerLemire = new float [Const::tsLength];
    auto* upperLemire = new float [Const::tsLength];
    SaxUtil::lower_upper_lemire(query, Const::tsLength, Const::dtw_window_size, lowerLemire, upperLemire);
    double *lowerPaa = SaxUtil::paaFromTs(lowerLemire, Const::tsLengthPerSegment, Const::segmentNum);
    double *upperPaa = SaxUtil::paaFromTs(upperLemire, Const::tsLengthPerSegment, Const::segmentNum);


    unordered_set<int>visited;
    set<PqItemTAR>pq;
    pq.insert(PqItemTAR(root, 0));

    PqItemTAR cur;
    while(!pq.empty()){
        cur = *pq.begin();
        if(cur.dist > bsf)  break;
        pq.erase(pq.begin());
        if(!cur.node->children.empty()){
            for(auto &iter:cur.node->children)
                // leaf node
                if(iter.second->children.empty()){
                    if(visited.contains(iter.second->pid) || iter.second->pid == target_pack)   continue;
                    pq.insert(PqItemTAR(iter.second, cur.dist));
                    visited.insert(iter.second->pid);
                } else{ // internal node
                    string invsax_str = iter.second->invSAX;
                    auto invsax = SaxUtil::str2Invsax(invsax_str);
                    auto sax = SaxUtil::invSax2Sax(invsax, iter.second->layer);
                    double lb_dist = SaxUtil::minidist_paa_to_isax_DTW(upperPaa,lowerPaa, sax, iter.second->layer);
                    if(lb_dist < bsf){
                        pq.insert(PqItemTAR(iter.second, lb_dist));
                    }
                    delete []sax;
                    delete[] invsax;
                }
        }else{
            string file_name = Const::tardisfn + to_string(cur.node->pid);
            TARLNode* l_root = TARLNode::loadFromDisk(file_name);

            exactSearchLocalDTW(l_root, heap, k ,query, index_dir + to_string(cur.node->pid));
            LOADED_NODE_CNT += l_root->getLeafNodeNbr();
            bsf = (*heap)[0]->dist;

            l_root->deleteDescendants();
            delete l_root;
        }
    }
    pq.clear();
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}

void TARSearcher::exactSearchLocal(TARLNode *l_root, vector<PqItemSeries *> *heap, int k, float *query) {
    if(l_root == nullptr)   return;
    if(l_root->rcdNbr <= Const::th){
        l_root->search(k, query, *heap);
        return;
    }
    for(auto &iter:l_root->descendants){
        exactSearchLocal(iter.second, heap, k, query);
    }
}

void TARSearcher::exactSearchLocalDTW(TARLNode* l_root, vector<PqItemSeries*>*heap, int k, float *query, const string &index_dir){
    if(l_root == nullptr)   return;
    if(l_root->rcdNbr <= Const::th){
        l_root->search_dtw(k, query, *heap, Const::dtw_window_size);
        return;
    }
    for(auto &iter:l_root->descendants){
        exactSearchLocalDTW(iter.second, heap, k, query, index_dir);
    }
}

void TARSearcher::incSearchLocal(TARLNode *l_root, vector<PqItemSeries *> *heap, int k, float *query, int *node_num,
                                 const string &query_invsax_str) {
    if(l_root == nullptr)   return;
    if(*node_num < 1)    return;
    int cur_leaf_node_num = l_root->getLeafNodeNbr();
    if(cur_leaf_node_num <= *node_num) {
        exactSearchLocal(l_root, heap, k, query);
        *node_num = *node_num - cur_leaf_node_num;
        return;
    }
    if(cur_leaf_node_num > *node_num){
        string key = SaxUtil::invSaxHeadKFromInvSax(query_invsax_str, l_root->layer + 1);
        if(l_root->descendants.contains(key)){
            incSearchLocal(l_root->descendants[key], heap, k, query, node_num, query_invsax_str);
        }
        if(*node_num > 0){
            for(auto &iter:l_root->descendants){
                if(iter.first != key){
                    incSearchLocal(iter.second, heap, k, query, node_num, query_invsax_str);
                }
                if(*node_num < 1)   return;
            }
        }
    }
}