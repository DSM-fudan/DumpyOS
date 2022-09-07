//
// Created by pengwang5 on 2022/1/16.
//

#include <set>
#include <unordered_set>
#include <chrono>
#include <cmath>
#include "../../include/Searchers/FADASSearcher.h"
#include "../../include/DataStructures/PqItemSeries.h"
#include "../../include/DataStructures/TimeSeries.h"
#include "../../include/Utils/SaxUtil.h"
#include "../../include/Utils/MathUtil.h"
#include "../../include/Const.h"

extern long LB_SERIES_TIME = 0, HEAP_TIME = 0,
    LB_NODE_TIME_STAT = 0, LB_NODE_CNT = 0, LOADED_NODE_CNT = 0;
extern double DIST_CALC_TIME = 0, READ_TIME = 0;

static FADASNode* targetNode;
vector<PqItemSeries *> * FADASSearcher::approxSearch(FADASNode *root, float *query, int k, vector<vector<int>> *g,
                                                     const string &index_dir) {
    auto* queryTs = new TimeSeries(query);
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    int head = SaxUtil::invSaxHeadFromSax(sax, Const::bitsCardinality, Const::segmentNum);
    FADASNode *cur = (root->children)[head];
    if(cur == nullptr){
        FADASNode *node = nullptr;
        for(int i=0;i<FADASNode::a + FADASNode::b + FADASNode::c;++i){
            if(root->children[(*g)[head][i]] != nullptr){
                node = root->children[(*g)[head][i]];
                break;
            }
        }
        assert(node!=nullptr);
        // we only concern whether the nearest node is a leaf or an internal node
        if(node->isInternalNode()){
            approxSearchInterNode(node, queryTs, sax, k, heap, index_dir);
        }else { node->search(k, queryTs, *heap, index_dir, nullptr, nullptr); targetNode = node;}
    }else if(!cur->isInternalNode()){
        { cur->search(k, queryTs, *heap, index_dir, nullptr, nullptr); targetNode = cur;}
    }else approxSearchInterNode(cur, queryTs, sax, k, heap, index_dir);
    
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void FADASSearcher::approxSearchInterNode(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                          vector<PqItemSeries *> *heap, const string &index_dir) {
    FADASNode *cur = root->route(sax);
    if(!cur->isInternalNode()){
        cur->search(k, queryTs, *heap, index_dir, nullptr, nullptr);
        targetNode = cur;
        return;
    }

    // below is only for a nullptr target leaf node, then we search the nearest sibling
    double min_dist = numeric_limits<double>::max(), max_size = 0;
    FADASNode *node;
    for(int i=0;i<cur->children.size();++i){
        if(cur->children[i] == nullptr)  continue;
        double dist;
        if(!cur->children[i]->isInternalNode())
            dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, cur->sax, cur->bits_cardinality, cur->children[i]->chosenSegments, i);
//        if(cur->children[i]->isLeafNode())  dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, cur->children[i]->sax, cur->children[i]->bits_cardinality);
        else dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, cur->children[i]->sax, cur->children[i]->bits_cardinality);
        if(dist < min_dist){
            min_dist = dist;
            max_size = cur->children[i]->size;
            node = cur->children[i];
        }else if(dist == min_dist && cur->children[i]->size > max_size){
            max_size = cur->children[i]->size;
            node = cur->children[i];
        }
    }

    // we only concern whether the nearest node is a leaf or an internal node
    if(node->isInternalNode()){
        approxSearchInterNode(node, queryTs, sax, k, heap, index_dir);
        return;
    }else { node->search(k, queryTs, *heap, index_dir, nullptr, nullptr); targetNode = node;}
}

vector<PqItemSeries *> * FADASSearcher::approxSearchDTW(FADASNode *root, float *query, int k, vector<vector<int>> *g,
                                                     const string &index_dir) {
    auto* queryTs = new TimeSeries(query);
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    int head = SaxUtil::invSaxHeadFromSax(sax, Const::bitsCardinality, Const::segmentNum);
    FADASNode *cur = (root->children)[head];
    if(cur == nullptr){
        FADASNode *node = nullptr;
        for(int i=0;i<FADASNode::a + FADASNode::b + FADASNode::c;++i){
            if(root->children[(*g)[head][i]] != nullptr){
                node = root->children[(*g)[head][i]];
                break;
            }
        }
        assert(node!=nullptr);
        // we only concern whether the nearest node is a leaf or an internal node
        if(node->isInternalNode()){
            approxSearchInterNodeDTW(node, queryTs, sax, k, heap, index_dir);
        }else { node->searchDTW(k, queryTs, *heap, index_dir); targetNode = node;}
    }else if(!cur->isInternalNode()){
        { cur->searchDTW(k, queryTs, *heap, index_dir); targetNode = cur;}
    }else approxSearchInterNodeDTW(cur, queryTs, sax, k, heap, index_dir);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void FADASSearcher::approxSearchInterNodeDTW(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                          vector<PqItemSeries *> *heap, const string &index_dir) {
    FADASNode *cur = root->route(sax);
    if(!cur->isInternalNode()){
        cur->searchDTW(k, queryTs, *heap, index_dir);
        targetNode = cur;
        return;
    }

    auto* lowerLemire = new float [Const::tsLength];
    auto* upperLemire = new float [Const::tsLength];
    SaxUtil::lower_upper_lemire(queryTs->ts, Const::tsLength, Const::dtw_window_size, lowerLemire, upperLemire);
    double *lowerPaa = SaxUtil::paaFromTs(lowerLemire, Const::tsLengthPerSegment, Const::segmentNum);
    double *upperPaa = SaxUtil::paaFromTs(upperLemire, Const::tsLengthPerSegment, Const::segmentNum);


    // below is only for a nullptr target leaf node, then we search the nearest sibling
    double min_dist = numeric_limits<double>::max(), max_size = 0;
    FADASNode *node;
    for(auto & i : cur->children){
        if(i == nullptr)  continue;
        double dist = SaxUtil::minidist_paa_to_isax_DTW(upperPaa, lowerPaa, i->sax, i->bits_cardinality);
        if(dist < min_dist){
            min_dist = dist;
            max_size = i->size;
            node = i;
        }else if(dist == min_dist && i->size > max_size){
            max_size = i->size;
            node = i;
        }
    }

    // we only concern whether the nearest node is a leaf or an internal node
    if(node->isInternalNode()){
        approxSearchInterNodeDTW(node, queryTs, sax, k, heap, index_dir);
        return;
    }else { node->searchDTW(k, queryTs, *heap, index_dir); targetNode = node;}
}

vector<PqItemSeries *> * FADASSearcher::approxSearchLessPack(FADASNode *root, float *query, int k, vector<vector<int>> *g,
                                                     const string &index_dir) {
    auto* queryTs = new TimeSeries(query);
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    int head = SaxUtil::invSaxHeadFromSax(sax, Const::bitsCardinality, Const::segmentNum);
    FADASNode *cur = (root->children)[head];
    if(cur == nullptr){
        FADASNode *node = nullptr;
        for(int i=0;i<FADASNode::a + FADASNode::b + FADASNode::c;++i){
            if(root->children[(*g)[head][i]] != nullptr){
                node = root->children[(*g)[head][i]];
                break;
            }
        }
        assert(node!=nullptr);
        // we only concern whether the nearest node is a leaf or an internal node
        if(node->isInternalNode()){
            approxSearchInterNodeLessPack(node, queryTs, sax, k, heap, index_dir);
        }else { node->searchLessPack(k, queryTs, *heap, index_dir); targetNode = node;}
    }else if(!cur->isInternalNode()){
        { cur->searchLessPack(k, queryTs, *heap, index_dir); targetNode = cur;}
    }else approxSearchInterNodeLessPack(cur, queryTs, sax, k, heap, index_dir);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void FADASSearcher::approxSearchInterNodeLessPack(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                          vector<PqItemSeries *> *heap, const string &index_dir) {
    FADASNode *cur = root->route(sax);
    if(!cur->isInternalNode()){
        cur->searchLessPack(k, queryTs, *heap, index_dir);
        targetNode = cur;
        return;
    }

    // below is only for a nullptr target leaf node, then we search the nearest sibling
    cout <<"nullptr target leaf node" <<endl;
    double min_dist = numeric_limits<double>::max(), max_size = 0;
    FADASNode *node;
    for(int i=0;i<cur->children.size();++i){
        if(cur->children[i] == nullptr)  continue;
        double dist;
        if(!cur->children[i]->isInternalNode())
            dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, cur->sax, cur->bits_cardinality, cur->children[i]->chosenSegments, i);
//        if(cur->children[i]->isLeafNode())  dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, cur->children[i]->sax, cur->children[i]->bits_cardinality);
        else dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, cur->children[i]->sax, cur->children[i]->bits_cardinality);
        if(dist < min_dist){
            min_dist = dist;
            max_size = cur->children[i]->size;
            node = cur->children[i];
        }else if(dist == min_dist && cur->children[i]->size > max_size){
            max_size = cur->children[i]->size;
            node = cur->children[i];
        }
    }

    // we only concern whether the nearest node is a leaf or an internal node
    if(node->isInternalNode()){
        approxSearchInterNode(node, queryTs, sax, k, heap, index_dir);
        return;
    }else { node->search(k, queryTs, *heap, index_dir, nullptr, nullptr); targetNode = node;}
}

vector<PqItemSeries *> * FADASSearcher::approxSearchPos(FADASNode *root, float *query, int k, vector<vector<int>> *g,
                                                     const string &index_dir) {
    auto* queryTs = new TimeSeries(query);
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    int head = SaxUtil::invSaxHeadFromSax(sax, Const::bitsCardinality, Const::segmentNum);
    FADASNode *cur = (root->children)[head];
    if(cur == nullptr){
        FADASNode *node = nullptr;
        for(int i=0;i<FADASNode::a + FADASNode::b + FADASNode::c;++i){
            if(root->children[(*g)[head][i]] != nullptr){
                node = root->children[(*g)[head][i]];
                break;
            }
        }
        assert(node!=nullptr);
        // we only concern whether the nearest node is a leaf or an internal node
        if(node->isInternalNode()){
            approxSearchInterNodePos(node, queryTs, sax, k, heap, index_dir);
        }else { node->search_offset(k, queryTs, *heap, index_dir); targetNode = node;}
    }else if(!cur->isInternalNode()){
        { cur->search_offset(k, queryTs, *heap, index_dir); targetNode = cur;}
    }else approxSearchInterNodePos(cur, queryTs, sax, k, heap, index_dir);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void FADASSearcher::approxSearchInterNodePos(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                          vector<PqItemSeries *> *heap, const string &index_dir) {
    FADASNode *cur = root->route(sax);
    if(!cur->isInternalNode()){
        cur->search_offset(k, queryTs, *heap, index_dir);
        targetNode = cur;
        return;
    }

    double min_dist = numeric_limits<double>::max(), max_size = 0;
    FADASNode *node;
    for(int i=0;i<cur->children.size();++i){
        if(cur->children[i] == nullptr)  continue;
        double dist;
        if(!cur->children[i]->isInternalNode())  dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, cur->sax, cur->bits_cardinality, cur->children[i]->chosenSegments, i);
        else dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, cur->children[i]->sax, cur->children[i]->bits_cardinality);
        if(dist < min_dist){
            min_dist = dist;
            max_size = cur->children[i]->size;
            node = cur->children[i];
        }else if(dist == min_dist && cur->children[i]->size > max_size){
            max_size = cur->children[i]->size;
            node = cur->children[i];
        }
    }

    // we only concern whether the nearest node is a leaf or an internal node
    if(node->isInternalNode()){
        approxSearchInterNodePos(node, queryTs, sax, k, heap, index_dir);
        return;
    }else { node->search_offset(k, queryTs, *heap, index_dir); targetNode = node;}
}

struct PqItemFadas{
    FADASNode* node{};
    double dist{};

    PqItemFadas(FADASNode* n, double d){node = n; dist = d;}
    PqItemFadas(){node = nullptr; dist = 0;}

    bool operator <(const PqItemFadas & pt) const{
        if(node == pt.node) return false;
        else {
            if(dist != pt.dist) return dist < pt.dist;
            if(node->isInternalNode())  return true;
            return false;
        }
    }
};

struct PqItemFadasId{
    FADASNode* parent;
    int id;
    double dist{};

    PqItemFadasId(FADASNode* _parent, int _id, double d){parent = _parent;id= _id; dist = d;}
    PqItemFadasId(){id=-1; dist = 0;parent = nullptr;}

    bool operator <(const PqItemFadasId & pt) const{
        if(dist != pt.dist)
            return dist < pt.dist;
        if(parent->layer != pt.parent->layer)
            return parent->layer < pt.parent->layer;
        return false;
    }
};

vector<PqItemSeries*>*FADASSearcher::exactSearch(FADASNode* root, float *query, int k, vector<vector<int>> *g){
    vector<PqItemSeries*>* heap = approxSearch(root, query, k, g, Const::fidxfn);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);

    set<PqItemFadas>pq;
    pq.insert(PqItemFadas(root, 0));

    PqItemFadas cur;
    while(!pq.empty()){
        cur = *pq.begin();
        if(cur.dist > bsf)  break;
        pq.erase(pq.begin());
        if(cur.node->isInternalNode()){
            unordered_set<FADASNode*>inserted;
            for(FADASNode* node:cur.node->children)
                if(node != nullptr && node != targetNode && inserted.find(node) == inserted.end()) {
//                    auto start = chrono::system_clock::now();
                    inserted.insert(node);
                    double lb_dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax, node->bits_cardinality);
//                    auto end = chrono::system_clock::now();
//                    LB_NODE_TIME_STAT += chrono::duration_cast<chrono::microseconds>(end - start).count();
                    if(lb_dist < bsf){
                        pq.insert(PqItemFadas(node, lb_dist));
                        LB_NODE_CNT++;
                    }
                }
            inserted.clear();
        }else{
            cur.node->search(k, queryTs, *heap, Const::fidxfn, nullptr, nullptr);
            ++LOADED_NODE_CNT;
            bsf = (*heap)[0]->dist;
        }
    }
    pq.clear();
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

vector<PqItemSeries*>*FADASSearcher::exactSearchDTW(FADASNode* root, float *query, int k, vector<vector<int>> *g){
    vector<PqItemSeries*>* heap = approxSearchDTW(root, query, k, g, Const::fidxfn);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);
    auto* lowerLemire = new float [Const::tsLength];
    auto* upperLemire = new float [Const::tsLength];
    SaxUtil::lower_upper_lemire(query, Const::tsLength, Const::dtw_window_size, lowerLemire, upperLemire);
    double *lowerPaa = SaxUtil::paaFromTs(lowerLemire, Const::tsLengthPerSegment, Const::segmentNum);
    double *upperPaa = SaxUtil::paaFromTs(upperLemire, Const::tsLengthPerSegment, Const::segmentNum);

    set<PqItemFadas>pq;
    pq.insert(PqItemFadas(root, 0));

    PqItemFadas cur;
    while(!pq.empty()){
        cur = *pq.begin();
        if(cur.dist > bsf)  break;
        pq.erase(pq.begin());
        if(cur.node->isInternalNode()){
            unordered_set<FADASNode*>inserted;
            for(FADASNode* node:cur.node->children)
                if(node != nullptr && node != targetNode && inserted.find(node) == inserted.end()) {
//                    auto start = chrono::system_clock::now();
                    inserted.insert(node);
                    double lb_dist = SaxUtil::minidist_paa_to_isax_DTW(upperPaa, lowerPaa, node->sax, node->bits_cardinality);
//                    auto end = chrono::system_clock::now();
//                    LB_NODE_TIME_STAT += chrono::duration_cast<chrono::microseconds>(end - start).count();
                    if(lb_dist < bsf){
                        pq.insert(PqItemFadas(node, lb_dist));
                        LB_NODE_CNT++;
                    }
                }
            inserted.clear();
        }else{
            cur.node->searchDTW(k, queryTs, *heap, Const::fidxfn);
            ++LOADED_NODE_CNT;
            bsf = (*heap)[0]->dist;
        }
    }
    pq.clear();
    delete queryTs;
    delete[] lowerPaa;
    delete[] upperPaa;
    delete[] lowerLemire;
    delete[] upperLemire;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

vector<PqItemSeries*>*FADASSearcher::exactSearchPos(FADASNode* root, float *query, int k, vector<vector<int>> *g){
    vector<PqItemSeries*>* heap = approxSearchPos(root, query, k, g, Const::posidxfn);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);

    set<PqItemFadas>pq;
    pq.insert(PqItemFadas(root, 0));

    PqItemFadas cur;
    while(!pq.empty()){
        cur = *pq.begin();
        if(cur.dist >= bsf)  break;
        pq.erase(pq.begin());
        if(cur.node->isInternalNode()){
            for(FADASNode* node:cur.node->children)
                if(node != nullptr && node != targetNode) {
//                    auto start = chrono::system_clock::now();
                    double lb_dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax, node->bits_cardinality);
//                    auto end = chrono::system_clock::now();
//                    LB_NODE_TIME_STAT += chrono::duration_cast<chrono::microseconds>(end - start).count();
                    if(lb_dist < bsf) {
                        pq.insert(PqItemFadas(node, lb_dist));
                        LB_NODE_CNT++;
                    }
                }
        }else{
            cur.node->search_offset(k, queryTs, *heap, Const::posidxfn);
            ++LOADED_NODE_CNT;
            bsf = (*heap)[0]->dist;
        }
    }
    pq.clear();
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

static float *t_paa;
static double *low_paa, *up_paa;
bool comp_fadas(const FADASNode* x, const FADASNode* y){
    if(x == nullptr)    return false;
    if(y == nullptr)    return true;
    return SaxUtil::LowerBound_Paa_iSax(t_paa, x->sax, x->layer) < SaxUtil::LowerBound_Paa_iSax(t_paa, y->sax, y->layer);
}

bool comp_fadas_dtw(const FADASNode* x, const FADASNode* y){
    if(x == nullptr)    return false;
    if(y == nullptr)    return true;
    return SaxUtil::minidist_paa_to_isax_DTW(up_paa, low_paa, x->sax, x->bits_cardinality) < SaxUtil::minidist_paa_to_isax_DTW(up_paa, low_paa, y->sax, y->bits_cardinality);
}


vector<PqItemSeries *> * FADASSearcher::approxIncSearch(FADASNode *root, float *query, int k, const string &index_dir,
                                                        int node_num, float * query_reordered, int*ordering) {
    auto* queryTs = new TimeSeries(query);
    t_paa = queryTs->paa;
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    approxIncSearchInterNode(root, queryTs, sax, k, heap, index_dir, node_num, query_reordered, ordering);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}

double getMinDist1stLayer(const float *paa, int id){
    double ret = 0;
    for(int i=Const::segmentNum-1;i>=0;--i){
        if(id %2 == 1){
            if(paa[i] <0)
                ret += (-paa[i]);
        }else{
            if(paa[i] > 0)
                ret += (paa[i]);
        }
        id >>=1;
    }
    return ret;
}

vector<PqItemSeries *> *FADASSearcher::ngSearch(FADASNode *root, float *query, float *query_reordered, int *ordering,
                                                            int k,
                                                            vector<vector<int>> *g, int nprobes){
    vector<PqItemSeries*>* heap = approxSearch(root, query, k, g, Const::fidxfn);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);

    set<PqItemFadasId>pq;
    for(int i =0;i<Const::vertexNum;++i){
        if(root->children[i] == nullptr)    continue;
        double dist  = getMinDist1stLayer(queryTs->paa , i);
        pq.insert(PqItemFadasId(root,i, dist));
    }
    unordered_set<FADASNode*>visited;
    int cur_probe = 1;
    while(!pq.empty() && cur_probe < nprobes){
        double top_dist;
        FADASNode* node;
        top_dist = pq.begin()->dist;
        if(top_dist > bsf)  break;
        node = pq.begin()->parent->children[pq.begin()->id];
        pq.erase(pq.begin());
        if(visited.count(node) > 0) continue;
        visited.insert(node);

        if(node->isInternalNode()){
            int len = (1 << (node->chosenSegments.size()));
            for(int i =0;i<len;++i){
                if(node->children[i] == nullptr || node->children[i] == targetNode)    continue;
                double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax, node->bits_cardinality, node->chosenSegments, i);
                if(dist < bsf){
                    pq.insert(PqItemFadasId(node, i, dist));
                }
            }
        }else{
            node->search(k, queryTs, *heap, Const::fidxfn, query_reordered, ordering);
            ++cur_probe;
            bsf = (*heap)[0]->dist;
        }
    }
    pq.clear();
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}


vector<PqItemSeries *> *FADASSearcher::ngSearchIdLevelNaive(FADASNode *root, float *query, float *query_reordered, int *ordering,
                                                int k,
                                                vector<vector<int>> *g, int nprobes){
    vector<PqItemSeries*>* heap = approxSearch(root, query, k, g, Const::fidxfn);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);

    set<PqItemFadasId>pq;
    for(int i =0;i<Const::vertexNum;++i){
        if(root->children[i] == nullptr)    continue;
        double dist  = getMinDist1stLayer(queryTs->paa , i);
        pq.insert(PqItemFadasId(root,i, dist));
    }
    unordered_set<FADASNode*>visited;
    int cur_probe = 1;
    while(!pq.empty() && cur_probe < nprobes){
        double top_dist;
        FADASNode* node;
        top_dist = pq.begin()->dist;
        if(top_dist > bsf)  break;
        node = pq.begin()->parent->children[pq.begin()->id];
        pq.erase(pq.begin());
        if(visited.count(node) > 0) continue;
        visited.insert(node);

        if(node->isInternalNode()){
            int len = (1 << (node->chosenSegments.size()));
            for(int i =0;i<len;++i){
                if(node->children[i] == nullptr || node->children[i] == targetNode)    continue;
                double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax, node->bits_cardinality, node->chosenSegments, i);
                if(dist < bsf){
                    pq.insert(PqItemFadasId(node, i, dist));
                }
            }
        }else{
            node->search(k, queryTs, *heap, Const::fidxfn, query_reordered, ordering);
            ++cur_probe;
            bsf = (*heap)[0]->dist;
        }
    }
    pq.clear();
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}

vector<PqItemSeries *> *FADASSearcher::ngSearchNaive(FADASNode *root, float *query, float *query_reordered, int *ordering,
                                                int k,
                                                vector<vector<int>> *g, int nprobes) {
    vector<PqItemSeries*>* heap = approxSearch(root, query, k, g, Const::fidxfn);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);

    set<PqItemFadas>pq;
    pq.insert(PqItemFadas(root, 0));

    PqItemFadas cur;
    int cur_probe = 1;
    while(!pq.empty() && cur_probe < nprobes){
        cur = *pq.begin();
        if(cur.dist > bsf)  break;
        pq.erase(pq.begin());
        if(cur.node->isInternalNode()){
            unordered_set<FADASNode*>inserted;
            for(FADASNode* node:cur.node->children)
                if(node != nullptr && node != targetNode && inserted.find(node) == inserted.end()) {
//                    auto start = chrono::system_clock::now();
                    inserted.insert(node);
                    double lb_dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax, node->bits_cardinality);
//                    auto end = chrono::system_clock::now();
//                    LB_NODE_TIME_STAT += chrono::duration_cast<chrono::microseconds>(end - start).count();
                    if(lb_dist < bsf){
                        pq.insert(PqItemFadas(node, lb_dist));
                        LB_NODE_CNT++;
                    }
                }
            inserted.clear();
        }else{
            cur.node->search(k, queryTs, *heap, Const::fidxfn, query_reordered, ordering);
//            ++LOADED_NODE_CNT;
            ++cur_probe;
            bsf = (*heap)[0]->dist;
        }
    }
    pq.clear();
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}



vector<PqItemSeries *> * FADASSearcher::approxIncSearchDTW(FADASNode *root, float *query, int k, const string &index_dir,
                                                        int node_num) {
    auto* queryTs = new TimeSeries(query);
    auto* lowerLemire = new float [Const::tsLength];
    auto* upperLemire = new float [Const::tsLength];
    SaxUtil::lower_upper_lemire(query, Const::tsLength, Const::dtw_window_size, lowerLemire, upperLemire);
    low_paa= SaxUtil::paaFromTs(lowerLemire, Const::tsLengthPerSegment, Const::segmentNum);
    up_paa = SaxUtil::paaFromTs(upperLemire, Const::tsLengthPerSegment, Const::segmentNum);

    t_paa = queryTs->paa;
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    approxIncSearchInterNodeDTW(root, queryTs, sax, k, heap, index_dir, node_num);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}

void FADASSearcher::approxIncSearchInterNodeDTW(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                             vector<PqItemSeries *> *heap, const string &index_dir,int &node_num) {
    if(!root->isInternalNode() || node_num <= 0)  return;
    FADASNode *cur = root->route1step(sax), *parent = root;
    while (cur!= nullptr && cur->isInternalNode() && cur->getLeafNodeNum() > node_num) {
        parent = cur;
        cur = cur->route(sax);
    }

    if(cur!= nullptr){
        if(!cur->isInternalNode()){
            cur->searchDTW(k, queryTs, *heap, index_dir);
            --node_num;
        }else{
            approxIncSearchInterNodeDTW(cur, queryTs, sax, k, heap, index_dir, node_num);
        }
    }

    if(node_num <=0)    return;

    vector<FADASNode*>candidates;
    unordered_set<FADASNode*>cands;
    for(FADASNode *node: parent->children)
        if(node != nullptr && node!=cur && cands.find(node) == cands.end()) {
            candidates.push_back(node);
            cands.insert(node);
        }
    cands.clear();
    sort(candidates.begin(), candidates.end(), comp_fadas_dtw);



    for(int i=0;i<candidates.size() && node_num > 0;++i){
        if(!candidates[i]->isInternalNode()) {
            candidates[i]->searchDTW(k, queryTs, *heap, index_dir);
            --node_num;
        }
        else {
            approxIncSearchInterNodeDTW(candidates[i], queryTs, sax, k, heap, index_dir, node_num);
        }
    }

}

int q_id;
bool comp_fadas_id(const int i, const int j){
    return MathUtil::bitDiffNum(q_id, i ) < MathUtil::bitDiffNum(q_id,j);
}

void FADASSearcher::approxIncSearchInterNode(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                             vector<PqItemSeries *> *heap, const string &index_dir,
                                             int &node_num, float *query_reordered, int *ordering) {
    if(!root->isInternalNode() || node_num <= 0)  return;
    FADASNode *cur = root->route1step(sax), *parent = root;
    while (cur!= nullptr && cur->isInternalNode() && cur->leaf_num > node_num) {
        parent = cur;
        cur = cur->route(sax);
    }

    if(cur!= nullptr){
        if(!cur->isInternalNode()){
            cur->search(k, queryTs, *heap, index_dir, query_reordered, ordering);
            --node_num;
        }else{
            approxIncSearchInterNode(cur, queryTs, sax, k, heap, index_dir, node_num, query_reordered, ordering);
        }
    }

    if(node_num <=0)    return;

//    double bsf = (*heap)[0]->dist;
//    vector<PqItemFadasId>candidates;
//    int len = (1 << (parent->chosenSegments.size()));
//    for(int i =0;i<len;++i){
//        if(parent->children[i] == nullptr || parent->children[i] == cur)    continue;
//        double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, parent->sax, parent->bits_cardinality, parent->chosenSegments, i);
//        if(dist < bsf){
//            candidates.emplace_back(parent, i, dist);
//        }
//    }
//    sort(candidates.begin(),  candidates.end());

    vector<PqItemFadas>candidates;
    unordered_set<FADASNode*>cands;
    for(FADASNode *node: parent->children)
        if(node != nullptr && node!=cur && cands.find(node) == cands.end()) {
            double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax, node->layer);
            cands.insert(node);
            candidates.emplace_back(node , dist);
        }
    cands.clear();
    sort(candidates.begin(), candidates.end(), comp_fadas);

//    q_id = SaxUtil::extendSax(sax, parent->bits_cardinality, parent->chosenSegments);
//    vector<int>ids(1 << parent->chosenSegments.size(), 0);
//    for(int i=0;i<(1<<parent->chosenSegments.size());++i)
//        ids[i] = i;
//    sort(ids.begin(),  ids.end(), comp_fadas_id);
//
//    unordered_set<FADASNode*>visited;
//    for(int i=0;i<ids.size() && node_num > 0;++i){
//        FADASNode *tmp = parent->children[ids[i]];
//        if(tmp != nullptr && tmp != cur && visited.find(tmp) == visited.end()){
//            if(tmp->isLeafNode()){
//                --node_num;
//                tmp->search(k, queryTs, *heap, index_dir);
//                visited.insert(tmp);
//            }else{
//                approxIncSearchInterNode(tmp, queryTs, sax, k, heap, index_dir, node_num);
//                visited.insert(tmp);
//            }
//        }
//    }

//    unordered_set<FADASNode*>visited;
//    for(int i=0;i<candidates.size() && node_num > 0;++i){
//        FADASNode* node = parent->children[candidates[i].id];
//        if(visited.count(node) > 0) continue;
//        visited.insert(node);
//        if(!node->isInternalNode()) {
//            node->search(k, queryTs, *heap, index_dir, query_reordered, ordering);
//            --node_num;
//        }
//        else {
//            approxIncSearchInterNode(node, queryTs, sax, k, heap, index_dir, node_num, query_reordered, ordering);
//        }
//    }

    for(int i=0;i<candidates.size() && node_num > 0;++i){
        FADASNode* node = candidates[i];
        if(!node->isInternalNode()) {
            node->search(k, queryTs, *heap, index_dir, query_reordered, ordering);
            --node_num;
        }
        else {
            approxIncSearchInterNode(node, queryTs, sax, k, heap, index_dir, node_num, query_reordered, ordering);
        }
    }


}


vector<PqItemSeries *> * FADASSearcher::approxIncSearchFuzzy(FADASNode *root, float *query, int k, const string &index_dir,
                                                        int node_num) {
    auto* queryTs = new TimeSeries(query);
    t_paa = queryTs->paa;
    auto*heap = new vector<PqItemSeries*>();
    auto*hash_set = new unordered_set<float*, createhash, isEqual>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    approxIncSearchInterNodeFuzzy(root, queryTs, sax, k, heap, index_dir, node_num, hash_set);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}


void FADASSearcher::approxIncSearchInterNodeFuzzy(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                             vector<PqItemSeries *> *heap, const string &index_dir,int &node_num,
                                             unordered_set<float*, createhash, isEqual>*hash_set) {
    if(!root->isInternalNode() || node_num <= 0)  return;
    FADASNode *cur = root->route1step(sax), *parent = root;
    while (cur!= nullptr && cur->isInternalNode() && cur->getLeafNodeNum() > node_num) {
        parent = cur;
        cur = cur->route(sax);
    }

    if(cur!= nullptr){
        if(!cur->isInternalNode()){
            cur->search(k, queryTs, *heap, index_dir, hash_set);
            --node_num;
        }else{
            approxIncSearchInterNodeFuzzy(cur, queryTs, sax, k, heap, index_dir, node_num, hash_set);
        }
    }

    if(node_num <=0)    return;

    vector<FADASNode*>candidates;
    unordered_set<FADASNode*>cands;
    for(FADASNode *node: parent->children)
        if(node != nullptr && cands.find(node) == cands.end()) {
            candidates.push_back(node);
            cands.insert(node);
        }
    cands.clear();

    sort(candidates.begin(), candidates.end(), comp_fadas);

    for(int i=0;i<candidates.size() && node_num > 0;++i){
        if(!candidates[i]->isInternalNode()) {
            candidates[i]->search(k, queryTs, *heap, index_dir, hash_set);
            --node_num;
        }
        else {
            approxIncSearchInterNodeFuzzy(candidates[i], queryTs, sax, k, heap, index_dir, node_num, hash_set);
        }
    }

}