//
// Created by pengwang5 on 2022/1/16.
//

#include <set>
#include <unordered_set>
#include <chrono>
#include <mutex>
#include <atomic>
#include <cmath>
#include <queue>
#include <io.h>
#include "liburing.h"
#include "../../include/Searchers/FADASSearcher.h"
#include "../../include/Utils/MathUtil.h"
#include "../../include/Utils/TimeSeriesUtil.h"
#include "../../include/DataStructures/SafePq.h"
#include "../../include/DataStructures/SafeHashMap.h"

extern long LB_SERIES_TIME = 0, HEAP_TIME = 0,
    LB_NODE_TIME_STAT = 0, LB_NODE_CNT = 0, LOADED_NODE_CNT = 0;
extern double DIST_CALC_TIME = 0, READ_TIME = 0, PREPARE_TIME = 0, SEARCH_TIME = 0;
extern int layer, _search_num;
static FADASNode* targetNode;
vector<PqItemSeries *> * FADASSearcher::approxSearch(FADASNode *root, float *query, int k, vector<vector<int>> *g,
                                                     const string &index_dir,
                                                     float *query_reordered, int *ordering) {
    auto* queryTs = new TimeSeries(query);
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    int head = SaxUtil::invSaxHeadFromSax(sax, Const::bitsCardinality, Const::segmentNum);
    FADASNode *cur = (root->children)[head];
    if(cur == nullptr){
        FADASNode *node = nullptr;
        double min_dist = numeric_limits<double>::max();
        for(int i =0;i<Const::vertexNum;++i){
            if(root->children[i] == nullptr)    continue;
            double dist  = SaxUtil::getMinDist1stLayer(queryTs->paa , i);
            if(dist < min_dist){
                min_dist = dist;
                node = root->children[i];
            }
        }

//        for(int i=0;i<FADASNode::a + FADASNode::b + FADASNode::c;++i){
//            if(root->children[(*g)[head][i]] != nullptr){
//                node = root->children[(*g)[head][i]];
//                break;
//            }
//        }
        assert(node!=nullptr);
        // we only concern whether the nearest node is a leaf or an internal node
        if(node->isInternalNode()){
            approxSearchInterNode(node, queryTs, sax, k, heap, index_dir, query_reordered, ordering);
        }else { node->search(k, queryTs, *heap, index_dir, query_reordered, ordering); targetNode = node;}
    }else if(!cur->isInternalNode()){
        { cur->search(k, queryTs, *heap, index_dir, query_reordered, ordering); targetNode = cur;}
    }else approxSearchInterNode(cur, queryTs, sax, k, heap, index_dir, query_reordered, ordering);
    
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void FADASSearcher::approxSearchInner4IncSearch(FADASNode *root, TimeSeries *queryTs, int k,
                                                vector<PqItemSeries *> *heap,
                                                vector<vector<int>> *g, unsigned short *sax,
                                                const string &index_dir,
                                                float *query_reordered, int *ordering) {
    int head = SaxUtil::invSaxHeadFromSax(sax, Const::bitsCardinality, Const::segmentNum);

    FADASNode *cur = (root->children)[head];
    if(cur == nullptr){
        FADASNode *node = nullptr;
        double min_dist = numeric_limits<double>::max();
        for(int i =0;i<Const::vertexNum;++i){
            if(root->children[i] == nullptr)    continue;
            double dist  = SaxUtil::getMinDist1stLayer(queryTs->paa , i);
            if(dist < min_dist){
                min_dist = dist;
                node = root->children[i];
            }
        }

        // we only concern whether the nearest node is a leaf or an internal node
        if(node->isInternalNode()){
            approxSearchInterNode(node, queryTs, sax, k, heap, index_dir, query_reordered, ordering);
        }else { node->search(k, queryTs, *heap, index_dir, query_reordered, ordering); targetNode = node;}

    }else if(!cur->isInternalNode()){
        { cur->search(k, queryTs, *heap, index_dir, query_reordered, ordering); targetNode = cur;
        }
    }else approxSearchInterNode(cur, queryTs, sax, k, heap, index_dir, query_reordered, ordering);

}

void FADASSearcher::approxSearchInterNode(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                          vector<PqItemSeries *> *heap, const string &index_dir,
                                          float *query_reordered, int *ordering) {
    FADASNode *cur = root->route(sax);

    if(!cur->isInternalNode()){
        cur->search(k, queryTs, *heap, index_dir, query_reordered, ordering);
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
            dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, cur->sax, cur->bits_cardinality, cur->chosenSegments, i);
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
        approxSearchInterNode(node, queryTs, sax, k, heap, index_dir, query_reordered, ordering);
        return;
    }else { node->search(k, queryTs, *heap, index_dir, query_reordered, ordering); targetNode = node;}
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
        approxSearchInterNode(node, queryTs, sax, k, heap, index_dir, nullptr, nullptr);
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
            return parent->layer > pt.parent->layer;
        return parent < pt.parent;
    }

    bool operator >(const PqItemFadasId& pt) const{
        if(dist != pt.dist)
            return dist > pt.dist;
        if(parent->layer != pt.parent->layer)
            return parent->layer < pt.parent->layer;
        return parent > pt.parent;
    }
};

struct cmp_PqItemDumpyId{
    bool operator()(const PqItemFadasId& a, const PqItemFadasId& pt) const{
        if(a.dist != pt.dist)
            return a.dist > pt.dist;
        if(a.parent->layer != pt.parent->layer)
            return a.parent->layer < pt.parent->layer;
        return a.parent > pt.parent;
    }
};

vector<PqItemSeries*>*FADASSearcher::exactSearch(FADASNode* root, float *query, int k, vector<vector<int>> *g){
    vector<PqItemSeries*>* heap = approxSearch(root, query, k, g, Const::fidxfn, nullptr, nullptr);
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

vector<PqItemSeries*>*FADASSearcher::exactSearchIdLevel(FADASNode* root, float *query, int k, vector<vector<int>> *g,
                                                        float *query_reordered, int *ordering){
    vector<PqItemSeries*>* heap = approxSearch(root, query, k, g, Const::fidxfn, query_reordered, ordering);
    unordered_set<FADASNode*>visited;
    visited.insert(targetNode);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);

    set<PqItemFadasId>pq;
    for(int i =0;i<Const::vertexNum;++i){
        if(root->children[i] == nullptr || visited.count(root->children[i]) > 0)    continue;
        double dist  = SaxUtil::getMinDist1stLayer(queryTs->paa , i);
        pq.insert(PqItemFadasId(root,i, dist));
    }

    PqItemFadas cur;
    double top_dist;
    FADASNode* node;
    int len;
    while(!pq.empty()){
        top_dist = pq.begin()->dist;
        if(top_dist > bsf)  break;
        node = pq.begin()->parent->children[pq.begin()->id];
        pq.erase(pq.begin());
        if(visited.count(node) > 0) continue;
        visited.insert(node);
        if(node->isInternalNode()){
            len = (1 << (node->chosenSegments.size()));
            for(int i =0;i<len;++i){
                if(node->children[i] == nullptr)    continue;
                double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax,
                                                           node->bits_cardinality, node->chosenSegments, i);
//                double dist_simd = SaxUtil::LowerBound_Paa_iSax_SIMD(queryTs->paa, node->sax,
//                                                                 node->bits_cardinality, node->chosenSegments, i);

                if(dist < bsf){
                    pq.insert(PqItemFadasId(node, i, dist));
                }
            }
        }else{
            ++LOADED_NODE_CNT;
//            node->search_SIMD_reordered(k, queryTs, *heap, Const::fidxfn, query_reordered, ordering);
            node->search_SIMD(k,queryTs,*heap,Const::fidxfn);
            for(auto item:*heap){
                if(fabs(item->dist - 11.36) < 1e-2)
                    cout << "here";
            }
            bsf = (*heap)[0]->dist;
        }
    }
    pq.clear();
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

atomic<int> load_node_cnt;
atomic<int> search_num;
mutex m, pq_mutex;
void process_pq(double bsf, TimeSeries* queryTs, int k, vector<PqItemSeries*>* heap, SafeHashMap<FADASNode*, char>*visited,
                PriorityBlockingCollection<PqItemFadasId>* pq){
    PqItemFadasId cur;
    double top_dist;
    FADASNode* node;
    int internal_node = 0, leaf_node = 0;
    char _ = 0;
    int len;
    while(!pq->is_empty()){
        pq->take_prio(cur);
        top_dist = cur.dist;
        if(top_dist >= bsf) break;
        node = cur.parent->children[cur.id];
        if(visited->find(node, _)) continue;
        visited->insert(node, _);
        if(node->isInternalNode()){
            ++internal_node;
            len = (1 << (node->chosenSegments.size()));
            for(int i =0;i<len;++i){
                if(node->children[i] == nullptr)    continue;
                double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax,
                                                           node->bits_cardinality, node->chosenSegments, i);
                if(dist < bsf)
                    pq->emplace(node, i, dist);
            }
        }else{
            load_node_cnt++;
            search_num += node->size;
            ++leaf_node;
//            node->search_SIMD_reordered(k, queryTs, *heap, Const::fidxfn, query_reordered, ordering);
            auto ret = node->search_SIMD(k,queryTs, Const::fidxfn, bsf);
            int i=0;
            {
                lock_guard<mutex> g(m);
                bsf = (*heap)[0]->dist;
                for(;i<ret->size() && (*ret)[i]->dist < bsf ;++i){
                    pop_heap(heap->begin(),  heap->end(), PqItemSeriesMaxHeap());
                    delete heap->back();
                    heap->pop_back();
                    heap->push_back((*ret)[i]);
                    push_heap(heap->begin(),  heap->end(), PqItemSeriesMaxHeap());
                    bsf = (*heap)[0]->dist;
                }

            }
            for(; i < ret->size(); ++i)
                delete (*ret)[i];
            delete ret;
        }
    }
//    cout << this_thread::get_id() << endl;
//    Const::logPrint( ": has finished after processing " +
//                            to_string(internal_node) + " internal nodes and " +
//                            to_string(leaf_node) + " leaf nodes");
}

void process_pq_mypq(double bsf, TimeSeries* queryTs, int k, vector<PqItemSeries*>* heap, SafeHashMap<FADASNode*, char>*visited,
                     priority_queue<PqItemFadasId, vector<PqItemFadasId>, cmp_PqItemDumpyId>* pq){
    PqItemFadasId cur;
    double top_dist;
    FADASNode* node;
    int internal_node = 0, leaf_node = 0;
    char _ = 0;
    int len;
    while(true){
        {
            lock_guard<mutex> gm(pq_mutex);
            if (pq->empty()) break;
            cur = pq->top();
            pq->pop();
        }
        top_dist = cur.dist;
        if(top_dist >= bsf) break;
        node = cur.parent->children[cur.id];
        if(visited->find(node, _)) continue;
        visited->insert(node, _);
        if(node->isInternalNode()){
            ++internal_node;
            len = (1 << (node->chosenSegments.size()));
            for(int i =0;i<len;++i){
                if(node->children[i] == nullptr)    continue;
                double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax,
                                                           node->bits_cardinality, node->chosenSegments, i);
                if(dist < bsf){
                    lock_guard<mutex>gm(pq_mutex);
                    pq->emplace(node, i, dist);
                }

            }
        }else{
            load_node_cnt++;
            search_num += node->size;
            ++leaf_node;
//            node->search_SIMD_reordered(k, queryTs, *heap, Const::fidxfn, query_reordered, ordering);
            auto ret = node->search_SIMD(k,queryTs, Const::fidxfn, bsf);
            int i=0;
            {
                lock_guard<mutex> g(m);
                bsf = (*heap)[0]->dist;
                for(;i<ret->size() && (*ret)[i]->dist < bsf ;++i){
                    pop_heap(heap->begin(),  heap->end(), PqItemSeriesMaxHeap());
                    delete heap->back();
                    heap->pop_back();
                    heap->push_back((*ret)[i]);
                    push_heap(heap->begin(),  heap->end(), PqItemSeriesMaxHeap());
                    bsf = (*heap)[0]->dist;
                }

            }
            for(; i < ret->size(); ++i)
                delete (*ret)[i];
            delete ret;
        }
    }
//    cout << this_thread::get_id() << endl;
//    Const::logPrint( ": has finished after processing " +
//                            to_string(internal_node) + " internal nodes and " +
//                            to_string(leaf_node) + " leaf nodes");
}

void process_pq_messi(double bsf, TimeSeries* queryTs, int k, vector<PqItemSeries*>* heap, SafeHashMap<FADASNode*, char>*visited,
                      priority_queue<PqItemFadasId, vector<PqItemFadasId>, cmp_PqItemDumpyId>* messi_pq, mutex* mu){
    PqItemFadasId cur;
    double top_dist;
    FADASNode* node;
    int internal_node = 0, leaf_node = 0;
    char _ = 0;
    int len;
    while (true){
        {
            lock_guard<mutex> gm(*mu);
            if(messi_pq->empty())   break;
            cur = messi_pq->top();
            messi_pq->pop();
        }
        top_dist = cur.dist;
        if(top_dist >= bsf) break;
        node = cur.parent->children[cur.id];
        if(visited->find(node, _)) continue;
        visited->insert(node, _);
        assert(!node->isInternalNode());
        load_node_cnt++;
        search_num += node->size;
        ++leaf_node;
//            node->search_SIMD_reordered(k, queryTs, *heap, Const::fidxfn, query_reordered, ordering);
        auto ret = node->search_SIMD(k,queryTs, Const::fidxfn, bsf);
        int i=0;
        {
            lock_guard<mutex> g(m);
            bsf = (*heap)[0]->dist;
            for(;i<ret->size() && (*ret)[i]->dist < bsf ;++i){
                pop_heap(heap->begin(),  heap->end(), PqItemSeriesMaxHeap());
                delete heap->back();
                heap->pop_back();
                heap->push_back((*ret)[i]);
                push_heap(heap->begin(),  heap->end(), PqItemSeriesMaxHeap());
                bsf = (*heap)[0]->dist;
            }

        }
        for(; i < ret->size(); ++i)
            delete (*ret)[i];
        delete ret;
    }
}

vector<PqItemSeries*>*FADASSearcher::Par_exactSearchIdLevel(FADASNode* root, float *query, int k, vector<vector<int>> *g,
                                                        float *query_reordered, int *ordering){
    vector<PqItemSeries*>* heap = approxSearch(root, query, k, g, Const::fidxfn, query_reordered, ordering);
    SafeHashMap<FADASNode*, char>visited;
    char _ = 0;
    visited.insert(targetNode, _);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);


    timeval io;
    Const::timer_start(&io);
    PriorityBlockingCollection<PqItemFadasId> pq;
    for(int i =0;i<Const::vertexNum;++i){
        if(root->children[i] == nullptr || visited.find(root->children[i], _))    continue;
        double dist  = SaxUtil::getMinDist1stLayer(queryTs->paa , i);
        pq.emplace(root, i, dist);
    }
    PREPARE_TIME += Const::timer_end(&io);

    Const::timer_start(&io);
    int pre_step = 20;
    PqItemFadasId cur;
    double top_dist;
    FADASNode* node;
    int len;
    while(!pq.is_empty() && LOADED_NODE_CNT < pre_step){
        pq.take_prio(cur);
        top_dist = cur.dist;
        if(top_dist >= bsf)  break;
        node =cur.parent->children[cur.id];
        if(visited.find(node, _)) continue;
        visited.insert(node, _);
        if(node->isInternalNode()){
            len = (1 << (node->chosenSegments.size()));
            for(int i =0;i<len;++i){
                if(node->children[i] == nullptr)    continue;
                double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax,
                                                           node->bits_cardinality, node->chosenSegments, i);

                if(dist < bsf){
                    pq.emplace(node, i, dist);
                }
            }
        }else{
            ++LOADED_NODE_CNT;
//            node->search_SIMD_reordered(k, queryTs, *heap, Const::fidxfn, query_reordered, ordering);
            node->search_SIMD(k,queryTs,*heap, Const::fidxfn);
            bsf = (*heap)[0]->dist;
        }
    }
    DIST_CALC_TIME += Const::timer_end(&io);

    Const::timer_start(&io);
    if(!pq.is_empty() && top_dist < bsf){
        // start a thread pool
        int th_num = Const::thread_num;
        load_node_cnt = LOADED_NODE_CNT;
        search_num = _search_num;
        vector<thread>thread_pool;
        for(int i =0; i < th_num; ++i)
            thread_pool.emplace_back(thread(process_pq,  bsf,queryTs, k, heap, &visited, &pq));
        for(int i = 0; i < th_num; ++i)
            thread_pool[i].join();
        LOADED_NODE_CNT = load_node_cnt;
        _search_num = search_num;
    }
    SEARCH_TIME += Const::timer_end(&io);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

vector<PqItemSeries*>*FADASSearcher::Par_exactSearchIdLevel_MyPq(FADASNode* root, float *query, int k, vector<vector<int>> *g,
                                                            float *query_reordered, int *ordering){
    vector<PqItemSeries*>* heap = approxSearch(root, query, k, g, Const::fidxfn, query_reordered, ordering);
    SafeHashMap<FADASNode*, char>visited;
    char _ = 0;
    visited.insert(targetNode, _);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);


    timeval io;
    Const::timer_start(&io);
    priority_queue<PqItemFadasId, vector<PqItemFadasId>, cmp_PqItemDumpyId> pq;
    for(int i =0;i<Const::vertexNum;++i){
        if(root->children[i] == nullptr || visited.find(root->children[i], _))    continue;
        double dist  = SaxUtil::getMinDist1stLayer(queryTs->paa , i);
        pq.emplace(root, i, dist);
    }
    PREPARE_TIME += Const::timer_end(&io);

    Const::timer_start(&io);
    int pre_step = 20;
    PqItemFadasId cur;
    double top_dist;
    FADASNode* node;
    int len;
    while(!pq.empty() && LOADED_NODE_CNT < pre_step){
        cur = pq.top();
        pq.pop();
        top_dist = cur.dist;
        if(top_dist >= bsf)  break;
        node =cur.parent->children[cur.id];
        if(visited.find(node, _)) continue;
        visited.insert(node, _);
        if(node->isInternalNode()){
            len = (1 << (node->chosenSegments.size()));
            for(int i =0;i<len;++i){
                if(node->children[i] == nullptr)    continue;
                double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax,
                                                           node->bits_cardinality, node->chosenSegments, i);
                if(dist < bsf)
                    pq.emplace(node, i, dist);
            }
        }else{
            ++LOADED_NODE_CNT;
//            node->search_SIMD_reordered(k, queryTs, *heap, Const::fidxfn, query_reordered, ordering);
            node->search_SIMD(k,queryTs,*heap, Const::fidxfn);
            bsf = (*heap)[0]->dist;
        }
    }
    DIST_CALC_TIME += Const::timer_end(&io);

    Const::timer_start(&io);
    if(!pq.empty() && top_dist < bsf){
        // start a thread pool
        int th_num = Const::thread_num;
        load_node_cnt = LOADED_NODE_CNT;
        search_num = _search_num;
        vector<thread>thread_pool;
        for(int i =0; i < th_num; ++i)
            thread_pool.emplace_back(thread(process_pq_mypq,  bsf, queryTs, k, heap, &visited, &pq));
        for(int i = 0; i < th_num; ++i)
            thread_pool[i].join();
        LOADED_NODE_CNT = load_node_cnt;
        _search_num = search_num;
    }
    SEARCH_TIME += Const::timer_end(&io);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

size_t hash_pointer(FADASNode* node){
    static const auto shift = (size_t)log2(1 + sizeof(FADASNode));
    return (size_t)(node) >> shift;
}

vector<PqItemSeries*>*FADASSearcher::Par_exactSearchIdLevel_MESSI(FADASNode* root, float *query, int k, vector<vector<int>> *g,
                                                                 float *query_reordered, int *ordering){
    vector<PqItemSeries*>* heap = approxSearch(root, query, k, g, Const::fidxfn, query_reordered, ordering);
    vector<SafeHashMap<FADASNode*, char>>visited(Const::messi_pq_num);
    char _ = 0;
    for(auto & vis: visited)    vis.insert(targetNode, _);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);

//    hash<FADASNode*> pointer_hash; // warning!!!  default is the pointer itself, so it must a multiple of 4.
    vector<priority_queue<PqItemFadasId, vector<PqItemFadasId>, cmp_PqItemDumpyId>> messi_pqs(Const::messi_pq_num);

    FADASNode* node;
    int len;
    double top_dist;

    timeval io;
    Const::timer_start(&io);
    queue<FADASNode*>nodes_list;
    nodes_list.push(root);
    while(!nodes_list.empty()){
        node = nodes_list.front();
        nodes_list.pop();
        assert(node->isInternalNode());
        len = (1 << (node->chosenSegments.size()));
        for(int i =0;i<len;++i){
            if(node->children[i] == nullptr || visited[0].find(node->children[i], _))    continue;
            top_dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax, node->bits_cardinality, node->chosenSegments, i);
            if(top_dist >= bsf)  continue;
            if(node->children[i]->isInternalNode()){
                nodes_list.push(node->children[i]);
            }else{
                // the same node must be in the same queue
                size_t pq_id = hash_pointer(node->children[i]) % Const::messi_pq_num;
                messi_pqs[pq_id].emplace(node, i , top_dist);
            }
        }
    }
    PREPARE_TIME += Const::timer_end(&io);

    Const::timer_start(&io);
    {
        // start a thread pool
        int th_num = Const::thread_num;
        load_node_cnt = LOADED_NODE_CNT;
        search_num = _search_num;
        vector<thread>thread_pool;
//        vector<int>pq_sizes;
//        pq_sizes.reserve(messi_pqs.size());
//        for(auto & _: messi_pqs)
//            pq_sizes.push_back(_.size());
        vector<mutex>pq_mutexs(Const::messi_pq_num);
        for(int i =0; i < th_num; ++i)
            thread_pool.emplace_back(thread(process_pq_messi,  bsf, queryTs, k, heap, &visited[i % Const::messi_pq_num],
                                            &messi_pqs[i % Const::messi_pq_num], &pq_mutexs[i % Const::messi_pq_num]));
        for(int i = 0; i < th_num; ++i)
            thread_pool[i].join();
        LOADED_NODE_CNT = load_node_cnt;
        _search_num = search_num;
    }
    SEARCH_TIME += Const::timer_end(&io);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

struct io_data{
    int node_size{};
    float lb_dist{};
    int fd{};
    float* tss{};
};

void raw_search_SIMD(float *tss, int size, float* query, vector<PqItemSeries *> &heap, int k, double& bsf){
    double dist;
    for(int i = 0; i < size; ++i){
        dist = TimeSeriesUtil::euclideanDist_SIMD(query, tss + i * Const::tsLength , Const::tsLength, bsf);
        if(heap.size() < k){
            heap.push_back(new PqItemSeries(tss + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(tss + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }
        if(heap.size() >= k)    bsf = heap[0]->dist;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
}

vector<PqItemSeries*>*FADASSearcher::Par_exactSearchIdLevel_SSD(FADASNode* root, float *query, int k, vector<vector<int>> *g,
                                                                  float *query_reordered, int *ordering){
    vector<PqItemSeries*>* heap = approxSearch(root, query, k, g, Const::fidxfn, query_reordered, ordering);
    SafeHashMap<FADASNode*, char>visited;
    char _ = 0;
    visited.insert(targetNode, _);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);

    priority_queue<PqItemFadasId, vector<PqItemFadasId>, cmp_PqItemDumpyId> pq;
    for(int i =0;i<Const::vertexNum;++i){
        if(root->children[i] == nullptr || visited.find(root->children[i], _))    continue;
        double dist  = SaxUtil::getMinDist1stLayer(queryTs->paa , i);
        pq.emplace(root, i, dist);
    }

    io_uring ring{};
    auto ret = io_uring_queue_init(Const::SSD_pq_num, &ring, 0);
    io_uring_sqe* sqe;
    io_uring_cqe* cqe;
    vector<vector<io_data>>io_buffer(2, vector<io_data>(Const::SSD_pq_num));
    for(int i = 0; i < 2 ; ++i){
        for(int j = 0; j < Const::SSD_pq_num; ++j){
            io_buffer[i][j].tss = new float [Const::th * Const::tsLength];
        }
    }
    int buffer_flag = 0, prev_buf_size =0, cur_buf_size = 0, err, len;
    bool first_time = true;
    PqItemFadasId cur;
    double top_dist;
    FADASNode* node; io_data* data;
    size_t fd;
    while(!pq.empty()){
        cur = pq.top();
        pq.pop();
        top_dist = cur.dist;
        if(top_dist >= bsf)  break;
        node = cur.parent->children[cur.id];
        if(node->isInternalNode()){
            len = (1 << (node->chosenSegments.size()));
            for(int i =0;i<len;++i){
                if(node->children[i] == nullptr)    continue;
                double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax,
                                                           node->bits_cardinality, node->chosenSegments, i);
                if(dist < bsf)
                    pq.emplace(node, i, dist);
            }
        }else{
            if(visited.find(node, _)) continue;
            visited.insert(node, _);
            ++LOADED_NODE_CNT;

            sqe = io_uring_get_sqe(&ring);
            assert(sqe);
            auto & tmp = io_buffer[buffer_flag][cur_buf_size];
            tmp.node_size = node->size;
            tmp.lb_dist = top_dist;
            string file = Const::fidxfn + node->getFileNameWrapper();
            tmp.fd = open((file).c_str(), O_RDONLY);
            assert(tmp.fd > 0);
            io_uring_prep_read(sqe, tmp.fd, tmp.tss, ((unsigned)node->size) * Const::tsLength * sizeof(float), 0);
            io_uring_sqe_set_data(sqe, &tmp);
            ++cur_buf_size;
            if(cur_buf_size == Const::SSD_pq_num){
                if(!first_time){
                    err = io_uring_wait_cqe_nr(&ring, &cqe, prev_buf_size);
                    assert(err == 0);
                    io_uring_cq_advance(&ring, prev_buf_size);
                    err = io_uring_submit(&ring);
                    assert(err == cur_buf_size);

                    for(int i = 0 ; i < prev_buf_size; ++i){
                        void * t = io_uring_cqe_get_data(cqe + i);
                        data = reinterpret_cast<io_data *>(t);
                        if(data->lb_dist < bsf){
                            _search_num += data->node_size;
                            raw_search_SIMD(data->tss, data->node_size, query, *heap, k, bsf);
                        }
                        assert(data->fd > 0);
                        close(data->fd);
                    }

                }else{
                    err = io_uring_submit(&ring);
                    assert(err == cur_buf_size);
                }
                first_time = false;
                prev_buf_size = cur_buf_size;
                cur_buf_size = 0;
                buffer_flag = 1 - buffer_flag;

            }
        }
    }

    err = io_uring_wait_cqe_nr(&ring, &cqe, prev_buf_size);
    assert(err == 0);
    io_uring_cq_advance(&ring, prev_buf_size);
    if(cur_buf_size > 0) {
        err = io_uring_submit(&ring);
        assert(err == cur_buf_size);
    }

    for(int i = 0 ; i < prev_buf_size; ++i){
        data = (io_data *)io_uring_cqe_get_data(cqe + i);
        if(data->lb_dist < bsf){
            _search_num += data->node_size;
            raw_search_SIMD(data->tss, data->node_size, query, *heap, k, bsf);
        }
        assert(data->fd > 0);
        close(data->fd);
    }

    if(cur_buf_size > 0){
        err = io_uring_wait_cqe_nr(&ring, &cqe, cur_buf_size);
        assert(err == 0);
        io_uring_cq_advance(&ring, cur_buf_size);

        for(int i = 0 ; i < cur_buf_size; ++i){
            data = (io_data *)io_uring_cqe_get_data(cqe + i);
            if(data->lb_dist < bsf){
                _search_num += data->node_size;
                raw_search_SIMD(data->tss, data->node_size, query, *heap, k, bsf);
            }
            assert(data->fd > 0);
            close(data->fd);
        }
    }

    delete queryTs;
    io_uring_queue_exit(&ring);
    for(int i = 0; i < 2 ; ++i)
        for(int j = 0; j < Const::SSD_pq_num; ++j)
            delete io_buffer[i][j].tss;

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


void searchSubTree(FADASNode *root, TimeSeries *queryTs, int k, vector<PqItemSeries *> *heap, const string &index_dir,int &node_num, float *query_reordered, int *ordering){
    unordered_set<FADASNode*>visited;
    if(!root->isInternalNode()){
        if(root != targetNode) {
            root->search(k, queryTs, *heap, index_dir, query_reordered, ordering);
            --node_num;
        }
        return;
    }
    for(auto child:root->children){
        if(child == nullptr || child == targetNode || visited.count(child) > 0)    continue;
        visited.insert(child);
        if(!child->isInternalNode()){
            child->search(k,queryTs,*heap, index_dir, query_reordered, ordering);
            --node_num;
        }else{
            searchSubTree(child, queryTs, k, heap, index_dir, node_num, query_reordered, ordering);
        }
    }
}

void searchSubTree(FADASNode *root, TimeSeries *queryTs, int k, vector<PqItemSeries *> *heap, const string &index_dir,
                   int &node_num, float *query_reordered, int *ordering, unordered_set<float*, createhash, isEqual>*hash_set){
    unordered_set<FADASNode*>visited;
    if(!root->isInternalNode()){
        if(root != targetNode) {
            root->search(k, queryTs, *heap, index_dir, hash_set, query_reordered, ordering);
            --node_num;
        }
        return;
    }
    for(auto child:root->children){
        if(child == nullptr || child == targetNode || visited.count(child) > 0)    continue;
        visited.insert(child);
        if(!child->isInternalNode()){
            child->search(k,queryTs,*heap, index_dir, hash_set, query_reordered, ordering);
            --node_num;
        }else{
            searchSubTree(child, queryTs, k, heap, index_dir, node_num, query_reordered, ordering, hash_set);
        }
    }
}

extern int nrest;
vector<PqItemSeries *> * FADASSearcher::approxIncSearch(FADASNode *root, float *query, int k, const string &index_dir,
                                                        int node_num,
                                                        float *query_reordered, int *ordering, vector<vector<int>> *g) {
    auto* queryTs = new TimeSeries(query);
    t_paa = queryTs->paa;
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    approxSearchInner4IncSearch(root, queryTs, k, heap, g, sax, index_dir, query_reordered, ordering);
    --node_num;
    if(node_num >= 0)
        approxIncSearchInterNode(root, queryTs, sax, k, heap, index_dir, node_num, query_reordered, ordering);

    nrest = node_num;
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}

void FADASSearcher::approxIncSearchInner(FADASNode *root, float *query, int k, const string &index_dir,
                                                        int node_num, vector<PqItemSeries*>*heap, TimeSeries* queryTs,unsigned short *sax,
                                                        float *query_reordered, int *ordering, vector<vector<int>> *g) {

    if(node_num >= 0)
        approxIncSearchInterNode(root, queryTs, sax, k, heap, index_dir, node_num, query_reordered, ordering);

}



vector<PqItemSeries *> *FADASSearcher::ngSearch(FADASNode *root, float *query, float *query_reordered, int *ordering,
                                                            int k,
                                                            vector<vector<int>> *g, int nprobes){
    auto *queryTs = new TimeSeries(query);
    auto heap = new vector<PqItemSeries*>();
    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];
//    approxSearchInner(root,queryTs, k,heap,g,sax, Const::fidxfn, query_reordered, ordering);
//    --nprobes;

    FADASNode* root_subtree = root->route1step(sax);
    unordered_set<FADASNode*>visited;
    if(root_subtree){
        if(!root_subtree->isInternalNode()){
            visited.insert(root_subtree);
            root_subtree->search(k, queryTs, *heap, Const::fidxfn, query_reordered, ordering);
            --nprobes;
        }else if(root_subtree->leaf_num <= Const::pre_read){
            // a small subtree
            if(nprobes >= root_subtree->leaf_num) {
                int _ =root_subtree->leaf_num;
                visited.insert(root_subtree);
                searchSubTree(root_subtree, queryTs, k, heap, Const::fidxfn, _, query_reordered, ordering);
                nprobes -= root_subtree->leaf_num;
            }else{
                int rest = nprobes;
                approxIncSearchInterNode(root_subtree, queryTs, sax, k, heap, Const::fidxfn, rest,  visited,
                                     query_reordered, ordering);
                nprobes = rest;
            }
        }else{
            // a big subtree
            int to_search = min(nprobes, Const::pre_read);
            int _ = to_search;
            approxIncSearchInterNode(root_subtree, queryTs, sax, k, heap, Const::fidxfn, to_search,  visited,
                                     query_reordered, ordering);
            nprobes = nprobes - _ + to_search;
        }
    }

    if(nprobes <= 0){
        delete queryTs;
        sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
        return heap;
    }

    double bsf = heap->size() < k ? numeric_limits<double>::max(): (*heap)[0]->dist;

    set<PqItemFadasId>pq;
    for(int i =0;i<Const::vertexNum;++i){
        if(root->children[i] == nullptr || visited.count(root->children[i]) > 0)    continue;
        double dist  = SaxUtil::getMinDist1stLayer(queryTs->paa , i);
        pq.insert(PqItemFadasId(root,i, dist));
    }
    int cur_probe = 0;
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
                if(node->children[i] == nullptr)    continue;
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
    vector<PqItemSeries*>* heap = approxSearch(root, query, k, g, Const::fidxfn, nullptr, nullptr);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);

    set<PqItemFadasId>pq;
    for(int i =0;i<Const::vertexNum;++i){
        if(root->children[i] == nullptr)    continue;
        double dist  = SaxUtil::getMinDist1stLayer(queryTs->paa , i);
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
    vector<PqItemSeries*>* heap = approxSearch(root, query, k, g, Const::fidxfn, nullptr, nullptr);
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


extern int __layer;
void FADASSearcher::approxIncSearchInterNode(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                             vector<PqItemSeries *> *heap, const string &index_dir,
                                             int &node_num, float *query_reordered, int *ordering) {
    if(!root->isInternalNode() || node_num <= 0)  return;
    FADASNode *cur = root->route1step(sax), *parent = root;
    while (cur!= nullptr && cur->isInternalNode() && cur->leaf_num > node_num) {
        parent = cur;
        cur = cur->route1step(sax);
    }

    if(root->layer == 0)
        __layer = cur->layer;
    if(cur!= nullptr && cur != targetNode){
        if(!cur->isInternalNode()){
            cur->search(k, queryTs, *heap, index_dir, query_reordered, ordering);
            --node_num;
        }else{
            searchSubTree(cur, queryTs,k, heap, index_dir, node_num, query_reordered, ordering);
//            approxIncSearchInterNode(cur, queryTs,sax,k, heap, index_dir, node_num, query_reordered, ordering);
        }
    }

    if(node_num <=0)    return;

    double bsf = (*heap)[0]->dist;
    vector<PqItemFadasId>candidates;
    int len = (1 << (parent->chosenSegments.size()));
    for(int i =0;i<len;++i){
        if(parent->children[i] == nullptr || parent->children[i] == cur)    continue;
        double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, parent->sax, parent->bits_cardinality, parent->chosenSegments, i);
        if(dist < bsf){
            candidates.emplace_back(parent, i, dist);
        }
    }
    sort(candidates.begin(),  candidates.end());

//    vector<FADASNode*>candidates;
//    unordered_set<FADASNode*>cands;
//    for(FADASNode *node: parent->children)
//        if(node != nullptr && node!=cur && cands.find(node) == cands.end()) {
//            double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax, node->layer);
//            cands.insert(node);
//            candidates.push_back(node);
//        }
//    cands.clear();
//    sort(candidates.begin(), candidates.end(), comp_fadas);

//    vector<PqItemFadas>candidates;
//    unordered_set<FADASNode*>cands;
//    double bsf = (*heap)[0]->dist;
//    for(FADASNode *node: parent->children)
//        if(node != nullptr && node!=cur && cands.find(node) == cands.end()) {
//            double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax, node->bits_cardinality);
//            if(dist > bsf)  continue;
//            cands.insert(node);
//            candidates.emplace_back(node, dist);
//        }
//    cands.clear();
//    sort(candidates.begin(), candidates.end());


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

    unordered_set<FADASNode*>visited;
    for(int i=0;i<candidates.size() && node_num > 0;++i){
        FADASNode* node = parent->children[candidates[i].id];
        if(candidates[i].dist > bsf)    break;
        if(visited.count(node) > 0) continue;
        visited.insert(node);
        if(!node->isInternalNode()) {
            node->search(k, queryTs, *heap, index_dir, query_reordered, ordering);
            --node_num;
        }
        else {
            if(node->leaf_num <= node_num){
                searchSubTree(node, queryTs,k, heap, index_dir, node_num, query_reordered, ordering);
            }else
                approxIncSearchInterNode(node, queryTs, sax, k, heap, index_dir, node_num, query_reordered, ordering);
        }
    }

//    for(int i=0;i<candidates.size() && node_num > 0;++i){
//        FADASNode* node = candidates[i];
//        if(!node->isInternalNode()) {
//            node->search(k, queryTs, *heap, index_dir, query_reordered, ordering);
//            --node_num;
//        }
//        else {
//            approxIncSearchInterNode(node, queryTs, sax, k, heap, index_dir, node_num, query_reordered, ordering);
//        }
//    }

//    for(int i=0;i<candidates.size() && node_num > 0;++i){
//        if(bsf < candidates[i].dist) break;
//        FADASNode* node = candidates[i].node;
//        if(!node->isInternalNode()) {
//            node->search(k, queryTs, *heap, index_dir, query_reordered, ordering);
//            --node_num;
//        }
//        else {
//            approxIncSearchInterNode(node, queryTs, sax, k, heap, index_dir, node_num, query_reordered, ordering);
//        }
//        bsf = (*heap)[0]->dist;
//    }



}

// for ng-search
void FADASSearcher::approxIncSearchInterNode(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                             vector<PqItemSeries *> *heap, const string &index_dir,
                                             int &node_num, unordered_set<FADASNode*>&visit, float *query_reordered, int *ordering) {
    if(node_num <= 0)  return;
    FADASNode *cur = root->route1step(sax), *parent = root;
    while (cur!= nullptr && cur->isInternalNode() && cur->leaf_num > node_num) {
        parent = cur;
        cur = cur->route1step(sax);
    }

    if(cur!= nullptr){
        visit.insert(cur);
        if(!cur->isInternalNode()){
            cur->search(k, queryTs, *heap, index_dir, query_reordered, ordering);
            --node_num;
        }else{
            searchSubTree(cur, queryTs,k, heap, index_dir, node_num, query_reordered, ordering);
        }
    }

    if(node_num <=0)    return;

    double bsf = (*heap)[0]->dist;
    vector<PqItemFadasId>candidates;
    int len = (1 << (parent->chosenSegments.size()));
    for(int i =0;i<len;++i){
        if(parent->children[i] == nullptr || parent->children[i] == cur)    continue;
        double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, parent->sax, parent->bits_cardinality, parent->chosenSegments, i);
        if(dist < bsf){
            candidates.emplace_back(parent, i, dist);
        }
    }
    sort(candidates.begin(),  candidates.end());

    for(int i=0;i<candidates.size() && node_num > 0;++i){
        FADASNode* node = parent->children[candidates[i].id];
        if(candidates[i].dist > bsf)    break;
        if(visit.count(node) > 0) continue;
        visit.insert(node);
        if(!node->isInternalNode()) {
            node->search(k, queryTs, *heap, index_dir, query_reordered, ordering);
            --node_num;
        }
        else {
            searchSubTree(node, queryTs,k, heap, index_dir, node_num, query_reordered, ordering);
        }
    }

}

// for ng-search fuzzy
void FADASSearcher::approxIncSearchInterNode(FADASNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                             vector<PqItemSeries *> *heap, const string &index_dir,
                                             int &node_num, unordered_set<FADASNode*>&visit, float *query_reordered, int *ordering,
                                             unordered_set<float*, createhash, isEqual>*hash_set) {
    if(node_num <= 0)  return;
    FADASNode *cur = root->route1step(sax), *parent = root;
    while (cur!= nullptr && cur->isInternalNode() && cur->leaf_num > node_num) {
        parent = cur;
        cur = cur->route1step(sax);
    }

    if(cur!= nullptr){
        visit.insert(cur);
        if(!cur->isInternalNode()){
            cur->search(k, queryTs, *heap, index_dir, hash_set, query_reordered, ordering);
            --node_num;
        }else{
            searchSubTree(cur, queryTs,k, heap, index_dir, node_num, query_reordered, ordering, hash_set);
        }
    }

    if(node_num <=0)    return;

    double bsf = (*heap)[0]->dist;
    vector<PqItemFadasId>candidates;
    int len = (1 << (parent->chosenSegments.size()));
    for(int i =0;i<len;++i){
        if(parent->children[i] == nullptr || parent->children[i] == cur)    continue;
        double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, parent->sax, parent->bits_cardinality, parent->chosenSegments, i);
        if(dist < bsf){
            candidates.emplace_back(parent, i, dist);
        }
    }
    sort(candidates.begin(),  candidates.end());

    for(int i=0;i<candidates.size() && node_num > 0;++i){
        FADASNode* node = parent->children[candidates[i].id];
        if(candidates[i].dist > bsf)    break;
        if(visit.count(node) > 0) continue;
        visit.insert(node);
        if(!node->isInternalNode()) {
            node->search(k, queryTs, *heap, index_dir, hash_set, query_reordered, ordering);
            --node_num;
        }
        else {
            searchSubTree(node, queryTs,k, heap, index_dir, node_num, query_reordered, ordering, hash_set);
        }
    }

}


vector<PqItemSeries *> *FADASSearcher::ngSearchFuzzy(FADASNode *root, float *query, float *query_reordered, int *ordering,
                                                int k, int nprobes){
    auto *queryTs = new TimeSeries(query);
    auto heap = new vector<PqItemSeries*>();
    auto*hash_set = new unordered_set<float*, createhash, isEqual>();
    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];
//    approxSearchInner(root,queryTs, k,heap,g,sax, Const::fidxfn, query_reordered, ordering);
//    --nprobes;

    FADASNode* root_subtree = root->route1step(sax);
    unordered_set<FADASNode*>visited;
    if(root_subtree){
        if(!root_subtree->isInternalNode()){
            visited.insert(root_subtree);
            root_subtree->search(k, queryTs, *heap, Const::fuzzyidxfn, hash_set, query_reordered, ordering);
            --nprobes;
        }else if(root_subtree->leaf_num <= Const::pre_read){
            // a small subtree
            if(nprobes >= root_subtree->leaf_num) {
                int _ = root_subtree->leaf_num;
                visited.insert(root_subtree);
                searchSubTree(root_subtree, queryTs, k, heap, Const::fuzzyidxfn, _, query_reordered, ordering, hash_set);
                nprobes -= root_subtree->leaf_num;
            }else{
                int rest = nprobes;
                approxIncSearchInterNode(root_subtree, queryTs, sax, k, heap, Const::fuzzyidxfn, rest,  visited,
                                         query_reordered, ordering, hash_set);
                nprobes = rest;
            }
        }else{
            // a big subtree
            int to_search = min(nprobes, Const::pre_read);
            int _ = to_search;
            approxIncSearchInterNode(root_subtree, queryTs, sax, k, heap, Const::fuzzyidxfn, to_search,  visited,
                                     query_reordered, ordering, hash_set);
            nprobes = nprobes - _ + to_search;
        }
    }

    if(nprobes <= 0){
        delete queryTs;
        sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
        return heap;
    }

    double bsf = heap->size() < k ? numeric_limits<double>::max(): (*heap)[0]->dist;

    set<PqItemFadasId>pq;
    for(int i =0;i<Const::vertexNum;++i){
        if(root->children[i] == nullptr || visited.count(root->children[i]) > 0)    continue;
        double dist  = SaxUtil::getMinDist1stLayer(queryTs->paa , i);
        pq.insert(PqItemFadasId(root,i, dist));
    }
    int cur_probe = 0;
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
                if(node->children[i] == nullptr)    continue;
                double dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax, node->bits_cardinality, node->chosenSegments, i);
                if(dist < bsf){
                    pq.insert(PqItemFadasId(node, i, dist));
                }
            }
        }else{
            node->search(k, queryTs, *heap, Const::fuzzyidxfn, hash_set, query_reordered, ordering);
            ++cur_probe;
            bsf = (*heap)[0]->dist;
        }
    }
    pq.clear();
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

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