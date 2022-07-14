//
// Created by wzy on 2021/8/9.
//

#include "../../include/Searchers/DSTreeExactSearcher.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Utils/TimeSeriesUtil.h"
#include "../../include/Const.h"
#include <cassert>

//extern long LBNodeTime, DSTreeProcessPqTime, distBigTime, DSTreeNodeApproxTime, IOBigTime, heapBigTime,
//LBSeriesTime2 , LBNodeTime2, DSTreeSearchTime2, DSTreePreparePqTime2, distBigTime2, DSTreeNodeApproxTime2, IOBigTime2, DSTreeProcessPqTime2;

vector<PqItemSeries *> & DSTreeExactSearcher::exactSearchKnn(DSTreeNode* root, float* queryTs, int k){
    auto q = new InsertedSeries(queryTs);
    DSTreeNode* bsfDSTreeNode = root->approximateSearch(q);
    vector<PqItemSeries*> &resultKnn = TimeSeriesUtil::knn(*bsfDSTreeNode, *q, k);
    //initialize priority queue;
    vector<PqItemNode*> pq;
    auto tempItem = new PqItemNode();
    tempItem->node = root;
    tempItem->dist = TimeSeriesUtil::minDistEstimation(root, q);
    pq.push_back(tempItem);

    //process the priority queue
    processPqInKnnSearch(q, k, resultKnn, bsfDSTreeNode, pq);
    sort_heap(resultKnn.begin(),  resultKnn.end(), PqItemSeriesMaxHeap());    //ascending order
    return resultKnn;
}

void DSTreeExactSearcher::exactSearchKnnWithBsf(DSTreeNode* root, float* queryTs, int k, vector<PqItemSeries *> &heap){
    auto q = new InsertedSeries(queryTs);
    DSTreeNode* bsfDSTreeNode = root->approximateSearch(q);
    TimeSeriesUtil::knnWithBsf(*bsfDSTreeNode,*q, k, heap);

    //initialize priority queue;
    vector<PqItemNode*> pq;
    //initialize the priority queue
    auto tempItem = new PqItemNode();
    tempItem->node = root;
    tempItem->dist = TimeSeriesUtil::minDistEstimation(root, q);
    pq.push_back(tempItem);

    //process the priority queue
    processPqInKnnSearch(q, k, heap, bsfDSTreeNode, pq);

    delete q;
}

void clearHeap(vector<PqItemNode*>& pq){
    for(auto * x: pq)
        delete x;
    pq.clear();
}

//void DSTreeExactSearcher::processPqInKnnSearch(InsertedSeries* q, int k, vector<PqItemSeriesVector *> &heap, DSTreeNode* bsfDSTreeNode, vector<PqItemNode*>& pq){
//    float* queryTs = q->ts;
//    PqItemNode* minPqItemNode;
//    PqItemNode* tempItem;
//    while (!pq.empty()) {
//        minPqItemNode = pq[0];
//        pop_heap(pq.begin(),  pq.end(), PqItemNodeComparator());
//        pq.pop_back();
//        if(minPqItemNode->node == bsfDSTreeNode){
//            delete minPqItemNode;
//            continue;
//        }
//        if(heap.size()>=k){
//            if(minPqItemNode->dist > heap[0]->dist) {
//                delete minPqItemNode;
//                break;
//            }
//        }
//
//        if (minPqItemNode->node->isLeafNode()) {
//            string fileName = minPqItemNode->node->getFileName();
//            FILE *f = fopen(fileName.c_str(), "rb");
//            assert(!heap.empty());
//            double maxPossibleDist = heap[0]->dist;
//            int cnt = 0;
//            while (cnt < minPqItemNode->node->size){
//                auto* ts = FileUtil::readSeriesVector(f);
//                double dist = TimeSeriesUtil::euclideanDist(ts, queryTs, Const::tsLength);
//
//                if(heap.size()<k){
//                    heap.push_back(new PqItemSeriesVector(ts, dist, true));
//                    push_heap(heap.begin(),  heap.end(), PqItemSeriesVectorMaxHeap());
//                }
//                else if(maxPossibleDist > dist){
//                    pop_heap(heap.begin(),  heap.end(), PqItemSeriesVectorMaxHeap());
//                    delete heap.back();
//                    heap.pop_back();
//                    heap.push_back(new PqItemSeriesVector(ts, dist, true));
//                    push_heap(heap.begin(),  heap.end(), PqItemSeriesVectorMaxHeap());
//                    maxPossibleDist = heap[0]->dist;
//                } else  delete ts;
//                ++cnt;
//            }
//            fclose(f);
//        } else {     //minPqItem is internal
//            //for left
//            tempItem = new PqItemNode();
//            tempItem->node = minPqItemNode->node->left;
//            tempItem->dist = TimeSeriesUtil::minDistEstimation(tempItem->node, q);
//
//            assert(!heap.empty());
//            if(heap[0]->dist > tempItem->dist)
//            {
//                pq.push_back(tempItem);
//                push_heap(pq.begin(),  pq.end(), PqItemNodeComparator());
//            } else delete tempItem;
//
//            //for right
//            tempItem = new PqItemNode();
//            tempItem->node = minPqItemNode->node->right;
//            tempItem->dist = TimeSeriesUtil::minDistEstimation(tempItem->node, q);
//
//            if(heap[0]->dist > tempItem->dist) {
//                pq.push_back(tempItem);
//                push_heap(pq.begin(),  pq.end(), PqItemNodeComparator());
//            } else delete tempItem;
//        }
//        delete minPqItemNode;
//    }
//    clearHeap(pq);
//}

void DSTreeExactSearcher::processPqInKnnSearch(InsertedSeries* q, int k, vector<PqItemSeries *> &heap, DSTreeNode* bsfDSTreeNode, vector<PqItemNode*>& pq){
    float* queryTs = q->ts;
    PqItemNode* minPqItemNode;
    PqItemNode* tempItem;
    while (!pq.empty()) {
        minPqItemNode = pq[0];
        pop_heap(pq.begin(),  pq.end(), PqItemNodeComparator());
        pq.pop_back();
        if(minPqItemNode->node == bsfDSTreeNode){
            delete minPqItemNode;
            continue;
        }
        if(heap.size()>=k){
            if(minPqItemNode->dist > heap[0]->dist) {
                delete minPqItemNode;
                break;
            }
        }

        if (minPqItemNode->node->isLeafNode()) {
            //read series in leaf node
            string fileName = minPqItemNode->node->getFileName();
            int size = minPqItemNode->node->size;
            FILE *f = fopen(fileName.c_str(), "rb");
            auto *tss = new float[size * Const::tsLength];
            fread(tss, sizeof(float ), size * Const::tsLength, f);
            fclose(f);

            assert(!heap.empty());
            double maxPossibleDist = heap[0]->dist;

            for(int i=0;i<size; ++i){
                double dist;
                if(heap.size()<k){
                    dist = TimeSeriesUtil::euclideanDist(tss + i * Const::tsLength, queryTs, Const::tsLength);
                    heap.push_back(new PqItemSeries(tss + i * Const::tsLength, dist, false, true));
                    push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
                    continue;
                }
                dist = TimeSeriesUtil::euclideanDist(tss + i * Const::tsLength, queryTs, Const::tsLength, maxPossibleDist);
                if(maxPossibleDist > dist){
                    pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
                    delete heap.back();
                    heap.pop_back();
                    heap.push_back(new PqItemSeries(tss + i * Const::tsLength, dist, false, true));
                    push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
                    maxPossibleDist = heap[0]->dist;
                }
            }

            TimeSeriesUtil::heap_data_copy(heap);
            delete[] tss;
        } else {     //minPqItem is internal
            //for left
            tempItem = new PqItemNode();
            tempItem->node = minPqItemNode->node->left;
            tempItem->dist = TimeSeriesUtil::minDistEstimation(tempItem->node, q);

            assert(!heap.empty());
            if(heap[0]->dist > tempItem->dist)
            {
                pq.push_back(tempItem);
                push_heap(pq.begin(),  pq.end(), PqItemNodeComparator());
            } else delete tempItem;

            //for right
            tempItem = new PqItemNode();
            tempItem->node = minPqItemNode->node->right;
            tempItem->dist = TimeSeriesUtil::minDistEstimation(tempItem->node, q);

            if(heap[0]->dist > tempItem->dist) {
                pq.push_back(tempItem);
                push_heap(pq.begin(),  pq.end(), PqItemNodeComparator());
            } else delete tempItem;
        }
        delete minPqItemNode;
    }
    clearHeap(pq);
}



void DSTreeExactSearcher::processPqInKnnSearchWithThreshold(InsertedSeries* q, int k, vector<PqItemSeries *> &heap,
                                                            DSTreeNode* bsfDSTreeNode, vector<PqItemNode*>& pq, int threshold){
    PqItemNode* minPqItemNode;
    PqItemNode* tempItem;
    float* queryTs = q->ts;
    int cnt = 0;
    while (!pq.empty() && cnt < threshold) {
        minPqItemNode = pq[0];
        pop_heap(pq.begin(),  pq.end(), PqItemNodeComparator());
        pq.pop_back();
        if(minPqItemNode->node == bsfDSTreeNode)   {
            delete minPqItemNode;
            continue;
        }
        if(heap.size()>=k){
            if(minPqItemNode->dist > heap[0]->dist) {
                delete minPqItemNode;
                break;
            }
        }

        if (minPqItemNode->node->isLeafNode()) {
            string fileName = minPqItemNode->node->getFileName();
            FILE *f = fopen(fileName.c_str(), "rb");
            int size= minPqItemNode->node->size;
            auto *tss = new float[size * Const::tsLength];
            fclose(f);

            assert(!heap.empty());
            double maxPossibleDist = heap[0]->dist;
            for(int i=0;i<size;++i){
                double dist = TimeSeriesUtil::euclideanDist(tss + i * Const::tsLength, queryTs, Const::tsLength, maxPossibleDist);
                if(heap.size()<k){
                    heap.push_back(new PqItemSeries(tss + i * Const::tsLength, dist, false, true));
                    push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
                }
                else if(maxPossibleDist > dist){
                    pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
                    delete heap.back();
                    heap.pop_back();
                    heap.push_back(new PqItemSeries(tss + i * Const::tsLength, dist, false, true));
                    push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
                    maxPossibleDist = heap[0]->dist;
                };
                if(++cnt >= threshold)  break;
            }
            TimeSeriesUtil::heap_data_copy(heap);
            delete[] tss;

        } else {     //minPqItem is internal
            //for left
            tempItem = new PqItemNode();
            tempItem->node = minPqItemNode->node->left;
            tempItem->dist = TimeSeriesUtil::minDistEstimation(tempItem->node, q);

            assert(!heap.empty());
            if(heap[0]->dist > tempItem->dist)
            {
                pq.push_back(tempItem);
                push_heap(pq.begin(),  pq.end(), PqItemNodeComparator());
            }   else delete tempItem;

            //for right
            tempItem = new PqItemNode();
            tempItem->node = minPqItemNode->node->right;
            tempItem->dist = TimeSeriesUtil::minDistEstimation(tempItem->node, q);

            if(heap[0]->dist > tempItem->dist) {
                pq.push_back(tempItem);
                push_heap(pq.begin(),  pq.end(), PqItemNodeComparator());
            }   else delete tempItem;

        }
        delete minPqItemNode;
    }
    clearHeap(pq);
}
//void DSTreeExactSearcher::processPqInKnnSearchWithThreshold(InsertedSeries* q, int k, vector<PqItemSeriesVector *> &heap,
//                                                            DSTreeNode* bsfDSTreeNode, vector<PqItemNode*>& pq, int threshold){
//    PqItemNode* minPqItemNode;
//    PqItemNode* tempItem;
//    float* queryTs = q->ts;
//    int cnt = 0;
//    while (!pq.empty() && cnt < threshold) {
//        minPqItemNode = pq[0];
//        pop_heap(pq.begin(),  pq.end(), PqItemNodeComparator());
//        pq.pop_back();
//        if(minPqItemNode->node == bsfDSTreeNode)   {
//            delete minPqItemNode;
//            continue;
//        }
//        if(heap.size()>=k){
//            if(minPqItemNode->dist > heap[0]->dist) {
//                delete minPqItemNode;
//                break;
//            }
//        }
//
//        if (minPqItemNode->node->isLeafNode()) {
//            string fileName = minPqItemNode->node->getFileName();
//            FILE *f = fopen(fileName.c_str(), "rb");
//            assert(!heap.empty());
//            double maxPossibleDist = heap[0]->dist;
//            int local_cnt = 0;
//            while (local_cnt < minPqItemNode->node->size){
//                auto* ts = FileUtil::readSeriesVector(f);
//                double dist = TimeSeriesUtil::euclideanDist(ts, queryTs, TimeSeries::tsLength);
//                if(heap.size()<k){
//                    heap.push_back(new PqItemSeriesVector(ts, dist, true));
//                    push_heap(heap.begin(),  heap.end(), PqItemSeriesVectorMaxHeap());
//                }
//                else if(maxPossibleDist > dist){
//                    pop_heap(heap.begin(),  heap.end(), PqItemSeriesVectorMaxHeap());
//                    delete heap.back();
//                    heap.pop_back();
//                    heap.push_back(new PqItemSeriesVector(ts, dist, true));
//                    push_heap(heap.begin(),  heap.end(), PqItemSeriesVectorMaxHeap());
//                    maxPossibleDist = heap[0]->dist;
//                } else delete ts;
//                if(++cnt >= threshold)  break;
//                ++local_cnt;
//            }
//            fclose(f);
//        } else {     //minPqItem is internal
//            //for left
//            tempItem = new PqItemNode();
//            tempItem->node = minPqItemNode->node->left;
//            tempItem->dist = TimeSeriesUtil::minDistEstimation(tempItem->node, q);
//
//            assert(!heap.empty());
//            if(heap[0]->dist > tempItem->dist)
//            {
//                pq.push_back(tempItem);
//                push_heap(pq.begin(),  pq.end(), PqItemNodeComparator());
//            }   else delete tempItem;
//
//            //for right
//            tempItem = new PqItemNode();
//            tempItem->node = minPqItemNode->node->right;
//            tempItem->dist = TimeSeriesUtil::minDistEstimation(tempItem->node, q);
//
//            if(heap[0]->dist > tempItem->dist) {
//                pq.push_back(tempItem);
//                push_heap(pq.begin(),  pq.end(), PqItemNodeComparator());
//            }   else delete tempItem;
//
//        }
//        delete minPqItemNode;
//    }
//    clearHeap(pq);
//}
//void DSTreeExactSearcher::exactSearchKnnWithBsfAndThreshold(DSTreeNode* root, float* queryTs, int k, vector<PqItemSeriesVector *> &heap, int threshold){
//    auto* q = new InsertedSeries(queryTs);
//    DSTreeNode* bsfDSTreeNode = root->approximateSearch(q);
//    TimeSeriesUtil::knnWithBsf(*bsfDSTreeNode,*q,k, heap);
//    //initialize priority queue;
//    vector<PqItemNode*> pq;
//    //initialize the priority queue
//    auto* tempItem = new PqItemNode();
//    tempItem->node = root;
//    tempItem->dist = TimeSeriesUtil::minDistEstimation(root, q);
//    pq.push_back(tempItem);
//
//    //process the priority queue
//    processPqInKnnSearchWithThreshold(q, k, heap, bsfDSTreeNode, pq, threshold);
//}

void DSTreeExactSearcher::exactSearchKnnWithBsfAndThreshold(DSTreeNode* root, float* queryTs, int k, vector<PqItemSeries *> &heap, int threshold){
    auto* q = new InsertedSeries(queryTs);
    DSTreeNode* bsfDSTreeNode = root->approximateSearch(q);
    TimeSeriesUtil::knnWithBsf(*bsfDSTreeNode,*q,k, heap);
    //initialize priority queue;
    vector<PqItemNode*> pq;
    //initialize the priority queue
    auto* tempItem = new PqItemNode();
    tempItem->node = root;
    tempItem->dist = TimeSeriesUtil::minDistEstimation(root, q);
    pq.push_back(tempItem);

    //process the priority queue
    processPqInKnnSearchWithThreshold(q, k, heap, bsfDSTreeNode, pq, threshold);
}

//void DSTreeExactSearcher::exactSearchKnnWithBsfAndThreshold2(DSTreeNode* root, float* queryTs, int k, vector<PqItemSeriesVector *> &heap, int threshold){
//    auto* q = new InsertedSeries(queryTs);
//
//    DSTreeNode* bsfDSTreeNode = root->approximateSearch(q);
//
//    //        TimeSeriesUtil.knnWithBsf(bsfDSTreeNode,q,k, heap);
//    //initialize priority queue;
//    vector<PqItemNode*> pq;
//    //initialize the priority queue
//    auto* tempItem = new PqItemNode();
//    tempItem->node = root;
//    tempItem->dist = TimeSeriesUtil::minDistEstimation(root, q);
//    pq.push_back(tempItem);
//    push_heap(pq.begin(),  pq.end(), PqItemNodeComparator());
//    pq.push_back(new PqItemNode(bsfDSTreeNode, -1));
//    push_heap(pq.begin(),  pq.end(), PqItemNodeComparator());
//
//    //process the priority queue
//    processPqInKnnSearchWithThreshold2(q, k, heap, bsfDSTreeNode, pq, threshold);
//}

//void DSTreeExactSearcher::processPqInKnnSearchWithThreshold2(InsertedSeries*  q, int k, vector<PqItemSeriesVector *> &heap,
//                                                             DSTreeNode* bsfDSTreeNode, vector<PqItemNode*>& pq, int threshold){
//    PqItemNode* minPqItemNode;
//    vector<PqItemNodeSeries*>seriesMaxHeap;
//    unordered_map<DSTreeNode*, vector<PqItemNodeSeries*>>nodeSeriesMap;
//    float* queryTs = q->ts;
//    int cnt = 0;
//    while (!pq.empty()) {
//        minPqItemNode = pq[0];
//        pop_heap(pq.begin(),  pq.end(), PqItemNodeComparator());
//        pq.pop_back();
//        auto node = minPqItemNode->node;
//        if(node == bsfDSTreeNode && cnt > 0)   continue;
//        ++cnt;
//
//        if(seriesMaxHeap.size()>=k){
//            if(minPqItemNode->dist >  seriesMaxHeap[0]->dist) break;
//        }
//
//        if (node->isLeafNode()) {
//            for(int id=0;id<node->means.size();++id){
//                double minDist = TimeSeriesUtil::minDistEstimation(*q, node->splitPoints, node->splitPointsLen, node->means[id], node->stdevs[id]);
//
//                if(seriesMaxHeap.size()<threshold){
//                    auto* item = new PqItemNodeSeries(node, id, minDist);
//                    seriesMaxHeap.push_back(item);
//                    push_heap(seriesMaxHeap.begin(),  seriesMaxHeap.end(), PqItemNodeSeriesMaxHeap());
//                    nodeSeriesMap[node].push_back(item);
//                    push_heap(nodeSeriesMap[node].begin(),  nodeSeriesMap[node].end(), PqItemNodeSeriesMaxHeap());
//                }
//                else {  // seriesMaxHeap is full
//                    if(minDist >= seriesMaxHeap[0]->dist)    continue;
//                    pop_heap(nodeSeriesMap[seriesMaxHeap[0]->node].begin(),  nodeSeriesMap[seriesMaxHeap[0]->node].end(), PqItemNodeSeriesMaxHeap());
//                    nodeSeriesMap[seriesMaxHeap[0]->node].pop_back();
//                    pop_heap(seriesMaxHeap.begin(),  seriesMaxHeap.end(), PqItemNodeSeriesMaxHeap());
//                    seriesMaxHeap.pop_back();
//
//                    auto* item = new PqItemNodeSeries(node, id, minDist);
//                    seriesMaxHeap.push_back(item);
//                    push_heap(seriesMaxHeap.begin(),  seriesMaxHeap.end(), PqItemNodeSeriesMaxHeap());
//                    nodeSeriesMap[node].push_back(item);
//                    push_heap(nodeSeriesMap[node].begin(),  nodeSeriesMap[node].end(), PqItemNodeSeriesMaxHeap());
//                }
//            }
//        } else {     //minPqItem is internal
//            //for left
//            pq.push_back(new PqItemNode(node->left, TimeSeriesUtil::minDistEstimation(node->left, q)));
//            push_heap(pq.begin(),  pq.end(), PqItemNodeComparator());
//            //for right
//            pq.push_back(new PqItemNode(node->right, TimeSeriesUtil::minDistEstimation(node->right, q)));
//            push_heap(pq.begin(),  pq.end(), PqItemNodeComparator());
//
//
//        }
//    }
//    clearHeap(pq);
//
//    for (auto cur : nodeSeriesMap) {
//        if (cur.second.empty()) continue;
//        FILE  * f= fopen(cur.first->getFileName().c_str(), "rb");
//        while (!cur.second.empty()) {
//            PqItemNodeSeries* curPq = cur.second[0];
//            pop_heap(cur.second.begin(),  cur.second.end(), PqItemNodeSeriesMaxHeap());
//            cur.second.pop_back();
//            fseek(f, curPq->offset * TimeSeries::tsLengthBytes, SEEK_SET);
//            auto* ts = FileUtil::readSeriesVector(f);
//            double dist = TimeSeriesUtil::euclideanDist(ts,queryTs, TimeSeries::tsLength);
//
//            if (heap.size() < k) {
//                heap.push_back(new PqItemSeriesVector(ts, dist, true));
//                push_heap(heap.begin(),  heap.end(), PqItemSeriesVectorMaxHeap());
//            }
//            else if (dist < heap[0]->dist) {
//                pop_heap(heap.begin(),  heap.end(), PqItemSeriesVectorMaxHeap());
//                delete heap.back();
//                heap.pop_back();
//                heap.push_back(new PqItemSeriesVector(ts, dist, true));
//                push_heap(heap.begin(),  heap.end(), PqItemSeriesVectorMaxHeap());
//            }
//        }
//        fclose(f);
//    }
//}
