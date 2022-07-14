//
// Created by wzy on 2021/8/9.
//

#include "../../include/Searchers/DSTreeApproxSearcher.h"
#include "../../include/Utils/TimeSeriesUtil.h"

void DSTreeApproxSearcher::approxKnnSearch(DSTreeNode *root, float* queryTs, int k, vector<PqItemSeries *> &heap){
    // find the target node
    auto* q = new InsertedSeries(queryTs);
    DSTreeNode* node = root->approximateSearch(q);
    //calc the knn
    TimeSeriesUtil::knnWithBsf(*node, *q, k, heap);
    delete q;
}

vector<PqItemSeries *> * DSTreeApproxSearcher::approxKnnSearchWithThreshold(DSTreeNode *root, float *queryTs, int k, int threshold) {
    auto *heap = new vector<PqItemSeries*>();
    auto* q = new InsertedSeries(queryTs);
    approxKnnSearchWithThresholdSub(root, k, threshold, *heap, q);
    delete q;
    return heap;
}

void DSTreeApproxSearcher::approxKnnSearchWithThresholdSub(DSTreeNode *root, int k,  int threshold,
                                                        vector<PqItemSeries *> &heap, InsertedSeries* q) {
    DSTreeNode* node = root->approximateSearch(q, threshold);
    if(node->size == threshold)
        TimeSeriesUtil::knnWithBsf(*node, *q, k, heap);
    else if(node->isLeafNode())
        node->knnRaw(*q, k, heap, threshold);
    else if(node->splitInfo->routeToLeft(q)){
        TimeSeriesUtil::knnWithBsf(*node->left, *q, k, heap);
        if(node->right->isLeafNode()){
            node->right->knnRaw( *q, k, heap, threshold - node->left->size);
        }else approxKnnSearchWithThresholdSub(node->right, k, threshold - node->left->size, heap, q);
    }else{
        TimeSeriesUtil::knnWithBsf(*node->right, *q, k, heap);
        if(node->left->isLeafNode()){
            node->left->knnRaw(*q, k, heap, threshold - node->right->size);
        }else approxKnnSearchWithThresholdSub(node->left, k, threshold - node->right->size, heap, q);
    }

}