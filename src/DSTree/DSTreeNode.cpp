//
// Created by wzy on 2021/8/7.
//

#include "../../include/DSTree/DSTreeNode.h"
#include "../../include/Utils/MathUtil.h"
#include "../../include/DSTree/DSTreeNodeConstruction.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Utils/TimeSeriesUtil.h"
#include "../../include/Const.h"
#include "../../include/DataStructures/PqItemIndex.h"
#include <cstdio>
#include <utility>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>

using namespace std;

INodeSegmentSplitPolicy* DSTreeNode::nodeSegmentSplitPolicies[2]{new MeanNodeSegmentSplitPolicy(), new StdevNodeSegmentSplitPolicy()};
unordered_map<const DSTreeNode*, float*> *DSTreeNode::rawdata = new unordered_map<const DSTreeNode*, float*>();



void DSTreeNode::disableTss(){
    for(InsertedSeries* i : tss)
        delete i;
    vector<InsertedSeries*>().swap(tss);
}

DSTreeNode::DSTreeNode(string _indexPath) {
    indexPath = std::move(_indexPath);
}

DSTreeNode::DSTreeNode(const DSTreeNode& parent) noexcept{
    indexPath = parent.indexPath;
    level = parent.level + 1;
    splitInfo = parent.splitInfo;
}

int DSTreeNode::getSegmentSize() const {
    return splitPointsLen;
}

int DSTreeNode::getSegmentStart(const vector<int> &points, int idx) {
    if (idx == 0)
        return 0;
    else
        return points[idx - 1];
}

int DSTreeNode::getSegmentEnd(const vector<int> &points, int idx) {
    return points[idx];
}

int DSTreeNode::getSegmentLength(int i) const {
    if (i == 0)
        return splitPoints[i];
    else
        return splitPoints[i] - splitPoints[i - 1];
}

int DSTreeNode::getSegmentLength(const vector<int> &points, int i) {
    if (i == 0)
        return points[i];
    else
        return points[i] - points[i - 1];
}

void DSTreeNode::initSegments(vector<int>& segmentPoints, int len) {
    splitPoints.resize(len);
    splitPointsLen = len;
    copy(segmentPoints.begin(), segmentPoints.begin() + len, splitPoints.begin());
    verticalSplitPointsLen = MathUtil::split(segmentPoints, (int) 1, len, verticalSplitPoints); //min length is 1

    //init nodeSegmentSketches and hsNodeSegmentSketches
    nodeSegmentSketches = new vector<Sketch>(len);
    for (int i = 0; i < len; i++)
        (*nodeSegmentSketches)[i] = *(new Sketch());
    verticalNodeSegmentSketches = new vector<Sketch>(verticalSplitPointsLen);
    for (int i = 0; i < verticalSplitPointsLen; i++)
        (*verticalNodeSegmentSketches)[i] = *(new Sketch());

}

int DSTreeNode::getVerticalSplitPoint(vector<int>& points, int from, int to) const {
    if (!binary_search(points.begin(), points.begin() + splitPointsLen, to)) {
        return to;
    } else
        return from;
}

bool DSTreeNode::isLeafNode() const {
    return (left == nullptr && right == nullptr);
}

void DSTreeNode::append(InsertedSeries* timeSeries){
    tss.push_back(timeSeries);
}

// 更新该节点中每个segment的均值、方差的极值，其中某个节点才会表示该节点的范围
void DSTreeNode::updateStatistics(InsertedSeries* timeSeries) {
    size++;
    for (int i = 0; i < splitPointsLen; i++)
        DSTreeNodeConstruction::updateSketch((*nodeSegmentSketches)[i], timeSeries, getSegmentStart(splitPoints, i),getSegmentEnd(splitPoints, i));
    for (int i = 0; i < verticalSplitPointsLen; i++)
        DSTreeNodeConstruction::updateSketch((*verticalNodeSegmentSketches)[i], timeSeries, getSegmentStart(verticalSplitPoints, i), getSegmentEnd(verticalSplitPoints, i));
}

double DSTreeNode::checkSplitPolicy(int i, double parentQoS, INodeSegmentSplitPolicy* nodeSegmentSplitPolicy, Sketch &nodeSegmentSketch, vector<int> &splitPoints) {
    auto* childNodeSegmentSketches = nodeSegmentSplitPolicy->split(nodeSegmentSketch);
    int len = getSegmentLength(splitPoints, i);
    double leftChildQoS = DSTreeNodeConstruction::calcQoS((*childNodeSegmentSketches)[0], len),
    rightChildQoS = DSTreeNodeConstruction::calcQoS((*childNodeSegmentSketches)[1], len);

    return parentQoS - (leftChildQoS + rightChildQoS) / 2;
}

void DSTreeNode::insert(InsertedSeries *timeSeries, unordered_set<DSTreeNode *> &leafnodes) {
    //update statistics dynamically for leaf and branch

    updateStatistics(timeSeries);

    if (isLeafNode()) {
        append(timeSeries);            //append to file first

        if (Const::th == size) {   //do split

//            splitInfo = new SplitInfo();
            //init the vars used in loop
            // max DiffValue is the BSF answer for target function to assess the split strategy
            double maxDiffValue = -numeric_limits<double>::max();
            int verticalSplitPoint = -1; //default not do vertical split

            //we want to test every horizontal split policy for each segment
            for (int i = 0; i < splitPointsLen; i++) {
                //for each segment
                // QoS of original node
                double parentQoS = DSTreeNodeConstruction::calcQoS((*nodeSegmentSketches)[i], getSegmentLength(splitPoints, i));

                //for every split policy
                for (INodeSegmentSplitPolicy* nodeSegmentSplitPolicy : nodeSegmentSplitPolicies) {
                    double diffValue = checkSplitPolicy(i, parentQoS, nodeSegmentSplitPolicy, (*nodeSegmentSketches)[i], splitPoints);
                    if (diffValue > maxDiffValue) {
                        maxDiffValue = diffValue;
//                        splitInfo->splitFrom = getSegmentStart(splitPoints, id);
//                        splitInfo->splitTo = getSegmentEnd(splitPoints, id);
                        INodeSegmentSplitPolicy *tmp;
                        if(nodeSegmentSplitPolicy->getIndicatorSplitIdx() == 0)
                            tmp = new MeanNodeSegmentSplitPolicy(nodeSegmentSplitPolicy->getIndicatorSplitValue());
                        else
                            tmp = new StdevNodeSegmentSplitPolicy(nodeSegmentSplitPolicy->getIndicatorSplitValue());
                        splitInfo = new SplitInfo(getSegmentStart(splitPoints, i), getSegmentEnd(splitPoints, i), tmp);
                    }
                }
            }

            //wy add trade off for horizontal split; bias for horizontal split for no need to copy series data
            maxDiffValue = maxDiffValue * hsTradeOffFactor;

            for (int i = 0; i < verticalSplitPointsLen; i++) {
                //for each segment
                double parentQoS = DSTreeNodeConstruction::calcQoS((*verticalNodeSegmentSketches)[i], getSegmentLength(verticalSplitPoints, i));

                //for every split policy
                for (INodeSegmentSplitPolicy* nodeSegmentSplitPolicy : nodeSegmentSplitPolicies) {
                    double diffValue = checkSplitPolicy(i, parentQoS, nodeSegmentSplitPolicy, (*verticalNodeSegmentSketches)[i], verticalSplitPoints);

                    if (diffValue > maxDiffValue) {
                        maxDiffValue = diffValue;
//                        splitInfo->splitFrom = getSegmentStart(verticalSplitPoints, id);
//                        splitInfo->splitTo = getSegmentEnd(verticalSplitPoints, id);
//                        if(nodeSegmentSplitPolicy->getIndicatorSplitIdx() == 0)
//                            splitInfo->nodeSegmentSplitPolicy = new MeanNodeSegmentSplitPolicy(nodeSegmentSplitPolicy->getIndicatorSplitValue());
//                        else
//                            splitInfo->nodeSegmentSplitPolicy = new StdevNodeSegmentSplitPolicy(nodeSegmentSplitPolicy->getIndicatorSplitValue());
                        INodeSegmentSplitPolicy *tmp;
                        if(nodeSegmentSplitPolicy->getIndicatorSplitIdx() == 0)
                            tmp = new MeanNodeSegmentSplitPolicy(nodeSegmentSplitPolicy->getIndicatorSplitValue());
                        else
                            tmp = new StdevNodeSegmentSplitPolicy(nodeSegmentSplitPolicy->getIndicatorSplitValue());
                        splitInfo = new SplitInfo(getSegmentStart(verticalSplitPoints, i), getSegmentEnd(verticalSplitPoints, i), tmp);

                        verticalSplitPoint = getVerticalSplitPoint(splitPoints, splitInfo->splitFrom, splitInfo->splitTo);  // find the new points
                    }
                }
            }

            vector<int> childNodePoint;
            int len;
            if (verticalSplitPoint < 0) //not vs
            {
                childNodePoint.resize(splitPointsLen);
                len = splitPointsLen;
                copy(splitPoints.begin(), splitPoints.begin()+splitPointsLen, childNodePoint.begin());
            }
            else {
                childNodePoint.resize(splitPointsLen + 1);
                len = splitPointsLen + 1;
                copy(splitPoints.begin(), splitPoints.begin() + splitPointsLen, childNodePoint.begin());
                childNodePoint[splitPointsLen] = verticalSplitPoint;
                sort(childNodePoint.begin(), childNodePoint.begin()+splitPointsLen+1);
            }


            //init children node
            left = new DSTreeNode(*this);
            left->initSegments(childNodePoint, len);
            left->isLeft = true;

            right = new DSTreeNode(*this);
            right->initSegments(childNodePoint, len);
            right->isLeft = false;
            //read the time series from file, all the time series in this file include the new insert one
            //using buffer
            for (InsertedSeries* ts : tss)
                if (splitInfo->routeToLeft(ts))
                    left->insert(ts, leafnodes);
                else
                    right->insert(ts, leafnodes);


            auto iter = leafnodes.find(this);
            if(iter != leafnodes.end())
                leafnodes.erase(iter);
            if(left->isLeafNode())
                leafnodes.insert(left);
            if(right->isLeafNode())
                leafnodes.insert(right);
            vector<InsertedSeries*>().swap(tss);

        }
    } else { //not terminal
        //find the children and insert recursively
        if (splitInfo->routeToLeft(timeSeries))
            left->insert(timeSeries, leafnodes);
        else
            right->insert(timeSeries, leafnodes);
    }
}

string DSTreeNode::getFileName() const {
    assert(isLeafNode());
    string ret = Const::dstreefn ;
//    if (indexPath[indexPath.length() - 1] != '/')
//        ret += '/';
    ret = ret + formatInt(getSegmentSize(), maxSegmentLength);

    if(level == 0)
        return  ret + "_root";

    if (isLeft)
        ret = ret + "_L";
    else
        ret = ret + "_R";

    //add parent split policy
    ret += "_" + splitInfo->nodeSegmentSplitPolicy->getName() + "_";
    ret += "(" + to_string(splitInfo->splitFrom) + "," + to_string(splitInfo->splitTo) + "," +
            formatDouble(splitInfo->nodeSegmentSplitPolicy->getIndicatorSplitValue(), maxValueLength) + ")";
    ret += "_" + to_string(level);
    return ret;
}

string DSTreeNode::formatInt(int value, int length) {
    string ret = to_string(value);
    if (ret.length() > length) {
        cout << "exceed length:" << length;
    }
    while (ret.length() < length) {
        ret.insert(0, "0");
    }
    return ret;
}

string DSTreeNode::formatDouble(double value, int length) {
    string ret = to_string(value);
    if (ret.length() > length) {
        ret = ret.substr(0, length - 1);
    }
    return ret;
}

void DSTreeNode::saveToFile(){
    string fn = getFileName();
    float *start;
    FILE *f = fopen(fn.c_str(), "wb");
    for(InsertedSeries* item:tss){
        start = item->ts;
        fwrite(start, sizeof(float ), Const::tsLength, f);
    }
    fclose(f);
}

void DSTreeNode::flushLeafNodes2Disk(unordered_set<DSTreeNode *>& leafNodes) {
    for(auto node:leafNodes) {
        node->statFormatLeafNode();
        node->saveToFile();
        node->disableTss();
    }
    leafNodes.clear();
}

void DSTreeNode::statFormatLeafNode(){
    assert(isLeafNode());
    int num = tss.size();
    means.resize(num);
    stdevs.resize(num);
    int starts[splitPointsLen] ;
    int ends[splitPointsLen];
    for(int i=0;i<splitPointsLen;++i){
        starts[i] = getSegmentStart(splitPoints, i);
        ends[i] = getSegmentEnd(splitPoints, i);
    }
    InsertedSeries* tmp;
    for(int i = 0; i < num; ++i){
        tmp = tss[i];
        means[i].resize(splitPointsLen);
        stdevs[i].resize(splitPointsLen);
        for(int j=0; j<splitPointsLen; ++j){
            means[i][j] = tmp->getMean(starts[j], ends[j]);
            stdevs[i][j] = tmp->getStdEv(starts[j], ends[j]);
        }
    }
}


float* DSTreeNode::loadTssRaw(int maxIndex) const{
    int limit = min(maxIndex + 1, size);
    auto * res = new float[limit * Const::tsLength];
    FILE *f = fopen(getFileName().c_str(), "rb");
    fread(res, sizeof(float), limit * Const::tsLength, f);
    fclose(f);
    return res;
}

void DSTreeNode::save2Disk() {
    string fileName = Const::dstreefn + "root.idx";
    cout << "Save DSTree Node Root to file: "<<fileName << endl;
    ofstream ofs(fileName, ios::binary);
    boost::archive::binary_oarchive oa(ofs);
    oa << (*this);
    ofs.close();
}

DSTreeNode * DSTreeNode::loadFromFile() {
    string fileName = Const::dstreefn + "root.idx";
    cout << "Load DSTree Node Root From file " << fileName << endl;
    ifstream ifs(fileName, ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    auto *g = new DSTreeNode();
    ia >> (*g);
    ifs.close();
    return g;
}


//vector<float*>* DSTreeNode::loadTss(int maxIndex) const{
//    auto * res = new vector<float*>(size, new float[TimeSeries::tsLength]);
//    FILE *f = fopen(getFileName().c_str(), "rb");
//    int limit = min(maxIndex, size);
//    for(int id=0;id<limit;++id)
//        fread(&((*res)[id][0]), sizeof(float), TimeSeries::tsLength, f);
//    fclose(f);
//    return res;
//}
//vector<vector<float>*>* DSTreeNode::loadTssVector(int maxIndex) const{
//
//    FILE *f = fopen(getFileName().c_str(), "rb");
//    int limit = min(maxIndex + 1, size);
//    auto * res = new vector<vector<float>*>(limit);
//    for(int id=0;id<limit;++id)
//        (*res)[id] = FileUtil::readSeriesVector(f);
//    fclose(f);
//    return res;
//}
//static DSTreeNode loadFromFile(string fileName){
//    FSTObjectInput in = new FSTObjectInput(new FileInputStream(fileName));
//    DSTreeNode node = null;
//    try {
//        node = (DSTreeNode) in.readObject(DSTreeNode.class);
//    } catch (Exception e) {
//        throw new IOException(e.getMessage());
//    }
//    in.close();
//    return node;
//}
//    TreeInfo printTreeInfo() {
//        //1.get total node
//        //2.get terminal node(empty and not empty)
//        //3.get average amount for none empty terminal node
//
//        //using dept-first
//        //init first level node
//        TreeInfo tInfo = new TreeInfo();
//        List<Node> list = new vector<Node>();
//        list.add(this);
//
//        int totalCount = 0;
//        int emptyNodeCount = 0;
//        int noneEmptyNodeCount = 0;
//        int tsCount = 0;
//        int sumofLevel = 0;
//
//        //stack
//        while (list.size() > 0) {
//            //pop
//            Node node = list.remove(list.size() - 1);
//
//            totalCount++;
//
//            if (node.isTerminal()) {
//                //for terminal node
//                if (node.size > 0) {
//                    noneEmptyNodeCount++;
//                    tsCount += node.size;
//                    sumofLevel += node.level;
//                } else
//                    emptyNodeCount++;
//            } else {
//                //for internal node
//                //push their children if it is internal node
//                list.add(node.left);
//                list.add(node.right);
//            }
//        }
//        double avgLevel = sumofLevel / (double) noneEmptyNodeCount;
//        double avgPerNode = tsCount / (noneEmptyNodeCount == 0 ? 1 : noneEmptyNodeCount);
//
//        tInfo.setTotalCount(totalCount);
//        tInfo.setTsCount(tsCount);
//        tInfo.setEmptyNodeCount(emptyNodeCount);
//        tInfo.setNoneEmptyNodeCount(noneEmptyNodeCount);
//        tInfo.setAvgPerNode(avgPerNode);
//        tInfo.setAvgLevel(avgLevel);
//        return tInfo;
//    }

//    void toXml(stringBuffer xml) {
//        for (int id = 0; id < level; id++) {
//            xml.append(" ");
//        }
//        xml.append("<node").append(" ");
//        xml.append("level=\"" + level).append("\" ");
//        xml.append("segmentSize=\"" + getSegmentSize()).append("\" ");
//        xml.append("isTerminal=\"" + isTerminal()).append("\" ");
//        xml.append("fileName=\"" + getFileName()).append("\" ");
//        xml.append("size=\"" + size).append("\" ");
//        if (!isTerminal()) {
//            xml.append("splitFrom=\"" + splitPolicy.splitFrom).append("\" ");
//            xml.append("splitTo=\"" + splitPolicy.splitTo).append("\" ");
//            xml.append("splitIndex=\"" + splitPolicy.indicatorIdx).append("\" ");
//            xml.append("splitValue=\"" + splitPolicy.indicatorSplitValue).append("\" ");
//            xml.append("splitPolicy=\"" + splitPolicy.getNodeSegmentSplitPolicy().getClass().getSimpleName()).append("\" ");
//
//        }
//        xml.append(" ");
//
//        xml.append(">").append("\n");
//
//        if (left != null)
//            left.toXml(xml);
//
//        if (right != null)
//            right.toXml(xml);
//
//        for (int id = 0; id < level; id++) {
//            xml.append(" ");
//        }
//        xml.append("</node>").append("\n");
//    }

DSTreeNode* DSTreeNode::approximateSearch(InsertedSeries* queryTs) {
    //        System.out.println("this.getFileName() = " + this.getFileName());
    if (isLeafNode())
        return this;
    else //internal node
    {
        if (splitInfo->routeToLeft(queryTs))
            return left->approximateSearch(queryTs);
        else
            return right->approximateSearch(queryTs);
    }
}

DSTreeNode* DSTreeNode::approximateSearch(InsertedSeries* queryTs,int threshold) {
    //        System.out.println("this.getFileName() = " + this.getFileName());
    if (isLeafNode() || (left->size < threshold && right->size < threshold))
        return this;
    else //internal node
    {
        if (splitInfo->routeToLeft(queryTs))
            return left->approximateSearch(queryTs, threshold);
        else
            return right->approximateSearch(queryTs, threshold);
    }
}

void DSTreeNode::load2Memory(DSTreeNode *root) {
    if(root == nullptr) return;
    if(root->isLeafNode()){
        (*rawdata)[root] = root->loadTssRaw(root->size);
    }else{
        load2Memory(root->left);
        load2Memory(root->right);
    }
}


void DSTreeNode::knnRaw(InsertedSeries &q, int k, vector<PqItemSeries *> &heap, int threshold) const{
    auto *raw_tss = (*rawdata)[this];
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    int bound;
    if(threshold == -1) bound = size;
    else bound = min(threshold, size);
    for(int i=0;i<bound;++i){
        double dist = TimeSeriesUtil::euclideanDist(q.ts, raw_tss + i * Const::tsLength, Const::tsLength);
        if(heap.size() < k){
            heap.push_back(new PqItemSeries(raw_tss + i * Const::tsLength, dist, false, false));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            heap.pop_back();
            heap.push_back(new PqItemSeries(raw_tss + i * Const::tsLength, dist, false, false));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }
}


