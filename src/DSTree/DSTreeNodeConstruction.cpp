//
// Created by wzy on 2021/8/7.
//

#include "../../include/DSTree/DSTreeNodeConstruction.h"
#include "../../include/Utils/TimeSeriesUtil.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Const.h"
#include <iostream>

//string DSTreeNodeConstruction::indexPath = "../data/expr_dstree_index/Index_%s_" + to_string(DSTreeNode::threshold) + "/";
//string DSTreeNodeConstruction::indexPath = "../data/index/Index_%s_" + to_string(Const::segmentNum)+ "_" + to_string(Const::th) + "/";


template<typename ... Args>
string str_format(const string &format, Args ... args)
{
    auto size_buf = std::snprintf(nullptr, 0, format.c_str(), args ...) + 1;
    std::unique_ptr<char[]> buf(new(std::nothrow) char[size_buf]);

    if (!buf)
        return string{};

    std::snprintf(buf.get(), size_buf, format.c_str(), args ...);
    return std::string(buf.get(), buf.get() + size_buf - 1);
}

//void DSTreeNodeConstruction::setIndexPath(const string& fileName){
//    // fileName for row data file name
//    int startIndex = fileName.find("Series"), endIndex = fileName.rfind('.');
//    indexPath = str_format(indexPath, fileName.substr(startIndex + 7, endIndex).c_str());
//    FileUtil::checkDirClean(indexPath.c_str());
//    indexPath += "Vertex_";
//}

DSTreeNode * DSTreeNodeConstruction::buildIndex(long series_num) {
    cout <<"Start build DS tree. ";
    // fileName for merged vertex file name
    int tsLength = Const::tsLength;
    if(series_num == -1){
        long size = FileUtil::getFileSize(Const::datafn.c_str());
        series_num = size / Const::tsLengthBytes;
    }

    cout<<"Series Number = " << series_num <<endl;
//    int startIndex = fileName.find("Vertex"),endIndex = fileName.rfind('.');
//    string curIndexPath = indexPath + fileName.substr(startIndex, endIndex);

//    int startIndex = fileName.find("Series"), endIndex = fileName.rfind(".bin");
//    string curIndexPath = str_format(indexPath, fileName.substr(startIndex + 7, endIndex).c_str());


    FileUtil::checkDirClean(Const::dstreefn.c_str());
    unordered_set<DSTreeNode*>leafNodes;
    auto* root = new DSTreeNode(Const::dstreefn);
    leafNodes.insert(root);
    //calc the split points by segmentSize
    auto* points = calcPoints(tsLength, initSegmentNum);
    root->initSegments(*points, initSegmentNum);

    FILE *f = fopen(Const::datafn.c_str(), "rb");
    auto *ts = new float[Const::tsLength * series_num];
    fread(ts, sizeof(float), Const::tsLength * series_num, f);
    for(long i=0;i<series_num;++i) {
        if(i % 1000000 == 0)    cout << i << endl;
        root->insert(new InsertedSeries(ts + i * Const::tsLength), leafNodes);
    }
    fclose(f);
    cout << "Start flush leaf nodes to disk.\t";
    DSTreeNode::flushLeafNodes2Disk(leafNodes);
    cout << "DS Tree Build Finished.\n";

    delete[] ts;
    return root;
}

//DSTreeNode * DSTreeNodeConstruction::buildIndex(const Vertex &v, const string &fileName, int parent_vertex_id) {
//
//    // fileName for merged vertex file name
//    int tsLength = TimeSeries::tsLength;
//    int seriesNum =v.cnt;
//    string curIndexPath;
//    if(parent_vertex_id == -1){
//        curIndexPath = indexPath + to_string(v.id);
//        cout <<"Start build DS tree. " << "Series Number = " << seriesNum <<endl;
//        cout << curIndexPath << endl;
//    }
//    else{
//        curIndexPath = indexPath + to_string(parent_vertex_id);
//        if(!FileUtil::checkFileExists(curIndexPath.c_str()))
//            FileUtil::createDir(curIndexPath.c_str());
//        curIndexPath += "/" + to_string(v.id);
//        cout <<"\tStart build DS tree. " << "Series Number = " << seriesNum <<endl;
//        cout << "\t"<<curIndexPath << endl;
//    }
//
//    FileUtil::checkDirClean(curIndexPath.c_str());
//
//    unordered_set<DSTreeNode*>leafNodes;
//    auto* root = new DSTreeNode(curIndexPath);
//    leafNodes.insert(root);
//    //calc the split points by segmentSize
//    auto* points = calcPoints(tsLength, initSegmentNum);
//    root->initSegments(*points, initSegmentNum);
//
//    FILE *f = fopen(fileName.c_str(), "rb");
//    auto *tss = new float[TimeSeries::tsLength * seriesNum];
//    int before = -1;
//    rewind(f);
//    long off;
//    for(int i=0;i<seriesNum;++i) {
//        off = (long)(v.tss_offsets[i] - before - 1) * TimeSeries::tsLengthBytes;
//        fseek(f, off, SEEK_CUR);
//        before = v.tss_offsets[i];
//        fread(tss + i*TimeSeries::tsLength, sizeof(float), TimeSeries::tsLength, f);
//        root->insert(new InsertedSeries(tss + i * TimeSeries::tsLength), leafNodes);
//    }
//    fclose(f);
//    if(parent_vertex_id != -1)  cout << "\t";
//    cout << "Start flush leaf nodes to disk.\t";
//    DSTreeNode::flushLeafNodes2Disk(leafNodes);
//    cout << "DS Tree Build Finished.\n";
//    delete[] tss;
//    return root;
//}
//
//DSTreeNode * DSTreeNodeConstruction::buildIndex(const Vertex &v, int parent_vertex_id, float *tss) {
//    // fileName for merged vertex file name
//    int tsLength = TimeSeries::tsLength;
//    int seriesNum =v.cnt;
//    string curIndexPath;
//    if(parent_vertex_id == -1){
//        curIndexPath = indexPath + to_string(v.id);
//        cout <<"Start build DS tree. " << "Series Number = " << seriesNum <<endl;
//        cout << curIndexPath << endl;
//    }
//    else{
//        curIndexPath = indexPath + to_string(parent_vertex_id);
//        if(!FileUtil::checkFileExists(curIndexPath.c_str()))
//            FileUtil::createDir(curIndexPath.c_str());
//        curIndexPath += "/" + to_string(v.id);
//        cout <<"\tStart build DS tree. " << "Series Number = " << seriesNum <<endl;
//        cout << "\t"<<curIndexPath << endl;
//    }
//
//    FileUtil::checkDirClean(curIndexPath.c_str());
//
//    unordered_set<DSTreeNode*>leafNodes;
//    auto* root = new DSTreeNode(curIndexPath);
//    leafNodes.insert(root);
//    //calc the split points by segmentSize
//    auto* points = calcPoints(tsLength, initSegmentNum);
//    root->initSegments(*points, initSegmentNum);
//
//    for(int i=0;i<seriesNum;++i)
//        root->insert(new InsertedSeries(tss + i * TimeSeries::tsLength), leafNodes);
//
//    if(parent_vertex_id != -1)  cout << "\t";
//    cout << "Start flush leaf nodes to disk.\t";
//    fflush(stdout);
//    DSTreeNode::flushLeafNodes2Disk(leafNodes);
//    cout << "DS Tree Build Finished.\n";
//    delete[] tss;
//    return root;
//}

vector<int> * DSTreeNodeConstruction::calcPoints(int tsLength, int segmentNo) {
    int avgLegnth = tsLength / segmentNo;
    auto* points = new vector<int>(segmentNo);
    for (int i = 0; i < segmentNo; i++) {
        (*points)[i] = (int) ((i + 1) * avgLegnth);
    }

    //set the last one
    (*points)[segmentNo - 1] = (int) tsLength;
    return points;
}

double DSTreeNodeConstruction::calcQoS(Sketch &sketch, int len){
    double mean_width = sketch.indicators[0] - sketch.indicators[1];
    double stdev_upper = sketch.indicators[2];
    return len * (mean_width * mean_width + stdev_upper * stdev_upper);
}

// This function is the bottleneck DSTree Construction, and most frequently called.
void DSTreeNodeConstruction::updateSketch(Sketch& nodeSegmentSketch, InsertedSeries* series, int fromIdx, int toIdx) {
    if (nodeSegmentSketch.indicators.size() != 4) //not initial
        {
        nodeSegmentSketch.indicators.resize(4);              //float
        nodeSegmentSketch.indicators[0] = -numeric_limits<double>::max(); //for max mean
        nodeSegmentSketch.indicators[1] = numeric_limits<double>::max(); //for min mean
        nodeSegmentSketch.indicators[2] = -numeric_limits<double>::max(); //for max stdev0
        nodeSegmentSketch.indicators[3] = numeric_limits<double>::max(); //for min stdev
        }

    double mean = series->getMean(fromIdx, toIdx);
    double stdev = series->getStdEv(fromIdx, toIdx);

    nodeSegmentSketch.indicators[0] = max(nodeSegmentSketch.indicators[0], mean);
    nodeSegmentSketch.indicators[1] = min(nodeSegmentSketch.indicators[1], mean);
    nodeSegmentSketch.indicators[2] = max(nodeSegmentSketch.indicators[2], stdev);
    nodeSegmentSketch.indicators[3] = min(nodeSegmentSketch.indicators[3], stdev);

}

