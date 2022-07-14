//
// Created by wzy on 2021/8/7.
//

#ifndef MULGIFT_DSTREENODE_H
#define MULGIFT_DSTREENODE_H
#include "InsertedSeries.h"
#include "MeanNodeSegmentSplitPolicy.h"
#include "StdevNodeSegmentSplitPolicy.h"
#include "SplitInfo.h"

#include <string>
#include <vector>
#include<unordered_map>
#include <unordered_set>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
using namespace std;


class DSTreeNode {
public:
    static INodeSegmentSplitPolicy* nodeSegmentSplitPolicies[2];
    const static int maxSegmentLength = 2;
    const static int maxValueLength = 10;
    static unordered_map<const DSTreeNode*, float*>*rawdata;
    constexpr const static double hsTradeOffFactor = 2;
    vector<int> splitPoints{}, verticalSplitPoints{};
    int splitPointsLen{}, verticalSplitPointsLen{};
    vector<Sketch>* nodeSegmentSketches{},*verticalNodeSegmentSketches{}; //for horizontal splitting

    DSTreeNode* left{}, *right{};
    int level = 0;      //start from 0
    int size = 0;
    bool isLeft{};
    string indexPath;
    SplitInfo* splitInfo{};
    vector<InsertedSeries*>tss;
    vector<vector<double> >means;
    vector<vector<double> >stdevs;

    friend class boost::serialization::access;
    template<class Archive>
            void serialize(Archive &ar, const unsigned int version){
                ar & splitPoints; ar & splitPointsLen;
                ar & nodeSegmentSketches;
                ar & left; ar & right;
                ar & level; ar & size; ar & isLeft; ar & indexPath;
                ar & splitInfo;
                ar & means; ar & stdevs;
            }

    explicit DSTreeNode(string _indexPath);
    DSTreeNode(const DSTreeNode& parent) noexcept;
    DSTreeNode(){;}

    static void load2Memory(DSTreeNode* root);

    string getFileName() const;
    bool isLeafNode() const;

    void updateStatistics(InsertedSeries* timeSeries);
    void append(InsertedSeries* timeSeries);
    void disableTss();
    static double checkSplitPolicy(int i, double parentQoS, INodeSegmentSplitPolicy *nodeSegmentSplitPolicy, Sketch &nodeSegmentSketch,
                                   vector<int> &splitPoints);
    void insert(InsertedSeries *timeSeries, unordered_set<DSTreeNode *> &leafnodes);

    void initSegments(vector<int>&segmentPoints, int len);
    int getSegmentSize() const;
    static int getSegmentStart(const vector<int> &points, int idx);
    static int getSegmentEnd(const vector<int> &points, int idx);
    int getSegmentLength(int i) const;
    static int getSegmentLength(const vector<int> &points, int i);
    int getVerticalSplitPoint(vector<int>&points, int from, int to) const;


    static string formatInt(int value, int length);

    static string formatDouble(double value, int length);

    static void flushLeafNodes2Disk(unordered_set<DSTreeNode *>& leafNodes);

    void statFormatLeafNode();

    void saveToFile();

    DSTreeNode *approximateSearch(InsertedSeries *queryTs);

//    vector<float*> *loadTss(int maxIndex) const;
//
//    vector<vector<float> *> *loadTssVector(int maxIndex) const;

    void save2Disk();

    static DSTreeNode *loadFromFile();

    float *loadTssRaw(int maxIndex) const;

    DSTreeNode *approximateSearch(InsertedSeries *queryTs, int threshold);

    void knnRaw(InsertedSeries &q, int k, vector<PqItemSeries *> &heap, int threshold) const;
};


#endif //MULGIFT_DSTREENODE_H
