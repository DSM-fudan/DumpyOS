//
// Created by wzy on 2021/11/8.
//

#ifndef MULGIFT_IPGNODE_H
#define MULGIFT_IPGNODE_H
#include <vector>
#include <unordered_map>
#include <fstream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>
#include "../../include/Const.h"
#include "PqItemSeries.h"
#include "IPGDataNode.h"
#include "TimeSeries.h"
#include "OffsetDist.h"

using namespace std;

class IPGNode {
    IPGNode();

    IPGNode(int _id, int _layer) {
        id = _id;
        layer = _layer;
    }

    vector<int> offsets{};

    static double scoreSegmentmean(IPGNode *node, int i);

    static vector<int> chooseSegment(IPGNode *node, int chosen_num, const string &method);

    void insert(int offset, bool isReal);

    void fakeInsert(int offset){
        offsets.push_back(offset);
        ++size;
    }

    void generateSaxAndCardinality(IPGNode *node, int new_id, vector<int> &chosen_segments);

    void statPaa();

    void growIPGIndex(const string &method);

    static int loadSax(const string &saxfn);

    void build1stLayer(long series_number) const;

    static void loadPaa(const string &paafn);

public:
    const static int power_2[];
    static unsigned short *saxes;
    static float *paas;
    static int ***combines;
    static string rowDataFileName;
    unsigned short sax[Const::segmentNum]{};
    float paa_mu[Const::segmentNum]{}, paa_max[Const::segmentNum]{}, paa_min[Const::segmentNum]{}, paa_sigma[Const::segmentNum]{},
            paa_dist[Const::segmentNum]{}, paa_up_min[Const::segmentNum]{}, paa_below_max[Const::segmentNum]{},
            paa_up_avg[Const::segmentNum]{}, paa_below_avg[Const::segmentNum]{},
            paa_up_sigma[Const::segmentNum]{}, paa_below_sigma[Const::segmentNum]{};
    int paa_up_size[Const::segmentNum]{}, paa_below_size[Const::segmentNum]{};
    vector<int> chosenSegments{};
    vector<double> split_line{};
    const static int mask_16[];
    int actual_size{};
    int sax_cnt = -1;
    IPGDataNode *dataNode{};

    static vector<int> (*method_list[6])(IPGNode *node, int chosen_num);


    static IPGNode *BuildIPG(string &saxfn, string &paafn, vector<vector<int>> *g, const string &method);

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & id;ar & offsets;ar & size;ar & layer;ar & children;ar & partitionId;ar & dataNode;
        ar & actual_size;
        ar & sax;   ar & bits_cardinality;
        // for dynamic split line
//        ar & paa_mu;
    }

    void save2Disk(const string &output) {
        ofstream ofs(output.c_str(), ios::binary);
        boost::archive::binary_oarchive oa(ofs);
        oa << (*this);
        ofs.close();
    }

    static IPGNode *loadFromDisk(const string &saxfn, const string &paafn, const string &idxfn, bool loadFull);

//    static double scoreGraph(unordered_map<int, IPGNode *> *nodes_map);

    int size = 0;

    bool isLeafNode() const;

    int id = -1, partitionId = -1;

    static void partition(IPGNode *root, vector<vector<int>> *g);

    int getHeight() const;

    int getFakeLeafNodeNumber() const;

    static long generateSaxAndPaaTbl() ;

    static IPGNode *loadFromDisk(const string &saxfn, const string &idxfn);

    void exactSearchKnn(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap) const;

    vector<OffsetDist *> *preparePqUsingSaxLbNew(double bsf, const float *queryPaa) const;

    int layer = 0;
    unordered_map<int, IPGNode *> *children{};

    IPGNode *route(const unsigned short *asax) const;

    int bits_cardinality[Const::segmentNum]{};

    void statPaaMean();

    void statPaaOne();

    static double scoreSegmentstdev(IPGNode *node, int i);

    static double scoreSegmentrange(IPGNode *node, int i);

    static double scoreSegmentdistance(IPGNode *node, int i);

    static vector<int> chooseSegmentsigmarange(IPGNode *node, int chosen_num);

    static IPGNode *BuildIPGFuzzy(string &saxfn, string &paafn, vector<vector<int>> *g, const string &method);

    void growIPGIndexFuzzy(vector<vector<int>> *g, const string &method);

    static void fuzzyFirstLayer(IPGNode *root);

    static void fuzzy(IPGNode *parent, vector<int> &chosen_segments);

    static vector<int> chooseSegmentdistance(IPGNode *node, int chosen_num);

    static vector<int> chooseSegmentdunn(IPGNode *node, int chosen_num);

    static double scoreSegmentdunn(IPGNode *node, int i);

    static double scoreSegment3sigma(IPGNode *node, int i);

    static vector<int> chooseSegment3sigma(IPGNode *node, int chosen_num);

    static IPGDataNode *buildDataNode(const vector<IPGNode *> &nodes);

    static IPGDataNode *buildDataNode(IPGNode *node);

    void fuzzyNodeInFirstLayer(IPGNode *node) const;

    void fuzzyNode(IPGNode *node, int chosen_num, vector<int> &chosen_segments) const;

    static double scoreSegment3sigmaNorm(IPGNode *node, int i);

    static vector<int> chooseSegment3sigmaNorm(IPGNode *node, int chosen_num);

    static vector<int> chooseSegmentStdev(IPGNode *node, int chosen_num);

    static vector<int> chooseSegmentRange(IPGNode *node, int chosen_num);

    void statPaa1stLayer();

    static int chooseOneSegment(IPGNode *node);

    static IPGNode *BuildDynamicIPG(string &saxfn, string &paafn, vector<vector<int>> *g);

    static double scoreSegmentDynamic(IPGNode *node, int i, vector<float *> &paa_list, double *line);

    static void chooseSegmentDynamic(IPGNode *node, int chosen_num);

    void growDynamicIPGIndex(vector<vector<int>> *g);

    int getDynamicId(const float *paa) const;

    IPGNode *route(const float *paa) const;

    void build1stLayerDynamic(int series_number);

    static void chooseSegmentDynamicRange(IPGNode *node, int chosen_num);

    void growDynamicIPGIndexRange(vector<vector<int>> *g);

    IPGNode *routeIn1stLayer(const float *paa) const;

    IPGNode *routeBelow1stLayer(const float *paa) const;

    void statPaaRange();

    void build1stLayerDynamicMean(int series_number);

    static IPGNode *BuildIPGGrid(string &saxfn, string &paafn, vector<vector<int>> *g);

    static int ***readCombines();

    void countDistinctSaxWords();

    int countDistinctSaxWords(int pid) const;

    void trySplitting2Children(const int *chosen_segments, int chosen_num);

    void growIPGIndexGrid(vector<vector<int>> *g);

    static IPGNode *BuildIPGCluster(string &saxfn, string &paafn, vector<vector<int>> *g);

    void growIPGIndexCluster();

    void growIPGIndexClusterDynamic();

    static IPGNode *BuildIPGClusterDynamic(string &saxfn, string &paafn, vector<vector<int>> *g);

    IPGNode *routePaaMu(const float *paa) const;

    void statPaaMaxMinMean();
};

#endif //MULGIFT_IPGNODE_H
