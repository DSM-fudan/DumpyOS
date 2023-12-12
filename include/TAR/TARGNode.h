//
// Created by wzy on 2022/7/1.
//

#ifndef FADAS_TARGNODE_H
#define FADAS_TARGNODE_H
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <fstream>
#include "TARLNode.h"
#include "../Const.h"

struct dataRDDItem{
    string invsax_full{};
    int freq;
    dataRDDItem(string& a, int b):invsax_full(a), freq(b){};
};


using namespace std;
class TARGNode {
    static int pid_factory;
    static int block_threshold;


    static unsigned int *loadSampling();

    static bool order(const TARGNode*a, const TARGNode *b){
        return a->size > b->size;
    }

    static void computeFreqTbl(unordered_map<string, int> *invsax_freq, unsigned int *sample_invSAX_tbl, long tbl_size);

    static vector<unordered_map<string, int> *> *buildTARGbySampling(unordered_map<string, int> *invsax_freq);

    static void computeFreqTblSpecLayer(unordered_map<string, int> *pairRDD, unordered_map<string, dataRDDItem> *dataRDD,
                                 unordered_map<string, int> *saxNbrPairRDD, int layer);

    void addTreeNodeHexRobust(TARGNode *start);

    static void loadCollectionSampling(vector<unordered_map<string, int> *> *node_list, TARGNode *root);

    static void assignPIdSampling(TARGNode *node);

    void assignLocalRootForAllLeaves(vector<TARLNode *> &local_roots);

    static int getLeafNodeNum();

public:
    // general
    int layer{};
    int size{};
    string invSAX;
    TARGNode* ancestor{};

    //only for internal node

    unordered_map<string,TARGNode*>children{};

    // only for leaf pack
    int pid{-1};

    static TARGNode *buildIndex();

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & layer;
        ar & size;  ar & ancestor;
        ar & invSAX;   ar & pid;
        ar & children;
//        ar & local_root;
    }
    void save2Disk(const string &output) {
        ofstream ofs(output.c_str(), ios::binary);
        boost::archive::binary_oarchive oa(ofs);
        oa << (*this);
        ofs.close();
    }
    static TARGNode *loadFromDisk(const string &idxfn){
        ifstream ifs(idxfn, ios::binary);
        boost::archive::binary_iarchive ia(ifs);
        auto *g = new TARGNode();
        ia >> (*g);
        ifs.close();
        return g;
    };

    TARGNode * route(string &invsax, int *ret_pid);

    static void partition(TARGNode *root, long series_num);

    void getIndexStats();

    static int getTotalSize();

    int getGlobalNodeNum();

    int getNodeNum();

    int getMaxHeightGlobal(vector<int> &local_heights);

    int getMaxHeight();

    long getAvgHeight();

    long getAvgHeightGlobal(vector<int>&local_heights, vector<int>&leaf_nbr);
};


#endif //FADAS_TARGNODE_H
