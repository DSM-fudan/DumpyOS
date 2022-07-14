//
// Created by wzy on 2021/9/24.
//

#ifndef MULGIFT_TARDISTREENODE_H
#define MULGIFT_TARDISTREENODE_H
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include "../Const.h"
#include "../Utils/SaxUtil.h"
#include "../DataStructures/OffsetDist.h"
#include "../DataStructures/TimeSeries.h"
#include "../DataStructures/PqItemSeries.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/unordered_map.hpp>

using namespace std;

class TardisTreeNode {

    TardisTreeNode(int _id, int _layer, const unsigned short *_sax, TardisTreeNode* parent): id(_id), layer(_layer){
        offsets.clear();
        for(int i=0;i<Const::segmentNum;++i)
            sax[i] = _sax[i];
        SaxUtil::id2Sax2(_id, sax, Const::segmentNum);
        if(_layer == 1)
            file_id = to_string(_id);
        else if (_layer == 2){
            file_id = parent->file_id;
        }else{
            if(!parent->flag){
                parent->file_id += to_string(parent->id);
                parent->flag = true;
            }
            file_id = parent->file_id;
        }
    }

    TardisTreeNode(){for(unsigned short & i : sax)   i = 0; offsets.clear();}

    [[nodiscard]] string getFileName() const{
        return to_string(layer) + "-" + file_id + "_" + to_string(partitionId);
    }

    void insert(int offset);

public:
    vector<int>offsets{};
    int id = -1;
    int partitionId = -1;
    void exactSearchKnn(int k, TimeSeries *queryTs, vector<PqItemSeries*>&heap) const;
    TardisTreeNode *routeToTarget(unsigned short *asax);
    static TardisTreeNode *BuildTardisTree(string &saxfn, bool on_disk);

    int size = 0;
    bool flag = false;
    unsigned short sax[Const::segmentNum]{};
    int layer = 0;
    string file_id;
    unordered_map<int, TardisTreeNode*>*children{};
    static unsigned short * saxes;
    static float *dataset;
    static int a,b,c;

    vector<OffsetDist *> *preparePqUsingSaxLbNew(double bsf, float *queryPaa) const;

    friend class boost::serialization::access;
    template<class Archive>
            void serialize(Archive &ar, const unsigned int version){
                ar & id; ar & offsets;
                ar & size;
                ar & sax;
                ar & layer;
                ar & children;
                ar & partitionId;
            }

    void save2Disk(const string& output){
        ofstream ofs(output.c_str(), ios::binary);
        boost::archive::binary_oarchive oa(ofs);
        oa << (*this);
        ofs.close();
    }

    static TardisTreeNode *loadFromDisk(const string &saxfn, const string &idxfn,bool on_disk);

    int getLeafNodeNumber(TardisTreeNode *root);

    int getHeight(TardisTreeNode *root);

    bool isLeafNode() const{return children == nullptr;}

    int getGraphNum(TardisTreeNode *root);

//    void partition(TardisTreeNode *root, vector<vector<int>> *g);
//
//    static double scoreGraph(unordered_map<int, TardisTreeNode *> *nodes_map);

    TardisTreeNode *route(unsigned short *asax) const;

    static bool order(const TardisTreeNode*a, const TardisTreeNode *b){
        return a->size > b->size;
    }

    static void mergingTardis(TardisTreeNode *root);

    static TardisTreeNode *loadFromDisk(const string &idxfn);

    void exactSearchKnnLeaf(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap) const;

    void exactSearchKnnInMemory(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap) const;

    void exactSearchKnnInMemory(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, int &threshold) const;

    int getTotalSize(TardisTreeNode *root);

    int getNodeNumber(TardisTreeNode *root);

    void getIndexStats();

    int getSumHeight(TardisTreeNode *root) const;

    void materialize1stLayer(const string &index_dir) const;

    void materializeInterNode(const string &index_dir);

    void materialize(const string &index_dir) const;

    static long generateSaxTbl();

    void exactSearchKnn(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, int &node_number) const;

    void getLeafNodeSize(TardisTreeNode *root, ofstream &f);
};


#endif //MULGIFT_TARDISTREENODE_H
