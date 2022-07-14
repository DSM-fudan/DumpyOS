//
// Created by wzy on 2021/11/20.
//

#ifndef MULGIFT_ISAXNODE_H
#define MULGIFT_ISAXNODE_H
#include "../Const.h"
#include "TimeSeries.h"
#include "PqItemSeries.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <vector>
#include <fstream>

class iSAXNode {
    int segment_id = -1;

public:
    int bits_cardinality[Const::segmentNum]{};

    vector<int>offsets;
    vector<int>fbl;
    int size= 0, actual_size = 0;

    iSAXNode(int _layer, unsigned short *_sax): layer(_layer){
        for(int i=0;i<Const::segmentNum;++i)
            sax[i] = _sax[i] >> (Const::bitsCardinality - 1);
        size = 0;
    }

    iSAXNode(iSAXNode* parent, bool isLeft){
        layer = parent->layer + 1;
        copy(parent->bits_cardinality, parent->bits_cardinality + Const::segmentNum, bits_cardinality);
        bits_cardinality[parent->segment_id] +=1;
        if(bits_cardinality[parent->segment_id] > Const::bitsCardinality)  exit(-1);
        sax[parent->segment_id] <<=1;
        if(!isLeft) sax[parent->segment_id] +=1;
    }

    iSAXNode(){;}

    bool route2Left(const unsigned short *_sax);

    bool isLeafNode(){return left == nullptr && right == nullptr;}

    iSAXNode *route2Target(const unsigned short *_sax, bool update);

    void append(int offset, bool isReal);

    void split(bool fuzzy);

    int chooseSegment(bool is1st);

    void insert(int offset, bool isReal, bool fuzzy);

    void exactSearchKnn(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap);

    void exactSearchKnnInMemoryLeafNode(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap);

    friend class boost::serialization::access;
    template<class Archive>
            void serialize(Archive &ar, const unsigned int version){
                ar & offsets;
                ar & size;
                ar & sax;   ar & bits_cardinality;
                ar & layer;
                ar & left;  ar & right;
                ar & segment_id;
            }

    static void fuzzyDuplicate(iSAXNode *parent);

    iSAXNode *route2Target(const unsigned short *_sax, int threshold);

    iSAXNode *right = nullptr;

    void exactSearchKnnInMemory(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap);

    iSAXNode *left = nullptr;

    void exactSearchKnnInMemory(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, int &threshold);

    unsigned short sax[Const::segmentNum]{};
    int layer = -1;

    void exactSearchKnnInMemoryEntry(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, int &rest);

    iSAXNode *route1Step(const unsigned short *_sax);
};

class iSAXRoot{
public:
    iSAXNode* children[1 << Const::segmentNum];

    static unsigned short* saxes;
    static float *paas;
    static float *dataset;
    iSAXNode *route2Target(unsigned short *_sax, bool update);
    static iSAXRoot *buildIndex(const string &saxfn, long series_num);

    iSAXRoot(){
        for(auto & i : children) i = nullptr;
    }

    friend class boost::serialization::access;
    template<class Archive>
            void serialize(Archive &ar, const unsigned int version){
                ar & children;
            }

    void save2Disk(const string& output){
        ofstream ofs(output.c_str(), ios::binary);
        boost::archive::binary_oarchive oa(ofs);
        oa << (*this);
        ofs.close();
    }

    static long loadSax(const string &saxfn, long series_num);

    static iSAXRoot *loadFromDisk(const string &saxfn, const string &idxfn, long series_num);

    iSAXNode *route2firstLayer(unsigned short *_sax);

    void insertPrev(int offset);

    void bulkInsert(bool fuzzy);

    static void loadPaa(const string &paafn);

    static iSAXRoot *loadFromDisk(const string &saxfn, const string &paafn, const string &idxfn);

    static iSAXRoot *buildIndexFuzzy(const string &saxfn, const string &paafn);

    iSAXNode *route2Target(unsigned short *_sax, int threshold);
};


#endif //MULGIFT_ISAXNODE_H
