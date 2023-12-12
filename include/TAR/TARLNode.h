//
// Created by wzy on 2022/7/1.
//

#ifndef FADAS_TARLNODE_H
#define FADAS_TARLNODE_H
#include <string>
#include <vector>
#include <unordered_map>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#include "../DataStructures/PqItemSeries.h"
using namespace std;

struct item{
    long id;
    float *ts;
    string invsax;
};

struct ser_item{
    long id{};
    vector<float>ts{};
    string invsax{};

    ser_item(){;}

    ser_item(const ser_item& other){
        id = other.id;
        ts = other.ts;
        invsax = other.invsax;
    }

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & id;    ar & ts;    ar & invsax;
    }
};


class TARLNode {
    string invSAX{};

    // only for leaf node
    long file_offset{-1};
    long file_length{};
    vector<ser_item>buffer{};

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & layer; ar & rcdNbr;
        ar & invSAX; ar & descendants;
        ar & file_offset; ar & file_length;
        ar & buffer;
    }

    void split(int pid);

    void addRecord(item *data, int pid);

    void collectAllLeaves(vector<TARLNode *> &leaves);

public:
    int layer{};
    long rcdNbr{};
    // only for internal node
    unordered_map<string, TARLNode*> descendants{};

    void save2Disk(const string &output);
    static TARLNode *loadFromDisk(const string &idxfn);

    void insertBatch(vector<item *> &datas, int pid);

    TARLNode *route2Leaf(string &invsax);

    void search(int k, const float *query, vector<PqItemSeries *> &heap) const;

    void search_dtw(int k, float *query, vector<PqItemSeries *> &heap, int window_size) const;

    int getLeafNodeNbr();

    int getNodeNbr();

    static TARLNode * buildLocalIndex(int pid);

    void addRecordORIGIN(const ser_item &data, int pid);

    static TARLNode *prueSingleElementLayer(TARLNode *root);

    void fetchAllRecords(vector<ser_item> &buf);

    void shrinkNodeLayer();

    void deleteDescendants();

    int getMaxHeight();

    long getSumHeight();
};


#endif //FADAS_TARLNODE_H
