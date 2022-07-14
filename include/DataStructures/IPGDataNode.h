//
// Created by pengwang5 on 2021/12/7.
//

#ifndef MULGIFT_IPGDATANODE_H
#define MULGIFT_IPGDATANODE_H
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <vector>
#include "../Const.h"
using namespace std;

// a data node must be a leaf node
class IPGDataNode {
public:
    int size{};

    IPGDataNode(){
        size = 0;
    }

    explicit IPGDataNode(int _size){
        size = _size;
    }

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version){
        ar & size;
    }
};


#endif //MULGIFT_IPGDATANODE_H
