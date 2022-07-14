//
// Created by wzy on 2021/8/7.
//

#ifndef MULGIFT_SKETCH_H
#define MULGIFT_SKETCH_H
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <vector>
using namespace std;

class Sketch {
public:
    vector<double> indicators;    // length=4, respectively max mean, min mean, max stdev, min stdev

    Sketch(){
        indicators.resize(4);
        indicators[0] = -numeric_limits<double>::max(); //for max mean
        indicators[1] = numeric_limits<double>::max(); //for min mean
        indicators[2] = -numeric_limits<double>::max(); //for max stdev
        indicators[3] = numeric_limits<double>::max(); //for min stdev
    }

    friend class boost::serialization::access;
    template<class Archive>
            void serialize(Archive &ar, const unsigned int version){
                ar & indicators;
            }
};


#endif //MULGIFT_SKETCH_H
