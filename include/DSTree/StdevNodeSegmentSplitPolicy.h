//
// Created by wzy on 2021/8/7.
//

#ifndef MULGIFT_STDEVNODESEGMENTSPLITPOLICY_H
#define MULGIFT_STDEVNODESEGMENTSPLITPOLICY_H
#include "Sketch.h"
#include "INodeSegmentSplitPolicy.h"
#include <string>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

using namespace std;


class StdevNodeSegmentSplitPolicy: public INodeSegmentSplitPolicy{
public:
    StdevNodeSegmentSplitPolicy();

    explicit StdevNodeSegmentSplitPolicy(double v);

    double indicatorSplitValue{};

    friend class boost::serialization::access;
    template<class Archive>
            void serialize(Archive &ar, const unsigned int version){
        boost::serialization::void_cast_register<StdevNodeSegmentSplitPolicy, INodeSegmentSplitPolicy>();
                ar & indicatorSplitValue;
            }

    vector<Sketch> * split(Sketch &nodeSegmentSketch) override;

    int getIndicatorSplitIdx() override;

    double getIndicatorSplitValue() override;

    string getName() override;
};


#endif //MULGIFT_STDEVNODESEGMENTSPLITPOLICY_H
