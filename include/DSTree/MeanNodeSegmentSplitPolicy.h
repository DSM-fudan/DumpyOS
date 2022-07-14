//
// Created by wzy on 2021/8/7.
//

#ifndef MULGIFT_MEANNODESEGMENTSPLITPOLICY_H
#define MULGIFT_MEANNODESEGMENTSPLITPOLICY_H
#include "INodeSegmentSplitPolicy.h"
#include <string>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

using namespace std;

class MeanNodeSegmentSplitPolicy : public INodeSegmentSplitPolicy{
public:
    MeanNodeSegmentSplitPolicy();

    explicit MeanNodeSegmentSplitPolicy(double v);

    double indicatorSplitValue{};

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version){
        boost::serialization::void_cast_register<MeanNodeSegmentSplitPolicy, INodeSegmentSplitPolicy>();
                ar & indicatorSplitValue;
            }

    int getIndicatorSplitIdx() override;

    double getIndicatorSplitValue() override;

    string getName() override;

    vector<Sketch> * split(Sketch &nodeSegmentSketch) override;
};


#endif //MULGIFT_MEANNODESEGMENTSPLITPOLICY_H
