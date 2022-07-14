//
// Created by wzy on 2021/8/7.
//

#ifndef MULGIFT_INODESEGMENTSPLITPOLICY_H
#define MULGIFT_INODESEGMENTSPLITPOLICY_H
#include <string>
#include "Sketch.h"


class INodeSegmentSplitPolicy {
public:
    virtual vector<Sketch> * split(Sketch &nodeSegmentSketch)=0;

    virtual int getIndicatorSplitIdx()=0;

    virtual double getIndicatorSplitValue()=0;

    virtual std::string getName() = 0;

};


#endif //MULGIFT_INODESEGMENTSPLITPOLICY_H
