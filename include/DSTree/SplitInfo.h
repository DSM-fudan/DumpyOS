//
// Created by wzy on 2021/8/7.
//

#ifndef MULGIFT_SPLITINFO_H
#define MULGIFT_SPLITINFO_H
#include "InsertedSeries.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>


class SplitInfo {
public:
    int splitFrom{};    // segment start point
    int splitTo{};      // segment end point
    INodeSegmentSplitPolicy* nodeSegmentSplitPolicy{};     // segment split policy

    friend class boost::serialization::access;
    template<class Archive>
            void serialize(Archive &ar, const unsigned int version){
                ar & splitFrom; ar & splitTo;
                ar & nodeSegmentSplitPolicy;
//                if(nodeSegmentSplitPolicy->getIndicatorSplitIdx() == 0)
//                    ar & (MeanNodeSegmentSplitPolicy*)nodeSegmentSplitPolicy;
//                else
//                    ar & (StdevNodeSegmentSplitPolicy*)nodeSegmentSplitPolicy;
            }

    bool routeToLeft(InsertedSeries* series) const {
        int idx = nodeSegmentSplitPolicy->getIndicatorSplitIdx();
        if(idx == 0)
            return series->getMean(splitFrom, splitTo) < nodeSegmentSplitPolicy->getIndicatorSplitValue();
        return series->getStdEv(splitFrom, splitTo) < nodeSegmentSplitPolicy->getIndicatorSplitValue();
    }

    SplitInfo(int from, int to, INodeSegmentSplitPolicy* policy){
        splitFrom = from;
        splitTo = to;
        nodeSegmentSplitPolicy = policy;
    }

    SplitInfo(){;}
};


#endif //MULGIFT_SPLITINFO_H
