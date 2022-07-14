//
// Created by wzy on 2021/8/7.
//

#include "../../include/DSTree/StdevNodeSegmentSplitPolicy.h"
#include <algorithm>

vector<Sketch> * StdevNodeSegmentSplitPolicy::split(Sketch &nodeSegmentSketch) {
    double max_stdev = nodeSegmentSketch.indicators[2];
    double min_stdev = nodeSegmentSketch.indicators[3];
    indicatorSplitValue = (float) ((max_stdev + min_stdev) / 2);  //the mean of stdeve value is split value

    auto* ret = new vector<Sketch>(2); //split into 2 node
    (*ret)[0].indicators.resize(4);
    (*ret)[1].indicators.resize(4);
    std::copy(nodeSegmentSketch.indicators.begin(), nodeSegmentSketch.indicators.begin() + 4, (*ret)[0].indicators.begin());
    std::copy(nodeSegmentSketch.indicators.begin(), nodeSegmentSketch.indicators.begin() + 4, (*ret)[1].indicators.begin());
    (*ret)[0].indicators[2] = indicatorSplitValue;
    (*ret)[1].indicators[3] = indicatorSplitValue;
    return ret;  //To change body of implemented methods use File | Settings | File Templates.
}

int StdevNodeSegmentSplitPolicy::getIndicatorSplitIdx() {
    return 1;
}   // 1 for standard deviation

double StdevNodeSegmentSplitPolicy::getIndicatorSplitValue() {
    return indicatorSplitValue;
}

StdevNodeSegmentSplitPolicy::StdevNodeSegmentSplitPolicy(double v){
    indicatorSplitValue = v;
}

StdevNodeSegmentSplitPolicy::StdevNodeSegmentSplitPolicy()= default;

string StdevNodeSegmentSplitPolicy::getName() {
    return "Stdev";
}
