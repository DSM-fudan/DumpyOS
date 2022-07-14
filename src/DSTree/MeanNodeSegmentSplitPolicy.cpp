//
// Created by wzy on 2021/8/7.
//

#include "../../include/DSTree/MeanNodeSegmentSplitPolicy.h"

// 新生成两个子节点的nodeSegmentSketch， 将均值范围对半切分
vector<Sketch> * MeanNodeSegmentSplitPolicy::split(Sketch &nodeSegmentSketch) {
    double max_mean = nodeSegmentSketch.indicators[0];
    double min_mean = nodeSegmentSketch.indicators[1];
    indicatorSplitValue = (max_mean + min_mean) / 2;  //the mean value is split value

    auto* ret = new vector<Sketch>(2);
    (*ret)[0].indicators.resize(4);
    (*ret)[1].indicators.resize(4);
    std::copy(nodeSegmentSketch.indicators.begin(), nodeSegmentSketch.indicators.begin() + 4, (*ret)[0].indicators.begin());
    std::copy(nodeSegmentSketch.indicators.begin(), nodeSegmentSketch.indicators.begin() + 4, (*ret)[1].indicators.begin());
    (*ret)[0].indicators[1] = indicatorSplitValue;
    (*ret)[1].indicators[0] = indicatorSplitValue;
    return ret;  //To change body of implemented methods use File | Settings | File Templates.
}

int MeanNodeSegmentSplitPolicy::getIndicatorSplitIdx() {
    return 0;
}   // 0 for mean

double MeanNodeSegmentSplitPolicy::getIndicatorSplitValue() {
    return indicatorSplitValue;
}

MeanNodeSegmentSplitPolicy::MeanNodeSegmentSplitPolicy(double v){
    indicatorSplitValue = v;
}

MeanNodeSegmentSplitPolicy::MeanNodeSegmentSplitPolicy()= default;

string MeanNodeSegmentSplitPolicy::getName() {
    return "Mean";
}