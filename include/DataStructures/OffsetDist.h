//
// Created by zywang on 2021/12/12.
//

#ifndef MULGIFT_OFFSETDIST_H
#define MULGIFT_OFFSETDIST_H
struct OffsetDist{
    long offset;
    double dist;

    OffsetDist(int o, double d){offset =  o; dist =d;}
};
struct OffsetDistMinHeap{
    bool operator()(const OffsetDist* x, const OffsetDist* y){
        return x->dist > y->dist;
    }
};
struct OffsetDistMaxHeap{
    bool operator()(const OffsetDist* x, const OffsetDist* y){
        return x->dist < y->dist;
    }
};
#endif //MULGIFT_OFFSETDIST_H
