//
// Created by wzy on 2021/8/7.
//

#ifndef MULGIFT_PQITEMINDEX_H
#define MULGIFT_PQITEMINDEX_H


class PqItemIndex {
public:
    int index;
    double dist;

    PqItemIndex(int _index, double _dist){
        index = _index;
        dist = _dist;
    }

    PqItemIndex(const PqItemIndex& a){
        index = a.index;
        dist = a.dist;
    }

    PqItemIndex() {
        index = dist = 0;
    }

    PqItemIndex& operator=(const PqItemIndex & a)= default;

//    @Override
//    public int compareTo(PqItemIndex o) {
//        return Double.compare(o.dist, dist);
//    }
};

struct PqItemIndexMaxHeap{
    bool operator()(const PqItemIndex* x, const PqItemIndex* y){
        return x->dist < y->dist;
    }
};

struct PqItemIndexComparator{
    bool operator()(const PqItemIndex* x, const PqItemIndex* y){
        return x->dist > y->dist;
    }
};


#endif //MULGIFT_PQITEMINDEX_H
