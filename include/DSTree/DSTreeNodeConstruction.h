//
// Created by wzy on 2021/8/7.
//

#ifndef MULGIFT_DSTREENODECONSTRUCTION_H
#define MULGIFT_DSTREENODECONSTRUCTION_H
#include "DSTreeNode.h"


class DSTreeNodeConstruction {
public:
    const static int initSegmentNum = 4;
    //    const static String indexPath = "data/expr_dstree_index/Index_%s_" + DSTreeNode.threshold + "/";
//    static string indexPath;

//    static void setIndexPath(const string &fileName);

    static double calcQoS(Sketch &sketch, int len);

    static vector<int> * calcPoints(int tsLength, int segmentNo);

    static void updateSketch(Sketch& nodeSegmentSketch, InsertedSeries* series, int fromIdx, int toIdx);

    static DSTreeNode *buildIndex(long series_num);

//    static DSTreeNode *buildIndex(const Vertex &v, const string &fileName, int parent_vertex_id = -1);
//
//    static DSTreeNode *buildIndex(const Vertex &v, int parent_vertex_id, float *tss);
};


#endif //MULGIFT_DSTREENODECONSTRUCTION_H
