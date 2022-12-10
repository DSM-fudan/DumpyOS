//
// Created by wzy on 2021/8/7.
//

#ifndef MULGIFT_TIMESERIESUTIL_H
#define MULGIFT_TIMESERIESUTIL_H
#include "../DataStructures/PqItemIndex.h"
#include "../DataStructures/PqItemSeries.h"
#include "../../include/DataStructures/IPGNode.h"
#include "../DSTree/DSTreeNode.h"
#include "../DataStructures/IPGNode.h"
#include "../Const.h"
#include <vector>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <bitset>


class TimeSeriesUtil {
public:
    static bool isSame(float* ts1, float* ts2, int, int);

    static bool isSame(const PqItemSeries& ts1, const PqItemSeries& ts2);

    static vector<PqItemSeries *> & intersectionTsSets(const vector<PqItemSeries *> *tsSet1, const vector<float *> *tsSet2);

    static int intersectionTsSetsCardinality(const vector<PqItemSeries *> &tsSet1, const vector<PqItemSeries *> &tsSet2);

    static vector<PqItemSeries> & intersectionTsSets(const vector<PqItemSeries>& tsSet1, const vector<PqItemSeries>& tsSet2);


    static string timeSeries2Line(float* timeSeries) {
        string sb;
        for (int i=0; i < Const::tsLength; ++i) {
            sb += formatFloatValue(timeSeries[i], 4) + " , ";
        }
        sb.append("\n");
        return sb;
    }

    static string timeSeries2Line(vector<float>* timeSeries) {
        string sb;
        for (int i=0; i < Const::tsLength; ++i) {
            sb += formatFloatValue((*timeSeries)[i], 4) + " , ";
        }
        sb.append("\n");
        return sb;
    }

    static string formatFloatValue(float val, int fixed) {
        ostringstream oss;
        oss << setprecision(fixed) << val;
        return oss.str();
    }


    static double* avgBySegments(float* timeSeries, const int* segments, int segment_number);

    static double* devBySegments(float* timeSeries, const int* segments, int segment_number);

    static double euclideanDist(const float* ts_1, const float* ts_2, int len);

    /**
     * using avg and std
     *
     * @param node
     * @param queryTs
     * @return
     * @throws IOException
     */

    template<typename ... Args>
    static string str_format(const string &format, Args... args);

    static void knnWithBsf(const DSTreeNode &node, InsertedSeries &queryTs, int k, vector<PqItemSeries *> &heap);

    static void knnWithBsfAndThreshold(const DSTreeNode &node, InsertedSeries &queryTs, int k, vector<PqItemSeries *> &heap,
                                       int threshold);

    static vector<PqItemSeries *> & knn(const DSTreeNode &node, InsertedSeries &q, int k);

    static double minDistEstimation(const DSTreeNode *node, InsertedSeries *queryTs);

    static double euclideanDist(const vector<float> *ts_1, const float *ts_2, int len);

    static double dtw(const float* A, const float* B, int len, int r, double bsf);

    static double euclideanDist(float *query_reordered, float *ts, int size, double bound, int *order);

    static float euclideanDist_SIMD(float * t, float * s, int size, float bound);

    static float euclideanDist_SIMD(float *query_reordered, float *ts, int size, double bound, int *order);

//    static void knnWithBsf(const DSTreeNode &node, InsertedSeries &queryTs, int k, vector<PqItemSeriesVector *> &heap);
//
//    static vector<PqItemSeriesVector *> & knnVector(const DSTreeNode &node, InsertedSeries &q, int k);
//
//    static int
//    intersectionTsSetsCardinality(const vector<PqItemSeriesVector *> *tsSet1, const vector<PqItemSeriesVector *> *tsSet2);
//
//    static vector<PqItemSeriesVector> &
//    intersectionTsSets(const vector<PqItemSeriesVector> &tsSet1, const vector<PqItemSeriesVector> &tsSet2);
//
//    static bool containsSeries(const vector<PqItemSeriesVector *> &tss, float *q);

    static double minDistEstimation(InsertedSeries &q, const vector<int> &points, int pointsLen, const vector<double> &mean,
                                    const vector<double> &stdev);

//    static int
//    intersectionTsSetsCardinality(const vector<PqItemSeriesVector> &tsSet1, const vector<PqItemSeriesVector> &tsSet2);
//
//    static int
//    intersectionTsSetsCardinality(const vector<PqItemSeriesVector *> &tsSet1,
//                                  const vector<PqItemSeriesVector *> &tsSet2);
//
//    static int intersectionTsSetsCardinality(const vector<PqItemSeriesVector *> *tsSet1, vector<float *> *tsSet2);

    static void heap_data_copy(vector<PqItemSeries *> &heap);

    static int intersectionTsSetsCardinality(const vector<PqItemSeries *> *tsSet1, vector<float *> *tsSet2);

    static bool containsSeries(const vector<PqItemSeries *> &tss, float *q);

    static double euclideanDist(const float *ts_1, const float *ts_2, int len, double bound);

    static double distanceEstimation(const IPGNode *n1, const IPGNode *n2);

    static void knnRow(const DSTreeNode &node, InsertedSeries &q, int k, vector<PqItemSeries *> &heap, int threshold);

    static double euclideanDist(const float *ts_1, const vector<float> &ts_2, int len, double bound);

    static double processSingleSeriesInKnn(float* queryTs, float*ts, int k, vector<PqItemSeries *> &heap, double bsfMin);



    //    static double maxDistEstimation(DSTreeNode node, float* queryTs) {
    //        double sum = 0;
    //        int* points = node.splitPoints;
    //        double[] avg = avgBySegments(queryTs, points);
    //        double[] stdDev = devBySegments(queryTs, points);
    //
    //        for (int id = 0; id < avg.length; id++) {
    //            double tempDist = 0;
    //            //using max std
    //            tempDist += Math.pow(stdDev[id] + node.nodeSegmentSketches[id].indicators[2], 2);
    //            //max of (avg to min mean and max mean)
    //            tempDist += Math.pow(Math.max(Math.abs(avg[id] - node.nodeSegmentSketches[id].indicators[0]), Math.abs(avg[id] - node.nodeSegmentSketches[id].indicators[1])), 2);
    //            sum += tempDist * node.getSegmentLength(id);
    //        }
    //        sum = Math.sqrt(sum);
    //        return sum;
    //    }
    //    static void timeSeries2Line(PqItemSeries[] tss){
    //        for(PqItemSeries ts:tss)
    //            System.out.println(timeSeries2Line(ts.ts));
    //    }


    static void
    prepareTsWithBoundsAndThreshold(InsertedSeries &queryTs, double bsfMin, const vector<int> &points, int pointsLen, int threshold,
                                    const vector<vector<double>> &means, const vector<vector<double>> &stdevs, vector<int> &res);

    static void prepareTsWithBounds(InsertedSeries &queryTs, double bsfMin, const vector<int> &points, int pointsLen,
                                    const vector<vector<double>> &means, const vector<vector<double>> &stdevs,
                                    vector<PqItemIndex *> &res);

    static void
    knnRow(const DSTreeNode &node, InsertedSeries &q, int k, vector<PqItemSeries *> &heap, int pointsLen, double bsfMin);

//    static void knnRow(const DSTreeNode &node, InsertedSeries &q, int k, vector<PqItemSeriesVector *> &heap, int pointsLen,
//                       double bsfMin);
//
//    static double
//    processSingleSeriesInKnn(float *queryTs, vector<float> *ts, int k, vector<PqItemSeriesVector *> &heap,
//                             double bsfMin);
//
//    static bool isSame(const PqItemSeriesVector &ts1, const PqItemSeriesVector &ts2);
//
//    static bool isSame(const PqItemSeriesVector *ts1, float *ts2);
//
//    static bool isSame(const PqItemSeriesVector *ts1, const PqItemSeriesVector *ts2);
//

    static bool isSame(const PqItemSeries *ts1, const float *ts2);


    static double dtw(const float *A, const vector<float> &B, int len, int r, double bsf);
};


#endif //MULGIFT_TIMESERIESUTIL_H
