//
// Created by wzy on 2021/8/7.
//

#ifndef MULGIFT_MATHUTIL_H
#define MULGIFT_MATHUTIL_H
#include <vector>
#include "../DataStructures/PqItemSeries.h"
using namespace std;
class MathUtil {
public:
    static double deviation(const float *ts, int start, int end);

    static double deviation(const float *ts, int start, int end, double dAve);

    static double avg(const double *timeSeries, int start, int end);

    static double avg(const float *timeSeries, int start, int end);

    static double avg(const float *timeSeries, int len);

    static double deviation(double *timeSeries, int start, int end);

    static double deviation(double *timeSeries, int len);

    static double deviation(const float *timeSeries, int len);

    static float *z_Normalize(const float *timeSeries, int len);

    static int split(const vector<int> &points, int minLength, int len, vector<int> &ret);

    static int bitDiffNum(int i, int j, int n);

    static int bitCount(int input);

    static int nChooseK(int n, int k);

    static void get1BitChangedNums(int x, int n, int *array, int start);

    static void get2BitsChangedNums(int x, int n, int *array, int start);

    static void get3BitsChangedNums(int x, int n, int *array, int start);

    static void get4BitsChangedNums(int x, int n, int *array, int start);

    static int bitDiffNum(int i, int j);

    static double errorRatio(vector<PqItemSeries *> &approx, vector<PqItemSeries *> &exact, int k);

    static double invertedErrorRatio(vector<PqItemSeries *> &approx, vector<PqItemSeries *> &exact, int k);

    static int *generateMask(int segno);

    static int generateMaskSettingKbits(const int *pos, int k, int len);
};


#endif //MULGIFT_MATHUTIL_H
