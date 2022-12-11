//
// Created by wzy on 2021/8/7.
//

#ifndef MULGIFT_SAXUTIL_H
#define MULGIFT_SAXUTIL_H
#include <string>
#include <vector>

using namespace std;

typedef struct dequeue
{   int *dq;
    int size,capacity;
    int f,r;
} dequeue;

class SaxUtil {
public:
    static double* breakpoints;
    static float* breakpoints_f;
    static string breakpointsFile;
    static string SEPARATOR;
    static double *bp8;

    static int * generatePrefixMask();

    static double* readFromstring(string str);

    static float* readFromstringFloat(string str);

    static double* readDoubleFromFileAtOnce(const string& fileName);

    static float* readFloatFromFileAtOnce(const string& fileName);

    static double * paaFromTs(const float* ts, int tsLengthPerSegment, int segmentNum);

    static vector<int> * saxFromTs(float*ts, int tsLengthPerSegment, int segmentNum, int cardinality);

    static int findFirstGE(const double* array, int start, int length, double target);// satisfy condition: array[?] >= target  and the first one

    static int* invSaxFromSax(std::vector<int> *sax, int bitsCardinality, int segmentNum);

    /**
     * The MINDIST_PAA_iSAX function is used in kNN Search and exact search to find a baseline
     * distance between a PAA representation and a SAX representation.
     */
    static double LowerBound_Paa_iSax(const float *paa, std::vector<int> *sax);

    /**
     This function prints a sax record.
     */
    static void saxPrint(int* sax, int bits_cardinality, int segment_num);

    static void printBinary(long n, int size);

    static vector<unsigned short> * saxFromPaa(float *paa, int segmentNum, int cardinality);

    static int invSaxHeadFromSax(vector<int> *sax, int bitsCardinality, int segmentNum);

    static int invSaxHeadFromTs(const float *ts, int tsLengthPerSegment, int segmentNum);

    static int invSaxHead2FromTs(const float *ts, int tsLengthPerSegment, int segmentNum);

    static void saxFromTs(const float *ts, unsigned short *sax, int tsLengthPerSegment, int segmentNum, int cardinality);

    static double LowerBound_Paa_iSax(const float *paa, const unsigned short *sax, int bits_cardinality);

    static double  LowerBound_Paa_iSax_SIMD(const float *paa, unsigned char *sax_short, unsigned char *bits_cardinality_short);

    static double LowerBound_Paa_iSax_SIMD(const float *paa, const unsigned short *sax, const int* bits_cardinality, vector<int>&chosen_segs, int new_i);

    static double LowerBound_Paa_iSax_SIMD(const float *paa, const unsigned short *sax, const int* bits_cardinality) ;

        static double  minidist_paa_to_isax_rawa_SIMD(const float *paa, unsigned char*sax, unsigned char*sax_cardinalities);

    static double getMinDist1stLayer(const float *paa, int id);

    static int invSaxHead2FromSax(vector<int> *sax, int bitsCardinality, int segmentNum);

    static void generateSaxFile(const string &fn, const string &output);

//    static int invSaxHeadFromSax(const int *sax, int bitsCardinality, int segmentNum);

    static int invSaxHead2FromSax(const int *sax, int bitsCardinality, int segmentNum);

    static void id2Sax(int id, int *sax, int segment_num);

    static void id2Sax2(int id, unsigned short *sax, int segment_num);

    static int invSaxHeadkFromSax(const unsigned short *sax, int bitsCardinality, int segmentNum, int k);

    static void paaFromTs(const float *ts, float *paa, int tsLengthPerSegment, int segmentNum);

    static void generatePaaFile(const string &fn, const string &output);

    static int findFirstGE(const int *array, int start, int length, int target);

    static void getValueRange(int sax_single, int bits_cardinality, double *lb, double *ub);

    static int extendSax(float *paa, const int *bits_cardinality, vector<int> &segments);

    static double LowerBound_Paa_iSax(const double *paa, const int *sax, const int *bits_cardinality);

    static int invSaxHeadFromPaa(const float *paa, int tsLengthPerSegment, int segmentNum);

    static double LowerBound_Paa_iSax(const double *paa, const int *sax, const int *bits_cardinality, vector<int> &chosen_segs,
                               int new_id);

    static void invSaxPrintDec(int *invsax, int bc);

    static int extendSax(const unsigned short *sax, const int *bits_cardinality);

    static int getNewId(const float *paa, const float *split_line);

    static int getNewId(const float *paa, const float *split_line, vector<int> &segments);

    static int extendSax(const unsigned short *sax, const int *bits_cardinality, vector<int> &segments);

    static int invSaxHeadFromSax(const unsigned short *sax, int bitsCardinality, int segmentNum);

    static double
    LowerBound_Paa_iSax(const float *paa, const unsigned short *sax, const int *bits_cardinality,
                        vector<int> &chosen_segs,
                        int new_id);

    static double LowerBound_Paa_iSax(const float *paa, const unsigned short *sax, const int *bits_cardinality);

    static int extendSax(const unsigned short *sax, const int *bits_cardinality, vector<int> &segments,
                  const unsigned short *parent_sax);

    static void paaAndSaxFromTs(const float *ts, float *paa, unsigned short *sax, int tsLengthPerSegment, int segmentNum,
                         int cardinality);

    static double LowerBound_Paa_iSax(const float *paa, const unsigned short *sax);

    static double getMidLineFromSaxSymbolbc8(unsigned short symbol);

    static double minidist_paa_to_isax_DTW(const double *paaU, const double *paaL, const unsigned short *sax,
                                    const int *sax_cardinalities);

    static void lower_upper_lemire(const float *t, int len, int r, float *l, float *u);

    static float *paaFromTsFloat(const float *ts, int tsLengthPerSegment, int segmentNum);

    static void invSaxFromSax(vector<int> *sax, unsigned int *invSAX, int bitsCardinality, int segmentNum);

    static string invSax2String(unsigned int *invsax);

    static unsigned int *str2Invsax(string str);

    static string invSaxHeadKFromInvSax(string invsax, int layer);

    static void string2invSAX(string invsax, unsigned int *ret);

    static unsigned short *invSax2Sax(const unsigned *invsax, int unit_num);

    static double
    minidist_paa_to_isax_DTW(const double *paaU, const double *paaL, const unsigned short *sax, int bitsCardinality);

    static void
    extendSax(const unsigned short *parent_sax, const int *parent_bits_cardinality, vector<int> &chosen_segs,
              int new_id,
              unsigned short *sax, int *bits_cardinality);
};



#endif //MULGIFT_SAXUTIL_H
