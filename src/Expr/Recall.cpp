//
// Created by wzy on 2021/8/11.
//

#include "../../include/Expr/Recall.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Utils/TimeSeriesUtil.h"
#include "../../include/Const.h"
#include "../../include/Tardis/TardisTreeNode.h"
#include "../../include/Searchers/TardisApproxSearch.h"
#include "../../include/Searchers/iSAXSearcher.h"
#include "../../include/TAR/TARSearcher.h"
#include "../../include/DataStructures/IPGNode.h"
#include "../../include/Searchers/IPGApproxSearcher.h"
#include "../../include/Searchers/FADASSearcher.h"
#include "../../include/Searchers/DSTreeApproxSearcher.h"
#include "../../include/Const.h"
#include "../../include/TAR/TARGNode.h"
#include <string>
#include <iostream>
#include <cstdio>
#include <vector>
#include <ctime>
#include <chrono>
using namespace std;
//extern long determineCandidatesTime = 0, searchTime = 0, dataCopyTime = 0,
//            LBSeriesTime = 0, LBNodeTime=0, DSTreeProcessPqTime =0, IOBigTime = 0, distBigTime = 0, DSTreeNodeApproxTime = 0, heapBigTime = 0,
//            LBSeriesTime2 = 0, LBNodeTime2=0, DSTreeSearchTime2 =0, IOBigTime2 = 0, distBigTime2 = 0, DSTreeNodeApproxTime2 = 0, heapBigTime2 = 0, DSTreePreparePqTime2 =0,DSTreeProcessPqTime2 =0,
//            approxSearchTime = 0, exactSearchTime = 0, wallClockTime = 0, hashTableTime = 0, LBSeriesCnt = 0;
extern long approxSearchTimeDetails[]{0,0,0,0}, approxSearchUnits[]{0,0,0,0};
extern long LB_SERIES_TIME , HEAP_TIME ,
        LB_NODE_TIME_STAT , LB_NODE_CNT, LOADED_NODE_CNT, LOADED_PACK_CNT=0;
extern double DIST_CALC_TIME ,READ_TIME, PREPARE_TIME, SEARCH_TIME;

//void free_heap(vector<PqItemSeriesVector*>*heap){
//    for(auto* x:*heap)
//        delete x;
//    delete heap;
//}

void free_heap(vector<PqItemSeries*>*heap){
    for(auto* x:*heap)
        delete x;
    delete heap;
}

/// Data structure for sorting the query.
typedef struct q_index
{   double value;
    int    index;
} q_index;

int znorm_comp(const void *a, const void* b)
{
    q_index* x = (q_index*)a;
    q_index* y = (q_index*)b;

    //    return abs(y->value) - abs(x->value);

    if (fabsf(y->value) > fabsf(x->value) )
        return 1;
    else if (fabsf(y->value) == fabsf(x->value))
        return 0;
    else
        return -1;

}

void reorder_query(float * query_ts, float * query_ts_reordered, int * query_order)
{

    auto *q_tmp = (q_index*)malloc(sizeof(q_index) * Const::tsLength);
    int i;

    if( q_tmp == NULL )
        return ;

    for( i = 0 ; i < Const::tsLength ; i++ )
    {
        q_tmp[i].value = query_ts[i];
        q_tmp[i].index = i;
    }

    qsort(q_tmp, Const::tsLength, sizeof(q_index),znorm_comp);

    for( i=0; i<Const::tsLength; i++)
    {

        query_ts_reordered[i] = q_tmp[i].value;
        query_order[i] = q_tmp[i].index;
    }
    free(q_tmp);

}

//void Recall::doExpr(Graph* g){
//    int maxExprRound = 50;
//    int lastIndex= g->rowDataFileName.rfind(".bin");
//    string fn = g->rowDataFileName.substr(0, lastIndex) + "_query.bin";
//    int ks[]{10 ,50 ,100, 250, 500};
//    int thresholds[]{1000, 5000, 10000, 20000, 30000, 50000};
//    float *query;
//    FILE *f = fopen(fn.c_str(), "rb");
//    for(int threshold:thresholds){
//        for(int k:ks){
////            cout<<"k: " << k << endl;
////            cout<<"threshold: " << threshold<< endl;
//            int recallNums[maxExprRound];
//            for(int curRound = 0; curRound < maxExprRound; ++curRound){
//                //                    cout<<"Round : " + (curRound + 1));
//                query = FileUtil::readSeries(f);
//                vector<PqItemSeries*> *approxKnn = ApproximateSearcher::approxKnnSearchNew(g, query, k, threshold);
//                //                    cout<<"Approximate Search Finished.\nExact Search Start!");
//                vector<PqItemSeries*> *exactKnn = ExactSearcher::exactKnnSearch(g, query, k, nullptr);
//                //                    PqItemSeries[] exactKnn = readExactResult(g, cnt, k);
//                //                    cout<<"Exact Search Finished");
//                //                    if(TimeSeriesUtil::containsSeries(approxKnn, query)) cout<<"Approximate Search find the target series!");
//                //                    else cout<<"WARNING!  Approximate Search does not find the target series! ");
//                //                    if(TimeSeriesUtil::containsSeries(exactKnn, query)) cout<<"EXACT Search find the target series!");
//                //                    else cout<<"ERROR!  EXACT Search does not find the target series! ");
//                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(*approxKnn, *exactKnn);
//                free_heap(approxKnn);
//                free_heap(exactKnn);
//                delete query;
//            }
//
//            cout<<"------------------Experiment Finished--------------------" << endl;
//                            cout<<"k: " << k<<endl;
//                            cout<<"threshold: " << threshold<<endl;
//            int totalRecallNum = 0;
//            for(int temp:recallNums){
//                cout << temp << ",";
//                totalRecallNum += temp;
//            }
//            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
//        }
//        rewind(f);
//    }
//    fclose(f);
//}

vector<float*>* Recall::getResult(const string& fn, int queryNo, int k){
    auto * res = new vector<float*>(k);
    FILE *f= fopen(fn.c_str(), "rb");
    long off = (long)queryNo * Const::tsLengthBytes * Const::maxK;
    fseek(f, off, SEEK_SET);
    for(int i=0;i<k;++i)
    {
        auto *ts = new float [Const::tsLength];
        fread(ts, sizeof(float ), Const::tsLength, f);
        (*res)[i] = ts;
    }
    fclose(f);
    return res;
}

//void Recall::doExprWithRes(Graph* g){
//    int maxExprRound = 50;
//    int lastIndex= g->rowDataFileName.rfind(".bin"), index = g->rowDataFileName.find("generator");
//    string fn = g->rowDataFileName.substr(0, lastIndex) + "_query.bin";
//    string resFile = g->rowDataFileName.substr(0, index) + "results" + g->rowDataFileName.substr(index + 9);
//    cout << "result file is " << resFile <<endl;
//    int ks[]{10,50,100,500};
//    int thresholds[]{10000,50000};
//    float *query;
//    FILE *f = fopen(fn.c_str(), "rb");
//    for(int threshold:thresholds){
//        int _k = 0;
//        for(int k:ks){
//            int recallNums[maxExprRound];
//            cout<<"------------------Experiment--------------------" << endl;
//            cout<<"k: " << k << endl;
//            cout<<"threshold: " << threshold<< endl;
//
//            for(int curRound = 0; curRound < maxExprRound; ++curRound){
//                //                    cout<<"Round : " + (curRound + 1));
//                query = FileUtil::readSeries(f);
//                vector<PqItemSeries*> *approxKnn = ApproximateSearcher::approxKnnSearchNew(g, query, k, threshold);
//                vector<float*>* exactKnn = getResult(resFile, _k * maxExprRound + curRound, k);
//                //                    cout<<"Approximate Search Finished.\nExact Search Start!");
//                //                    PqItemSeries[] exactKnn = readExactResult(g, cnt, k);
//                //                    cout<<"Exact Search Finished");
//                //                    if(TimeSeriesUtil::containsSeries(approxKnn, query)) cout<<"Approximate Search find the target series!");
//                //                    else cout<<"WARNING!  Approximate Search does not find the target series! ");
//                //                    if(TimeSeriesUtil::containsSeries(exactKnn, query)) cout<<"EXACT Search find the target series!");
//                //                    else cout<<"ERROR!  EXACT Search does not find the target series! ");
//                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
//                cout << recallNums[curRound] << ",";
//                fflush(stdout);
//                free_heap(approxKnn);
//                for(int i=0;i<k;++i)
//                    delete[] (*exactKnn)[i];
//                delete exactKnn;
//                delete[] query;
//            }
//            ++_k;
////            cout<<"k: " << k<<endl;
////            cout<<"threshold: " << threshold<<endl;
//            int totalRecallNum = 0;
//            for(int temp:recallNums){
//                totalRecallNum += temp;
//            }
//            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
//        }
//        rewind(f);
//    }
//
//    cout << "small: " << approxSearchUnits[0] << " series, cost " << approxSearchTimeDetails[0] << " us. "<< endl;
//    cout << "middle: " << approxSearchUnits[1] << " series, cost " << approxSearchTimeDetails[1] << " us. "<< endl;
//    cout << "big: " << approxSearchUnits[2] << " series, cost " << approxSearchTimeDetails[2] << " us. "<< endl;
//    cout << "huge: " << approxSearchUnits[3] << " series, cost " << approxSearchTimeDetails[3] << " us. "<< endl;
//    fclose(f);
//}

extern int layer=0;

void analyzePrintSaxFADAS(vector<PqItemSeries*> *approxKnn, vector<float*>* exactKnn, float *query){
    cout << "Query:" << endl;
    vector<int>* sax = SaxUtil::saxFromTs(query, Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
    int *invsax = SaxUtil::invSaxFromSax(sax, Const::bitsCardinality, Const::segmentNum);
    SaxUtil::invSaxPrintDec(invsax, Const::bitsCardinality);

    cout << "Exact:" << endl;
    int i=1;
    for(float *ts:*exactKnn){
        cout << i++ << ":" << endl;
        sax = SaxUtil::saxFromTs(ts, Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
        invsax = SaxUtil::invSaxFromSax(sax, Const::bitsCardinality, Const::segmentNum);
        SaxUtil::invSaxPrintDec(invsax, Const::bitsCardinality);
    }

    cout << "Approx:" << endl;
    i=1;
    for(PqItemSeries *ts:*approxKnn){
        cout << i++ << ":" << endl;
        sax = SaxUtil::saxFromTs(ts->ts, Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
        invsax = SaxUtil::invSaxFromSax(sax, Const::bitsCardinality, Const::segmentNum);
        SaxUtil::invSaxPrintDec(invsax, Const::bitsCardinality);
    }
}

void analyzePrintSax(vector<PqItemSeries*> *approxKnn, vector<float*>* exactKnn, float *query){
    cout << "Query:" << endl;
    vector<int>* sax = SaxUtil::saxFromTs(query, Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
    int *invsax = SaxUtil::invSaxFromSax(sax, Const::bitsCardinality, Const::segmentNum);
    SaxUtil::saxPrint(invsax, Const::segmentNum, layer);

    cout << "Exact:" << endl;
    int i=1;
    for(float *ts:*exactKnn){
        cout << i++ << ":" << endl;
        sax = SaxUtil::saxFromTs(ts, Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
        invsax = SaxUtil::invSaxFromSax(sax, Const::bitsCardinality, Const::segmentNum);
        SaxUtil::saxPrint(invsax, Const::segmentNum, layer);
    }

    cout << "Approx:" << endl;
    i=1;
    for(PqItemSeries *ts:*approxKnn){
        cout << i++ << ":" << endl;
        sax = SaxUtil::saxFromTs(ts->ts, Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
        invsax = SaxUtil::invSaxFromSax(sax, Const::bitsCardinality, Const::segmentNum);
        SaxUtil::saxPrint(invsax, Const::segmentNum, layer);
    }
}

void Recall::doExprWithResInMemory(TardisTreeNode *root, const string &resFile, const string &queryFile) {
    int maxExprRound = 200;
    cout << "result file is " << resFile <<endl;
    int ks[]{1,3,5,10,25,50};
    int thresholds[]{10000};
    float *query;
    FILE *f = fopen(queryFile.c_str(), "rb");
    for(int threshold:thresholds){
        int _k = 0;
        for(int k:ks){
            int recallNums[maxExprRound];
            double error_ratio[maxExprRound];
            double inv_error_ratio[maxExprRound];
            cout<<"------------------Experiment--------------------" << endl;
            cout<<"k: " << k << endl;
            cout<<"threshold: " << threshold<< endl;

            for(int curRound = 0; curRound < maxExprRound; ++curRound){
                //                    cout<<"Round : " + (curRound + 1));
                query = FileUtil::readSeries(f);
                vector<PqItemSeries*> *approxKnn = TardisApproxSearch::tardisApproxKnnSearch(root, query, k, threshold);
                vector<float*>* exactKnn = getResult(resFile, curRound, k);
                vector<PqItemSeries*> exactKnn2;
                for(float *t: *exactKnn)
                    exactKnn2.push_back(new PqItemSeries(t, query));
                //                analyzePrintSax(approxKnn, exactKnn, query);
                //                vector<PqItemSeries*>exact;
                //                exact.reserve(k);
                //                for(float *ts:*exactKnn){
                //                    exact.push_back(new PqItemSeries(ts, query));
                //                }
                //                vector<PqItemSeries*> gt = ExactSearcher::groundTruthKnnSearch(TardisTreeNode::rowDataFileName, query, k);
                //                    cout<<"Approximate Search Finished.\nExact Search Start!");
                //                    PqItemSeries[] exactKnn = readExactResult(g, cnt, k);
                //                    cout<<"Exact Search Finished");
                //                    if(TimeSeriesUtil::containsSeries(approxKnn, query)) cout<<"Approximate Search find the target series!");
                //                    else cout<<"WARNING!  Approximate Search does not find the target series! ");
                //                    if(TimeSeriesUtil::containsSeries(exactKnn, query)) cout<<"EXACT Search find the target series!");
                //                    else cout<<"ERROR!  EXACT Search does not find the target series! ");
                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
                cout << recallNums[curRound] << ",";
                error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
                inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
                fflush(stdout);
                free_heap(approxKnn);
                for(int i=0;i<k;++i)
                    delete[] (*exactKnn)[i];
                delete exactKnn;
                delete[] query;
            }
            ++_k;
            rewind(f);
            //            cout<<"k: " << k<<endl;
            //            cout<<"threshold: " << threshold<<endl;
            int totalRecallNum = 0;
            for(int temp:recallNums){
                totalRecallNum += temp;
            }
            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
            double totalErrorRatio = 0;
            for(double _:error_ratio)   totalErrorRatio += _;
            cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
            double totalinvErrorRatio = 0;
            for(double _:inv_error_ratio)   totalinvErrorRatio += _;
            cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;

        }
        rewind(f);
    }

    //    cout << "small: " << approxSearchUnits[0] << " series, cost " << approxSearchTimeDetails[0] << " us. "<< endl;
    //    cout << "middle: " << approxSearchUnits[1] << " series, cost " << approxSearchTimeDetails[1] << " us. "<< endl;
    //    cout << "big: " << approxSearchUnits[2] << " series, cost " << approxSearchTimeDetails[2] << " us. "<< endl;
    //    cout << "huge: " << approxSearchUnits[3] << " series, cost " << approxSearchTimeDetails[3] << " us. "<< endl;
    fclose(f);
}

void Recall::doExprWithResTardis(TardisTreeNode *root, const string &resFile, const string &queryFile) {
    int maxExprRound = Const::query_num;
    cout << "result file is " << resFile <<endl;
    int ks[]{1,3,5,10,25,50};
    int thresholds[]{10000};
    float *query;
    FILE *f = fopen(queryFile.c_str(), "rb");
    for(int threshold:thresholds){
        int _k = 0;
        for(int k:ks){
            int recallNums[maxExprRound];
            double error_ratio[maxExprRound];
            double inv_error_ratio[maxExprRound];
            cout<<"------------------Experiment--------------------" << endl;
            cout<<"k: " << k << endl;
            cout<<"threshold: " << threshold<< endl;

            for(int curRound = 0; curRound < maxExprRound; ++curRound){
                //                    cout<<"Round : " + (curRound + 1));
                query = FileUtil::readSeries(f);
                vector<PqItemSeries*> *approxKnn = TardisApproxSearch::tardisApproxKnnIncSearchModel(root, query, k,1);
                vector<float*>* exactKnn = getResult(resFile, curRound, k);
                vector<PqItemSeries*> exactKnn2;
                for(float *t: *exactKnn)
                    exactKnn2.push_back(new PqItemSeries(t, query));
//                analyzePrintSax(approxKnn, exactKnn, query);
//                vector<PqItemSeries*>exact;
//                exact.reserve(k);
//                for(float *ts:*exactKnn){
//                    exact.push_back(new PqItemSeries(ts, query));
//                }
//                vector<PqItemSeries*> gt = ExactSearcher::groundTruthKnnSearch(TardisTreeNode::rowDataFileName, query, k);
                //                    cout<<"Approximate Search Finished.\nExact Search Start!");
                //                    PqItemSeries[] exactKnn = readExactResult(g, cnt, k);
                //                    cout<<"Exact Search Finished");
                //                    if(TimeSeriesUtil::containsSeries(approxKnn, query)) cout<<"Approximate Search find the target series!");
                //                    else cout<<"WARNING!  Approximate Search does not find the target series! ");
                //                    if(TimeSeriesUtil::containsSeries(exactKnn, query)) cout<<"EXACT Search find the target series!");
                //                    else cout<<"ERROR!  EXACT Search does not find the target series! ");
                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
                cout << recallNums[curRound] << ",";
                error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
                inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
                fflush(stdout);
                free_heap(approxKnn);
                for(int i=0;i<k;++i)
                    delete[] (*exactKnn)[i];
                delete exactKnn;
                delete[] query;
            }
            ++_k;
            //            cout<<"k: " << k<<endl;
            //            cout<<"threshold: " << threshold<<endl;
            int totalRecallNum = 0;
            for(int temp:recallNums){
                totalRecallNum += temp;
            }
            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
            double totalErrorRatio = 0;
            for(double _:error_ratio)   totalErrorRatio += _;
            cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
            double totalinvErrorRatio = 0;
            for(double _:inv_error_ratio)   totalinvErrorRatio += _;
            cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
            rewind(f);
        }

    }

//    cout << "small: " << approxSearchUnits[0] << " series, cost " << approxSearchTimeDetails[0] << " us. "<< endl;
//    cout << "middle: " << approxSearchUnits[1] << " series, cost " << approxSearchTimeDetails[1] << " us. "<< endl;
//    cout << "big: " << approxSearchUnits[2] << " series, cost " << approxSearchTimeDetails[2] << " us. "<< endl;
//    cout << "huge: " << approxSearchUnits[3] << " series, cost " << approxSearchTimeDetails[3] << " us. "<< endl;
    fclose(f);
}

extern int actual_node_number =0;
void Recall::doExprWithResIncTardis(TardisTreeNode *root, const string &resFile, const string &queryFile) {
    int maxExprRound = Const::query_num;
    cout << "result file is " << resFile <<endl;
    int k = 1;
    int node_nums[]{1,2,3,4,5,10,25};
    float *query;
    FILE *f = fopen(queryFile.c_str(), "rb");
    for(int node_num:node_nums){
        int recallNums[maxExprRound];
        int actual_rest_node_numbers[maxExprRound];
        double error_ratio[maxExprRound];
        double inv_error_ratio[maxExprRound];
        cout<<"------------------Experiment--------------------" << endl;
        cout<<"k: " << k << endl;
        cout<<"node number: " << node_num<< endl;

        for(int curRound = 0; curRound < maxExprRound; ++curRound){
            //                    cout<<"Round : " + (curRound + 1));
            query = FileUtil::readSeries(f);
            actual_node_number = node_num;
            vector<PqItemSeries*> *approxKnn = TardisApproxSearch::tardisApproxKnnIncSearchModel(root, query, k, node_num);
            vector<float*>* exactKnn = getResult(resFile, curRound, k);
            vector<PqItemSeries*> exactKnn2;
            for(float *t: *exactKnn)
                exactKnn2.push_back(new PqItemSeries(t, query));
            recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
            cout << recallNums[curRound] << ",";
            actual_rest_node_numbers[curRound] = actual_node_number;
            error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
            inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
            fflush(stdout);
            free_heap(approxKnn);
            for(int i=0;i<k;++i)
                delete[] (*exactKnn)[i];
            delete exactKnn;
            delete[] query;
        }
        //            cout<<"k: " << k<<endl;
        //            cout<<"threshold: " << threshold<<endl;
        int totalRecallNum = 0;
        for(int temp:recallNums){
            totalRecallNum += temp;
        }
        cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
        double totalErrorRatio = 0;
        for(double _:error_ratio)   totalErrorRatio += _;
        cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
        double totalinvErrorRatio = 0;
        for(double _:inv_error_ratio)   totalinvErrorRatio += _;
        cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
        for(int _:actual_rest_node_numbers) cout <<_<<",";
        cout <<endl;
        rewind(f);

    }

//    cout << "small: " << approxSearchUnits[0] << " series, cost " << approxSearchTimeDetails[0] << " us. "<< endl;
//    cout << "middle: " << approxSearchUnits[1] << " series, cost " << approxSearchTimeDetails[1] << " us. "<< endl;
//    cout << "big: " << approxSearchUnits[2] << " series, cost " << approxSearchTimeDetails[2] << " us. "<< endl;
//    cout << "huge: " << approxSearchUnits[3] << " series, cost " << approxSearchTimeDetails[3] << " us. "<< endl;
    fclose(f);
}

void Recall::progressiveSearchInMeomoryFADAS(TardisTreeNode *root, vector<vector<int>> *g) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);
    int k = Const::k;
    int series_num = Const::series_num;
    float thresholds[]{0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05};
//    float thresholds[]{ 0.002,0.003,0.004, 0.007,0.013,0.015,0.017, 0.02,0.03,0.04, 0.07,0.08};
//    float thresholds[]{ 0.0001};
    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    long offset = 0;
    fseek(f, offset * Const::tsLengthBytes, SEEK_SET);

    for(float threshold:thresholds){
        int recallNums[maxExprRound];
        int search_number[maxExprRound];
        int layers[maxExprRound];
        long duration[maxExprRound];
        double error_ratio[maxExprRound];
        double inv_error_ratio[maxExprRound];
        cout<<"------------------Experiment--------------------" << endl;
        cout<<"k: " << k << endl;
        int s_num = threshold * series_num;
        cout<<"threshold: " << s_num<< endl;

        for(int curRound = 0; curRound < maxExprRound; ++curRound){
            query = FileUtil::readSeries(f);
            actual_node_number = s_num;
            auto start = chrono::system_clock::now();
            vector<PqItemSeries*> *approxKnn = TardisApproxSearch::tardisApproxKnnSearch(root, query, k, s_num);
            auto end = chrono::system_clock::now();
            vector<float*>* exactKnn = getResult(Const::resfn, offset + curRound, k);
            vector<PqItemSeries*> exactKnn2;
            for(float *t: *exactKnn)
                exactKnn2.push_back(new PqItemSeries(t, query));

            duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count();
//            layers[curRound] = layer;
            recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
            search_number[curRound] = actual_node_number;
            error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
            inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
//                cout << curRound << ":"<<recallNums[curRound] << endl;
            cout << recallNums[curRound] << ",";
            fflush(stdout);
            free_heap(approxKnn);
            for(int i=0;i<k;++i)
                delete[] (*exactKnn)[i];
            delete exactKnn;
            delete[] query;
        }

        cout << fixed  << endl;
//        for(int _:layers)   cout << _ << ",";
//        cout << endl;
//            for(double _:error_ratio)   cout << _ << ",";
        int totalRecallNum = 0;
        for(int temp:recallNums)
            totalRecallNum += temp;
        cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
        double totalErrorRatio = 0;
        for(double _:error_ratio)   totalErrorRatio += _;
        cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
        double totalinvErrorRatio = 0;
        for(double _:inv_error_ratio)   totalinvErrorRatio += _;
        cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
        double total_duration = 0;
        for(long _:duration)    total_duration += _;
        total_duration /= (double ) maxExprRound;
        cout<<"Average duration is : " << total_duration << "us. "
            <<"And QPS = "<< 1000000.0 / total_duration <<endl;
        for(int _:search_number)    cout << _ <<",";
        cout << endl;
        rewind(f);

    }

    fclose(f);
}

void Recall::progressiveSearchInMeomoryDSTree(DSTreeNode *root) {
    int maxExprRound = Const::query_num;
    int k = Const::k;
    int series_num = Const::series_num;
    //    float thresholds[]{0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05};
    float thresholds[]{0.1};
    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    long offset = 0;
    fseek(f, offset * Const::tsLengthBytes, SEEK_SET);

    Const::logPrint("loading index.");
    DSTreeNode::load2Memory(root);

    for(float threshold:thresholds){
        int recallNums[maxExprRound];
        int search_number[maxExprRound];
        int layers[maxExprRound];
        long duration[maxExprRound];
        double error_ratio[maxExprRound];
        double inv_error_ratio[maxExprRound];
        cout<<"------------------Experiment--------------------" << endl;
        cout<<"k: " << k << endl;
        int s_num = threshold * series_num;
        cout<<"threshold: " << s_num<< endl;

        for(int curRound = 0; curRound < maxExprRound; ++curRound){
            query = FileUtil::readSeries(f);
            auto start = chrono::system_clock::now();
            vector<PqItemSeries*> *approxKnn = DSTreeApproxSearcher::approxKnnSearchWithThreshold(root, query, k, s_num);
            auto end = chrono::system_clock::now();
            vector<float*>* exactKnn = getResult(Const::resfn, offset + curRound, k);
            vector<PqItemSeries*> exactKnn2;
            for(float *t: *exactKnn)
                exactKnn2.push_back(new PqItemSeries(t, query));

            duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count();
//            layers[curRound] = layer;
            recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
            error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
            inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
//                cout << curRound << ":"<<recallNums[curRound] << endl;
            cout << recallNums[curRound] << ",";
            fflush(stdout);
            free_heap(approxKnn);
            for(int i=0;i<k;++i)
                delete[] (*exactKnn)[i];
            delete exactKnn;
            delete[] query;
        }

        cout << fixed  << endl;
//        for(int _:layers)   cout << _ << ",";
//        cout << endl;
//            for(double _:error_ratio)   cout << _ << ",";
        int totalRecallNum = 0;
        for(int temp:recallNums)
            totalRecallNum += temp;
        cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
        double totalErrorRatio = 0;
        for(double _:error_ratio)   totalErrorRatio += _;
        cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
        double totalinvErrorRatio = 0;
        for(double _:inv_error_ratio)   totalinvErrorRatio += _;
        cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
        double total_duration = 0;
        for(long _:duration)    total_duration += _;
        total_duration /= (double ) maxExprRound;
        cout<<"Average duration is : " << total_duration << "us. "
            <<"And QPS = "<< 1000000.0 / total_duration <<endl;
//        for(int _:search_number)    cout << _ <<",";
        cout << endl;
        rewind(f);

    }

    fclose(f);
}

void Recall::progressiveSearchInMeomoryiSAX(iSAXRoot *root) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);
    int k = Const::k;
    int series_num = Const::series_num;
//    float thresholds[]{0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1};
    float thresholds[]{ 0.002,0.003,0.004, 0.007,0.013,0.015,0.017, 0.02,0.03,0.04, 0.07,0.08};
//    float thresholds[]{};
    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    long offset = 0;

    for(float threshold:thresholds){
        fseek(f, offset * Const::tsLengthBytes, SEEK_SET);
        int recallNums[maxExprRound];
        int search_number[maxExprRound];
        int layers[maxExprRound];
        long duration[maxExprRound];
        vector<int>candidates;
        double error_ratio[maxExprRound];
        double inv_error_ratio[maxExprRound];
        cout<<"------------------Experiment--------------------" << endl;
        cout<<"k: " << k << endl;
        int s_num = threshold * series_num;
        cout<<"threshold: " << s_num<< endl;

        for(int curRound = 0; curRound < maxExprRound; ++curRound){
            query = FileUtil::readSeries(f);
            auto start = chrono::system_clock::now();
            vector<PqItemSeries*> *approxKnn = iSAXSearcher::iSAXApproxKnnSearch(root, query, k, s_num);
            auto end = chrono::system_clock::now();
            vector<float*>* exactKnn = getResult(Const::resfn, offset + curRound, k);
            vector<PqItemSeries*> exactKnn2;
            for(float *t: *exactKnn)
                exactKnn2.push_back(new PqItemSeries(t, query));

            duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count();
            layers[curRound] = layer;
            if(layer != 1)  candidates.push_back(curRound);
            recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
            error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
            inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
//                cout << curRound << ":"<<recallNums[curRound] << endl;
            cout << recallNums[curRound] << ",";
            fflush(stdout);
            free_heap(approxKnn);
            for(int i=0;i<k;++i)
                delete[] (*exactKnn)[i];
            delete exactKnn;
            delete[] query;
        }

        cout << fixed  << endl;
        for(int _:layers)   cout << _ << ",";
//        cout << endl;
//            for(double _:error_ratio)   cout << _ << ",";
        int totalRecallNum = 0;
        for(int temp:recallNums)
            totalRecallNum += temp;
        cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
        double totalErrorRatio = 0;
        for(double _:error_ratio)   totalErrorRatio += _;
        cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
        double totalinvErrorRatio = 0;
        for(double _:inv_error_ratio)   totalinvErrorRatio += _;
        cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
        double total_duration = 0;
        for(long _:duration)    total_duration += _;
        total_duration /= (double ) maxExprRound;
        cout<<"Average duration is : " << total_duration << "us. "
            <<"And QPS = "<< 1000000.0 / total_duration <<endl;
//        for(int _:search_number)    cout << _ <<",";
        cout << endl;
//        for(auto _: candidates) cout << _ <<",";
//        cout << endl;
        rewind(f);

    }

    fclose(f);
}

extern vector<IPGNode*>c_nodes(0);

void analyzePrintSaxIPG(vector<PqItemSeries*> *approxKnn, vector<float*>* exactKnn, float *query){
//    layer = 0;
//    for(IPGNode* node:c_nodes)
//        for(int _:node->bits_cardinality)   layer = max(layer, _);
    cout << "Query:" << endl;
    unsigned short q_sax[Const::segmentNum];
    float paa[Const::segmentNum];
    SaxUtil::saxFromTs(query, q_sax, Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
    SaxUtil::paaFromTs(query, paa, Const::tsLengthPerSegment, Const::segmentNum);
//    SaxUtil::saxPrint(q_sax, layer, Const::segmentNum);
    for(unsigned short & j : q_sax){
        j >>= (Const::bitsCardinality - layer);
    }
    for(int j=0;j<Const::segmentNum;++j){
        cout<<j << "\t";
        SaxUtil::printBinary(q_sax[j], layer);
        cout << "\t"<<paa[j] << "\t";
        cout << endl;
    }

    unsigned short sax[Const::segmentNum];
    cout << "Exact:" << fixed <<setprecision(3)<< endl;
    int i=1;
    for(float *ts:*exactKnn){
        cout << i++ << ":" << endl;
        SaxUtil::saxFromTs(ts, sax, Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
        SaxUtil::paaFromTs(ts, paa, Const::tsLengthPerSegment, Const::segmentNum);
//        SaxUtil::saxPrint(sax, layer, Const::segmentNum);
        for(unsigned short & j : sax){
            j >>= (Const::bitsCardinality - layer);
        }
        for(int j=0;j<Const::segmentNum;++j){
            cout<<j << "\t";
            SaxUtil::printBinary(sax[j], layer);
            cout << "\t"<<paa[j] << "\t";
            if(sax[j] != q_sax[j]) { cout << "\t --- "; SaxUtil::printBinary(q_sax[j], layer);}
            cout << endl;
        }
    }

//    cout << "Approx:" << endl;
//    cout << "Node:" <<endl;
//    for(IPGNode* node:c_nodes){
//        for(int j=0;j<Const::segmentNum;++j){
//            cout<<j << "\t";
//            SaxUtil::printBinary(node->sax[j], node->bits_cardinality[j]);
//            cout<<endl;
//        }
//    }
//    i=1;
//    for(PqItemSeries *ts:*approxKnn){
//        cout << i++ << ":" << endl;
//        SaxUtil::saxFromTs(ts->ts, sax, TimeSeries::tsLengthPerSegment, Const::segmentNum, TimeSeries::cardinality);
//        for(int & j : sax){
//            j >>= (TimeSeries::bitsCardinality - layer);
//        }
//        for(int j=0;j<Const::segmentNum;++j){
//            cout<<j << "\t";
//            SaxUtil::printBinary(sax[j], layer);
//            if(sax[j] != q_sax[j]) { cout << "\t --- "; SaxUtil::printBinary(q_sax[j], layer);}
//            cout << endl;
//        }
//    }
}

extern int _search_num = 0;
void Recall::doExprWithResFADASNonMat(IPGNode *root, const string &queryFile, const string &resFile) {
    int maxExprRound = 200;
    cout << "result file is " << resFile <<endl;
    int ks[]{1,3,5,10,25,50};
//    int ks[]{1};
    int thresholds[]{10000};
    float *query;
    FILE *f = fopen(queryFile.c_str(), "rb");
    for(int threshold:thresholds){
        int _k = 0;
        for(int k:ks){
            int recallNums[maxExprRound];
            int search_number[maxExprRound];
            int layers[maxExprRound];
            double error_ratio[maxExprRound];
            double inv_error_ratio[maxExprRound];
            cout<<"------------------Experiment--------------------" << endl;
            cout<<"k: " << k << endl;
            cout<<"threshold: " << threshold<< endl;

            for(int curRound = 0; curRound < maxExprRound; ++curRound){
                //                    cout<<"Round : " + (curRound + 1));
                c_nodes.clear();
                _search_num = 0;
                query = FileUtil::readSeries(f);
                vector<PqItemSeries*> *approxKnn = IPGApproxSearcher::approxKnnSearchModel(root, query, k, threshold);
                vector<float*>* exactKnn = getResult(resFile, curRound, k);
                vector<PqItemSeries*> exactKnn2;
                for(float *t: *exactKnn) {
//                    cout << TimeSeriesUtil::timeSeries2Line(t)<<endl;
                    exactKnn2.push_back(new PqItemSeries(t, query));
                }

                layers[curRound] = layer;
                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
                search_number[curRound] = _search_num;
                error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
                inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
                cout << recallNums[curRound] << ",";
                fflush(stdout);
//                analyzePrintSaxIPG(approxKnn, exactKnn, query);
                free_heap(approxKnn);
                for(int i=0;i<k;++i)
                    delete[] (*exactKnn)[i];
                delete exactKnn;
                delete[] query;
            }
            ++_k;
            cout << fixed  << endl;
            for(int _:layers)   cout << _ << ",";
            cout << endl;
//            for(double _:error_ratio)   cout << _ << ",";
            int totalRecallNum = 0;
            for(int temp:recallNums)
                totalRecallNum += temp;
            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
            double totalErrorRatio = 0;
            for(double _:error_ratio)   totalErrorRatio += _;
            cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
            double totalinvErrorRatio = 0;
            for(double _:inv_error_ratio)   totalinvErrorRatio += _;
            cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
//            for(int _:search_number)    cout << _ <<",";
            cout << endl;
            rewind(f);
        }

    }

    //    cout << "small: " << approxSearchUnits[0] << " series, cost " << approxSearchTimeDetails[0] << " us. "<< endl;
    //    cout << "middle: " << approxSearchUnits[1] << " series, cost " << approxSearchTimeDetails[1] << " us. "<< endl;
    //    cout << "big: " << approxSearchUnits[2] << " series, cost " << approxSearchTimeDetails[2] << " us. "<< endl;
    //    cout << "huge: " << approxSearchUnits[3] << " series, cost " << approxSearchTimeDetails[3] << " us. "<< endl;
    fclose(f);
}

void Recall::doExprWithResIPGDynamicPaaMu(IPGNode *root, const string &queryFile, const string &resFile) {
    int maxExprRound = 200;
    cout << "result file is " << resFile <<endl;
    int ks[]{1,3,5,10,25,50};
//    int ks[]{1};
    int thresholds[]{10000};
    float *query;
    FILE *f = fopen(queryFile.c_str(), "rb");
    for(int threshold:thresholds){
        int _k = 0;
        for(int k:ks){
            int recallNums[maxExprRound];
            int search_number[maxExprRound];
            int layers[maxExprRound];
            double error_ratio[maxExprRound];
            double inv_error_ratio[maxExprRound];
            cout<<"------------------Experiment--------------------" << endl;
            cout<<"k: " << k << endl;
            cout<<"threshold: " << threshold<< endl;

            for(int curRound = 0; curRound < maxExprRound; ++curRound){
                //                    cout<<"Round : " + (curRound + 1));
                c_nodes.clear();
                _search_num = 0;
                query = FileUtil::readSeries(f);
                vector<PqItemSeries*> *approxKnn = IPGApproxSearcher::approxKnnSearchModelPaaMuAsSplit(root, query, k);
                vector<float*>* exactKnn = getResult(resFile, curRound, k);
                vector<PqItemSeries*> exactKnn2;
                for(float *t: *exactKnn) {
//                    cout << TimeSeriesUtil::timeSeries2Line(t)<<endl;
                    exactKnn2.push_back(new PqItemSeries(t, query));
                }

                layers[curRound] = layer;
                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
                search_number[curRound] = _search_num;
                error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
                inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
                cout << recallNums[curRound] << ",";
                fflush(stdout);
//                analyzePrintSaxIPG(approxKnn, exactKnn, query);
                free_heap(approxKnn);
                for(int i=0;i<k;++i)
                    delete[] (*exactKnn)[i];
                delete exactKnn;
                delete[] query;
            }
            ++_k;
            cout << fixed  << endl;
            for(int _:layers)   cout << _ << ",";
            cout << endl;
//            for(double _:error_ratio)   cout << _ << ",";
            int totalRecallNum = 0;
            for(int temp:recallNums)
                totalRecallNum += temp;
            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
            double totalErrorRatio = 0;
            for(double _:error_ratio)   totalErrorRatio += _;
            cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
            double totalinvErrorRatio = 0;
            for(double _:inv_error_ratio)   totalinvErrorRatio += _;
            cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
//            for(int _:search_number)    cout << _ <<",";
            cout << endl;
            rewind(f);
        }

    }

    //    cout << "small: " << approxSearchUnits[0] << " series, cost " << approxSearchTimeDetails[0] << " us. "<< endl;
    //    cout << "middle: " << approxSearchUnits[1] << " series, cost " << approxSearchTimeDetails[1] << " us. "<< endl;
    //    cout << "big: " << approxSearchUnits[2] << " series, cost " << approxSearchTimeDetails[2] << " us. "<< endl;
    //    cout << "huge: " << approxSearchUnits[3] << " series, cost " << approxSearchTimeDetails[3] << " us. "<< endl;
    fclose(f);
}

void Recall::doExprWithResFADAS(FADASNode *root, vector<vector<int>> *g, const string &index_dir) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);
    int ks[]{1,3,5,10,25,50};
    float query_reordered[Const::tsLength];
    int ordering[Const::tsLength];

//    int ks[]{25};
    int thresholds[]{10000};
    cout << fixed  << setprecision(3) << endl;
    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    long offset = 0;
    fseek(f, offset * Const::tsLengthBytes, SEEK_SET);
    for(int threshold:thresholds){
        int _k = 0;
        for(int k:ks){
            int recallNums[maxExprRound];
            int search_number[maxExprRound];
            int layers[maxExprRound];
            long duration[maxExprRound];
            double error_ratio[maxExprRound];
            double inv_error_ratio[maxExprRound];
            cout<<"------------------Experiment--------------------" << endl;
            cout<<"k: " << k << endl;
            cout<<"threshold: " << threshold<< endl;

            for(int curRound = 0; curRound < maxExprRound; ++curRound){
                //                    cout<<"Round : " + (curRound + 1));
                c_nodes.clear();
                _search_num = 0;
                query = FileUtil::readSeries(f);
                reorder_query(query, query_reordered, ordering);
                auto start = chrono::system_clock::now();
                vector<PqItemSeries*> *approxKnn = FADASSearcher::approxSearch(root, query, k, g, index_dir, query_reordered,
                                                                               ordering);
                auto end = chrono::system_clock::now();
//                for(int i=0;i<256;++i)
//                    cout << (*approxKnn)[0]->ts[i] <<",";
                vector<float*>* exactKnn = getResult(Const::resfn, offset + curRound, k);
                vector<PqItemSeries*> exactKnn2;
                for(float *t: *exactKnn)
                    exactKnn2.push_back(new PqItemSeries(t, query));

                duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count();
                layers[curRound] = layer;
                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
                search_number[curRound] = _search_num;
                error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
                inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
//                cout << curRound << ":"<<recallNums[curRound] << endl;
                cout << recallNums[curRound] << "," ;
                fflush(stdout);
                free_heap(approxKnn);
                for(int i=0;i<k;++i)
                    delete[] (*exactKnn)[i];
                delete exactKnn;
                delete[] query;
            }
            ++_k;
            cout << endl;
            for(int _:layers)   cout << _ << ",";
            cout << endl;
//            for(double _:error_ratio)   cout << _ << ",";
            int totalRecallNum = 0;
            for(int temp:recallNums)
                totalRecallNum += temp;
            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
            double totalErrorRatio = 0;
            for(double _:error_ratio)   totalErrorRatio += _;
            cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
//            double totalinvErrorRatio = 0;
//            for(double _:inv_error_ratio)   totalinvErrorRatio += _;
//            cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
            double total_duration = 0;
            for(long _:duration)    total_duration += _;
            total_duration /= (double ) maxExprRound;
            cout<<"Average duration is : " << total_duration << "us. "
              <<"And QPS = "<< 1000000.0 / total_duration <<endl;
            for(int _:search_number)    cout << _ <<",";
            cout << endl;
            rewind(f);
        }
    }

    fclose(f);
}

void Recall::ngSearchDumpy(FADASNode *root, vector<vector<int>> *g) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);
    int nprobes[]{1,3,5,10, 20, 35, 50 ,75, 100, 130};
//    int nprobes[]{50};
//    int ks[]{25};
    int thresholds[]{10000};
    root->assignLeafNum();
    float *query;
    float query_ts_reordered[Const::tsLength];
    int ordering[Const::tsLength];
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    long offset = 0;
    fseek(f, offset * Const::tsLengthBytes, SEEK_SET);
    for(int threshold:thresholds){
        int _k = 0;
        for(int probe:nprobes){
            int recallNums[maxExprRound];
            int search_number[maxExprRound];
            int layers[maxExprRound];
            long duration[maxExprRound];
            double error_ratio[maxExprRound];
            double inv_error_ratio[maxExprRound];
            cout<<"------------------Experiment--------------------" << endl;
            cout<<"nprobe: " << probe << endl;
            cout<<"threshold: " << threshold<< endl;

            for(int curRound = 0; curRound < maxExprRound; ++curRound){
                //                    cout<<"Round : " + (curRound + 1));
                c_nodes.clear();
                _search_num = 0;
                query = FileUtil::readSeries(f);
                reorder_query(query, query_ts_reordered, ordering);
                auto start = chrono::system_clock::now();
                vector<PqItemSeries*> *approxKnn = FADASSearcher::ngSearch(root, query, query_ts_reordered, ordering, Const::k,g, probe);
                auto end = chrono::system_clock::now();
//                for(int i=0;i<256;++i)
//                    cout << (*approxKnn)[0]->ts[i] <<",";
                vector<float*>* exactKnn = getResult(Const::resfn, offset + curRound, Const::k);
                vector<PqItemSeries*> exactKnn2;
                for(float *t: *exactKnn)
                    exactKnn2.push_back(new PqItemSeries(t, query));

//                layer = 4;
//                analyzePrintSax(approxKnn,exactKnn, query);
                duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count();
                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
                search_number[curRound] = _search_num;
                error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, Const::k);
                inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, Const::k);
//                cout << curRound << ":"<<recallNums[curRound] << endl;
                cout << recallNums[curRound] << ",";
                fflush(stdout);
                free_heap(approxKnn);
                for(int i=0;i<Const::k;++i)
                    delete[] (*exactKnn)[i];
                delete exactKnn;
                delete[] query;
            }
            ++_k;
            cout << fixed  << endl;
//            for(int _:layers)   cout << _ << ",";
//            cout << endl;
//            for(double _:error_ratio)   cout << _ << ",";
            int totalRecallNum = 0;
            for(int temp:recallNums)
                totalRecallNum += temp;
            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * Const::k)<< endl;
            double totalErrorRatio = 0;
            for(double _:error_ratio)   totalErrorRatio += _;
            cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
//            double totalinvErrorRatio = 0;
//            for(double _:inv_error_ratio)   totalinvErrorRatio += _;
//            cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
            double total_duration = 0;
            for(long _:duration)    total_duration += _;
            total_duration /= (double ) maxExprRound;
            cout<<"Average duration is : " << total_duration/1000.0 << "ms. "
                <<"And QPM = "<< 60000000.0 / total_duration <<endl;
            for(int _:search_number)    cout << _ <<",";
            cout << endl;
            rewind(f);
        }
    }

    fclose(f);
}

void Recall::ngSearchDumpyFuzzy(FADASNode *root) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);
    int nprobes[]{1,3,5,10, 20, 35, 50 ,75, 100, 150, 200};
//    int nprobes[]{50};
//    int ks[]{25};
    int thresholds[]{10000};
    root->assignLeafNum();
    float *query;
    float query_ts_reordered[Const::tsLength];
    int ordering[Const::tsLength];
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    long offset = 0;
    fseek(f, offset * Const::tsLengthBytes, SEEK_SET);
    for(int threshold:thresholds){
        int _k = 0;
        for(int probe:nprobes){
            int recallNums[maxExprRound];
            int search_number[maxExprRound];
            int layers[maxExprRound];
            long duration[maxExprRound];
            double error_ratio[maxExprRound];
            double inv_error_ratio[maxExprRound];
            cout<<"------------------Experiment--------------------" << endl;
            cout<<"nprobe: " << probe << endl;
            cout<<"threshold: " << threshold<< endl;

            for(int curRound = 0; curRound < maxExprRound; ++curRound){
                //                    cout<<"Round : " + (curRound + 1));
                c_nodes.clear();
                _search_num = 0;
                query = FileUtil::readSeries(f);
                reorder_query(query, query_ts_reordered, ordering);
                auto start = chrono::system_clock::now();
                vector<PqItemSeries*> *approxKnn = FADASSearcher::ngSearchFuzzy(root, query, query_ts_reordered, ordering, Const::k, probe);
                auto end = chrono::system_clock::now();
//                for(int i=0;i<256;++i)
//                    cout << (*approxKnn)[0]->ts[i] <<",";
                vector<float*>* exactKnn = getResult(Const::resfn, offset + curRound, Const::k);
                vector<PqItemSeries*> exactKnn2;
                for(float *t: *exactKnn)
                    exactKnn2.push_back(new PqItemSeries(t, query));

//                layer = 4;
//                analyzePrintSax(approxKnn,exactKnn, query);
                duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count();
                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
                search_number[curRound] = _search_num;
                error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, Const::k);
                inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, Const::k);
//                cout << curRound << ":"<<recallNums[curRound] << endl;
                cout << recallNums[curRound] << ",";
                fflush(stdout);
                free_heap(approxKnn);
                for(int i=0;i<Const::k;++i)
                    delete[] (*exactKnn)[i];
                delete exactKnn;
                delete[] query;
            }
            ++_k;
            cout << fixed  << endl;
//            for(int _:layers)   cout << _ << ",";
//            cout << endl;
//            for(double _:error_ratio)   cout << _ << ",";
            int totalRecallNum = 0;
            for(int temp:recallNums)
                totalRecallNum += temp;
            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * Const::k)<< endl;
            double totalErrorRatio = 0;
            for(double _:error_ratio)   totalErrorRatio += _;
            cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
//            double totalinvErrorRatio = 0;
//            for(double _:inv_error_ratio)   totalinvErrorRatio += _;
//            cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
            double total_duration = 0;
            for(long _:duration)    total_duration += _;
            total_duration /= (double ) maxExprRound;
            cout<<"Average duration is : " << total_duration/1000.0 << "ms. "
                <<"And QPM = "<< 60000000.0 / total_duration <<endl;
            for(int _:search_number)    cout << _ <<",";
            cout << endl;
            rewind(f);
        }
    }

    fclose(f);
}

extern int __layer=0, nrest = 0;
void Recall::doExprWithResIncFADAS(FADASNode *root, vector<vector<int>> *g, const string &index_dir) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);
    int k = Const::k;
//    int ks[]{10};
    int node_nums[]{1,2,3,4,5, 10, 25};
//    int node_nums[]{50};
    float *query;
    float query_reordered[Const::tsLength];
    int ordering[Const::tsLength];
    root->assignLeafNum();
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    long offset = 0;
    fseek(f, offset * Const::tsLengthBytes, SEEK_SET);
    for(int node_num: node_nums){
        int recallNums[maxExprRound];
        int search_number[maxExprRound];
        int layers[maxExprRound], rest_nodes[maxExprRound];
        long duration[maxExprRound];
        double error_ratio[maxExprRound];
        double inv_error_ratio[maxExprRound];
        cout<<"------------------Experiment--------------------" << endl;
        cout<<"k: " << k << endl;
        cout<<"node number: " << node_num<< endl;

        for(int curRound = 0; curRound < maxExprRound; ++curRound){
            //                    cout<<"Round : " + (curRound + 1));
            c_nodes.clear();
            _search_num = 0;
            query = FileUtil::readSeries(f);
            reorder_query(query, query_reordered, ordering);
            auto start = chrono::system_clock::now();
            vector<PqItemSeries*> *approxKnn = FADASSearcher::approxIncSearch(root, query, k, index_dir, node_num,
                                                                              query_reordered, ordering, g);
            auto end = chrono::system_clock::now();
//                for(int i=0;i<256;++i)
//                    cout << (*approxKnn)[0]->ts[i] <<",";
            vector<float*>* exactKnn = getResult(Const::resfn, offset + curRound, k);
            vector<PqItemSeries*> exactKnn2;
            for(float *t: *exactKnn)
                exactKnn2.push_back(new PqItemSeries(t, query));

            duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count();
            layers[curRound] = __layer;
            rest_nodes[curRound] = nrest;
            recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
            search_number[curRound] = _search_num;
            error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
            inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
//                cout << curRound << ":"<<recallNums[curRound] << endl;
            cout << recallNums[curRound] << ",";
            fflush(stdout);
            free_heap(approxKnn);
            for(int i=0;i<k;++i)
                delete[] (*exactKnn)[i];
            delete exactKnn;
            delete[] query;
        }
        cout << fixed  << endl;
        for(int _:layers)   cout << _ << ",";
        cout << endl;
        for(int _:rest_nodes)   cout << _ << ",";
        cout << endl;
//            for(double _:error_ratio)   cout << _ << ",";
        int totalRecallNum = 0;
        for(int temp:recallNums)
            totalRecallNum += temp;
        cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
        double totalErrorRatio = 0;
        for(double _:error_ratio)   totalErrorRatio += _;
        cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
        double totalinvErrorRatio = 0;
        for(double _:inv_error_ratio)   totalinvErrorRatio += _;
        cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
        double total_duration = 0;
        for(long _:duration)    total_duration += _;
        total_duration /= (double ) maxExprRound;
        cout<<"Average duration is : " << total_duration / 1000.0 << "us. "
            <<"And QPM = "<< 60000000.0 / total_duration <<endl;
        for(int _:search_number)    cout << _ <<",";
        cout << endl;
        rewind(f);

    }

    fclose(f);
}

void Recall::doExprWithResFADASDTW(FADASNode *root, vector<vector<int>> *g, const string &index_dir) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);
    int ks[]{1,3,5,10,25,50};
//    int ks[]{10};
    int thresholds[]{10000};
    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    long offset = 0;
    fseek(f, offset * Const::tsLengthBytes, SEEK_SET);
    for(int threshold:thresholds){
        int _k = 0;
        for(int k:ks){
            int recallNums[maxExprRound];
            int search_number[maxExprRound];
            int layers[maxExprRound];
            long duration[maxExprRound];
            double error_ratio[maxExprRound];
            double inv_error_ratio[maxExprRound];
            cout<<"------------------Experiment--------------------" << endl;
            cout<<"k: " << k << endl;
            cout<<"threshold: " << threshold<< endl;

            for(int curRound = 0; curRound < maxExprRound; ++curRound){
                //                    cout<<"Round : " + (curRound + 1));
                c_nodes.clear();
                _search_num = 0;
                query = FileUtil::readSeries(f);
                auto start = chrono::system_clock::now();
                vector<PqItemSeries*> *approxKnn = FADASSearcher::approxSearchDTW(root, query, k, g, index_dir);
                auto end = chrono::system_clock::now();
//                for(int i=0;i<256;++i)
//                    cout << (*approxKnn)[0]->ts[i] <<",";
                vector<float*>* exactKnn = getResult(Const::resfn, offset + curRound, k);
                vector<PqItemSeries*> exactKnn2;
                for(float *t: *exactKnn) {
                    double dist = TimeSeriesUtil::dtw(t, query, Const::tsLength, Const::dtw_window_size, numeric_limits<double>::max());
                    exactKnn2.push_back(new PqItemSeries(t, dist));
                }

                duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count();
                layers[curRound] = layer;
                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
                search_number[curRound] = _search_num;
                error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
                inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
//                cout << curRound << ":"<<recallNums[curRound] << endl;
                cout << recallNums[curRound] << ",";
                fflush(stdout);
                free_heap(approxKnn);
                for(int i=0;i<k;++i)
                    delete[] (*exactKnn)[i];
                delete exactKnn;
                delete[] query;
            }
            ++_k;
            cout << fixed  << endl;
            for(int _:layers)   cout << _ << ",";
            cout << endl;
//            for(double _:error_ratio)   cout << _ << ",";
            int totalRecallNum = 0;
            for(int temp:recallNums)
                totalRecallNum += temp;
            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
            double totalErrorRatio = 0;
            for(double _:error_ratio)   totalErrorRatio += _;
            cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
            double totalinvErrorRatio = 0;
            for(double _:inv_error_ratio)   totalinvErrorRatio += _;
            cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
            double total_duration = 0;
            for(long _:duration)    total_duration += _;
            total_duration /= (double ) maxExprRound;
            cout<<"Average duration is : " << total_duration << "us. "
                <<"And QPS = "<< 1000000.0 / total_duration <<endl;
            for(int _:search_number)    cout << _ <<",";
            cout << endl;
            rewind(f);
        }
    }

    fclose(f);
}

void Recall::doExprWithResIncFADASDTW(FADASNode *root, vector<vector<int>> *g, const string &index_dir) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);
    int k = 1;
//    int ks[]{10};
    int node_nums[]{1,2,3,4,5, 10, 25};
    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    long offset = 0;
    fseek(f, offset * Const::tsLengthBytes, SEEK_SET);
    for(int node_num: node_nums){
        int recallNums[maxExprRound];
        int search_number[maxExprRound];
        int layers[maxExprRound];
        long duration[maxExprRound];
        double error_ratio[maxExprRound];
        double inv_error_ratio[maxExprRound];
        cout<<"------------------Experiment--------------------" << endl;
        cout<<"k: " << k << endl;
        cout<<"node number: " << node_num<< endl;

        for(int curRound = 0; curRound < maxExprRound; ++curRound){
            //                    cout<<"Round : " + (curRound + 1));
            c_nodes.clear();
            _search_num = 0;
            query = FileUtil::readSeries(f);
            auto start = chrono::system_clock::now();
            vector<PqItemSeries*> *approxKnn = FADASSearcher::approxIncSearchDTW(root, query, k, index_dir, node_num);
            auto end = chrono::system_clock::now();
//                for(int i=0;i<256;++i)
//                    cout << (*approxKnn)[0]->ts[i] <<",";
            vector<float*>* exactKnn = getResult(Const::dtwresfn, offset + curRound, k);
            vector<PqItemSeries*> exactKnn2;
            for(float *t: *exactKnn) {
                double dist = TimeSeriesUtil::dtw(t, query, Const::tsLength, Const::dtw_window_size, numeric_limits<double>::max());
                exactKnn2.push_back(new PqItemSeries(t, dist));
            }

            duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count();
            layers[curRound] = layer;
            recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
            search_number[curRound] = _search_num;
            error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
            inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
//                cout << curRound << ":"<<recallNums[curRound] << endl;
            cout << recallNums[curRound] << ",";
            fflush(stdout);
            free_heap(approxKnn);
            for(int i=0;i<k;++i)
                delete[] (*exactKnn)[i];
            delete exactKnn;
            delete[] query;
        }
        cout << fixed  << endl;
        for(int _:layers)   cout << _ << ",";
        cout << endl;
//            for(double _:error_ratio)   cout << _ << ",";
        int totalRecallNum = 0;
        for(int temp:recallNums)
            totalRecallNum += temp;
        cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
        double totalErrorRatio = 0;
        for(double _:error_ratio)   totalErrorRatio += _;
        cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
        double totalinvErrorRatio = 0;
        for(double _:inv_error_ratio)   totalinvErrorRatio += _;
        cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
        double total_duration = 0;
        for(long _:duration)    total_duration += _;
        total_duration /= (double ) maxExprRound;
        cout<<"Average duration is : " << total_duration << "us. "
            <<"And QPS = "<< 1000000.0 / total_duration <<endl;
        for(int _:search_number)    cout << _ <<",";
        cout << endl;
        rewind(f);

    }

    fclose(f);
}

void Recall::doExprWithResIncFADASFuzzy(FADASNode *root, vector<vector<int>> *g, const string &index_dir) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);
    int k = 1;
//    int ks[]{10};
    int node_nums[]{1,2,3,4,5,10,25};
//    int node_nums[]{10};
    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    long offset = 0;
    fseek(f, offset * Const::tsLengthBytes, SEEK_SET);
    for(int node_num: node_nums){
        int recallNums[maxExprRound];
        int search_number[maxExprRound];
        int layers[maxExprRound];
        long duration[maxExprRound];
        double error_ratio[maxExprRound];
        double inv_error_ratio[maxExprRound];
        cout<<"------------------Experiment--------------------" << endl;
        cout<<"k: " << k << endl;
        cout<<"node number: " << node_num<< endl;

        for(int curRound = 0; curRound < maxExprRound; ++curRound){
            //                    cout<<"Round : " + (curRound + 1));
            c_nodes.clear();
            _search_num = 0;
            query = FileUtil::readSeries(f);
            auto start = chrono::system_clock::now();
            vector<PqItemSeries*> *approxKnn = FADASSearcher::approxIncSearchFuzzy(root, query, k, index_dir, node_num);
            auto end = chrono::system_clock::now();
//                for(int i=0;i<256;++i)
//                    cout << (*approxKnn)[0]->ts[i] <<",";
            vector<float*>* exactKnn = getResult(Const::resfn, offset + curRound, k);
            vector<PqItemSeries*> exactKnn2;
            for(float *t: *exactKnn)
                exactKnn2.push_back(new PqItemSeries(t, query));

            duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count();
            layers[curRound] = layer;
            recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
            search_number[curRound] = _search_num;
            error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
            inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
//                cout << curRound << ":"<<recallNums[curRound] << endl;
            cout << recallNums[curRound] << ",";
            fflush(stdout);
            free_heap(approxKnn);
            for(int i=0;i<k;++i)
                delete[] (*exactKnn)[i];
            delete exactKnn;
            delete[] query;
        }
        cout << fixed  << endl;
        for(int _:layers)   cout << _ << ",";
        cout << endl;
//            for(double _:error_ratio)   cout << _ << ",";
        int totalRecallNum = 0;
        for(int temp:recallNums)
            totalRecallNum += temp;
        cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
        double totalErrorRatio = 0;
        for(double _:error_ratio)   totalErrorRatio += _;
        cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
        double totalinvErrorRatio = 0;
        for(double _:inv_error_ratio)   totalinvErrorRatio += _;
        cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
        double total_duration = 0;
        for(long _:duration)    total_duration += _;
        total_duration /= (double ) maxExprRound;
        cout<<"Average duration is : " << total_duration << "us. "
            <<"And QPS = "<< 1000000.0 / total_duration <<endl;
        for(int _:search_number)    cout << _ <<",";
        cout << endl;
        rewind(f);

    }

    fclose(f);
}

void Recall::exactSearchFADAS(FADASNode *root, vector<vector<int>>*g) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);

    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    int k = Const::k;
    int recallNums[maxExprRound];
    LB_NODE_CNT = 0;
    int search_number[maxExprRound];
    float query_reordered[Const::tsLength];
    int ordering[Const::tsLength];
    double duration[maxExprRound];  // ms
    double read_time[maxExprRound], prepare_time[maxExprRound], search_time[maxExprRound], dist_calc_time[maxExprRound];
    int load_node_cnt[maxExprRound];
    cout<<"------------------Experiment--------------------" << endl;
    cout<<"k: " << k << endl;
    cout << fixed << setprecision(3);

    for(int curRound = 0; curRound < maxExprRound; ++curRound){
        //                    cout<<"Round : " + (curRound + 1));
        c_nodes.clear();
        _search_num = 0;
        READ_TIME = 0;  PREPARE_TIME = 0; SEARCH_TIME = 0; DIST_CALC_TIME = 0;
        LOADED_NODE_CNT = 0;
        query = FileUtil::readSeries(f);
//        if(curRound < 25)   continue;

        reorder_query(query, query_reordered, ordering);
        auto start = chrono::system_clock::now();
        vector<PqItemSeries*> *approxKnn = FADASSearcher::Par_exactSearchIdLevel_SSD(root, query, k, g, query_reordered, ordering);
        auto end = chrono::system_clock::now();
        duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;
//                vector<float*>* exactKnn = getResult(Const::resfn, _k*maxExprRound + curRound, k);
        vector<float*>* exactKnn = getResult(Const::resfn, curRound, k);
        vector<PqItemSeries*> exactKnn2;
        for(float *t: *exactKnn)
            exactKnn2.push_back(new PqItemSeries(t, query));

        load_node_cnt[curRound] = LOADED_NODE_CNT;
        read_time[curRound] = READ_TIME;
        prepare_time[curRound] = PREPARE_TIME;
        search_time[curRound] = SEARCH_TIME;
        dist_calc_time[curRound] = DIST_CALC_TIME;
        recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
        search_number[curRound] = _search_num;
        cout << recallNums[curRound] << ";"  << (double)_search_num / root->size << ", ";
        fflush(stdout);
//                analyzePrintSaxFADAS(approxKnn, exactKnn, query);
        free_heap(approxKnn);
        for(int i=0;i<k;++i)
            delete[] (*exactKnn)[i];
        delete exactKnn;
        delete[] query;
    }
    cout << endl;

    double totalDuration = 0;
    for(double _:duration)  {totalDuration += _; cout<< _ << ",";}
    cout << endl;
    double total_read_time = 0;
    for(double _:read_time) {total_read_time += _; cout<<_/1000.0<<",";}
    double total_prepare_time = 0;
    for(double _:prepare_time)  total_prepare_time += _;
    double total_search_time = 0;
    for(auto _:search_time) total_search_time += _;
    double total_dist_calc_time = 0;
    for(auto _:dist_calc_time)  total_dist_calc_time += _;
    cout << endl;
    long total_node_cnt = 0;
    for(int _:load_node_cnt) { total_node_cnt += _; cout << (_) << ","; }
    cout << endl;
    long totalSearchNum = 0;
    for(int _:search_number) { totalSearchNum += _; cout << (_ / (double)root->size) << ","; }
    cout << endl;
    int totalRecallNum = 0;
    for(int temp:recallNums)
        totalRecallNum += temp;
    cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
    cout << "Average Search rate is " <<(double)totalSearchNum / ((double)maxExprRound * root->size)<<endl;
    cout <<"Average Duration cost is " << totalDuration / maxExprRound << " ms" << endl;
//            cout << "Avg. Distance Calculation time is " << DIST_CALC_TIME / 1000.0 / maxExprRound <<"ms."<<endl;
    cout << "AVG. I/O read time is " << total_read_time / 1000.0 / maxExprRound << "ms."<<endl;
    cout << "AVG. prepare time is " << total_prepare_time / 1000.0 / maxExprRound << "ms."<<endl;
    cout << "AVG. pre calc  time is " << total_dist_calc_time / 1000.0 / maxExprRound << "ms."<<endl;
    cout << "AVG. parallel time is " << total_search_time / 1000.0 / maxExprRound << "ms."<<endl;
    cout << "Avg. Priority Queue nodes count = " << LB_NODE_CNT / maxExprRound <<endl;
    cout << "Avg. Loaded nodes count = " << total_node_cnt / maxExprRound << endl;
    cout << endl;
    rewind(f);


    fclose(f);
}

void Recall::exactSearchFADASDTW(FADASNode *root, vector<vector<int>>*g) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);

    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    int k = Const::k;
    int recallNums[maxExprRound];
    LB_NODE_CNT = 0;
    int search_number[maxExprRound];
    double duration[maxExprRound];  // ms
    double read_time[maxExprRound];
    int load_node_cnt[maxExprRound];
    cout<<"------------------Experiment--------------------" << endl;
    cout<<"k: " << k << endl;

    for(int curRound = 0; curRound < maxExprRound; ++curRound){
        //                    cout<<"Round : " + (curRound + 1));
        c_nodes.clear();
        _search_num = 0;
        READ_TIME = 0;
        LOADED_NODE_CNT = 0;
        query = FileUtil::readSeries(f);
        auto start = chrono::system_clock::now();
        vector<PqItemSeries*> *approxKnn = FADASSearcher::exactSearchDTW(root, query, k, g);
        auto end = chrono::system_clock::now();
        duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;
//                vector<float*>* exactKnn = getResult(Const::resfn, _k*maxExprRound + curRound, k);
        vector<float*>* exactKnn = getResult(Const::dtwresfn, curRound, k);
        vector<PqItemSeries*> exactKnn2;
        for(float *t: *exactKnn) {
            double dist = TimeSeriesUtil::dtw(t, query, Const::tsLength, Const::dtw_window_size, numeric_limits<double>::max());
            exactKnn2.push_back(new PqItemSeries(t, dist));
        }
        load_node_cnt[curRound] = LOADED_NODE_CNT;
        read_time[curRound] = READ_TIME;
        recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
        search_number[curRound] = _search_num;
        cout << recallNums[curRound] << ",";
        fflush(stdout);
//                analyzePrintSaxFADAS(approxKnn, exactKnn, query);
        free_heap(approxKnn);
        for(int i=0;i<k;++i)
            delete[] (*exactKnn)[i];
        delete exactKnn;
        delete[] query;
    }
    cout << fixed  << endl;

    double totalDuration = 0;
    for(double _:duration)  {totalDuration += _; cout<< _ << ",";}
    cout << endl;
    double total_read_time = 0;
    for(double _:read_time) {total_read_time += _; cout<<_/1000.0<<",";}
    cout << endl;
    long total_node_cnt = 0;
    for(int _:load_node_cnt) { total_node_cnt += _; cout << (_) << ","; }
    cout << endl;
    long totalSearchNum = 0;
    for(int _:search_number) { totalSearchNum += _; cout << (_ / (double)root->size) << ","; }
    cout << endl;
    int totalRecallNum = 0;
    for(int temp:recallNums)
        totalRecallNum += temp;
    cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
    cout << "Average Search rate is " <<(double)totalSearchNum / ((double)maxExprRound * root->size)<<endl;
    cout <<"Average Duration cost is " << totalDuration / maxExprRound << " ms" << endl;
//            cout << "Avg. Distance Calculation time is " << DIST_CALC_TIME / 1000.0 / maxExprRound <<"ms."<<endl;
    cout << "AVG. I/O read time is " << total_read_time / 1000.0 / maxExprRound << "ms."<<endl;
    cout << "Avg. Priority Queue nodes count = " << LB_NODE_CNT / maxExprRound <<endl;
    cout << "Avg. Loaded nodes count = " << total_node_cnt / maxExprRound << endl;
    cout << endl;
    rewind(f);


    fclose(f);
}

void Recall::exactSearchFADASNoExpr(FADASNode *root, vector<vector<int>>*g) {
//    int series_num = FileUtil::getFileSize(Const::datafn.c_str()) / Const::tsLengthBytes;
    int series_num = root->size;
    int maxExprRound = Const::query_num;
    int k = 50;
    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    int search_number[maxExprRound];
    double duration[maxExprRound];  // ms
    cout<<"------------------Experiment--------------------" << endl;


    for(int curRound = 0; curRound < maxExprRound; ++curRound){
        //                    cout<<"Round : " + (curRound + 1));
        c_nodes.clear();
        _search_num = 0;
        query = FileUtil::readSeries(f);
        auto start = chrono::system_clock::now();
        vector<PqItemSeries*> *approxKnn = FADASSearcher::exactSearch(root, query, k, g);
        auto end = chrono::system_clock::now();

        duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;
        search_number[curRound] = _search_num;
        cout << duration[curRound] << ";" << 1 - (_search_num / (double)series_num) << " , ";
        fflush(stdout);
        FILE *outf = fopen(Const::resfn.c_str(), "ab");
        for(PqItemSeries* series : *approxKnn){
            float * ts = series->ts;
            fwrite(ts, sizeof(float), Const::tsLength, outf);
        }
        fclose(outf);

        free_heap(approxKnn);
        delete[] query;
    }
    cout << fixed  << endl;
    long totalSearchNum = 0;
    for(int _:search_number) { totalSearchNum += _; cout << 1.0 - (_ / (double )series_num) << ","; }
    double totalDuration = 0;
    for(double _:duration)  {totalDuration += _; }
    cout << endl;
    cout << "Average Prune rate is " << 1.0 - (double)totalSearchNum / (maxExprRound * (double )series_num)<<endl;
    cout <<"Average Duration cost is " << totalDuration / maxExprRound << " ms" << endl;
    cout << endl;
    fclose(f);

}

void Recall::exactSearchFADASPos(FADASNode *root, vector<vector<int>>*g) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);
//    int ks[]{1,5,10,50};
//    int ks[]{1};
//    int thresholds[]{10000};
    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    int k = Const::k;
    int recallNums[maxExprRound];
    int search_number[maxExprRound];
    double duration[maxExprRound];  // ms
    double read_time[maxExprRound];
    int load_node_cnt[maxExprRound];
    cout<<"------------------Experiment--------------------" << endl;
    cout<<"k: " << k << endl;

    for(int curRound = 0; curRound < maxExprRound; ++curRound){
        //                    cout<<"Round : " + (curRound + 1));
        c_nodes.clear();
        _search_num = 0;
        READ_TIME = 0;
        LOADED_NODE_CNT = 0;
        query = FileUtil::readSeries(f);
        auto start = chrono::system_clock::now();
        vector<PqItemSeries*> *approxKnn = FADASSearcher::exactSearchPos(root, query, k, g);
        auto end = chrono::system_clock::now();
        duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;
//                vector<float*>* exactKnn = getResult(Const::resfn, _k*maxExprRound + curRound, k);
        vector<float*>* exactKnn = getResult(Const::resfn, curRound, k);
        vector<PqItemSeries*> exactKnn2;
        for(float *t: *exactKnn)
            exactKnn2.push_back(new PqItemSeries(t, query));

        load_node_cnt[curRound] = LOADED_NODE_CNT;
        read_time[curRound] = READ_TIME;
        recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
        search_number[curRound] = _search_num;
        cout << recallNums[curRound] << ",";
        fflush(stdout);
//                analyzePrintSaxFADAS(approxKnn, exactKnn, query);
        free_heap(approxKnn);
        for(int i=0;i<k;++i)
            delete[] (*exactKnn)[i];
        delete exactKnn;
        delete[] query;
    }
    cout << fixed  << endl;

    double totalDuration = 0;
    for(double _:duration)  {totalDuration += _; cout<< _ << ",";}
    cout << endl;
    double total_read_time = 0;
    for(double _:read_time) {total_read_time += _; cout<<_/1000.0<<",";}
    cout << endl;
    long total_node_cnt = 0;
    for(int _:load_node_cnt) { total_node_cnt += _; cout << (_) << ","; }
    cout << endl;
    long totalSearchNum = 0;
    for(int _:search_number) { totalSearchNum += _; cout << (_ / (double)root->size) << ","; }
    cout << endl;
    int totalRecallNum = 0;
    for(int temp:recallNums)
        totalRecallNum += temp;
    cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
    cout << "Average Search rate is " <<(double)totalSearchNum / ((double)maxExprRound * root->size)<<endl;
    cout <<"Average Duration cost is " << totalDuration / maxExprRound << " ms" << endl;
//            cout << "Avg. Distance Calculation time is " << DIST_CALC_TIME / 1000.0 / maxExprRound <<"ms."<<endl;
    cout << "AVG. I/O read time is " << total_read_time / 1000.0 / maxExprRound << "ms."<<endl;
    cout << "Avg. Priority Queue nodes count = " << LB_NODE_CNT / maxExprRound <<endl;
    cout << "Avg. Loaded nodes count = " << total_node_cnt / maxExprRound << endl;
    cout << endl;
    rewind(f);

    fclose(f);
}

void Recall::exactSearchFADASPosNoExpr(FADASNode *root, vector<vector<int>>*g) {
//    int series_num = FileUtil::getFileSize(Const::datafn.c_str()) / Const::tsLengthBytes;
    int series_num = root->size;
    int maxExprRound = Const::query_num;
    int k = 1;
    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    int search_number[maxExprRound];
    double duration[maxExprRound];  // ms
    cout<<"------------------Experiment--------------------" << endl;
//    fseek(f, 39 * Const::tsLengthBytes, SEEK_SET);
    for(int curRound = 0; curRound < maxExprRound; ++curRound){
        //                    cout<<"Round : " + (curRound + 1));
        c_nodes.clear();
        _search_num = 0;
        query = FileUtil::readSeries(f);
        auto start = chrono::system_clock::now();
        vector<PqItemSeries*> *approxKnn = FADASSearcher::exactSearchPos(root, query, k, g);
        auto end = chrono::system_clock::now();

        duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;
        search_number[curRound] = _search_num;
        cout << duration[curRound] << ";" << 1 - (_search_num / (double)series_num) << " , ";
        fflush(stdout);
        FILE *outf = fopen(Const::resfn.c_str(), "ab");
        for(PqItemSeries* series : *approxKnn){
            float * ts = series->ts;
            fwrite(ts, sizeof(float), Const::tsLength, outf);
        }
        fclose(outf);

        free_heap(approxKnn);
        delete[] query;
    }
    cout << fixed  << endl;
    long totalSearchNum = 0;
    for(int _:search_number) { totalSearchNum += _; cout << 1.0 - (_ / (double )series_num) << ","; }
    double totalDuration = 0;
    for(double _:duration)  {totalDuration += _; }
    cout << endl;
    cout << "Average Prune rate is " << 1.0 - (double)totalSearchNum / (maxExprRound * (double )series_num)<<endl;
    cout <<"Average Duration cost is " << totalDuration / maxExprRound << " ms" << endl;
    cout << endl;
    fclose(f);

}

void Recall::approxTARDISORIGIN(TARGNode *root) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);
    int ks[]{1,3,5,10,25,50};
//    int ks[]{25};
    int thresholds[]{10000};
    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    long offset = 0;
    fseek(f, offset * Const::tsLengthBytes, SEEK_SET);
    for(int threshold:thresholds){
        int _k = 0;
        for(int k:ks){
            int recallNums[maxExprRound];
            int search_number[maxExprRound];
            int layers[maxExprRound];
            long duration[maxExprRound];
            double error_ratio[maxExprRound];
            double inv_error_ratio[maxExprRound];
            cout<<"------------------Experiment--------------------" << endl;
            cout<<"k: " << k << endl;
            cout<<"threshold: " << threshold<< endl;

            for(int curRound = 0; curRound < maxExprRound; ++curRound){
                //                    cout<<"Round : " + (curRound + 1));
                c_nodes.clear();
                _search_num = 0;
                query = FileUtil::readSeries(f);
                auto start = chrono::system_clock::now();
                vector<PqItemSeries*> *approxKnn = TARSearcher::approxSearch(root, query, k, Const::tardisfn);
                auto end = chrono::system_clock::now();
//                for(int i=0;i<256;++i)
//                    cout << (*approxKnn)[0]->ts[i] <<",";
                vector<float*>* exactKnn = getResult(Const::resfn, offset + curRound, k);
                vector<PqItemSeries*> exactKnn2;
                for(float *t: *exactKnn)
                    exactKnn2.push_back(new PqItemSeries(t, query));

                duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count();
                layers[curRound] = layer;
                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
                search_number[curRound] = _search_num;
                error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
                inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
//                cout << curRound << ":"<<recallNums[curRound] << endl;
                cout << recallNums[curRound] << ",";
                fflush(stdout);
                free_heap(approxKnn);
                for(int i=0;i<k;++i)
                    delete[] (*exactKnn)[i];
                delete exactKnn;
                delete[] query;
            }
            ++_k;
            cout << fixed  << endl;
            for(int _:layers)   cout << _ << ",";
            cout << endl;
//            for(double _:error_ratio)   cout << _ << ",";
            int totalRecallNum = 0;
            for(int temp:recallNums)
                totalRecallNum += temp;
            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
            double totalErrorRatio = 0;
            for(double _:error_ratio)   totalErrorRatio += _;
            cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
//            double totalinvErrorRatio = 0;
//            for(double _:inv_error_ratio)   totalinvErrorRatio += _;
//            cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
            double total_duration = 0;
            for(long _:duration)    total_duration += _;
            total_duration /= (double ) maxExprRound;
            cout<<"Average duration is : " << total_duration << "us. "
                <<"And QPS = "<< 1000000.0 / total_duration <<endl;
            for(int _:search_number)    cout << _ <<",";
            cout << endl;
            rewind(f);
        }
    }

    fclose(f);
}

void Recall::approxIncSearchTARDISORIGIN(TARGNode *root) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);
    int k = 1;
//    int ks[]{10};
    int node_nums[]{100000};
//    int node_nums[]{2,3,4,5,10,25};
//    int node_nums[]{25};
    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    long offset = 0;
    fseek(f, offset * Const::tsLengthBytes, SEEK_SET);
    for(int node_num: node_nums){
        int recallNums[maxExprRound];
        int search_number[maxExprRound];
        int layers[maxExprRound];
        int search_pack[maxExprRound];
        long duration[maxExprRound];
        double error_ratio[maxExprRound];
        double inv_error_ratio[maxExprRound];
        cout<<"------------------Experiment--------------------" << endl;
        cout<<"k: " << k << endl;
        cout<<"node number: " << node_num<< endl;

        for(int curRound = 0; curRound < maxExprRound; ++curRound){
            //                    cout<<"Round : " + (curRound + 1));
            c_nodes.clear();
            _search_num = 0;
            query = FileUtil::readSeries(f);
            LOADED_PACK_CNT = 0;
//            if(curRound < 133)  continue;
            auto start = chrono::system_clock::now();
            int _ = node_num;
            vector<PqItemSeries*> *approxKnn = TARSearcher::approxIncSearch(root, query, k, Const::tardisfn, &_);
            auto end = chrono::system_clock::now();
//                for(int i=0;i<256;++i)
//                    cout << (*approxKnn)[0]->ts[i] <<",";
            vector<float*>* exactKnn = getResult(Const::resfn, offset + curRound, k);
            vector<PqItemSeries*> exactKnn2;
            for(float *t: *exactKnn)
                exactKnn2.push_back(new PqItemSeries(t, query));

            duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count();
            layers[curRound] = layer;
            recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
            search_number[curRound] = _;
            search_pack[curRound] = LOADED_PACK_CNT;
            error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
            inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
//                cout << curRound << ":"<<recallNums[curRound] << endl;
            cout << recallNums[curRound] << ",";
            fflush(stdout);
            free_heap(approxKnn);
            for(int i=0;i<k;++i)
                delete[] (*exactKnn)[i];
            delete exactKnn;
            delete[] query;
        }
        cout << fixed  << endl;
        for(int _:layers)   cout << _ << ",";
        cout << endl;
        for(int _:search_number)   cout << _ << ",";
        cout << endl;
        for(int _:search_pack)   cout << _ << ",";
        cout << endl;


//            for(double _:error_ratio)   cout << _ << ",";
        int totalRecallNum = 0;
        for(int temp:recallNums)
            totalRecallNum += temp;
        cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
        double totalErrorRatio = 0;
        for(double _:error_ratio)   totalErrorRatio += _;
        cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
//        double totalinvErrorRatio = 0;
//        for(double _:inv_error_ratio)   totalinvErrorRatio += _;
//        cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
        double total_duration = 0;
        for(long _:duration)    total_duration += _;
        total_duration /= (double ) maxExprRound;
        cout<<"Average duration is : " << total_duration << "us. "
            <<"And QPM = "<< 60000000.0 / total_duration <<endl;
        rewind(f);

    }

    fclose(f);
}

void Recall::approxDTWTARDISORIGIN(TARGNode *root) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);
    int ks[]{1,3,5,10,25,50};
//    int ks[]{10};
    int thresholds[]{10000};
    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    long offset = 0;
    fseek(f, offset * Const::tsLengthBytes, SEEK_SET);
    for(int threshold:thresholds){
        int _k = 0;
        for(int k:ks){
            int recallNums[maxExprRound];
            int search_number[maxExprRound];
            int layers[maxExprRound];
            long duration[maxExprRound];
            double error_ratio[maxExprRound];
            double inv_error_ratio[maxExprRound];
            cout<<"------------------Experiment--------------------" << endl;
            cout<<"k: " << k << endl;
            cout<<"threshold: " << threshold<< endl;

            for(int curRound = 0; curRound < maxExprRound; ++curRound){
                //                    cout<<"Round : " + (curRound + 1));
                c_nodes.clear();
                _search_num = 0;
                query = FileUtil::readSeries(f);
                auto start = chrono::system_clock::now();
                vector<PqItemSeries*> *approxKnn =  TARSearcher::approxSearchDTW(root ,query, k, Const::tardisfn, Const::dtw_window_size);
                auto end = chrono::system_clock::now();
//                for(int i=0;i<256;++i)
//                    cout << (*approxKnn)[0]->ts[i] <<",";
                vector<float*>* exactKnn = getResult(Const::dtwresfn, offset + curRound, k);
                vector<PqItemSeries*> exactKnn2;
                for(float *t: *exactKnn) {
                    double dist = TimeSeriesUtil::dtw(t, query, Const::tsLength, Const::dtw_window_size, numeric_limits<double>::max());
                    exactKnn2.push_back(new PqItemSeries(t, dist));
                }

                duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count();
                layers[curRound] = layer;
                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
                search_number[curRound] = _search_num;
                error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
                inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
//                cout << curRound << ":"<<recallNums[curRound] << endl;
                cout << recallNums[curRound] << ",";
                fflush(stdout);
                free_heap(approxKnn);
                for(int i=0;i<k;++i)
                    delete[] (*exactKnn)[i];
                delete exactKnn;
                delete[] query;
            }
            ++_k;
            cout << fixed  << endl;
            for(int _:layers)   cout << _ << ",";
            cout << endl;
//            for(double _:error_ratio)   cout << _ << ",";
            int totalRecallNum = 0;
            for(int temp:recallNums)
                totalRecallNum += temp;
            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
            double totalErrorRatio = 0;
            for(double _:error_ratio)   totalErrorRatio += _;
            cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
//            double totalinvErrorRatio = 0;
//            for(double _:inv_error_ratio)   totalinvErrorRatio += _;
//            cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
            double total_duration = 0;
            for(long _:duration)    total_duration += _;
            total_duration /= (double ) maxExprRound;
            cout<<"Average duration is : " << total_duration << "us. "
                <<"And QPS = "<< 1000000.0 / total_duration <<endl;
            for(int _:search_number)    cout << _ <<",";
            cout << endl;
            rewind(f);
        }
    }

    fclose(f);
}

void Recall::exactSearchTARDISORIGIN(TARGNode *root) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);

    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    int k = Const::k;
    int recallNums[maxExprRound];
    LB_NODE_CNT = 0;
    int search_number[maxExprRound];
    double duration[maxExprRound];  // ms
    double read_time[maxExprRound];
    int load_node_cnt[maxExprRound];
    cout<<"------------------Experiment--------------------" << endl;
    cout<<"k: " << k << endl;

    for(int curRound = 0; curRound < maxExprRound; ++curRound){
        //                    cout<<"Round : " + (curRound + 1));
        c_nodes.clear();
        _search_num = 0;
        READ_TIME = 0;
        LOADED_NODE_CNT = 0;
        LOADED_PACK_CNT = 0;
        query = FileUtil::readSeries(f);
        auto start = chrono::system_clock::now();
        vector<PqItemSeries*> *approxKnn = TARSearcher::exactSearch(root, query, k, Const::tardisfn);
        auto end = chrono::system_clock::now();
        duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;
//                vector<float*>* exactKnn = getResult(Const::resfn, _k*maxExprRound + curRound, k);
        vector<float*>* exactKnn = getResult(Const::resfn, curRound, k);
        vector<PqItemSeries*> exactKnn2;
        for(float *t: *exactKnn)
            exactKnn2.push_back(new PqItemSeries(t, query));

        load_node_cnt[curRound] = LOADED_NODE_CNT;
        read_time[curRound] = READ_TIME;
        recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
        search_number[curRound] = _search_num;
        cout << recallNums[curRound] << ",";
        fflush(stdout);
//                analyzePrintSaxFADAS(approxKnn, exactKnn, query);
        free_heap(approxKnn);
        for(int i=0;i<k;++i)
            delete[] (*exactKnn)[i];
        delete exactKnn;
        delete[] query;
    }
    cout << fixed  << endl;
    long series_num = FileUtil::getFileSize(Const::datafn.c_str()) / Const::tsLengthBytes;
    cout << " series number = " << series_num << endl;
    double totalDuration = 0;
    for(double _:duration)  {totalDuration += _; cout<< _ << ",";}
    cout << endl;
    double total_read_time = 0;
    for(double _:read_time) {total_read_time += _; cout<<_/1000.0<<",";}
    cout << endl;
    long total_node_cnt = 0;
    for(int _:load_node_cnt) { total_node_cnt += _; cout << (_) << ","; }
    cout << endl;
    long totalSearchNum = 0;
    for(int _:search_number) { totalSearchNum += _; cout << (_ / (double)series_num) << ","; }
    cout << endl;
    int totalRecallNum = 0;
    for(int temp:recallNums)
        totalRecallNum += temp;

    cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
    cout << "Average Search rate is " <<(double)totalSearchNum / ((double)maxExprRound * series_num)<<endl;
    cout <<"Average Duration cost is " << totalDuration / maxExprRound << " ms" << endl;
//            cout << "Avg. Distance Calculation time is " << DIST_CALC_TIME / 1000.0 / maxExprRound <<"ms."<<endl;
    cout << "AVG. I/O read time is " << total_read_time / 1000.0 / maxExprRound << "ms."<<endl;
    cout << "Avg. Priority Queue nodes count = " << LB_NODE_CNT / maxExprRound <<endl;
    cout << "Avg. Loaded nodes count = " << total_node_cnt / maxExprRound << endl;
    cout << endl;
    rewind(f);


    fclose(f);
}

void Recall::ngSearchTARDISORIGIN(TARGNode *root) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);
    int nprobes[]{50,100};
//    int nprobes[]{50};
//    int ks[]{25};
    int thresholds[]{10000};
    float *query;
    float query_ts_reordered[Const::tsLength];
    int ordering[Const::tsLength];
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    long offset = 0;
    fseek(f, offset * Const::tsLengthBytes, SEEK_SET);
    for(int threshold:thresholds){
        int _k = 0;
        for(int probe:nprobes){
            int recallNums[maxExprRound];
            int search_number[maxExprRound];
            int layers[maxExprRound];
            long duration[maxExprRound];
            double error_ratio[maxExprRound];
            double inv_error_ratio[maxExprRound];
            cout<<"------------------Experiment--------------------" << endl;
            cout<<"nprobe: " << probe << endl;
            cout<<"threshold: " << threshold<< endl;

            for(int curRound = 0; curRound < maxExprRound; ++curRound){
                //                    cout<<"Round : " + (curRound + 1));
                c_nodes.clear();
                _search_num = 0;
                query = FileUtil::readSeries(f);
                reorder_query(query, query_ts_reordered, ordering);
                auto start = chrono::system_clock::now();
                vector<PqItemSeries*> *approxKnn = TARSearcher::ngSearch(root, query, Const::k, Const::tardisfn, probe);
                auto end = chrono::system_clock::now();
//                for(int i=0;i<256;++i)
//                    cout << (*approxKnn)[0]->ts[i] <<",";
                vector<float*>* exactKnn = getResult(Const::resfn, offset + curRound, Const::k);
                vector<PqItemSeries*> exactKnn2;
                for(float *t: *exactKnn)
                    exactKnn2.push_back(new PqItemSeries(t, query));

//                layer = 4;
//                analyzePrintSax(approxKnn,exactKnn, query);
                duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count();
                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
                search_number[curRound] = _search_num;
                error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, Const::k);
                inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, Const::k);
//                cout << curRound << ":"<<recallNums[curRound] << endl;
                cout << recallNums[curRound] << ",";
                fflush(stdout);
                free_heap(approxKnn);
                for(int i=0;i<Const::k;++i)
                    delete[] (*exactKnn)[i];
                delete exactKnn;
                delete[] query;
            }
            ++_k;
            cout << fixed  << endl;
//            for(int _:layers)   cout << _ << ",";
//            cout << endl;
//            for(double _:error_ratio)   cout << _ << ",";
            int totalRecallNum = 0;
            for(int temp:recallNums)
                totalRecallNum += temp;
            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * Const::k)<< endl;
            double totalErrorRatio = 0;
            for(double _:error_ratio)   totalErrorRatio += _;
            cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
//            double totalinvErrorRatio = 0;
//            for(double _:inv_error_ratio)   totalinvErrorRatio += _;
//            cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
            double total_duration = 0;
            for(long _:duration)    total_duration += _;
            total_duration /= (double ) maxExprRound;
            cout<<"Average duration is : " << total_duration/1000.0 << "ms. "
                <<"And QPM = "<< 60000000.0 / total_duration <<endl;
            for(int _:search_number)    cout << _ <<",";
            cout << endl;
            rewind(f);
        }
    }

    fclose(f);
}

void Recall::exactSearchTARDISORIGINDTW(TARGNode *root) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);

    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    int k = Const::k;
    int recallNums[maxExprRound];
    LB_NODE_CNT = 0;
    int search_number[maxExprRound];
    double duration[maxExprRound];  // ms
    double read_time[maxExprRound];
    int load_node_cnt[maxExprRound];
    cout<<"------------------Experiment--------------------" << endl;
    cout<<"k: " << k << endl;

    for(int curRound = 0; curRound < maxExprRound; ++curRound){
        //                    cout<<"Round : " + (curRound + 1));
        c_nodes.clear();
        _search_num = 0;
        READ_TIME = 0;
        LOADED_NODE_CNT = 0;
        query = FileUtil::readSeries(f);
        auto start = chrono::system_clock::now();
        vector<PqItemSeries*> *approxKnn = TARSearcher::exactSearchDTW(root, query, k, Const::tardisfn);
        auto end = chrono::system_clock::now();
        duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0;
//                vector<float*>* exactKnn = getResult(Const::resfn, _k*maxExprRound + curRound, k);
        vector<float*>* exactKnn = getResult(Const::dtwresfn, curRound, k);
        vector<PqItemSeries*> exactKnn2;
        for(float *t: *exactKnn) {
            double dist = TimeSeriesUtil::dtw(t, query, Const::tsLength, Const::dtw_window_size, numeric_limits<double>::max());
            exactKnn2.push_back(new PqItemSeries(t, dist));
        }
        load_node_cnt[curRound] = LOADED_NODE_CNT;
        read_time[curRound] = READ_TIME;
        recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
        search_number[curRound] = _search_num;
        cout << recallNums[curRound] << "," <<read_time[curRound]<<","<<load_node_cnt[curRound]<<","<<search_number[curRound]<<","<< duration[curRound]<<endl;
        fflush(stdout);
//                analyzePrintSaxFADAS(approxKnn, exactKnn, query);
        free_heap(approxKnn);
        for(int i=0;i<k;++i)
            delete[] (*exactKnn)[i];
        delete exactKnn;
        delete[] query;
    }
    cout << fixed  << endl;

    long series_num = FileUtil::getFileSize(Const::datafn.c_str()) / Const::tsLengthBytes;
    cout << " series number = " << series_num << endl;

    double totalDuration = 0;
    for(double _:duration)  {totalDuration += _; cout<< _ << ",";}
    cout << endl;
    double total_read_time = 0;
    for(double _:read_time) {total_read_time += _; cout<<_/1000.0<<",";}
    cout << endl;
    long total_node_cnt = 0;
    for(int _:load_node_cnt) { total_node_cnt += _; cout << (_) << ","; }
    cout << endl;
    long totalSearchNum = 0;
    for(int _:search_number) { totalSearchNum += _; cout << (_ / (double)series_num) << ","; }
    cout << endl;
    int totalRecallNum = 0;
    for(int temp:recallNums)
        totalRecallNum += temp;
    cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
    cout << "Average Search rate is " <<(double)totalSearchNum / ((double)maxExprRound * series_num)<<endl;
    cout <<"Average Duration cost is " << totalDuration / maxExprRound << " ms" << endl;
//            cout << "Avg. Distance Calculation time is " << DIST_CALC_TIME / 1000.0 / maxExprRound <<"ms."<<endl;
    cout << "AVG. I/O read time is " << total_read_time / 1000.0 / maxExprRound << "ms."<<endl;
    cout << "Avg. Priority Queue nodes count = " << LB_NODE_CNT / maxExprRound <<endl;
    cout << "Avg. Loaded nodes count = " << total_node_cnt / maxExprRound << endl;
    cout << endl;
    rewind(f);


    fclose(f);
}

void Recall::doExprWithResFADASPos(FADASNode *root, vector<vector<int>> *g, const string &index_dir) {
    int maxExprRound = Const::query_num;
    Const::logPrint( "result file is " + Const::resfn);
    int ks[]{1,3,5,10,25,50};
//    int ks[]{10};
    int thresholds[]{10000};
    float *query;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");
    long offset = 0;
    fseek(f, offset * Const::tsLengthBytes, SEEK_SET);
    for(int threshold:thresholds){
        int _k = 0;
        for(int k:ks){
            int recallNums[maxExprRound];
            int search_number[maxExprRound];
            int layers[maxExprRound];
            long duration[maxExprRound];
            double error_ratio[maxExprRound];
            double inv_error_ratio[maxExprRound];
            cout<<"------------------Experiment--------------------" << endl;
            cout<<"k: " << k << endl;
            cout<<"threshold: " << threshold<< endl;

            for(int curRound = 0; curRound < maxExprRound; ++curRound){
                //                    cout<<"Round : " + (curRound + 1));
                c_nodes.clear();
                _search_num = 0;
                query = FileUtil::readSeries(f);
                auto start = chrono::system_clock::now();
                vector<PqItemSeries*> *approxKnn = FADASSearcher::approxSearchPos(root, query, k, g, index_dir);
                auto end = chrono::system_clock::now();
//                for(int i=0;i<256;++i)
//                    cout << (*approxKnn)[0]->ts[i] <<",";
                vector<float*>* exactKnn = getResult(Const::resfn, offset + curRound, k);
                vector<PqItemSeries*> exactKnn2;
                for(float *t: *exactKnn)
                    exactKnn2.push_back(new PqItemSeries(t, query));

                duration[curRound] = chrono::duration_cast<chrono::microseconds>(end - start).count();
                layers[curRound] = layer;
                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
                search_number[curRound] = _search_num;
                error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
                inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
//                cout << curRound << ":"<<recallNums[curRound] << endl;
                cout << recallNums[curRound] << ",";
                fflush(stdout);
                free_heap(approxKnn);
                for(int i=0;i<k;++i)
                    delete[] (*exactKnn)[i];
                delete exactKnn;
                delete[] query;
            }
            ++_k;
            cout << fixed  << endl;
            for(int _:layers)   cout << _ << ",";
            cout << endl;
//            for(double _:error_ratio)   cout << _ << ",";
            int totalRecallNum = 0;
            for(int temp:recallNums)
                totalRecallNum += temp;
            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
            double totalErrorRatio = 0;
            for(double _:error_ratio)   totalErrorRatio += _;
            cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
            double totalinvErrorRatio = 0;
            for(double _:inv_error_ratio)   totalinvErrorRatio += _;
            cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
            double total_duration = 0;
            for(long _:duration)    total_duration += _;
            total_duration /= (double ) maxExprRound;
            cout<<"Average duration is : " << total_duration << "us. "
                <<"And QPS = "<< 1000000.0 / total_duration <<endl;
            for(int _:search_number)    cout << _ <<",";
            cout << endl;
            rewind(f);
        }
    }

    fclose(f);
}

void Recall::doExprWithResDynamicIPG(IPGNode *root, const string &queryFile, const string &resFile) {
    int maxExprRound = 200;
    cout << "result file is " << resFile <<endl;
    int ks[]{1,3 ,5,10,25,50};
//    int ks[]{1};
    int thresholds[]{10000};
    float *query;
    FILE *f = fopen(queryFile.c_str(), "rb");
    for(int threshold:thresholds){
        int _k = 0;
        for(int k:ks){
            int recallNums[maxExprRound];
            int search_number[maxExprRound];
            int layers[maxExprRound];
            double error_ratio[maxExprRound];
            double inv_error_ratio[maxExprRound];
            cout<<"------------------Experiment--------------------" << endl;
            cout<<"k: " << k << endl;
            cout<<"threshold: " << threshold<< endl;

            for(int curRound = 0; curRound < maxExprRound; ++curRound){
                //                    cout<<"Round : " + (curRound + 1));
                c_nodes.clear();
                _search_num = 0;
                query = FileUtil::readSeries(f);
                vector<PqItemSeries*> *approxKnn = IPGApproxSearcher::approxKnnSearchDynamicModelRange(root, query, k, threshold);
                vector<float*>* exactKnn = getResult(resFile, curRound, k);
                vector<PqItemSeries*> exactKnn2;
                for(float *t: *exactKnn)
                    exactKnn2.push_back(new PqItemSeries(t, query));

                layers[curRound] = layer;
                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
                search_number[curRound] = _search_num;
                error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
                inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
                cout << recallNums[curRound] << ",";
                fflush(stdout);
                free_heap(approxKnn);
                for(int i=0;i<k;++i)
                    delete[] (*exactKnn)[i];
                delete exactKnn;
                delete[] query;
            }
            ++_k;
            cout << fixed  << endl;
            for(int _:layers)   cout << _ << ",";
            cout << endl;
//            for(double _:error_ratio)   cout << _ << ",";
            int totalRecallNum = 0;
            for(int temp:recallNums)
                totalRecallNum += temp;
            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
            double totalErrorRatio = 0;
            for(double _:error_ratio)   totalErrorRatio += _;
            cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
            double totalinvErrorRatio = 0;
            for(double _:inv_error_ratio)   totalinvErrorRatio += _;
            cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
//            for(int _:search_number)    cout << _ <<",";
//            cout << endl;
            rewind(f);
        }

    }

    //    cout << "small: " << approxSearchUnits[0] << " series, cost " << approxSearchTimeDetails[0] << " us. "<< endl;
    //    cout << "middle: " << approxSearchUnits[1] << " series, cost " << approxSearchTimeDetails[1] << " us. "<< endl;
    //    cout << "big: " << approxSearchUnits[2] << " series, cost " << approxSearchTimeDetails[2] << " us. "<< endl;
    //    cout << "huge: " << approxSearchUnits[3] << " series, cost " << approxSearchTimeDetails[3] << " us. "<< endl;
    fclose(f);
}

void Recall::doExprWithResIPGWO1st(IPGNode *root, const string &queryFile, const string &resFile) {
    int maxExprRound = 200;
    cout << "result file is " << resFile <<endl;
    int ks[]{1,3 ,5,10,25,50};
//    int ks[]{1};
    int thresholds[]{10000};
    float *query;
    FILE *f = fopen(queryFile.c_str(), "rb");
    for(int threshold:thresholds){
        int _k = 0;
        for(int k:ks){
            int recallNums[maxExprRound];
            int search_number[maxExprRound];
            int layers[maxExprRound];
            double error_ratio[maxExprRound];
            double inv_error_ratio[maxExprRound];
//            bool flag[maxExprRound];
//            for(bool &_:flag)   _ = false;
            cout<<"------------------Experiment--------------------" << endl;
            cout<<"k: " << k << endl;
            cout<<"threshold: " << threshold<< endl;

            for(int curRound = 0; curRound < maxExprRound; ++curRound){
                //                    cout<<"Round : " + (curRound + 1));
                c_nodes.clear();
                _search_num = 0;
                query = FileUtil::readSeries(f);
                vector<PqItemSeries*> *approxKnn = IPGApproxSearcher::approxKnnSearchModel(root, query, k, threshold);
//                if(layer <= 1)  continue;
                vector<float*>* exactKnn = getResult(resFile, curRound, k);
                vector<PqItemSeries*> exactKnn2;
                for(float *t: *exactKnn)
                    exactKnn2.push_back(new PqItemSeries(t, query));
//                analyzePrintSaxIPG(approxKnn, exactKnn, query);
                //                vector<PqItemSeries*>exact;
                //                exact.reserve(k);
                //                for(float *ts:*exactKnn){
                //                    exact.push_back(new PqItemSeries(ts, query));
                //                }
                //                vector<PqItemSeries*> gt = ExactSearcher::groundTruthKnnSearch(IPGNode::rowDataFileName, query, k);
                //                    cout<<"Approximate Search Finished.\nExact Search Start!");
                //                    PqItemSeries[] exactKnn = readExactResult(g, cnt, k);
                //                    cout<<"Exact Search Finished");
                //                    if(TimeSeriesUtil::containsSeries(approxKnn, query)) cout<<"Approximate Search find the target series!");
                //                    else cout<<"WARNING!  Approximate Search does not find the target series! ");
                //                    if(TimeSeriesUtil::containsSeries(exactKnn, query)) cout<<"EXACT Search find the target series!");
                //                    else cout<<"ERROR!  EXACT Search does not find the target series! ");
//                flag[curRound] = true;
                layers[curRound] = layer;
                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
                search_number[curRound] = _search_num;
                error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
                inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
                cout << recallNums[curRound] << ",";
                fflush(stdout);
                free_heap(approxKnn);
                for(int i=0;i<k;++i)
                    delete[] (*exactKnn)[i];
                delete exactKnn;
                delete[] query;
            }

            //            cout<<"k: " << k<<endl;
            //            cout<<"threshold: " << threshold<<endl;

//            cout << fixed  << endl;
//            for(int i=0;i< maxExprRound;++i)
//                if(flag[i]) cout << layers[i] << ",";
//            cout << endl;
//            for(int i=0;i< maxExprRound;++i)
//                if(flag[i]) cout << search_number[i] << ",";
//            cout << endl;
//            int actualRound = 0;
//            for(int i=0;i< maxExprRound;++i) {
//                if (flag[i]) {
//                    actualRound++;
////                    cout << _k * maxExprRound + i << ",";
////                    cout << i << ",";
//
//                }
//            }
//            int totalRecallNum = 0;
//            for(int i=0;i< maxExprRound;++i)
//                if(flag[i]) totalRecallNum += recallNums[i];
//            cout<<"\nAverage Recall rate is : " << (double)totalRecallNum / (float) (actualRound * k)<< endl;
////            double totalErrorRatio = 0;
////            for(double _:error_ratio)   totalErrorRatio += _;
////            cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
//            double totalinvErrorRatio = 0;
//            for(int i=0;i< maxExprRound;++i)
//                if(flag[i]) totalinvErrorRatio += inv_error_ratio[i];
//            cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (actualRound)<< endl;
//            cout << endl;


            cout << fixed << setprecision(3) << endl;
            for(double _:error_ratio)   cout << _ << ",";
            int totalRecallNum = 0;
            cout << setprecision(6) << endl;
            for(int temp:recallNums)
                totalRecallNum += temp;
            double totalErrorRatio = 0;
            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
            for(double _:error_ratio)   totalErrorRatio += _;
            cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
            double totalinvErrorRatio = 0;
            for(double _:inv_error_ratio)   totalinvErrorRatio += _;
            cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
            for(int _:search_number)    cout << _ <<",";
            cout << endl;

            ++_k;
            rewind(f);
        }

    }

    //    cout << "small: " << approxSearchUnits[0] << " series, cost " << approxSearchTimeDetails[0] << " us. "<< endl;
    //    cout << "middle: " << approxSearchUnits[1] << " series, cost " << approxSearchTimeDetails[1] << " us. "<< endl;
    //    cout << "big: " << approxSearchUnits[2] << " series, cost " << approxSearchTimeDetails[2] << " us. "<< endl;
    //    cout << "huge: " << approxSearchUnits[3] << " series, cost " << approxSearchTimeDetails[3] << " us. "<< endl;
    fclose(f);
}

void Recall::doExprWithResiSAX(iSAXRoot *root, const string &resFile, const string &queryFile) {
    int maxExprRound = 40;
    cout << "result file is " << resFile <<endl;
    int ks[]{1,3,5,10,25,50};
    int thresholds[]{10000};
    float *query;
    FILE *f = fopen(queryFile.c_str(), "rb");
    for(int threshold:thresholds){
        int _k = 0;
        for(int k:ks){
            int recallNums[maxExprRound];
            int search_number[maxExprRound];
            double error_ratio[maxExprRound];
            double inv_error_ratio[maxExprRound];

            cout<<"------------------Experiment--------------------" << endl;
            cout<<"k: " << k << endl;
            cout<<"threshold: " << threshold<< endl;

            for(int curRound = 0; curRound < maxExprRound; ++curRound){
                //                    cout<<"Round : " + (curRound + 1));
                c_nodes.clear();
                _search_num = 0;
                query = FileUtil::readSeries(f);
                vector<PqItemSeries*> *approxKnn = iSAXSearcher::approxKnnSearch(root, query, k);
                vector<float*>* exactKnn = getResult(resFile, _k * maxExprRound + curRound, k);
                vector<PqItemSeries*> exactKnn2;
                for(float *t: *exactKnn)
                    exactKnn2.push_back(new PqItemSeries(t, query));
                //                analyzePrintSaxIPG(approxKnn, exactKnn, query);
                //                vector<PqItemSeries*>exact;
                //                exact.reserve(k);
                //                for(float *ts:*exactKnn){
                //                    exact.push_back(new PqItemSeries(ts, query));
                //                }
                //                vector<PqItemSeries*> gt = ExactSearcher::groundTruthKnnSearch(IPGNode::rowDataFileName, query, k);
                //                    cout<<"Approximate Search Finished.\nExact Search Start!");
                //                    PqItemSeries[] exactKnn = readExactResult(g, cnt, k);
                //                    cout<<"Exact Search Finished");
                //                    if(TimeSeriesUtil::containsSeries(approxKnn, query)) cout<<"Approximate Search find the target series!");
                //                    else cout<<"WARNING!  Approximate Search does not find the target series! ");
                //                    if(TimeSeriesUtil::containsSeries(exactKnn, query)) cout<<"EXACT Search find the target series!");
                //                    else cout<<"ERROR!  EXACT Search does not find the target series! ");
                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
                search_number[curRound] = _search_num;
                cout << recallNums[curRound] << ",";
                fflush(stdout);
                error_ratio[curRound] = MathUtil::errorRatio(*approxKnn, exactKnn2, k);
                inv_error_ratio[curRound] = MathUtil::invertedErrorRatio(*approxKnn, exactKnn2, k);
                free_heap(approxKnn);
                for(int i=0;i<k;++i)
                    delete[] (*exactKnn)[i];
                delete exactKnn;
                delete[] query;
            }
            ++_k;
            //            cout<<"k: " << k<<endl;
            //            cout<<"threshold: " << threshold<<endl;
            cout << fixed << setprecision(3) << endl;
            for(double _:error_ratio)   cout << _ << ",";
            int totalRecallNum = 0;
            cout << setprecision(6) << endl;
            for(int temp:recallNums)
                totalRecallNum += temp;
            double totalErrorRatio = 0;
            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
            for(double _:error_ratio)   totalErrorRatio += _;
            cout<<"Average Error ratio is : " << totalErrorRatio / (float) (maxExprRound)<< endl;
            double totalinvErrorRatio = 0;
            for(double _:inv_error_ratio)   totalinvErrorRatio += _;
            cout<<"Average inv Error ratio is : " << totalinvErrorRatio / (float) (maxExprRound)<< endl;
            for(int _:search_number)    cout << _ <<",";
            cout << endl;
        }
        rewind(f);
    }

    //    cout << "small: " << approxSearchUnits[0] << " series, cost " << approxSearchTimeDetails[0] << " us. "<< endl;
    //    cout << "middle: " << approxSearchUnits[1] << " series, cost " << approxSearchTimeDetails[1] << " us. "<< endl;
    //    cout << "big: " << approxSearchUnits[2] << " series, cost " << approxSearchTimeDetails[2] << " us. "<< endl;
    //    cout << "huge: " << approxSearchUnits[3] << " series, cost " << approxSearchTimeDetails[3] << " us. "<< endl;
    fclose(f);
}

void Recall::completeWorkload(){
//    auto* root  = FADASNode::BuildIndex(Const::datafn, Const::saxfn);
    auto root = FADASNode::loadFromDisk(Const::saxfn,Const::fidxfn + "root.idx",false);
    FILE *f = fopen(Const::datafn.c_str(), "rb");
    fseek(f, (long)Const::series_num*Const::tsLengthBytes, SEEK_SET);
    auto tss = new float [(long)Const::batch_size * Const::tsLength];
    for(int i=0;i<Const::batch_num;++i){
        Const::logPrint("Batch: " + to_string(i));
        fread(tss, sizeof(float),(long)Const::batch_size * Const::tsLength, f);
        root->insertBatch(tss, Const::batch_size);
    }
    Const::logPrint("Finished.");
    fclose(f);

}

int find_series(const string &fileName, float *ts){
    FILE *f = fopen(fileName.c_str(), "rb");
    long size = FileUtil::getFileSize(fileName.c_str());
    int num  = size / Const::tsLengthBytes;
    float t[Const::tsLength];
    for(int i=0;i<num;++i){
        fread(t, sizeof(float ), Const::tsLength, f);
        if(TimeSeriesUtil::isSame(t, ts, Const::tsLength, Const::tsLength))
            return i;
    }
    return  -1;
}

int count(const int x[], int a, int len){
    int res = 0;
    for(int i=0;i<len;++i)
        if(x[i] == a)  res++;
    return res;
}

int count(const int x[], int a, int len, const int y[], int head){
    int res = 0;
    for(int i=0;i<len;++i)
        if(x[i] == a && y[i] == head)  res++;
    return res;
}

//int reportKit(vector<Vertex *> & nnList, int startIndex, int endIndex, int approxWords[], int exactWords[], int accum , int dest, int app_size, int exact_size, int head, int approx1[], int exact1[]);
//void analyzeResult3(Graph* g, vector<PqItemSeries*> &approx, vector<float*> &exact, int k, float*query, int head, int approx1[], int exact1[]){
//    cout << "\n\t Now report deep level graph search result."<<endl;
//    assert(exact.size() == k);
//    TimeSeries ts(query);
//    int invSax = SaxUtil::invSaxHead2FromSax(ts.sax, TimeSeries::bitsCardinality, TimeSeries::segmentNum);
//    int approxWords[approx.size()];
//    int exactWords[k];
//    int dest = k + approx.size();
//
//    // compute invsax word head 2
//    for(int i=0;i<k;++i) {
//        auto *tsss = new float[TimeSeries::tsLength];
//        for(int j=0;j<TimeSeries::tsLength;++j)
//            tsss[j] = (exact[i])[j];
//        exactWords[i] = SaxUtil::invSaxHead2FromSax((new TimeSeries(tsss))->sax, TimeSeries::bitsCardinality, TimeSeries::segmentNum);
//        if(i < approx.size()) {
//            auto *tss = new float[TimeSeries::tsLength];
//            for(int j=0;j<TimeSeries::tsLength;++j)
//                tss[j] = approx[i]->ts[j];
//            approxWords[i] = SaxUtil::invSaxHead2FromSax(
//                    (new TimeSeries(tss))->sax, TimeSeries::bitsCardinality, TimeSeries::segmentNum);
//        }
//    }
//    vector<Vertex *> &nnList = g->NeighborsList[invSax];
//    int accum = 0;
//    cout<<"\n==================Layer 0:===================="<<endl;
//    int arn = count(approxWords, invSax, approx.size(), approx1, head), ern = count(exactWords, invSax, k, exact1, head);
//    SaxUtil::printBinary(invSax, TimeSeries::segmentNum);
//    cout<<": Count = " << g->VERTEXES[invSax]->cnt << " ; Approx Result Number = " << arn << " ; Exact Result Number = " << ern;
//    accum += arn + ern;
//    int pos = 0;
//    for(int i=1; i<= Const::bitsReserve; ++i){
//        cout<<"\n==================Layer " << i << " :========================"<<endl;
//        accum = reportKit(nnList, pos, pos + MathUtil::nChooseK(TimeSeries::segmentNum, i), approxWords, exactWords, accum, dest, approx.size(), k, head, approx1, exact1);
//        if(accum >= dest) return;
//        pos += MathUtil::nChooseK(TimeSeries::segmentNum, i);
//    }
//    cout << "\n\t Deep level graph search report finished.\n"<<endl;
//}
//int reportKit(vector<Vertex *> & nnList, int startIndex, int endIndex, int approxWords[], int exactWords[], int accum , int dest, int app_size, int exact_size, int head, int approx1[], int exact1[]){
//    for(int i=startIndex; i< endIndex; ++i){
//        int arn = count(approxWords, nnList[i]->id, app_size, approx1, head), ern = count(exactWords, nnList[i]->id, exact_size, exact1, head);
//        SaxUtil::printBinary(nnList[i]->id, TimeSeries::segmentNum);
//        cout<<": Count = " << nnList[i]->cnt << " ; Approx Result Number = " << arn <<" ; Exact Result Number = " << ern << endl;
//        accum += arn + ern;
//        if(accum >= dest)
//            return accum;
//    }
//    return accum;
//}
//void analyzeResult(Graph* g, vector<PqItemSeries*> &approx, vector<float*> &exact, int k, float*query){
//    assert(exact.size() == k);
//    TimeSeries ts(query);
//    int invSax = SaxUtil::invSaxHeadFromSax(ts.sax, TimeSeries::bitsCardinality, TimeSeries::segmentNum);
//    int approxWords[approx.size()];
//    int exactWords[k];
//    int dest = k + approx.size();
//
//    // compute approx and exact invsax word head
//    for(int i=0;i<k;++i) {
//        auto *tsss = new float[TimeSeries::tsLength];
//        for(int j=0;j<TimeSeries::tsLength;++j)
//            tsss[j] = exact[i][j];
//        exactWords[i] = SaxUtil::invSaxHeadFromSax((new TimeSeries(tsss))->sax, TimeSeries::bitsCardinality, TimeSeries::segmentNum);
//        if(i < approx.size()) {
//            auto *tss = new float[TimeSeries::tsLength];
//            for(int j=0;j<TimeSeries::tsLength;++j)
//                tss[j] = approx[i]->ts[j];
//            approxWords[i] = SaxUtil::invSaxHeadFromSax(
//                    (new TimeSeries(tss))->sax, TimeSeries::bitsCardinality, TimeSeries::segmentNum);
//        }
//    }
//
//
//    vector<Vertex *> &nnList = g->NeighborsList[invSax];
//    int accum = 0;
//    cout<<"\n==================Layer 0:===================="<<endl;
//    int arn = count(approxWords, invSax, approx.size()), ern = count(exactWords, invSax, k);
//    SaxUtil::printBinary(invSax, TimeSeries::segmentNum);
//    cout<<": Count = " << g->VERTEXES[invSax]->cnt << " ; Approx Result Number = " << arn << " ; Exact Result Number = " << ern;
//    if(g->VERTEXES[invSax]->level == huge) analyzeResult3(g->VERTEXES[invSax]->g, approx, exact, k, query, invSax, approxWords, exactWords);
//    accum += arn + ern;
//    int pos = 0;
//    for(int i=1; i<= Const::bitsReserve; ++i){
//        cout<<"\n==================Layer " << i << " :========================"<<endl;
//        for(int j=pos; j < pos + MathUtil::nChooseK(TimeSeries::segmentNum, i); ++j){
//            int arn2 = count(approxWords, nnList[j]->id, approx.size()), ern2 = count(exactWords, nnList[j]->id, k);
//            SaxUtil::printBinary(nnList[j]->id, TimeSeries::segmentNum);
//            cout << ": Count = " << nnList[j]->cnt << " ; Approx Result Number = " << arn2 << " ; Exact Result Number = " << ern2 << endl;
//            if(nnList[j]->level == huge) analyzeResult3(nnList[j]->g, approx, exact, k, query, nnList[j]->id, approxWords, exactWords);
//            accum += arn2 + ern2;
//            if(accum >= dest)
//                break;
//        }
//        if(accum >= dest) return;
//        pos += MathUtil::nChooseK(TimeSeries::segmentNum, i);
//    }
//}

void printResult(vector<PqItemSeries*> &res){
    int cnt = 0;
    for(const auto& ts:res) {
        cout << ++cnt << " : [ " << ts->dist << " ] ";
        for (int i=0; i < Const::tsLength; ++i)
            cout << ts->ts[i] << ",";
        cout << endl;
    }
    cout << endl;
}

void printResult2(vector<float*> &res, float *query){
    int cnt = 0;
    for(const auto& ts:res) {
        cout << ++cnt << " : [ " << TimeSeriesUtil::euclideanDist(ts, query, Const::tsLength) << " ] ";
        for (int i=0; i < Const::tsLength; ++i)
            cout << ts[i] << ",";
        cout << endl;
    }
    cout << endl;
}

//void Recall::doExprDebug(Graph* g) {
//    int maxExprRound = 1;
//    int lastIndex= g->rowDataFileName.rfind(".bin"), index = g->rowDataFileName.find("generator");
//    string resFile = g->rowDataFileName.substr(0, index) + "results" + g->rowDataFileName.substr(index + 9);
//    cout << "result file is " << resFile <<endl;
//    string fn = g->rowDataFileName.substr(0, lastIndex) + "_query.bin";
//    vector<int> ks{50};
//    vector<int> thresholds{50000};
//    float *query;
//    FILE *f = fopen(fn.c_str(), "rb");
//    for(int i=0;i<81;++i)
//        query = FileUtil::readSeries(f);
//    for(int threshold:thresholds){
//        for(int k:ks){
//            cout<<"k: " << k << endl;
//            cout<<"threshold: " << threshold<<endl;
//            int recallNums[maxExprRound];
//            for(int curRound = 0; curRound < maxExprRound; ++curRound){
////                query = FileUtil::readSeries(f);
//                vector<PqItemSeries*> *approxKnn = ApproximateSearcher::approxKnnSearchNew(g, query, k, threshold);
//                vector<float*>* gt = getResult(resFile, 80, k);
////                cout<<"Approximate Search Finished.\nExact Search Start!" << endl;
////                vector<PqItemSeries*> *exactKnn = ExactSearcher::exactKnnSearch(g, query, k, nullptr);
////                vector<PqItemSeries*> &gt = ExactSearcher::groundTruthKnnSearch(g->rowDataFileName, query, k);
//
//                cout << " ground truth result distribution: " << endl;
//                int i=0;
//                for(float * t: *gt){
//                    auto sax = SaxUtil::saxFromTs(t, TimeSeries::tsLengthPerSegment, TimeSeries::segmentNum, TimeSeries::cardinality);
//                    int head = SaxUtil::invSaxHeadFromSax(sax, TimeSeries::bitsCardinality, TimeSeries::segmentNum);
//                    int head2 = SaxUtil::invSaxHead2FromSax(sax, TimeSeries::bitsCardinality, TimeSeries::segmentNum);
////                    int off = find_series(g->rowDataFileName ,t);
////                    cout << ++id << " : " << head <<" , " << head2 << " , " <<off<<endl;
//                    cout << ++i << " : " << head <<" , " << head2 <<endl;
//                }
//
//                cout << " approx result distribution: " << endl;
//                i=0;
//                for(PqItemSeries* t: *approxKnn){
//                    auto sax = SaxUtil::saxFromTs(t->ts, TimeSeries::tsLengthPerSegment, TimeSeries::segmentNum, TimeSeries::cardinality);
//                    int head = SaxUtil::invSaxHeadFromSax(sax, TimeSeries::bitsCardinality, TimeSeries::segmentNum);
//                    int head2 = SaxUtil::invSaxHead2FromSax(sax, TimeSeries::bitsCardinality, TimeSeries::segmentNum);
////                    int off = find_series(g->rowDataFileName ,t->ts);
//                    cout << ++i << " : " << head <<" , " << head2 <<endl;
//                }
//
//                //                    PqItemSeries[] exactKnn = readExactResult(g, cnt, k);
////                cout<<"Exact Search Finished" << endl;
//                if(!TimeSeriesUtil::containsSeries(*approxKnn, query)) cout<<"WARNING!  Approximate Search does not find the target series! "<<endl;
////                if(!TimeSeriesUtil::containsSeries(*exactKnn, query)) cout<<"ERROR!  EXACT Search does not find the target series! "<<endl;
//                vector<PqItemSeries*> &intersection = TimeSeriesUtil::intersectionTsSets(approxKnn, gt);
////                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(*approxKnn, *exactKnn);
////                cout<<"Recall Number = " << recallNums[curRound] << endl;
////                cout << "Exact result hits ground Truth:" << TimeSeriesUtil::intersectionTsSetsCardinality(*exactKnn, gt)<<endl;
//                cout << "Approximate result hist ground Truth:" << TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, gt) << endl;
//
//
////                for(int id=0;id<recallNums[curRound];++id)
////                    approxKnn[id] = intersection[id];
//                analyzeResult(g, intersection, *gt, k, query);
//                cout << "****************Ground Truth is followed : *********************"<<endl;
//                printResult2(*gt, query);
//                cout << "****************Approx Result is followed : *********************"<<endl;
//                printResult(*approxKnn);
//            }
//
//            cout<<"------------------Experiment Finished--------------------" << endl;
//            //                cout<<"k: " + k);
//            //                cout<<"threshold: " + threshold);
//            int totalRecallNum = 0;
//            for(int temp:recallNums){
//                cout << temp  <<  "," ;
//                totalRecallNum += temp;
//            }
//            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<<endl;
//        }
//
//    }
//    fclose(f);
//}
//
//void Recall::doExprPerformance(Graph* g){
//    int maxExprRound = 2000;
//    int lastIndex= g->rowDataFileName.rfind(".bin");
//    string fn = g->rowDataFileName.substr(0, lastIndex) + "_query.bin";
//    int ks[]{10,50,100,250,500};
//    int thresholds[]{1000,5000,10000,25000,50000};
//    float *query;
//    FILE *f = fopen(fn.c_str(), "rb");
//    long sum = 0 ;
//
//    for(int threshold:thresholds){
//        for(int k:ks){
//            //            cout<<"k: " << k << endl;
//            //            cout<<"threshold: " << threshold<< endl;
//            //            int recallNums[maxExprRound];
//            for(int curRound = 0; curRound < maxExprRound; ++curRound){
//                query = FileUtil::readSeries(f);
//
//                auto t1 = chrono::system_clock::now();
//                vector<PqItemSeries*> *approxKnn = ApproximateSearcher::approxKnnSearchNew(g, query, k, threshold);
//                auto t2 = chrono::system_clock::now();
//                sum += chrono::duration_cast<chrono::microseconds>(t2 - t1).count();
//                free_heap(approxKnn);
//                delete[] query;
//
//                //                start = clock();
//                //                vector<PqItemSeriesVector> &exactKnn = ExactSearcher::exactKnnSearch(g, query, k);
//                //                exactSearchTime += (clock() - start);
//
//                //                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
//            }
//
//            //            cout<<"------------------Experiment Finished--------------------" << endl;
//            //                cout<<"k: " + k);
//            //                cout<<"threshold: " + threshold);
//            //            int totalRecallNum = 0;
//            //            for(int temp:recallNums){
//            //                cout << temp << ",";
//            //                totalRecallNum += temp;
//        }
//        rewind(f);
//        //            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
//        //            rewind(f);
//    }
//    cout << sum<<endl;
//    fclose(f);
//
//
//}
//void Recall::doExprPerformance(Graph* g){
//    int maxExprRound = 1000;
//    int lastIndex= g->rawDataFileName.rfind(".bin");
//    string fn = g->rawDataFileName.substr(0, lastIndex) + "_query.bin";
//    long wcStart = clock();
//    int ks[]{10 };
//    int thresholds[]{1000};
//    float *query;
//    FILE *f = fopen(fn.c_str(), "rb");
//
//    for(int threshold:thresholds){
//        for(int k:ks){
//            //            cout<<"k: " << k << endl;
//            //            cout<<"threshold: " << threshold<< endl;
//            //            int recallNums[maxExprRound];
//            for(int curRound = 0; curRound < maxExprRound; ++curRound){
//                query = FileUtil::readSeries(f);
//
//                long start = clock();
//                vector<PqItemSeriesVector> &approxKnn = ApproximateSearcher::approxKnnSearch(g, query, k, threshold);
//                approxSearchTime += (clock() - start);
//
//                //                start = clock();
//                //                vector<PqItemSeriesVector> &exactKnn = ExactSearcher::exactKnnSearch(g, query, k);
//                //                exactSearchTime += (clock() - start);
//
//                //                recallNums[curRound] = TimeSeriesUtil::intersectionTsSetsCardinality(approxKnn, exactKnn);
//            }
//            wallClockTime += (clock() - wcStart);
//            //            cout<<"------------------Experiment Finished--------------------" << endl;
//            //                cout<<"k: " + k);
//            //                cout<<"threshold: " + threshold);
//            //            int totalRecallNum = 0;
//            //            for(int temp:recallNums){
//            //                cout << temp << ",";
//            //                totalRecallNum += temp;
//        }
//        //            cout<<"\nAverage Recall rate is : " << (float)totalRecallNum / (float) (maxExprRound * k)<< endl;
//        //            rewind(f);
//    }
//
//    fclose(f);
//
//    cout << "Wall Clock Search Time is " << (double)wallClockTime / CLOCKS_PER_SEC << " s." << endl;
//    cout << "Approximate Search Time is " << (double)approxSearchTime / CLOCKS_PER_SEC << " s." << endl;
//    //            cout << "Exact Search Time is " << exactSearchTime << " s." << endl;
//    cout << endl << "In Approximate Search Process : "<< endl;
//    cout << "1. The process that determines which vertexes(candidates) to search with/out threshold spent " << (double )determineCandidatesTime / CLOCKS_PER_SEC<< " s."<<endl;
//    cout << "2. The operations that copy data into different STLs in order to satisfy the result format spent " <<(double)dataCopyTime / CLOCKS_PER_SEC<< " s."<<endl;
//    cout << "3. Searching candidate vertexes process(including multiple level vertexes) spent " << (double)searchTime / CLOCKS_PER_SEC << " s."<<endl;
//    cout << "\t In the process of searching candidate vertexes: " << endl;
//    cout << "\t 3.1. Small Vertexes: Search Number = " << approxSearchUnits[0] << " , "<< "Search Thresholds Sum = "
//    << approxSearchUnits[1] <<" , Search Time = " << (double )approxSearchTimeDetails[0] / CLOCKS_PER_SEC << " s."<<endl;
//    cout << "\t 3.2. Middle Vertexes: Search Number = " << approxSearchUnits[2] << " , "<< "Search Thresholds Sum = "
//    << approxSearchUnits[3] <<" , Search Time = " << (double )approxSearchTimeDetails[1]/ CLOCKS_PER_SEC << " s."<<endl;
//    cout << "\t 3.3. Big Vertexes: Search Number = " << approxSearchUnits[4] << " , "<< "Search Thresholds Sum = "
//    << approxSearchUnits[5] <<" , Search Time = " << (double )approxSearchTimeDetails[2]/ CLOCKS_PER_SEC << " s."<<endl;
//    cout<<"\t\t 3.3.1 With threshold: ApproxSearchTime = " << (double )DSTreeNodeApproxTime2 / CLOCKS_PER_SEC << " s, Exact Search Time = " << (double )DSTreeSearchTime2 / CLOCKS_PER_SEC <<" s." <<endl
//    << "\t\t\t 3.3.1.1 During Exact Search, Prepare PQ Time = " << (double )DSTreePreparePqTime2 / CLOCKS_PER_SEC << "s, "
//    << "Process Priority Queue Time = " << (double )DSTreeProcessPqTime2 / CLOCKS_PER_SEC << "s. "<<endl;
//    cout<<"\t\t\t\t 3.3.1.1.1 During Prepare PQ, calculating LB of node costs " << (double )LBNodeTime2/CLOCKS_PER_SEC << "s, "
//    <<"LB of series costs "<<(double )LBSeriesTime2 / CLOCKS_PER_SEC <<" s."<<endl;
//    cout <<"\t\t\t\t 3.3.1.1.2 During Process PQ, IO time(read series) = " << (double ) IOBigTime2 / CLOCKS_PER_SEC << " s, "
//    <<"Distance Computation costs " << (double )distBigTime2 / CLOCKS_PER_SEC << " s."<<endl;
//
//    cout<<"\t\t 3.3.2 Without threshold, Approximate Search Time = " << (double )DSTreeNodeApproxTime / CLOCKS_PER_SEC << " s, "
//    <<"Exact Search Time = " << (double )DSTreeProcessPqTime / CLOCKS_PER_SEC << " s." <<endl;
//    cout <<"\t\t\t 3.3.2.1 During Exact Search, Calculating LB of node costs " << (double )LBNodeTime / CLOCKS_PER_SEC << " s,"
//    <<"Adjust heap time = " << (double )heapBigTime / CLOCKS_PER_SEC << "s, Distance Computation costs " << (double )distBigTime / CLOCKS_PER_SEC << " s."<<endl;
//
//    cout << "\t\t 3.3.3 Hash Table Search Time = " << (double )hashTableTime / CLOCKS_PER_SEC << "s. LB Series Count = " << LBSeriesCnt <<endl;
//
//}
//void analyzeResult2(Graph* g, vector<PqItemSeries*> &approx, vector<PqItemSeries*> &exact, int k, float*query){
//    cout << "\n\t Now report deep level graph search result."<<endl;
//    assert(exact.size() == k);
//    TimeSeries ts(query);
//    int invSax = SaxUtil::invSaxHead2FromSax(ts.sax, TimeSeries::bitsCardinality, TimeSeries::segmentNum);
//    int approxWords[approx.size()];
//    int exactWords[k];
//    int dest = k + approx.size();
//    for(int id=0;id<k;++id) {
//        auto *tsss = new float[TimeSeries::tsLength];
//        for(int j=0;j<TimeSeries::tsLength;++j)
//            tsss[j] = (exact[id]->ts)[j];
//        exactWords[id] = SaxUtil::invSaxHead2FromSax((new TimeSeries(tsss))->sax, TimeSeries::bitsCardinality, TimeSeries::segmentNum);
//        if(id < approx.size()) {
//            auto *tss = new float[TimeSeries::tsLength];
//            for(int j=0;j<TimeSeries::tsLength;++j)
//                tss[j] = approx[id]->ts[j];
//            approxWords[id] = SaxUtil::invSaxHead2FromSax(
//                    (new TimeSeries(tss))->sax, TimeSeries::bitsCardinality, TimeSeries::segmentNum);
//        }
//    }
//    vector<Vertex *> &nnList = g->NeighborsList[invSax];
//    int accum = 0;
//    cout<<"\n==================Layer 0:===================="<<endl;
//    int arn = count(approxWords, invSax, approx.size()), ern = count(exactWords, invSax, k);
//    SaxUtil::printBinary(invSax, TimeSeries::segmentNum);
//    cout<<": Count = " << g->VERTEXES[invSax]->cnt << " ; Approx Result Number = " << arn << " ; Exact Result Number = " << ern;
//    accum += arn + ern;
//    int pos = 0;
//    for(int id=1; id<= Const::bitsReserve; ++id){
//        cout<<"\n==================Layer " << id << " :========================"<<endl;
//        accum = reportKit(nnList, pos, pos + MathUtil::nChooseK(TimeSeries::segmentNum, id), approxWords, exactWords, accum, dest, approx.size(), k);
//        if(accum >= dest) return;
//        pos += MathUtil::nChooseK(TimeSeries::segmentNum, id);
//    }
//    cout << "\n\t Deep level graph search report finished.\n"<<endl;
//}
