//
// Created by wzy on 2021/8/25.
//

#ifndef MULGIFT_CONST_H
#define MULGIFT_CONST_H
#include <string>
#include <iostream>
#include <sys/time.h>
#include "../include/Utils/INIReader.h"




using namespace std;
class Const {

public:
//    const static int segmentNum = 16, vertexNum = 1 << segmentNum, bitsReserve = 3;
//
//    constexpr const static char graphFileName[]{"../data/RowGraph_16_3.bin"};
//    const static int th = 10000, max_radius = 6, maxK = 100, max_replica = 4;
//    constexpr const static double boundary_1st = 0.2, boundary = 0.15, filling_factor_1st = 0.8, filling_factor = 0.5;

    // sec:expr
    static string dataset, method;
    static int tsLength, maxK, index, ops, materialized, method_code, query_num, series_num, k, dtw_window_size,
    batch_size, batch_num, pre_read;
    static double dtw_window_percent;

    //sec: parameter
    // segment number is special, it needs to be specified here
    const static int segmentNum = 16;
    const static long small_file_threshold = 32*1024;   // 32KB
    static int th, bitsCardinality, max_replica, max_radius, fbl_size, max_diff, fbl_series_num;
    static double imbalance, boundary_1st, boundary, filling_factor_1st, filling_factor, small_perc, f_low, f_high,
    alpha, max_mask_bit_percentage, tardis_sample_percent;

    //sec: others
    static string graphfn;
    static int bitsReserve;

    //sec: dataset
    static string paafn, saxfn, idxfn, fidxfn, posidxfn, fuzzyidxfn, datafn, queryfn, resfn, dtwresfn, dstreefn, tardisfn;

    // 2-nd level parameter

    static int tsLengthPerSegment, cardinality, tsLengthBytes, vertexNum, neighborNum;
    static long offset;

    static void readConfig(){
        INIReader reader("../config.ini");

        if (reader.ParseError() < 0) {
            cout << "Can't load '.ini'\n";
            return;
        }
        dataset = reader.Get("expr", "dataset","");
        cout << "dataset: " << dataset<<endl;
        if(dataset.empty())   exit(-1);

        method = reader.Get("expr", "method","");
        cout << "method: " << method<<endl;

        index = reader.GetInteger("expr", "index",-1);
        cout << "index: " << index<<endl;

        materialized = reader.GetInteger("expr", "materialized",-1);
        cout << "materialized: " << materialized<<endl;

        ops = reader.GetInteger("expr", "ops",-1);
        cout << "ops: " << ops<<endl;

        maxK = reader.GetInteger("expr", "maxK", -1);
        cout << "maxK: " << maxK << endl;
        if(maxK == -1)  exit(-1);

        query_num = reader.GetInteger("expr", "query_num", -1);
        cout << "query_num: " << query_num << endl;

        series_num = reader.GetInteger("expr", "series_num", -1);
        cout << "series_num: " << series_num << endl;

        k = reader.GetInteger("expr", "k", -1);
        cout << "k: " << k << endl;

        batch_size = reader.GetInteger("expr", "batch_size", -1);
        cout << "batch_size: " << batch_size << endl;

        batch_num = reader.GetInteger("expr", "batch_num", -1);
        cout << "batch_num: " << batch_num << endl;

        pre_read = reader.GetInteger("expr", "pre_read", -1);
        cout << "pre_read: " << pre_read << endl;

        dtw_window_percent = reader.GetReal("expr", "dtw_window_percent", -1);
        cout << "dtw_window_percent: " << dtw_window_percent << endl;

        th = reader.GetInteger("parameter", "th", 10000);
        cout << "th: " << th << endl;

        fbl_size = reader.GetInteger("parameter", "fbl_size", 1024 * 4);
        cout << "fbl_size: " << fbl_size << endl;

        max_diff = reader.GetInteger("parameter", "max_diff", 2);
        cout << "max_diff: " << max_diff << endl;

//        segmentNum = reader.GetInteger("parameter", "segmentNum", 16);
//        cout << "segmentNum: " << segmentNum << endl;

        bitsCardinality = reader.GetInteger("parameter", "bitsCardinality", 8);
        cout << "bitsCardinality: " << bitsCardinality << endl;

        boundary_1st = reader.GetReal("parameter", "boundary_1st", 0.2);
        cout << "boundary_1st: " << boundary_1st << endl;

        boundary = reader.GetReal("parameter", "boundary", 0.3);
        cout << "boundary: " << boundary << endl;

        imbalance = reader.GetReal("parameter", "imbalance", 0.3);
        cout << "imbalance: " << imbalance << endl;

        small_perc = reader.GetReal("parameter", "small_perc", 0.3);
        cout << "small_perc: " << small_perc << endl;

        f_low = reader.GetReal("parameter", "f_low", 0);
        cout << "f_low: " << f_low << endl;

        f_high = reader.GetReal("parameter", "f_high", 0);
        cout << "f_high: " << f_high << endl;

        alpha = reader.GetReal("parameter", "alpha", 0);
        cout << "alpha: " << alpha << endl;

        max_mask_bit_percentage = reader.GetReal("parameter", "max_mask_bit_percentage", 0);
        cout << "max_mask_bit_percentage: " << max_mask_bit_percentage << endl;

        tardis_sample_percent = reader.GetReal("parameter", "tardis_sample_percent", 0);
        cout << "tardis_sample_percent: " << tardis_sample_percent << endl;

        max_replica = reader.GetInteger("parameter", "max_replica", 10000);
        cout << "max_replica: " << max_replica << endl;

        filling_factor_1st = reader.GetReal("parameter", "filling_factor_1st", 0.8);
        cout << "filling_factor_1st: " << filling_factor_1st << endl;

        filling_factor = reader.GetReal("parameter", "filling_factor", 0.5);
        cout << "filling_factor: " << filling_factor << endl;

        max_radius = reader.GetInteger("parameter", "max_radius", 6);
        cout << "max_radius: " << max_radius << endl;

        graphfn = reader.Get("other", "graphfn","");
        cout << "graphfn: " << graphfn <<endl;

        bitsReserve = reader.GetInteger("other", "bitsReserve", 6);
        cout << "bitsReserve: " << bitsReserve << endl;

        tsLength = reader.GetInteger(dataset, "tsLength", -1);
        cout << "tsLength: " << tsLength << endl;
        if(tsLength == -1)  exit(-1);

        paafn = reader.Get(dataset, "paafn","");
        cout << "paafn: " << paafn <<endl;

        saxfn = reader.Get(dataset, "saxfn","");
        cout << "saxfn: " << saxfn <<endl;

        idxfn = reader.Get(dataset, "idxfn","");
        cout << "idxfn: " << idxfn <<endl;

        fidxfn = reader.Get(dataset, "fidxfn","");
        cout << "fidxfn: " << fidxfn <<endl;

        posidxfn = reader.Get(dataset, "posidxfn","");
        cout << "posidxfn: " << posidxfn <<endl;

        fuzzyidxfn = reader.Get(dataset, "fuzzyidxfn","");
        cout << "fuzzyidxfn: " << fuzzyidxfn <<endl;

        datafn = reader.Get(dataset, "datafn","");
        cout << "datafn: " << datafn <<endl;

        queryfn = reader.Get(dataset, "queryfn","");
        cout << "queryfn: " << queryfn <<endl;

        resfn = reader.Get(dataset, "resfn","");
        cout << "resfn: " << resfn <<endl;

        dtwresfn = reader.Get(dataset, "dtwresfn","");
        cout << "dtwresfn: " << dtwresfn <<endl;

        dstreefn = reader.Get(dataset, "dstreefn","");
        cout << "dstreefn: " << dstreefn <<endl;

        tardisfn = reader.Get(dataset, "tardisfn","");
        cout << "tardisfn: " << dstreefn <<endl;

        tsLengthPerSegment = tsLength / segmentNum;
        cardinality = 1<<bitsCardinality;
        tsLengthBytes = tsLength * 4;
        offset = ((long )(cardinality - 1) * (cardinality - 2)) / 2;
        vertexNum = 1 << segmentNum;
        dtw_window_size = dtw_window_percent * Const::tsLength;
        fbl_series_num = (fbl_size * 1024/ tsLengthBytes) * 1024 ;
        int max_series_num = numeric_limits<unsigned>::max() / Const::tsLength;
        fbl_series_num = fbl_series_num > max_series_num ? max_series_num: fbl_series_num;
        if(Const::method == "range")    method_code = 0;
        else if(Const::method == "stdev")   method_code =1;
        else if(Const::method == "dist")    method_code = 2;

        cout << "buffer_series_num: " << fbl_series_num << endl;
        neighborNum = 0;
        for(int i=1;i<=bitsReserve;++i){
            neighborNum += nChooseK(Const::segmentNum, i);
        }
        cout << "Neighbor number : " << neighborNum << endl;

    }

    static int nChooseK(int n, int r) {
        r = min(r, (n - r));
        if (r <= 1) {  // C(n, 0) = 1, C(n, 1) = n
            return r == 0 ? 1 : n;
        }
        int limit = numeric_limits<int>::max() >> (31 - n);
        int cnk = 0;
        for (int i = 3; i < limit; i++) {
            if (bitCount(i) == r) {
                cnk++;
            }
        }
        return cnk;
    }

    static int bitCount(int input){
        int res = 0;
        while (input!=0) {
            res += (input % 2);
            input /= 2;
        }
        return res;
    }

    static string now() {
        time_t t = time(nullptr);
        char buffer[9] = {0};
        strftime(buffer, 9, "%H:%M:%S", localtime(&t));
        return (buffer);
    }

    static void logPrint(const string& content){
        cout << now() << ": " << content <<endl;
    }

    static void timer_start(struct timeval *timer)
    {
        (void)gettimeofday(timer, nullptr);
    }

    static double timer_end(const struct timeval *start)
    {
        struct timeval now{};
        (void)gettimeofday(&now, NULL);
        return (now.tv_sec - start->tv_sec)*1000000 + ((double)now.tv_usec - start->tv_usec) ;
    }

};


#endif //MULGIFT_CONST_H
