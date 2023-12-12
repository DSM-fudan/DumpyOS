//
// Created by wzy on 2021/12/12.
//
#include "../include/Const.h"

// sec:expr
string Const::dataset = "",Const:: method = "";
int Const::tsLength = -1, Const::maxK = -1, Const::index = -1, Const::ops = -1, Const::materialized = -1,
Const::method_code = -1, Const::query_num = -1, Const::series_num = -1, Const::k = -1, Const::dtw_window_size = -1,
Const::batch_size = -1, Const::batch_num = -1, Const::pre_read = -1, Const::thread_num = -1, Const::messi_pq_num = -1,
Const::SSD_pq_num = -1;
double Const::dtw_window_percent = -1;

//sec: parameter
// segment number is special, it needs to be specified here
int Const::th = -1, Const::bitsCardinality = -1, Const::max_replica = -1, Const::max_radius = -1, Const::fbl_size = -1,Const::max_diff = -1, Const::fbl_series_num = -1;
double Const::imbalance = -1, Const::boundary_1st = -1, Const::boundary = -1, Const::filling_factor_1st = -1,
Const::filling_factor = -1, Const::small_perc = -1, Const::f_low = -1, Const::f_high = -1, Const::alpha = -1,
Const::max_mask_bit_percentage = -1, Const::tardis_sample_percent = -1;

//sec: others
string Const::graphfn = "";
int Const::bitsReserve = -1;

//sec: dataset
string Const::paafn = "", Const::saxfn = "", Const::idxfn = "",Const::fidxfn = "",Const::posidxfn = "",
Const::fuzzyidxfn, Const::datafn = "", Const::queryfn = "", Const::resfn = "",Const::dtwresfn = "", Const::dstreefn = "", Const::tardisfn = "";

// 2-nd level parameter

int Const::tsLengthPerSegment = -1, Const::cardinality = -1, Const::tsLengthBytes = -1, Const::vertexNum = -1, Const::neighborNum = 0;
long Const::offset = -1;