//
// Created by wzy on 2022/1/13.
//

#include <thread>
#include <cassert>
#include <cmath>
#include <bitset>
#include <unordered_map>
#include "../../include/DataStructures/FADASNode.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Utils/MathUtil.h"
#include "../../include/Utils/TimeSeriesUtil.h"
#include "../../include/Utils/SaxUtil.h"

unsigned short *FADASNode::saxes = nullptr;
float *FADASNode::paas = nullptr;
int *** FADASNode::combines = nullptr;
int * FADASNode::combine_num = nullptr;
int FADASNode::a = MathUtil::nChooseK(Const::segmentNum, 1), FADASNode::b = MathUtil::nChooseK(Const::segmentNum, 2), FADASNode::c = MathUtil::nChooseK(Const::segmentNum, 3);
const int FADASNode::power_2[]{1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536};
extern long SAX_PAA_CPU_TIME , SAX_PAA_READ_TIME ;
extern long LB_SERIES_TIME, HEAP_TIME;
extern double DIST_CALC_TIME,READ_TIME;
static long MAT1_CPU_TIME_STAT = 0, MAT1_CPU_TIME_COPY = 0, MAT1_WRITE_TIME = 0, MAT1_TOTAL_TIME = 0, MAT1_READ_TIME = 0,
        MAT2_CPU_TIME = 0, MAT2_WRITE_TIME = 0, MAT2_TOTAL_TIME = 0, SAX_PAA_TOTAL_TIME = 0, MAT2_READ_TIME = 0,
        GROW_CPU_TIME = 0,  GROW_CPU_TIME_1st = 0, GROW_TOTAL_TIME = 0, SMALL_FILES_BYTES_WRITE = 0, SMALL_FILES_BYTES_READ = 0,
        RAND_READ_CNT = 0, RAND_WRITE_CNT = 0, SEQ_READ_CNT = 0, SEQ_WRITE_CNT = 0,
        SAX_WRITE_TIME = 0;

void FADASNode::loadCombines(){
    string base = "../combines/" + to_string(Const::segmentNum) + "-";
    auto ret = new int**[Const::segmentNum + 1];
    combine_num = new int[Const::segmentNum];
    ifstream ff("../combines/cnum-"+ to_string(Const::segmentNum) + ".txt", ios::in);
    for(int i=0;i<Const::segmentNum;++i){
        ff >> combine_num[i];
    }
    ff.close();

    for(int i=1;i<=Const::segmentNum - 1;++i){
        ret[i] = new int*[combine_num[i]];
        ifstream f(base + to_string(i) + ".txt", ios::in);
        for(int j=0;j<combine_num[i];++j){
            ret[i][j] = new int[i];
            for(int k=0;k<i;++k) {
                f >> ret[i][j][k];
//                cout << ret[i][j][k] << ",";
            }
//            cout << endl;
        }
//        cout << endl;
        f.close();
    }
    ret[Const::segmentNum] = new int*[1];
    ret[Const::segmentNum][0] = new int[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i){
        ret[Const::segmentNum][0][i] = i;
    }
    combines = ret;

}

FADASNode* FADASNode::route(const unsigned short *_sax){
    if(!isInternalNode())
        return this;
    int nav_id = SaxUtil::extendSax(_sax, bits_cardinality, chosenSegments);
    if(children[nav_id] == nullptr) return this;
    return children[nav_id]->route(_sax);
}

void FADASNode::routeDuringInsertion(const unsigned short *_sax, int pos){
    ++size;
    if(size > Const::f_high * Const::th * (1 << chosenSegments.size())){
        pos_cache.push_back(pos);
        return;
    }
    if(children.empty()){
        pos_cache.push_back(pos);
        return;
    }
    int nav_id = SaxUtil::extendSax(_sax, bits_cardinality, chosenSegments);
    if(children[nav_id] == nullptr){
        children[nav_id] = new FADASNode(this, nav_id, true);
        generateSaxAndCardinality(children[nav_id], nav_id);
        children[nav_id]->pos_cache.push_back(pos);
    }
    return children[nav_id]->routeDuringInsertion(_sax, pos);
}

FADASNode* FADASNode::route1step(const unsigned short *_sax){
    assert(!isLeafNode() && !isLeafPack());
    int nav_id;
    if(layer >= 1)
        nav_id = SaxUtil::extendSax(_sax, bits_cardinality, chosenSegments);
    else
        nav_id = SaxUtil::invSaxHeadFromSax(_sax, Const::bitsCardinality, Const::segmentNum);
    return children[nav_id];
}

extern int _search_num, layer;
void FADASNode::search(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir,
                       float *query_reordered, int *ordering) const{
    assert(!isInternalNode());
    if(query_reordered == nullptr || ordering == nullptr) {
        search(k, queryTs, heap, index_dir);
        return;
    }
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    _search_num += size;
    string fn = index_dir+ getFileName();
    if(partition_id == -1)  fn += "_L";
//    long fs = FileUtil::getFileSize(fn.c_str());
//    int series_num = fs / Const::tsLengthBytes;
//    assert(series_num == size);

    FILE *f = fopen(fn.c_str(), "rb");
    struct timeval io{};
    Const::timer_start(&io);
    auto *ts = new float[size * Const::tsLength];
//    for(int i=0;i<size;++i)
//        fread(ts + i * Const::tsLength, sizeof(float), Const::tsLength, f);
    fread(ts, sizeof(float), size * Const::tsLength, f);
    READ_TIME += Const::timer_end(&io);

    _search_num += size; ::layer = FADASNode::layer;
    for(int i=0;i<size;++i){
//        struct timeval start{};
//        Const::timer_start(&start);
//        double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts + i * Const::tsLength, Const::tsLength, bsf);
        double dist = TimeSeriesUtil::euclideanDist(query_reordered, ts + i * Const::tsLength, Const::tsLength, bsf, ordering);
//        DIST_CALC_TIME += Const::timer_end(&start);

        if(heap.size() < k){
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
    delete[]ts;
    fclose(f);
}

//void FADASNode::search_SIMD(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir,
//                       float *query_reordered, int *ordering) const{
//    assert(!isInternalNode());
//    if(query_reordered == nullptr || ordering == nullptr) {
//        search(k, queryTs, heap, index_dir);
//        return;
//    }
//    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
//    _search_num += size;
//    string fn = index_dir+ getFileName();
//    if(partition_id == -1)  fn += "_L";
////    long fs = FileUtil::getFileSize(fn.c_str());
////    int series_num = fs / Const::tsLengthBytes;
////    assert(series_num == size);
//
//    FILE *f = fopen(fn.c_str(), "rb");
//    struct timeval io{};
//    Const::timer_start(&io);
//    auto *ts = new float[size * Const::tsLength];
////    for(int i=0;i<size;++i)
////        fread(ts + i * Const::tsLength, sizeof(float), Const::tsLength, f);
//    fread(ts, sizeof(float), size * Const::tsLength, f);
//    READ_TIME += Const::timer_end(&io);
//
//    _search_num += size; ::layer = FADASNode::layer;
//    for(int i=0;i<size;++i){
////        struct timeval start{};
////        Const::timer_start(&start);
////        double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts + i * Const::tsLength, Const::tsLength, bsf);
//        double dist = TimeSeriesUtil::euclideanDist_SIMD(query_reordered, ts + i * Const::tsLength, Const::tsLength, bsf, ordering);
////        DIST_CALC_TIME += Const::timer_end(&start);
//
//        if(heap.size() < k){
//            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
//            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
//        }else if(dist < bsf){
//            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
//            delete heap.back();
//            heap.pop_back();
//            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
//            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
//        }
//
//        if(heap.size() >= k)    bsf = heap[0]->dist;
//    }
//
//    for(PqItemSeries*s: heap){
//        if(s->needDeepCopy) s->copyData();
//    }
//    delete[]ts;
//    fclose(f);
//}

void FADASNode::search(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir) const{
    assert(!isInternalNode());
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    string fn = index_dir+ getFileName();
    if(partition_id == -1)  fn += "_L";
//    long fs = FileUtil::getFileSize(fn.c_str());
//    int series_num = fs / Const::tsLengthBytes;
//    assert(series_num == size);

    FILE *f = fopen(fn.c_str(), "rb");
    struct timeval io{};
    Const::timer_start(&io);
    auto *ts = new float[size * Const::tsLength];
//    for(int i=0;i<size;++i)
//        fread(ts + i * Const::tsLength, sizeof(float), Const::tsLength, f);
    fread(ts, sizeof(float), size * Const::tsLength, f);
    READ_TIME += Const::timer_end(&io);

    _search_num += size; ::layer = FADASNode::layer;
    for(int i=0;i<size;++i){
//        struct timeval start{};
//        Const::timer_start(&start);
        double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts + i * Const::tsLength, Const::tsLength, bsf);
//        double dist = TimeSeriesUtil::euclideanDist(query_reordered, ts + i * Const::tsLength, Const::tsLength, bsf, ordering);
//        DIST_CALC_TIME += Const::timer_end(&start);

        if(heap.size() < k){
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
    delete[]ts;
    fclose(f);
}


void FADASNode::searchDTW(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir) const{
    assert(!isInternalNode());
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    string fn = index_dir+ getFileName();

//    long fs = FileUtil::getFileSize(fn.c_str());
//    int series_num = fs / Const::tsLengthBytes;
//    assert(series_num == size);

    FILE *f = fopen(fn.c_str(), "rb");
    struct timeval io{};
    Const::timer_start(&io);
    auto *ts = new float[size * Const::tsLength];
//    for(int i=0;i<size;++i)
//        fread(ts + i * Const::tsLength, sizeof(float), Const::tsLength, f);
    fread(ts, sizeof(float), size * Const::tsLength, f);
    READ_TIME += Const::timer_end(&io);

    _search_num += size; ::layer = FADASNode::layer;
    for(int i=0;i<size;++i){
//        struct timeval start{};
//        Const::timer_start(&start);
        double dist = TimeSeriesUtil::dtw(queryTs->ts, ts + i * Const::tsLength, Const::tsLength,Const::dtw_window_size, bsf);
//        DIST_CALC_TIME += Const::timer_end(&start);

        if(heap.size() < k){
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
    delete[]ts;
    fclose(f);
}

void FADASNode::searchLessPack(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir) const{
    assert(!isInternalNode());
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    string fn = index_dir+ getFileNamePack();

//    long fs = FileUtil::getFileSize(fn.c_str());
//    int series_num = fs / Const::tsLengthBytes;
//    assert(series_num == size);

    FILE *f = fopen(fn.c_str(), "rb");
    struct timeval io{};
    Const::timer_start(&io);
    auto *ts = new float[size * Const::tsLength];
    for(int i=0;i<size;++i)
        fread(ts + i * Const::tsLength, sizeof(float), Const::tsLength, f);
//    fread(ts, sizeof(float), size * Const::tsLength, f);
    READ_TIME += Const::timer_end(&io);

    _search_num += size; ::layer = FADASNode::layer;
    for(int i=0;i<size;++i){
//        struct timeval start{};
//        Const::timer_start(&start);
        double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts + i * Const::tsLength, Const::tsLength, bsf);
//        DIST_CALC_TIME += Const::timer_end(&start);

        if(heap.size() < k){
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
    delete[]ts;
    fclose(f);
}

void FADASNode::search(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir,unordered_set<float*, createhash, isEqual>*hash_set) const{
    assert(!isInternalNode());
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    string fn = index_dir+ getFileName();

    FILE *f = fopen(fn.c_str(), "rb");
    auto *ts = new float[size * Const::tsLength];
    fread(ts, sizeof(float), size * Const::tsLength, f);

    _search_num += size; ::layer = FADASNode::layer;
    for(int i=0;i<size;++i){
        if(hash_set->find(ts + i * Const::tsLength) != hash_set->end()) continue;
        double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts + i * Const::tsLength, Const::tsLength, bsf);

        if(heap.size() < k){
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            hash_set->insert(ts + i * Const::tsLength);
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            hash_set->erase(heap.back()->ts);
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            hash_set->insert(ts + i * Const::tsLength);
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) {
            hash_set->erase(s->ts);
            s->copyData();
            hash_set->insert(s->ts);
        }
    }
    delete[]ts;
    fclose(f);
}

void FADASNode::search(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir,unordered_set<float*, createhash, isEqual>*hash_set,
                       float *query_reordered, int *ordering) const{
    assert(!isInternalNode());
    if(query_reordered == nullptr || ordering == nullptr) {
        search(k, queryTs, heap, index_dir, hash_set);
        return;
    }

    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    string fn = index_dir+ getFileName();

    FILE *f = fopen(fn.c_str(), "rb");
    auto *ts = new float[size * Const::tsLength];
    fread(ts, sizeof(float), size * Const::tsLength, f);

    _search_num += size; ::layer = FADASNode::layer;
    for(int i=0;i<size;++i){
        if(hash_set->find(ts + i * Const::tsLength) != hash_set->end()) continue;
        double dist = TimeSeriesUtil::euclideanDist(query_reordered, ts + i * Const::tsLength, Const::tsLength, bsf, ordering);

        if(heap.size() < k){
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            hash_set->insert(ts + i * Const::tsLength);
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            hash_set->erase(heap.back()->ts);
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(ts + i * Const::tsLength, dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            hash_set->insert(ts + i * Const::tsLength);
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) {
            hash_set->erase(s->ts);
            s->copyData();
            hash_set->insert(s->ts);
        }
    }
    delete[]ts;
    fclose(f);
}

void FADASNode::search_offset(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir) const{
    assert(!isInternalNode());
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    string fn = index_dir+ getFileName();
    FILE *f = fopen(fn.c_str(), "rb");
    auto *tss = new float[size * (Const::tsLength + 1)];
    struct timeval io{};
    Const::timer_start(&io);
//    for(int i=0;i<size;++i)
//        fread(tss + i * (Const::tsLength + 1), sizeof(float), Const::tsLength + 1, f);
    fread(tss, sizeof(float), size * (Const::tsLength + 1), f);
    READ_TIME += Const::timer_end(&io);
    ::layer = FADASNode::layer;

    for(int i=0;i<size;++i){
        if(heap.size() >= k) {
            int offset = *((int *) (tss + (i + 1) * (Const::tsLength + 1) - 1));
//            auto st = chrono::system_clock::now();
            double lb = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, saxes + offset * Const::segmentNum);
//            auto en = chrono::system_clock::now();
//            LB_SERIES_TIME += chrono::duration_cast<chrono::microseconds>(en - st).count();
            if (lb >= bsf) continue;
        }
        _search_num++;
        struct timeval start{};
        Const::timer_start(&start);
        double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, tss + i * (Const::tsLength + 1), Const::tsLength, bsf);
        DIST_CALC_TIME += Const::timer_end(&start);

        if(heap.size() < k){
            heap.push_back(new PqItemSeries(tss + i * (Const::tsLength + 1), dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(tss + i * (Const::tsLength + 1), dist, false, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
//        s = chrono::system_clock::now();
//        HEAP_TIME += chrono::duration_cast<chrono::microseconds>(s - e).count();;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
    delete[]tss;
    fclose(f);

}

// put actual series into disk file of nodes in 1st layer
void materialize1stLayer(string datafn, FADASNode* root, int *navids, string index_dir){
    auto start_t = chrono::system_clock::now();
    Const::logPrint("Start move data to disk file in 1st layer.");
    FILE *f = fopen(datafn.c_str(), "r");
    long rest = root->size, total = root->size, cur = 0;
    unordered_map<FADASNode*, FBL_UNIT>fbl;

    RAND_READ_CNT++;
    SEQ_READ_CNT += rest;

    // There is another implementation method that fbl in each node stores a pointer vector where each pointer points to a series.
    // This may incur many write calls.
    while(rest > 0){
        fbl.clear();
        long num;
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto end = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto start = chrono::system_clock::now();
        MAT1_READ_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

        // statistic the size of each node fbl size, and allocate memory for them
        for(long i=cur;i<cur+num;++i)
            fbl[root->children[navids[i]]].size++;
        for(auto &iter:fbl)
            iter.second.buffer = new float [(long)iter.second.size * Const::tsLength];
        end = chrono::system_clock::now();
        MAT1_CPU_TIME_STAT += chrono::duration_cast<chrono::microseconds>(end - start).count();

        //copy series to node to ensure write is only called once upon a node fbl
        for(long i = cur; i<cur+num;++i){
            FBL_UNIT* fbl_node = &fbl[root->children[navids[i]]];
            copy(tss + (i-cur) * Const::tsLength, tss + (i+1-cur)* Const::tsLength, fbl_node->buffer + (long)fbl_node->pos++ * Const::tsLength);
        }
        start = chrono::system_clock::now();
        MAT1_CPU_TIME_COPY += chrono::duration_cast<chrono::microseconds>(start - end).count();

        start = chrono::system_clock::now();
        // write series in order to node file from node fbl
        for(auto & iter:fbl){
            string outfile = index_dir ;
            int id= iter.first->id;
            if(iter.first->partition_id == -1)
                outfile += "U_" + to_string(iter.first->id);
            else
                outfile += to_string(iter.first->layer) + "_" + to_string(iter.first->partition_id);
            FILE *outf = fopen(outfile.c_str(), "a");

            RAND_WRITE_CNT++;
            SEQ_WRITE_CNT += iter.second.size;
//            long bytes = iter.second.size * Const::tsLengthBytes;
//            if(bytes >= Const::small_file_threshold)  SMALL_FILES_BYTES_WRITE += bytes;

            fwrite(iter.second.buffer, sizeof(float), iter.second.size * Const::tsLength, outf);
            fclose(outf);
            delete[]iter.second.buffer;
        }
        end = chrono::system_clock::now();
        MAT1_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();
        delete[] tss;

        rest-=num;
        cur += num;
        Const::logPrint("Now in 1st layer " + to_string((double)cur / (double)total * 100) + "% series have been written to disk.");

    }

    fclose(f);
    delete[] navids;
    auto end_t = chrono::system_clock::now();
    MAT1_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(end_t - start_t).count();
}

// put actual series into disk file of nodes in 1st layer
void materialize1stLayerWithSax(string datafn, FADASNode* root, int *navids, string index_dir, unsigned short*sax_tbl){
    auto start_sax = chrono::system_clock::now();
    Const::logPrint("Start move sax to disk file in 1st layer.");
    unordered_set<FADASNode*>visited;
    for(auto &child:root->children){
        if(child == nullptr || child->size > Const::th || visited.contains(child))  continue;
        visited.insert(child);
        long symbol_num = (long)child->size * Const::segmentNum;
        auto to_write = new unsigned short[symbol_num];

        for(int i=0;i< child->offsets.size();++i){
            long offset =child->offsets[i];
            copy(sax_tbl + offset * Const::segmentNum, sax_tbl + (offset + 1) * Const::segmentNum, to_write + i * Const::segmentNum);
        }
        string outfile = index_dir + "1_";
        if(child->partition_id == -1)
            outfile +=  to_string(child->id) + "_sax_L";
        else
            outfile += to_string(child->partition_id) + "_sax";
        FILE *outf = fopen(outfile.c_str(), "wb");
        fwrite(to_write, sizeof(unsigned short), symbol_num,outf);
        fclose(outf);
        delete[] to_write;
        vector<int>().swap(child->offsets);
    }
    unordered_set<FADASNode*>().swap(visited);
    auto end_sax = chrono::system_clock::now();
    SAX_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(end_sax - start_sax).count();

    auto start_t = chrono::system_clock::now();
    Const::logPrint("Start move data to disk file in 1st layer.");
    FILE *f = fopen(datafn.c_str(), "r");
    long rest = root->size, total = root->size, cur = 0;
    unordered_map<FADASNode*, FBL_UNIT>fbl;

    RAND_READ_CNT++;
    SEQ_READ_CNT += rest;

    // There is another implementation method that fbl in each node stores a pointer vector where each pointer points to a series.
    // This may incur many write calls.
    while(rest > 0){
        fbl.clear();
        long num;
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto end = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto start = chrono::system_clock::now();
        MAT1_READ_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

        // statistic the size of each node fbl size, and allocate memory for them
        for(long i=cur;i<cur+num;++i)
            fbl[root->children[navids[i]]].size++;
        for(auto &iter:fbl)
            iter.second.buffer = new float [(long)iter.second.size * Const::tsLength];
        end = chrono::system_clock::now();
        MAT1_CPU_TIME_STAT += chrono::duration_cast<chrono::microseconds>(end - start).count();

        //copy series to node to ensure write is only called once upon a node fbl
        for(long i = cur; i<cur+num;++i){
            FBL_UNIT* fbl_node = &fbl[root->children[navids[i]]];
            copy(tss + (i-cur) * Const::tsLength, tss + (i+1-cur)* Const::tsLength, fbl_node->buffer + (long)fbl_node->pos++ * Const::tsLength);
        }
        start = chrono::system_clock::now();
        MAT1_CPU_TIME_COPY += chrono::duration_cast<chrono::microseconds>(start - end).count();

        start = chrono::system_clock::now();
        // write series in order to node file from node fbl
        for(auto & iter:fbl){
            string outfile = index_dir ;
            int id= iter.first->id;
            if(iter.first->size > Const::th)
                outfile += "U_" + to_string(iter.first->id);
            else if(iter.first->partition_id != -1)
                outfile += to_string(iter.first->layer) + "_" + to_string(iter.first->partition_id);
            else
                outfile += to_string(iter.first->layer) + "_" + to_string(iter.first->id) + "_L";
            FILE *outf = fopen(outfile.c_str(), "a");

            RAND_WRITE_CNT++;
            SEQ_WRITE_CNT += iter.second.size;
//            long bytes = iter.second.size * Const::tsLengthBytes;
//            if(bytes >= Const::small_file_threshold)  SMALL_FILES_BYTES_WRITE += bytes;

            fwrite(iter.second.buffer, sizeof(float), iter.second.size * Const::tsLength, outf);
            fclose(outf);
            delete[]iter.second.buffer;
        }
        end = chrono::system_clock::now();
        MAT1_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();
        delete[] tss;

        rest-=num;
        cur += num;
        Const::logPrint("Now in 1st layer " + to_string((double)cur / (double)total * 100) + "% series have been written to disk.");

    }

    fclose(f);
    delete[] navids;
    auto end_t = chrono::system_clock::now();
    MAT1_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(end_t - start_t).count();
}

void materialize1stLayerWithSaxOnlyLeaf(string datafn, FADASNode* root, int *navids, string index_dir, unsigned short*sax_tbl){
    auto start_sax = chrono::system_clock::now();
    Const::logPrint("Start move sax to disk file in 1st layer.");
    unordered_set<FADASNode*>visited;
    for(auto &child:root->children){
        if(child == nullptr || child->size > Const::th || visited.contains(child))  continue;
        visited.insert(child);
        long symbol_num = (long)child->size * Const::segmentNum;
        auto to_write = new unsigned short[symbol_num];

        for(int i=0;i< child->offsets.size();++i){
            long offset =child->offsets[i];
            copy(sax_tbl + offset * Const::segmentNum, sax_tbl + (offset + 1) * Const::segmentNum, to_write + i * Const::segmentNum);
        }
        string outfile = index_dir + "1_";
        if(child->partition_id == -1)
            outfile +=  to_string(child->id) + "_sax_L";
        else
            outfile += to_string(child->partition_id) + "_sax";
        FILE *outf = fopen(outfile.c_str(), "wb");
        fwrite(to_write, sizeof(unsigned short), symbol_num,outf);
        fclose(outf);
        delete[] to_write;
        vector<int>().swap(child->offsets);
    }
    unordered_set<FADASNode*>().swap(visited);
    auto end_sax = chrono::system_clock::now();
    SAX_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(end_sax - start_sax).count();

    auto start_t = chrono::system_clock::now();
    Const::logPrint("Start move data to disk file in 1st layer.");
    FILE *f = fopen(datafn.c_str(), "r");
    long rest = root->size, total = root->size, cur = 0;
    unordered_map<FADASNode*, FBL_UNIT>fbl;

    RAND_READ_CNT++;
    SEQ_READ_CNT += rest;

    // There is another implementation method that fbl in each node stores a pointer vector where each pointer points to a series.
    // This may incur many write calls.
    while(rest > 0){
        fbl.clear();
        long num;
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto end = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto start = chrono::system_clock::now();
        MAT1_READ_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

        // statistic the size of each node fbl size, and allocate memory for them
        for(long i=cur;i<cur+num;++i)
            if(root->children[navids[i]]->size <= Const::th)
                fbl[root->children[navids[i]]].size++;
        for(auto &iter:fbl)
            iter.second.buffer = new float [(long)iter.second.size * Const::tsLength];
        end = chrono::system_clock::now();
        MAT1_CPU_TIME_STAT += chrono::duration_cast<chrono::microseconds>(end - start).count();

        //copy series to node to ensure write is only called once upon a node fbl
        for(long i = cur; i<cur+num;++i){
            if(root->children[navids[i]]->size > Const::th) continue;
            FBL_UNIT* fbl_node = &fbl[root->children[navids[i]]];
            copy(tss + (i-cur) * Const::tsLength, tss + (i+1-cur)* Const::tsLength, fbl_node->buffer + (long)fbl_node->pos++ * Const::tsLength);
        }
        start = chrono::system_clock::now();
        MAT1_CPU_TIME_COPY += chrono::duration_cast<chrono::microseconds>(start - end).count();

        start = chrono::system_clock::now();
        // write series in order to node file from node fbl
        for(auto & iter:fbl){
            string outfile = index_dir ;
            int id= iter.first->id;
            assert(iter.first->size <= Const::th);
            if(iter.first->partition_id != -1)
                outfile += to_string(iter.first->layer) + "_" + to_string(iter.first->partition_id);
            else
                outfile += to_string(iter.first->layer) + "_" + to_string(iter.first->id) + "_L";
            FILE *outf = fopen(outfile.c_str(), "a");

            RAND_WRITE_CNT++;
            SEQ_WRITE_CNT += iter.second.size;
//            long bytes = iter.second.size * Const::tsLengthBytes;
//            if(bytes >= Const::small_file_threshold)  SMALL_FILES_BYTES_WRITE += bytes;

            fwrite(iter.second.buffer, sizeof(float), iter.second.size * Const::tsLength, outf);
            fclose(outf);
            delete[]iter.second.buffer;
        }
        end = chrono::system_clock::now();
        MAT1_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();
        delete[] tss;

        rest-=num;
        cur += num;
        Const::logPrint("Now Materialize leaves in the 1st layer. Progress: " + to_string((double)cur / (double)total * 100) + "%");

    }

    fclose(f);
//    delete[] navids;
    auto end_t = chrono::system_clock::now();
    MAT1_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(end_t - start_t).count();
}

void materializeAllLeavesWithSax(string datafn, FADASNode* root, int *navids, string index_dir, unsigned short*sax_tbl){
    auto start_sax = chrono::system_clock::now();
    Const::logPrint("Start move sax to disk file in 1st layer.");

    unordered_map<FADASNode*, vector<unsigned short *>>sax_buffer;
    for(int i=0;i<root->size;++i){
        if(root->children[navids[i]]->size <= Const::th)    continue;
        auto * sax = sax_tbl + i * Const::segmentNum;
        FADASNode* node = root->route(sax);
        sax_buffer[node].push_back(sax);
    }
    for(auto &[node, buffer]:sax_buffer){
        string outfile = Const::fidxfn + node->getFileName() + "_sax";
        if(node->partition_id == -1)    outfile += "_L";
        FILE *outf = fopen(outfile.c_str(), "a");
        for(auto sax:buffer)
            fwrite(sax, sizeof(unsigned short ), Const::segmentNum, outf);
        fclose(outf);
    }
    sax_buffer.clear();
    auto end_sax = chrono::system_clock::now();
    SAX_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(end_sax - start_sax).count();

    auto start_t = chrono::system_clock::now();
    Const::logPrint("Start move data to disk file in 1st layer.");
    FILE *f = fopen(datafn.c_str(), "r");
    long rest = root->size, total = root->size, cur = 0;
    unordered_map<FADASNode*, LBL_UNIT>lbl;

    // There is another implementation method that fbl in each node stores a pointer vector where each pointer points to a series.
    // This may incur many write calls.
    while(rest > 0){
        lbl.clear();
        long num;
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto end = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto start = chrono::system_clock::now();
        MAT2_READ_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

        //copy series to node to ensure write is only called once upon a node fbl
        for(long i = cur; i<cur+num;++i){
            if(root->children[navids[i]]->size <= Const::th) continue;
            FADASNode* node = root->route(sax_tbl + i * Const::segmentNum);
            lbl[node].buffer.push_back(tss + (i-cur) * Const::tsLength);
        }
        end = chrono::system_clock::now();
        MAT2_CPU_TIME += chrono::duration_cast<chrono::microseconds>(end-start).count();

        // write series in order to node file from node fbl
        for(auto & [node,lbl_unit]:lbl){
            string outfile = Const::fidxfn + node->getFileName();
            if(node->partition_id == -1)  outfile += "_L";
            FILE *outf = fopen(outfile.c_str(), "a");

            RAND_WRITE_CNT++;
            for(float *dat:lbl_unit.buffer)
                fwrite(dat, sizeof(float), Const::tsLength, outf);
            fclose(outf);
        }
        start = chrono::system_clock::now();
        MAT2_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();
        delete[] tss;

        rest-=num;
        cur += num;
        Const::logPrint("Now Materialize all leaves below the 1st layer. Progress: " + to_string((double)cur / (double)total * 100) + "%");

    }

    fclose(f);
    delete[] navids;
    auto end_t = chrono::system_clock::now();
    MAT2_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(end_t - start_t).count();
}

// put actual series into disk file of nodes in 1st layer
void materialize1stLayerLessPack(string datafn, FADASNode* root, int *navids, string index_dir){
    auto start_t = chrono::system_clock::now();
    Const::logPrint("Start move data to disk file in 1st layer.");
    FILE *f = fopen(datafn.c_str(), "r");
    long rest = root->size, total = root->size, cur = 0;
    unordered_map<FADASNode*, FBL_UNIT>fbl;

    RAND_READ_CNT++;
    SEQ_READ_CNT += rest;

    // There is another implementation method that fbl in each node stores a pointer vector where each pointer points to a series.
    // This may incur many write calls.
    while(rest > 0){
        fbl.clear();
        long num;
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto end = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto start = chrono::system_clock::now();
        MAT1_READ_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

        // statistic the size of each node fbl size, and allocate memory for them
        for(long i=cur;i<cur+num;++i)
            fbl[root->children[navids[i]]].size++;
        for(auto &iter:fbl)
            iter.second.buffer = new float [(long)iter.second.size * Const::tsLength];
        end = chrono::system_clock::now();
        MAT1_CPU_TIME_STAT += chrono::duration_cast<chrono::microseconds>(end - start).count();

        //copy series to node to ensure write is only called once upon a node fbl
        for(long i = cur; i<cur+num;++i){
            FBL_UNIT* fbl_node = &fbl[root->children[navids[i]]];
            copy(tss + (i-cur) * Const::tsLength, tss + (i+1-cur)* Const::tsLength, fbl_node->buffer + (long)fbl_node->pos++ * Const::tsLength);
        }
        start = chrono::system_clock::now();
        MAT1_CPU_TIME_COPY += chrono::duration_cast<chrono::microseconds>(start - end).count();

        start = chrono::system_clock::now();
        // write series in order to node file from node fbl
        for(auto & iter:fbl){
            string outfile = index_dir ;
            int id= iter.first->id;
            if(iter.first->isInternalNode())
                outfile += "U_" + to_string(iter.first->id);
            else if(iter.first->isLeafPack())
                outfile += to_string(iter.first->layer) + "_P_" + to_string(iter.first->partition_id);
            else if(iter.first->isLeafNode())
                outfile += to_string(iter.first->layer) + "_" + to_string(iter.first->id);
            else
                Const::logPrint(" Some error in materialzing.");
            FILE *outf = fopen(outfile.c_str(), "a");

            RAND_WRITE_CNT++;
            SEQ_WRITE_CNT += iter.second.size;
//            long bytes = iter.second.size * Const::tsLengthBytes;
//            if(bytes >= Const::small_file_threshold)  SMALL_FILES_BYTES_WRITE += bytes;

            fwrite(iter.second.buffer, sizeof(float), iter.second.size * Const::tsLength, outf);
            fclose(outf);
            delete[]iter.second.buffer;
        }
        end = chrono::system_clock::now();
        MAT1_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();
        delete[] tss;

        rest-=num;
        cur += num;
        Const::logPrint("Now in 1st layer " + to_string((double)cur / (double)total * 100) + "% series have been written to disk.");

    }

    fclose(f);
    delete[] navids;
    auto end_t = chrono::system_clock::now();
    MAT1_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(end_t - start_t).count();
}

// put actual series into disk file of nodes in 1st layer
void materialize1stLayerPos(string datafn, FADASNode* root, int *navids, string index_dir){
    auto start_t = chrono::system_clock::now();
    Const::logPrint("Start move data to disk file in 1st layer.");
    FILE *f = fopen(datafn.c_str(), "r");
    int rest = root->size, total = root->size, cur = 0;
    unordered_map<FADASNode*, LBL_UNIT>fbl;

    RAND_READ_CNT++;
    SEQ_READ_CNT += rest;

    // There is another implementation method that fbl in each node stores a pointer vector where each pointer points to a series.
    // This may incur many write calls.
    while(rest > 0){
        fbl.clear();
        int num;
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto end = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto start = chrono::system_clock::now();
        MAT1_READ_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

        // statistic the size of each node fbl size, and allocate memory for them
//        for(int i=cur;i<cur+num;++i)
//            fbl[root->children[navids[i]]].size++;
//        for(auto &iter:fbl)
//            iter.second.buffer = new float [iter.second.size * Const::tsLength];
        for(int i=cur;i<cur+num;++i) {
            fbl[root->children[navids[i]]].offsets.push_back(i);
        }
        end = chrono::system_clock::now();
        MAT1_CPU_TIME_STAT += chrono::duration_cast<chrono::microseconds>(end - start).count();

        //copy series to node to ensure write is only called once upon a node fbl
//        for(int i = cur; i<cur+num;++i){
//            FBL_UNIT* fbl_node = &fbl[root->children[navids[i]]];
//            copy(tss + (i-cur) * Const::tsLength, tss + (i+1-cur)* Const::tsLength, fbl_node->buffer + fbl_node->pos++ * Const::tsLength);
//        }
//        start = chrono::system_clock::now();
//        MAT1_CPU_TIME_COPY += chrono::duration_cast<chrono::microseconds>(start - end).count();

        start = chrono::system_clock::now();
        // write series in order to node file from node fbl
        for(auto & iter:fbl){
            string outfile = index_dir ;
            if(iter.first->partition_id == -1)
                outfile += "U_" + to_string(iter.first->id);
            else outfile += to_string(iter.first->layer) + "_" + to_string(iter.first->partition_id);
            FILE *outf = fopen(outfile.c_str(), "a");
            RAND_WRITE_CNT++;

            for(int offset:iter.second.offsets){
                fwrite(tss + (offset - cur) * Const::tsLength, sizeof(float), Const::tsLength, outf);
                SEQ_WRITE_CNT++;
                if(iter.first->partition_id != -1) {
                    fwrite(&offset, sizeof(int), 1, outf);
                    SEQ_WRITE_CNT++;
                }
            }
//            fwrite(iter.second.buffer, sizeof(float), iter.second.size * Const::tsLength, outf);
            fclose(outf);
//            delete[]iter.second.buffer;
        }
        end = chrono::system_clock::now();
        MAT1_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();
        delete[] tss;

        rest-=num;
        cur += num;
        Const::logPrint("Now in 1st layer " + to_string((double)cur / (double)total * 100) + "% series have been written to disk.");

    }

    fclose(f);
    delete[] navids;
    auto end_t = chrono::system_clock::now();
    MAT1_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(end_t - start_t).count();
}

// put actual series into disk file of nodes below 1st layer from file in 1st layer
void materializeInterNode(FADASNode *node, unsigned short *saxes) {
    auto start_t = chrono::system_clock::now();

    FILE *f = fopen((Const::fidxfn + "U_" + to_string(node->id)).c_str(), "r");
    long rest = node->size, cur = 0, num;
    unordered_map<FADASNode*, LBL_UNIT>lbl;

//    long bytes = rest * Const::tsLengthBytes;
//    if(bytes >= Const::small_file_threshold)    SMALL_FILES_BYTES_READ += bytes;

    RAND_READ_CNT++;
    SEQ_READ_CNT+=rest;

    while(rest > 0){
        lbl.clear();
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto end = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto start = chrono::system_clock::now();
        MAT2_READ_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

        for(long i = cur; i < cur + num; ++i){
            FADASNode* target = node->route(saxes + (long)node->offsets[i] * Const::segmentNum);
            lbl[target].buffer.push_back(tss + (i - cur) * Const::tsLength);
        }

        end = chrono::system_clock::now();
        MAT2_CPU_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();


        for(auto &iter:lbl){
            string outfile = Const::fidxfn + iter.first->getFileName();
            FILE *outf = fopen(outfile.c_str(), "a");
//            SEQ_WRITE_CNT += iter.second.buffer.size();
//            bytes = iter.second.buffer.size() * Const::tsLengthBytes;
//            if(bytes >= Const::small_file_threshold)    SMALL_FILES_BYTES_WRITE += bytes;
            RAND_WRITE_CNT++;
            for(float *dat:iter.second.buffer)
                fwrite(dat, sizeof(float), Const::tsLength, outf);
            fclose(outf);
        }
        start = chrono::system_clock::now();
        MAT2_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

        delete[]tss;
        rest-=num;
        cur += num;
    }
    fclose(f);
    FileUtil::FileRemove((Const::fidxfn + "U_" + to_string(node->id)).c_str());
    auto end_t = chrono::system_clock::now();
    MAT2_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(end_t - start_t).count();
}

// put actual series into disk file of nodes below 1st layer from file in 1st layer
void materializeInterNodeWithSax(FADASNode *node, unsigned short *saxes) {
    auto start_sax = chrono::system_clock::now();
    // route each sax words to the buffer of a leaf pack
    unordered_map<FADASNode*, SAX_BUF_UNIT>sax_buffer;
    for(int i=0;i<node->size;++i){
        int offset = node->offsets[i];
        unsigned short *sax = saxes + offset * Const::segmentNum;
        FADASNode*target = node->route(sax);
        if(sax_buffer[target].size == 0){
            sax_buffer[target].buffer = new unsigned short [(long)target->size*Const::segmentNum];
        }
        copy(sax, sax+Const::segmentNum, sax_buffer[target].buffer + (long)sax_buffer[target].size * Const::segmentNum);
        sax_buffer[target].size++;
    }

    // write each buffer of the hash table
    for(auto &iter: sax_buffer){
        string outfile = Const::fidxfn + iter.first->getFileName() + "_sax";
        if(iter.first->partition_id == -1)  outfile += "_L";
        FILE *outf = fopen(outfile.c_str(), "a");
        fwrite(iter.second.buffer, sizeof(unsigned short ), (long)iter.first->size * Const::segmentNum, outf);
        fclose(outf);
        delete[] iter.second.buffer;
    }
    unordered_map<FADASNode*, SAX_BUF_UNIT>().swap(sax_buffer);
    auto end_sax = chrono::system_clock::now();
    SAX_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(end_sax - start_sax).count();


    auto start_t = chrono::system_clock::now();
    FILE *f = fopen((Const::fidxfn + "U_" + to_string(node->id)).c_str(), "r");
    long rest = node->size, cur = 0, num;
    unordered_map<FADASNode*, LBL_UNIT>lbl;

//    long bytes = rest * Const::tsLengthBytes;
//    if(bytes >= Const::small_file_threshold)    SMALL_FILES_BYTES_READ += bytes;

    RAND_READ_CNT++;
    SEQ_READ_CNT+=rest;

    while(rest > 0){
        lbl.clear();
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto end = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto start = chrono::system_clock::now();
        MAT2_READ_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

        for(long i = cur; i < cur + num; ++i){
            FADASNode* target = node->route(saxes + (long)node->offsets[i] * Const::segmentNum);
            lbl[target].buffer.push_back(tss + (i - cur) * Const::tsLength);
        }

        end = chrono::system_clock::now();
        MAT2_CPU_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();


        for(auto &iter:lbl){
            string outfile = Const::fidxfn + iter.first->getFileName();
            if(iter.first->partition_id == -1)  outfile += "_L";
            FILE *outf = fopen(outfile.c_str(), "a");
//            SEQ_WRITE_CNT += iter.second.buffer.size();
//            bytes = iter.second.buffer.size() * Const::tsLengthBytes;
//            if(bytes >= Const::small_file_threshold)    SMALL_FILES_BYTES_WRITE += bytes;
            RAND_WRITE_CNT++;
            for(float *dat:iter.second.buffer)
                fwrite(dat, sizeof(float), Const::tsLength, outf);
            fclose(outf);
        }
        start = chrono::system_clock::now();
        MAT2_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

        delete[]tss;
        rest-=num;
        cur += num;
    }
    fclose(f);
    FileUtil::FileRemove((Const::fidxfn + "U_" + to_string(node->id)).c_str());
    auto end_t = chrono::system_clock::now();
    MAT2_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(end_t - start_t).count();
}

// put actual series into disk file of nodes below 1st layer from raw dataset
void materializeOnePassWithSax(FADASNode *node, unsigned short *saxes, float *series) {
    // route each sax words to the buffer of a leaf pack
    unordered_map<FADASNode*, SAX_BUF_UNIT>sax_buffer;
    unordered_map<FADASNode*, LBL_UNIT>lbl;
    for(long i=0;i<node->size;++i){
        unsigned short *sax = saxes + (long)node->offsets[i] * Const::segmentNum;
        FADASNode*target = node->route(sax);
        if(sax_buffer[target].size == 0){
            sax_buffer[target].buffer = new unsigned short [(long)target->size*Const::segmentNum];
        }
        lbl[target].buffer.push_back(series + (long)node->offsets[i] * Const::tsLength);
        copy(sax, sax+Const::segmentNum, sax_buffer[target].buffer + (long)sax_buffer[target].size * Const::segmentNum);
        sax_buffer[target].size++;
    }

    // write each buffer of the hash table
    for(auto &iter: sax_buffer){
        string outfile = Const::fidxfn + iter.first->getFileName() + "_sax";
        if(iter.first->partition_id == -1)  outfile += "_L";
        FILE *outf = fopen(outfile.c_str(), "a");
        fwrite(iter.second.buffer, sizeof(unsigned short ), (long)iter.first->size * Const::segmentNum, outf);
        fclose(outf);
        delete[] iter.second.buffer;
    }
    unordered_map<FADASNode*, SAX_BUF_UNIT>().swap(sax_buffer);

    // write series
    for(auto &iter:lbl){
        string outfile = Const::fidxfn + iter.first->getFileName();
        if(iter.first->partition_id == -1)  outfile += "_L";
        FILE *outf = fopen(outfile.c_str(), "a");
        RAND_WRITE_CNT++;
        for(float *dat:iter.second.buffer)
            fwrite(dat, sizeof(float), Const::tsLength, outf);
        fclose(outf);
    }

    vector<int>().swap(node->offsets);
}

// put actual series into disk file of nodes below 1st layer from raw dataset
void materializeOnePassWithSax(FADASNode *node, unsigned short *saxes, vector<string>&leaf_files, vector<string>&sax_files) {
    // route each sax words to the buffer of a leaf pack
    unordered_map<FADASNode*, SAX_BUF_UNIT>sax_buffer;
    unordered_map<FADASNode*, LBL_UNIT>lbl;
    for(long i=0;i<node->size;++i){
        unsigned short *sax = saxes + (long)node->offsets[i] * Const::segmentNum;
        FADASNode*target = node->route(sax);
        if(sax_buffer[target].size == 0){
            sax_buffer[target].buffer = new unsigned short [(long)target->size*Const::segmentNum];
        }
        copy(sax, sax+Const::segmentNum, sax_buffer[target].buffer + (long)sax_buffer[target].size * Const::segmentNum);
        sax_buffer[target].size++;
    }

    // write each buffer of the hash table
    for(auto &iter: sax_buffer){
        string outfile = Const::fidxfn + iter.first->getFileName() + "_sax";
        if(iter.first->partition_id == -1)  outfile += "_L";
        FILE *outf = fopen(outfile.c_str(), "a");
        fwrite(iter.second.buffer, sizeof(unsigned short ), (long)iter.first->size * Const::segmentNum, outf);
        fclose(outf);
        delete[] iter.second.buffer;
    }
    unordered_map<FADASNode*, SAX_BUF_UNIT>().swap(sax_buffer);

    int cur_num = 0;
    auto buf = new float [(long)Const::fbl_series_num * Const::tsLength];
    auto buf_sax = new unsigned short [(long)Const::fbl_series_num * Const::segmentNum];
    for(int i =0;i<leaf_files.size();++i){
        string &leaf_file = leaf_files[i];
        string &sax_file = sax_files[i];
        int series_num = FileUtil::getFileSize(leaf_file.c_str()) / sizeof(float) / Const::tsLength;
        if(cur_num + series_num > Const::fbl_series_num){
            for(int j=0;j<cur_num;++j){
                FADASNode*target = node->route(buf_sax + (long)j * Const::segmentNum);
                lbl[target].buffer.push_back(buf + (long)j*Const::tsLength);
            }

            // write series
            for(auto &iter:lbl){
                string outfile = Const::fidxfn + iter.first->getFileName();
                if(iter.first->partition_id == -1)  outfile += "_L";
                FILE *outf = fopen(outfile.c_str(), "a");
                RAND_WRITE_CNT++;
                for(float *dat:iter.second.buffer)
                    fwrite(dat, sizeof(float), Const::tsLength, outf);
                fclose(outf);
            }

            cur_num = 0;
            lbl.clear();
        }

        FILE *f = fopen(leaf_file.c_str(), "rb");
        FILE *saxf = fopen(sax_file.c_str(), "rb");
        fread(buf + (long)cur_num * Const::tsLength, sizeof(float ), (long)series_num * Const::tsLength, f);
        fread(buf_sax + (long)cur_num * Const::segmentNum, sizeof(unsigned short ), (long)series_num * Const::segmentNum, saxf);
        cur_num += series_num;
        fclose(f);
        fclose(saxf);
        FileUtil::FileRemove(leaf_file.c_str());
        FileUtil::FileRemove(sax_file.c_str());
    }

    if(cur_num > 0) {
        for (int j = 0; j < cur_num; ++j) {
            FADASNode *target = node->route(buf_sax + (long) j * Const::segmentNum);
            lbl[target].buffer.push_back(buf + (long) j * Const::tsLength);
        }

        // write series
        for (auto &iter: lbl) {
            string outfile = Const::fidxfn + iter.first->getFileName();
            if (iter.first->partition_id == -1) outfile += "_L";
            FILE *outf = fopen(outfile.c_str(), "a");
            RAND_WRITE_CNT++;
            for (float *dat: iter.second.buffer)
                fwrite(dat, sizeof(float), Const::tsLength, outf);
            fclose(outf);
        }

        lbl.clear();
    }

    delete[] buf;
    delete[] buf_sax;

    vector<int>().swap(node->offsets);
}


// put actual series into disk file of nodes below 1st layer from file in 1st layer
void materializeInterNodeLessPack(FADASNode *node, unsigned short *saxes) {
    auto start_t = chrono::system_clock::now();

    FILE *f = fopen((Const::fidxfn + "U_" + to_string(node->id)).c_str(), "r");
    long rest = node->size, cur = 0, num;
    unordered_map<FADASNode*, LBL_UNIT>lbl;

//    long bytes = rest * Const::tsLengthBytes;
//    if(bytes >= Const::small_file_threshold)    SMALL_FILES_BYTES_READ += bytes;

    RAND_READ_CNT++;
    SEQ_READ_CNT+=rest;

    while(rest > 0){
        lbl.clear();
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto end = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto start = chrono::system_clock::now();
        MAT2_READ_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

        for(long i = cur; i < cur + num; ++i){
            FADASNode* target = node->route(saxes + (long)node->offsets[i] * Const::segmentNum);
            lbl[target].buffer.push_back(tss + (i - cur) * Const::tsLength);
        }

        end = chrono::system_clock::now();
        MAT2_CPU_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();


        for(auto &iter:lbl){
            assert(iter.first->layer > 1);
            string outfile = Const::fidxfn;
            if(iter.first->isLeafPack())
                outfile += "P-" + iter.first->getFileName();
            else
                outfile += iter.first->getFileName();
            FILE *outf = fopen(outfile.c_str(), "a");
//            SEQ_WRITE_CNT += iter.second.buffer.size();
//            bytes = iter.second.buffer.size() * Const::tsLengthBytes;
//            if(bytes >= Const::small_file_threshold)    SMALL_FILES_BYTES_WRITE += bytes;
            RAND_WRITE_CNT++;
            for(float *dat:iter.second.buffer)
                fwrite(dat, sizeof(float), Const::tsLength, outf);
            fclose(outf);
        }
        start = chrono::system_clock::now();
        MAT2_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

        delete[]tss;
        rest-=num;
        cur += num;
    }
    fclose(f);
    FileUtil::FileRemove((Const::fidxfn + "U_" + to_string(node->id)).c_str());
    auto end_t = chrono::system_clock::now();
    MAT2_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(end_t - start_t).count();
}

// put actual series into disk file of nodes below 1st layer from file in 1st layer
void materializeInterNodePos(FADASNode *node, unsigned short *saxes) {
    auto start_t = chrono::system_clock::now();

    FILE *f = fopen((Const::posidxfn + "U_" + to_string(node->id)).c_str(), "r");
    int rest = node->size, cur = 0, num;
    unordered_map<FADASNode*, LBL_UNIT>lbl;

    RAND_READ_CNT++;
    SEQ_READ_CNT+=rest;

    while(rest > 0){
        lbl.clear();
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto end = chrono::system_clock::now();
        fread(tss, sizeof(float ), num * Const::tsLength,  f);
        auto start = chrono::system_clock::now();
        MAT2_READ_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

//        for(int i = cur; i < cur + num; ++i){
//            FADASNode* target = node->route(saxes + node->offsets[i] * Const::segmentNum);
//            lbl[target].buffer.push_back(tss + (i - cur) * Const::tsLength);
//        }
        for(int i = cur; i < cur + num; ++i){
            FADASNode* target = node->route(saxes + node->offsets[i] * Const::segmentNum);
            lbl[target].offsets.push_back(node->offsets[i]);
            lbl[target].buffer.push_back(tss + (i - cur) * Const::tsLength);
        }

        end = chrono::system_clock::now();
        MAT2_CPU_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();


        for(auto &iter:lbl){
            string outfile = Const::fidxfn + iter.first->getFileName();
            FILE *outf = fopen(outfile.c_str(), "a");
            LBL_UNIT &tmp = iter.second;
            RAND_WRITE_CNT++;
            for(int i=0;i<tmp.offsets.size();++i){
                fwrite(tmp.buffer[i], sizeof(float), Const::tsLength, outf);
                fwrite(&tmp.offsets[i], sizeof(int), 1, outf);
                SEQ_WRITE_CNT+=2;
            }
//            for(float *dat:iter.second.buffer)
//                fwrite(dat, sizeof(float), Const::tsLength, outf);
            fclose(outf);
        }
        start = chrono::system_clock::now();
        MAT2_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

        delete[]tss;
        rest-=num;
        cur += num;
    }
    fclose(f);
    FileUtil::FileRemove((Const::posidxfn + "U_" + to_string(node->id)).c_str());
    auto end_t = chrono::system_clock::now();
    MAT2_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(end_t - start_t).count();
}

FADASNode *FADASNode::BuildIndex(string &datafn, string &saxfn) {
    Const::logPrint("Start building index.");
    auto start_t = chrono::system_clock::now();
    loadCombines();
    FileUtil::checkDirClean(Const::fidxfn.c_str());
    auto end = chrono::system_clock::now();
//    long series_num = generateSaxTbl();
    long series_num = loadSax(saxfn);
//    loadPaa(paafn);
    auto start = chrono::system_clock::now();
    Const::logPrint("Finish building sax table.");
    SAX_PAA_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

    auto* root = new FADASNode();
    root->size = series_num;
    for(int &i:root->bits_cardinality)  i=0;
    for(int i=0;i<Const::segmentNum;++i)    root->chosenSegments.push_back(i);
    vector<partUnit> nodeIn1stLayer(Const::vertexNum);
    int *navids = new int[series_num];
    for(int i=0;i<Const::vertexNum;++i)
        nodeIn1stLayer[i].id = i, nodeIn1stLayer[i].size=0, nodeIn1stLayer[i].pid = -1;

    // get 1st layer node size
    for(long i=0;i<series_num;++i){
        unsigned short *asax = saxes + i * Const::segmentNum;
        int nav_id = SaxUtil::invSaxHeadFromSax(asax, Const::bitsCardinality, Const::segmentNum);
        navids[i] = nav_id;
        nodeIn1stLayer[nav_id].size++;
    }

    Const::logPrint("Finish statistic size of nodes in the 1st layer.");

    // partition 1st layer
//    int partNum = partition1stLayer(nodeIn1stLayer, g, Const::filling_factor_1st);
    int partNum = partitionNew(nodeIn1stLayer, Const::segmentNum);
    Const::logPrint("Finish partition");
    // build rest data node if any
//    for(auto &node:nodeIn1stLayer)
//        if(node.size <= Const::th && node.pid == -1)
//            node.pid = partNum++;

    FADASNode* childrenList[partNum];
    for(int i=0;i<partNum;++i)  childrenList[i] = new FADASNode(1, i);
    root->children.resize(Const::vertexNum);
    for(int i=0;i<Const::vertexNum;++i){
        if(nodeIn1stLayer[i].size <= 0) continue;
        if(nodeIn1stLayer[i].size > Const::th) {
//            assert(nodeIn1stLayer[i].size > Const::th);
            root->children[i] = new FADASNode(1, nodeIn1stLayer[i].size, i);
            root->children[i]->generateSaxAndCardIn1stLayer(i);
        }else if(nodeIn1stLayer[i].pid == -1){
            root->children[i] = new FADASNode(1, nodeIn1stLayer[i].size, i);
            root->children[i]->generateSaxAndCardIn1stLayer(i);
        }
        else{
            int pid = nodeIn1stLayer[i].pid;
            root->children[i] = childrenList[pid];
            childrenList[pid]->size += nodeIn1stLayer[i].size;
            childrenList[pid]->generateSaxAndCardIn1stLayer4LeafNode(i);
        }
    }
    Const::logPrint("Finish build index structure 1st layer.");

    // put data offsets to internal nodes in 1st layer
    for(int i=0;i<Const::vertexNum;++i)
        if(root->children[i]!= nullptr)
            root->children[i]->offsets.reserve(nodeIn1stLayer[i].size);
    for(int i=0;i<series_num;++i){
        int nav_id = navids[i];
        root->children[nav_id]->offsets.push_back(i);
    }
    Const::logPrint("data offsets have been put into nodes in the 1st layer.");
    end = chrono::system_clock::now();
    GROW_CPU_TIME_1st += chrono::duration_cast<chrono::microseconds>(end - start).count();

    thread IO(materialize1stLayerWithSaxOnlyLeaf, datafn, root, navids, Const::fidxfn, saxes);

    int j = 0;
    int milestone = 0.1 * Const::vertexNum;
    Const::logPrint("start grow the index structure");
    for(int i=0;i<Const::vertexNum;++i){
        if(nodeIn1stLayer[i].size > Const::th) {
            root->children[i]->growIndex();
        }
        if(++j%milestone == 0)
            Const::logPrint(to_string(j) + " nodes in the 1st layer has been processed.");
    }
    start = chrono::system_clock::now();
    GROW_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

    Const::logPrint("build index skeleton finished.");

    IO.join();
    Const::logPrint("Start materialize internal nodes in the 1st layer");
    materializeAllLeavesWithSax(datafn, root, navids, Const::fidxfn, saxes);
//    for(int i=0;i<Const::vertexNum;++i)
//        if(nodeIn1stLayer[i].size > Const::th)
//            materializeInterNodeWithSax(root->children[i], saxes);
    Const::logPrint("build index successfully!");
    delete[] saxes;
    auto end_t = chrono::system_clock::now();
    cout << "Total time is " << chrono::duration_cast<chrono::microseconds>(end_t - start_t).count() / 1000 << "ms."<<endl;
    cout << "Building sax and paa total time is "<<SAX_PAA_TOTAL_TIME / 1000 <<"ms, cpu time is "
    << SAX_PAA_CPU_TIME / 1000 <<"ms, I/O read time is " << SAX_PAA_READ_TIME / 1000 << "ms."<<endl;

    cout << "During the process of building 1st layer index structure, CPU time is "<< GROW_CPU_TIME_1st / 1000 <<"ms, "<<endl;

    cout << "During the process of materializing 1st layer nodes, total time is "<< MAT1_TOTAL_TIME / 1000
         <<"ms, I/O read time is "<< MAT1_READ_TIME / 1000 <<"ms, CPU statistic time  is " << MAT1_CPU_TIME_STAT / 1000
         << "ms, CPU copy Time is " << MAT1_CPU_TIME_COPY / 1000 << "ms, while I/O write time  is " << MAT1_WRITE_TIME / 1000 << "ms. "<< endl;

    cout << "During the process of growing index structure, total time is "<< GROW_TOTAL_TIME / 1000
        <<"ms, other CPU time is "<< GROW_CPU_TIME / 1000 <<"ms."<<endl;

    cout << "During the process of materializing internal nodes, total time is "<<MAT2_TOTAL_TIME / 1000
         <<"ms, I/O read time is "<<MAT2_READ_TIME / 1000 <<"ms, CPU time is " << MAT2_CPU_TIME / 1000
         << "ms, while I/O write time is " << MAT2_WRITE_TIME / 1000 << "ms." << endl;

    cout << "write sax needs " << SAX_WRITE_TIME / 1000 << "ms" << endl;

    cout << "Random read count = " << RAND_READ_CNT << endl
        << "Random write count = " << RAND_WRITE_CNT << endl
        <<"Sequence read count = " << SEQ_READ_CNT << endl
        <<"Sequence write count = " << SEQ_WRITE_CNT << endl
        <<"Small file read bytes = "<<SMALL_FILES_BYTES_READ << endl
        <<"Small file write bytes = " << SMALL_FILES_BYTES_WRITE << endl;

    double output_time = (MAT1_WRITE_TIME + MAT2_WRITE_TIME + SAX_WRITE_TIME) / (1000.0 * 1000.0 * 60);
    double input_time = (SAX_PAA_READ_TIME + MAT1_READ_TIME + MAT2_READ_TIME) /(1000.0 * 1000.0 * 60);
    double cpu_time = (SAX_PAA_CPU_TIME + GROW_CPU_TIME_1st + MAT1_CPU_TIME_STAT + MAT1_CPU_TIME_COPY + GROW_TOTAL_TIME + MAT2_CPU_TIME) / (1000.0 * 1000.0 * 60);
    cout << "Output needs" << output_time << " mins" <<endl
    <<"Input needs" << input_time << " mins"<<endl
    <<"CPU needs" << cpu_time << "mins"<<endl;

    root->getIndexStats();

    return root;
}

void FADASNode::insertBatch(float *tss, int batch_size){
    generateSaxTbl(tss, batch_size);
    for(int i=0;i<batch_size;++i){
        int head = SaxUtil::invSaxHeadFromSax(saxes + i * Const::segmentNum, Const::bitsCardinality, Const::segmentNum);
        if(children[head] == nullptr){
            children[head] = new FADASNode(1,1,head);
            children[head]->generateSaxAndCardIn1stLayer(head);
            children[head]->pos_cache.push_back(i);
        }else
            children[head]->routeDuringInsertion(saxes + i * Const::segmentNum, i);
    }

//    for(FADASNode* child:children){
//        if(child!= nullptr && child->partition_id == 63)
//            cout <<"here" <<endl;
//    }
    reorganize(tss, nullptr);

    delete[] saxes;
}

void FADASNode::reorganize(float * tss, FADASNode* parent){
    if(pos_cache.empty()){
        if(!isInternalNode())   return;
        unordered_set<FADASNode*>visited;
        for(FADASNode* child: children){
            if(child == nullptr)    continue;
            if(visited.contains(child)) continue;
            child->reorganize(tss, this);
            visited.insert(child);
        }
        return;
    }

    if(size <= Const::th){
        // leaf pack or leaf node
        string sax_file, data_file;
        getFileNameInsert(Const::fidxfn, sax_file, data_file);
        FILE* sax_f = fopen(sax_file.c_str(), "ab");
        for(int pos: pos_cache){
            fwrite(saxes + pos * Const::segmentNum, sizeof(unsigned short), Const::segmentNum, sax_f);
        }
        fclose(sax_f);
        FILE* data_f = fopen(data_file.c_str(), "ab");
        for(int pos: pos_cache){
            fwrite(tss + pos * Const::tsLength, sizeof(float ), Const::tsLength, data_f);
        }
        fclose(data_f);
        vector<int>().swap(pos_cache);
        return;
    }
    // we need to do split (grow index) when it exceeds th
    if(children.empty()){
        // leaf node or pack
        string sax_file, data_file;
        getFileNameInsert(Const::fidxfn, sax_file, data_file);
        auto *node_saxes = new unsigned short[size * Const::segmentNum];
        auto *node_series = new float [size * Const::tsLength];
        for(int i=0;i<pos_cache.size();++i){
            copy(saxes + (long)pos_cache[i] * Const::segmentNum, saxes + ((long)pos_cache[i]+1) * Const::segmentNum, node_saxes + (long)i * Const::segmentNum);
            copy(tss+ (long)pos_cache[i] * Const::tsLength, tss + ((long)pos_cache[i]+1) * Const::tsLength, node_series + (long)i * Const::tsLength);
        }
        long on_disk_nbr = FileUtil::getFileSize(sax_file.c_str()) / sizeof(unsigned short ) / Const::segmentNum;
        FILE* sax_f = fopen(sax_file.c_str(), "rb");
        FILE* data_f = fopen(data_file.c_str(), "rb");
        for(long i = pos_cache.size(); i< pos_cache.size() + on_disk_nbr;++i){
            fread(node_saxes + i * Const::segmentNum, sizeof(unsigned short ), Const::segmentNum, sax_f);
            fread(node_series + i * Const::tsLength, sizeof(float ), Const::tsLength, data_f);
        }
        fclose(sax_f);
        fclose(data_f);
//        FileUtil::renameFile(sax_file, sax_file + "tmp");
//        FileUtil::renameFile(data_file, data_file + "tmp");

        if(partition_id == -1){ // leaf node
            offsets.resize(size);
            for(int i=0;i<size;++i) offsets[i] = i;
            growIndex(node_saxes, false);
            materializeOnePassWithSax(this, node_saxes, node_series);
            vector<int>().swap(pos_cache);
        }else{//leaf pack
            int chosen_num;
            if(layer == 1)
                chosen_num = Const::segmentNum;
            else
                chosen_num = parent->chosenSegments.size();
            // statistic children information in order to partition
            partUnit nodes[1<<chosen_num];
            for(int i=0;i<(1<<chosen_num);++i)
                nodes[i].id = i, nodes[i].size=0, nodes[i].pid = -1;
            vector<vector<int>>node_offsets(1<<chosen_num, vector<int>());

            for(int i=0;i<size;++i){
                int new_id;
                if(layer == 1)
                    new_id = SaxUtil::invSaxHeadFromSax(node_saxes + (long)i * (Const::segmentNum), Const::bitsCardinality, Const::segmentNum);
                else
                    new_id = SaxUtil::extendSax(node_saxes + (long)i * (Const::segmentNum), parent->bits_cardinality,parent->chosenSegments);
                nodes[new_id].size++;
                node_offsets[new_id].push_back(i);
            }

//            FADASNode * pack_node;
//            if(layer == 1){
//                pack_node = new FADASNode(1, partition_id);
//            }else{
//                pack_node = new FADASNode(parent, partition_id);
//            }

            for(int i=0;i<(1<<chosen_num);++i){
                if(nodes[i].size <= 0)  continue;
                auto new_node = new FADASNode(parent, nodes[i].size, i);

                if(nodes[i].size > Const::th){
                    new_node->offsets.resize(nodes[i].size);
                    copy(node_offsets[i].begin(), node_offsets[i].end(), new_node->offsets.begin());
                    new_node->growIndex(node_saxes, false);
                    materializeOnePassWithSax(new_node, node_saxes, node_series);
                    parent->children[i] = new_node;
                    vector<int>().swap(node_offsets[i]);
                }
//                else if(nodes[i].size >= Const::th * Const::small_perc)
                else{
                    // directly write series into a leaf node
                    string out_sax_file, out_data_file;
                    new_node->getFileNameInsert(Const::fidxfn, out_sax_file, out_data_file);
                    FILE* out_sax_f = fopen(out_sax_file.c_str(), "wb");
                    FILE* out_data_f = fopen(out_data_file.c_str(), "wb");
                    for(int offset: node_offsets[i]){
                        fwrite(node_saxes + (long)offset * Const::segmentNum, sizeof(unsigned short ), Const::segmentNum, out_sax_f);
                        fwrite(node_series + (long)offset * Const::tsLength, sizeof(float ), Const::tsLength, out_data_f);
                    }
                    vector<int>().swap(node_offsets[i]);
                    fclose(out_data_f);
                    fclose(out_sax_f);
                    parent->children[i] = new_node;
                }

//                else{
//                    delete new_node;
//                    int _pid = nodes[i].pid;
//                    pack_node->size += nodes[i].size;
//                    if(layer == 1)
//                        pack_node->generateSaxAndCardIn1stLayer4LeafNode(i);
//                    else
//                        parent->generateSaxAndCardinality4LeafNode(pack_node, i);
//                    for(int & j : node_offsets[i])
//                        pack_node->offsets.push_back(j);
//                    vector<int>().swap(node_offsets[i]);
//                    parent->children[i] = pack_node;
//                }
            }

//            string out_sax_file, out_data_file;
//            pack_node->getFileNameInsert(Const::fidxfn, out_sax_file, out_data_file);
//            FILE* out_sax_f = fopen(out_sax_file.c_str(), "wb");
//            FILE* out_data_f = fopen(out_data_file.c_str(), "wb");
//            for(int offset: pack_node->offsets){
//                fwrite(node_saxes + (long)offset * Const::segmentNum, sizeof(unsigned short ), Const::segmentNum, out_sax_f);
//                fwrite(node_series + (long)offset * Const::tsLength, sizeof(float ), Const::tsLength, out_data_f);
//            }
//            vector<int>().swap(pack_node->offsets);
//            fclose(out_data_f);
//            fclose(out_sax_f);
            vector<vector<int>>().swap(node_offsets);
            delete this;

        }

        delete[] node_saxes;
        delete[] node_series;
        FileUtil::FileRemove((sax_file).c_str());
        FileUtil::FileRemove((data_file).c_str());

    }else{
        // internal node
        auto node_saxes = new unsigned short [size * Const::segmentNum];
        vector<string>leaf_files, sax_files;
        int cur = 0;
        collectSAXwords(node_saxes, &cur, leaf_files, sax_files);
        assert(cur == size);
        offsets.resize(size);
        for(int i=0;i<size;++i) offsets[i] = i;
        for(FADASNode* child: children){
            if(child != nullptr)
                child->deleteSubtree();
        }
        chosenSegments.clear();
        children.clear();
        growIndex(node_saxes, false);
        materializeOnePassWithSax(this, node_saxes, leaf_files, sax_files);
        vector<string>().swap(leaf_files);
        vector<string>().swap(sax_files);
        vector<int>().swap(pos_cache);
        delete[] node_saxes;
    }
}

void
FADASNode::collectSAXwords(unsigned short *node_saxes, int *cur, vector<string> &leaf_files, vector<string> &sax_files) {
    if(!pos_cache.empty()){
        for(int offset:pos_cache){
            copy(saxes + (long)offset * Const::segmentNum, saxes + ((long)offset + 1) * Const::segmentNum, node_saxes + (long)(*cur) * Const::segmentNum);
            (*cur) = (*cur) + 1;
        }
    }
    string sax_file, data_file;
    getFileNameInsert(Const::fidxfn, sax_file, data_file);
    if(FileUtil::checkFileExists(sax_file.c_str())){
        long series_num = FileUtil::getFileSize(sax_file.c_str()) / sizeof(unsigned short ) / Const::segmentNum;
        FILE *sax_f = fopen(sax_file.c_str(), "rb");
        fread(node_saxes + (long)(*cur) * Const::segmentNum, sizeof(unsigned short), Const::segmentNum * series_num, sax_f);
        fclose(sax_f);
        (*cur) = (*cur) + series_num;
        FileUtil::renameFile(data_file,data_file + "tmp");
        FileUtil::renameFile(sax_file, sax_file + "tmp");
        leaf_files.push_back(data_file + "tmp");
        sax_files.push_back(sax_file + "tmp");
    }
    if(!children.empty()){
        for(FADASNode*child: children){
            if(child!= nullptr){
                child->collectSAXwords(node_saxes, cur, leaf_files, sax_files);
            }
        }
    }
}

void FADASNode::deleteSubtree(){
    if(children.empty())
        delete this;
    else{
        unordered_set<FADASNode*>childs;
        for(FADASNode* child: children){
            if(child != nullptr)
                childs.insert(child);
        }
        for(FADASNode* child:childs)
            child->deleteSubtree();
        delete this;
    }
}

void FADASNode::growIndex(unsigned short *node_saxes, bool need_free) {
    if(size <= Const::th)   return;
    auto start = chrono::system_clock::now();
//    int chosen_num = SaxUtil::findFirstGE(power_2, 1, Const::segmentNum + 1, size / Const::th + 1);
//    SAX_INFO* sax_info = statSAX();
//    chooseSegment(sax_info, chosen_num);
    determineSegments(node_saxes);
    int chosen_num = chosenSegments.size();
    // statistic children information in order to partition
    partUnit nodes[1<<chosen_num];
    for(int i=0;i<(1<<chosen_num);++i)
        nodes[i].id = i, nodes[i].size=0, nodes[i].pid = -1;
    vector<vector<int>>node_offsets(1<<chosen_num, vector<int>());

    for(int i=0;i<size;++i){
        int new_id = SaxUtil::extendSax(node_saxes + (long)offsets[i] * (Const::segmentNum), bits_cardinality,chosenSegments);
        nodes[new_id].size++;
        node_offsets[new_id].push_back(offsets[i]);
    }

    if(need_free) vector<int>().swap(offsets);

    int partNum = partitionNew(nodes, chosen_num);
    // build rest data node if any
//    for(auto &node:nodes)
//        if(node.size <= Const::th && node.pid == -1)
//            node.pid = partNum++;

    FADASNode* childrenList[partNum];
    for(int i=0;i<partNum;++i)  childrenList[i] = new FADASNode(this, i);
    children.resize(1 << chosen_num);
    for(int i=0;i<(1 << chosen_num);++i){
        if(nodes[i].size <= 0)  continue;
        else if(nodes[i].size > Const::th) {
            children[i] = new FADASNode(this, nodes[i].size, i);
            generateSaxAndCardinality(children[i], i);
            children[i]->offsets.resize(nodes[i].size);
            copy(node_offsets[i].begin(),  node_offsets[i].end(), children[i]->offsets.begin());
            vector<int>().swap(node_offsets[i]);
        }else if(partition_id == -1){
            children[i] = new FADASNode(this, nodes[i].size, i);
            generateSaxAndCardinality(children[i], i);
            vector<int>().swap(node_offsets[i]);
        }
        else{
            int _pid = nodes[i].pid;
            children[i] = childrenList[_pid];
            childrenList[_pid]->size += nodes[i].size;
            generateSaxAndCardinality4LeafNode(children[i], i);
            vector<int>().swap(node_offsets[i]);
        }
    }

    vector<vector<int>>().swap(node_offsets);
    auto end = chrono::system_clock::now();
    GROW_CPU_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();

    for(auto &child: children){
        if(child!= nullptr && child->size > Const::th){
//              cout << file_id << "+" << child->id<<endl;
            child->growIndex(node_saxes, true);
        }
    }

}


FADASNode *FADASNode::BuildIndexLessPack(string &datafn, string &saxfn, string &paafn, vector<vector<int>> *g) {
    Const::logPrint("Start building index.");
    auto start_t = chrono::system_clock::now();
    FileUtil::checkDirClean(Const::fidxfn.c_str());
    auto end = chrono::system_clock::now();
//    long series_num = generateSaxAndPaaTbl();
    long series_num = loadSax(saxfn);
//    loadPaa(paafn);
    auto start = chrono::system_clock::now();
    Const::logPrint("Finish building sax and paa table.");
    SAX_PAA_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

    auto* root = new FADASNode();
    root->size = series_num;
    for(int &i:root->bits_cardinality)  i=0;
    partUnit nodeIn1stLayer[Const::vertexNum];  // all the nodes in the 1st layer
    int *navids = new int[series_num];  // invSAX head for the whole dataset
    for(int i=0;i<Const::vertexNum;++i)
        nodeIn1stLayer[i].id = i, nodeIn1stLayer[i].size=0, nodeIn1stLayer[i].pid = -1;

    // get 1st layer node size
    for(long i=0;i<series_num;++i){
        unsigned short *asax = saxes + i * Const::segmentNum;
        int nav_id = SaxUtil::invSaxHeadFromSax(asax, Const::bitsCardinality, Const::segmentNum);
        navids[i] = nav_id;
        nodeIn1stLayer[nav_id].size++;
    }

    Const::logPrint("Finish statistic size of nodes in the 1st layer.");

    // partition 1st layer
    int partNum = partitionLessPack(nodeIn1stLayer, Const::segmentNum);
    Const::logPrint("Finish partition");
    vector<FADASNode*>childrenList;
    childrenList.reserve(partNum);
    for(int i=0;i<partNum;++i)
        childrenList.push_back(new FADASNode(1, i));

    root->children.resize(Const::vertexNum);
    for(int i=0;i<Const::vertexNum;++i){
        if(nodeIn1stLayer[i].pid == -1) {
            root->children[i] = new FADASNode(1, nodeIn1stLayer[i].size, i);
            root->children[i]->generateSaxAndCardIn1stLayer(i);
        }
        else{
            int pid = nodeIn1stLayer[i].pid;
            root->children[i] = childrenList[pid];
            childrenList[pid]->size += nodeIn1stLayer[i].size;
            childrenList[pid]->generateSaxAndCardIn1stLayer4LeafNode(i);
        }
    }
    Const::logPrint("Finish build index structure 1st layer.");

    thread IO(materialize1stLayerLessPack, datafn, root, navids, Const::fidxfn);

    // put data offsets to internal nodes in 1st layer
    for(int i=0;i<Const::vertexNum;++i)
        if(nodeIn1stLayer[i].size > Const::th)
            root->children[i]->offsets.reserve(nodeIn1stLayer[i].size);
    for(int i=0;i<series_num;++i){
        int nav_id = navids[i];
        if(nodeIn1stLayer[nav_id].size > Const::th) {
            root->children[nav_id]->offsets.push_back(i);
        }
    }
    Const::logPrint("data offsets have been put into nodes in the 1st layer.");
    end = chrono::system_clock::now();
    GROW_CPU_TIME_1st += chrono::duration_cast<chrono::microseconds>(end - start).count();

    int j = 0;
    Const::logPrint("start grow the index structure");
    for(int i=0;i<Const::vertexNum;++i){
        if(nodeIn1stLayer[i].size > Const::th)
            root->children[i]->growIndexLessPack();
        if(++j%10000 == 0)  Const::logPrint(to_string(j) + " nodes in the 1st layer has been processed.");
    }
    start = chrono::system_clock::now();
    GROW_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

    Const::logPrint("build index skeleton finished.");

    IO.join();
    Const::logPrint("Start materialize internal nodes in the 1st layer");
    for(int i=0;i<Const::vertexNum;++i)
        if(nodeIn1stLayer[i].size > Const::th)
            materializeInterNodeLessPack(root->children[i], saxes);
    Const::logPrint("build index successfully!");
    auto end_t = chrono::system_clock::now();
    cout << "Total time is " << chrono::duration_cast<chrono::microseconds>(end_t - start_t).count() / 1000 << "ms."<<endl;
    cout << "Building sax and paa total time is "<<SAX_PAA_TOTAL_TIME / 1000 <<"ms, cpu time is "
         << SAX_PAA_CPU_TIME / 1000 <<"ms, I/O read time is " << SAX_PAA_READ_TIME / 1000 << "ms."<<endl;

    cout << "During the process of building 1st layer index structure, CPU time is "<< GROW_CPU_TIME_1st / 1000 <<"ms, "<<endl;

    cout << "During the process of materializing 1st layer nodes, total time is "<< MAT1_TOTAL_TIME / 1000
         <<"ms, I/O read time is "<< MAT1_READ_TIME / 1000 <<"ms, CPU statistic time  is " << MAT1_CPU_TIME_STAT / 1000
         << "ms, CPU copy Time is " << MAT1_CPU_TIME_COPY / 1000 << "ms, while I/O write time  is " << MAT1_WRITE_TIME / 1000 << "ms. "<< endl;

    cout << "During the process of growing index structure, total time is "<< GROW_TOTAL_TIME / 1000
         <<"ms, other CPU time is "<< GROW_CPU_TIME / 1000 <<"ms."<<endl;

    cout << "During the process of materializing internal nodes, total time is "<<MAT2_TOTAL_TIME / 1000
         <<"ms, I/O read time is "<<MAT2_READ_TIME / 1000 <<"ms, CPU time is " << MAT2_CPU_TIME / 1000
         << "ms, while I/O write time is " << MAT2_WRITE_TIME / 1000 << "ms." << endl;

    cout << "Random read count = " << RAND_READ_CNT << endl
         << "Random write count = " << RAND_WRITE_CNT << endl
         <<"Sequence read count = " << SEQ_READ_CNT << endl
         <<"Sequence write count = " << SEQ_WRITE_CNT << endl
         <<"Small file read bytes = "<<SMALL_FILES_BYTES_READ << endl
         <<"Small file write bytes = " << SMALL_FILES_BYTES_WRITE << endl;

    return root;
}

FADASNode *FADASNode::BuildIndexPos(string &datafn, string &saxfn, string &paafn, vector<vector<int>> *g) {
    Const::logPrint("Start building index.");
    auto start_t = chrono::system_clock::now();
    FileUtil::checkDirClean(Const::fidxfn.c_str());
    auto end = chrono::system_clock::now();
    int series_num = generateSaxAndPaaTbl();
//    int series_num = loadSax(saxfn);
//    loadPaa(paafn);
    auto start = chrono::system_clock::now();
    Const::logPrint("Finish building sax and paa table.");
    SAX_PAA_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

    auto* root = new FADASNode();
    root->size = series_num;
    for(int &i:root->bits_cardinality)  i=0;
    partUnit nodeIn1stLayer[Const::vertexNum];
    int *navids = new int[series_num];
    for(int i=0;i<Const::vertexNum;++i)
        nodeIn1stLayer[i].id = i, nodeIn1stLayer[i].size=0, nodeIn1stLayer[i].pid = -1;

    // get 1st layer node size
    for(int i=0;i<series_num;++i){
        unsigned short *asax = saxes + i * Const::segmentNum;
        int nav_id = SaxUtil::invSaxHeadFromSax(asax, Const::bitsCardinality, Const::segmentNum);
        navids[i] = nav_id;
        nodeIn1stLayer[nav_id].size++;
    }

    Const::logPrint("Finish statistic size of nodes in the 1st layer.");

    // partition 1st layer
    int partNum = partition1stLayer(nodeIn1stLayer, g, Const::filling_factor_1st);
    Const::logPrint("Finish partition");
    FADASNode* childrenList[partNum];
    for(int i=0;i<partNum;++i)  childrenList[i] = new FADASNode(1, i);
    root->children.resize(Const::vertexNum);
    for(int i=0;i<Const::vertexNum;++i){
        if(nodeIn1stLayer[i].pid == -1) {
            root->children[i] = new FADASNode(1, nodeIn1stLayer[i].size, i);
            root->children[i]->generateSaxAndCardIn1stLayer(i);
        }
        else{
            int pid = nodeIn1stLayer[i].pid;
            root->children[i] = childrenList[pid];
            childrenList[pid]->size += nodeIn1stLayer[i].size;
            childrenList[pid]->generateSaxAndCardIn1stLayer4LeafNode(i);
        }
    }
    Const::logPrint("Finish build index structure 1st layer.");

    thread IO(materialize1stLayerPos, datafn, root, navids, Const::posidxfn);

    // put data offsets to internal nodes in 1st layer
    for(int i=0;i<Const::vertexNum;++i)
        if(nodeIn1stLayer[i].size > Const::th)
            root->children[i]->offsets.reserve(nodeIn1stLayer[i].size);
    for(int i=0;i<series_num;++i){
        int nav_id = navids[i];
        if(nodeIn1stLayer[nav_id].size > Const::th) {
            root->children[nav_id]->offsets.push_back(i);
        }
    }
    Const::logPrint("data offsets have been put into nodes in the 1st layer.");
    end = chrono::system_clock::now();
    GROW_CPU_TIME_1st += chrono::duration_cast<chrono::microseconds>(end - start).count();

    int j = 0;
    Const::logPrint("start grow the index structure");
    for(int i=0;i<Const::vertexNum;++i){
        if(nodeIn1stLayer[i].size > Const::th)
            root->children[i]->growIndex();
        if(++j%10000 == 0)  Const::logPrint(to_string(j) + " nodes in the 1st layer has been processed.");
    }
    start = chrono::system_clock::now();
    GROW_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

    Const::logPrint("build index skeleton finished.");

    IO.join();
    Const::logPrint("Start materialize internal nodes in the 1st layer");
    for(int i=0;i<Const::vertexNum;++i)
        if(nodeIn1stLayer[i].size > Const::th)
            materializeInterNodePos(root->children[i], saxes);
    Const::logPrint("build index successfully!");
    auto end_t = chrono::system_clock::now();
    cout << "Total time is " << chrono::duration_cast<chrono::microseconds>(end_t - start_t).count() / 1000 << "ms."<<endl;
    cout << "Building sax and paa total time is "<<SAX_PAA_TOTAL_TIME / 1000 <<"ms, cpu time is "
         << SAX_PAA_CPU_TIME / 1000 <<"ms, I/O read time is " << SAX_PAA_READ_TIME / 1000 << "ms."<<endl;

    cout << "During the process of building 1st layer index structure, CPU time is "<< GROW_CPU_TIME_1st / 1000 <<"ms, "<<endl;

    cout << "During the process of materializing 1st layer nodes, total time is "<< MAT1_TOTAL_TIME / 1000
         <<"ms, I/O read time is "<< MAT1_READ_TIME / 1000 <<"ms, CPU statistic time  is " << MAT1_CPU_TIME_STAT / 1000
         << "ms, while I/O write time  is " << MAT1_WRITE_TIME / 1000 << "ms. "<< endl;

    cout << "During the process of growing index structure, total time is "<< GROW_TOTAL_TIME / 1000
         <<"ms, other CPU time is "<< GROW_CPU_TIME / 1000 <<"ms."<<endl;

    cout << "During the process of materializing internal nodes, total time is "<<MAT2_TOTAL_TIME / 1000
         <<"ms, I/O read time is "<<MAT2_READ_TIME / 1000 <<"ms, CPU time  is " << MAT2_CPU_TIME / 1000
         << "ms, while I/O write time is " << MAT2_WRITE_TIME / 1000 << "ms." << endl;

    cout << "Random read count = " << RAND_READ_CNT << endl
         << "Random write count = " << RAND_WRITE_CNT << endl
         <<"Sequence read count = " << SEQ_READ_CNT << endl
         <<"Sequence write count = " << SEQ_WRITE_CNT << endl;

    return root;
}

//void FADASNode::growIndex() {
//    if(size <= Const::th)   return;
//    auto start = chrono::system_clock::now();
//    int chosen_num = SaxUtil::findFirstGE(power_2, 1, Const::segmentNum + 1, size / Const::th + 1);
//    PAA_INFO* paa = statPaa();
//    chooseSegment(paa, chosen_num);
//
//    // statistic children information in order to partition
//    partUnit nodes[1<<chosen_num];
//    for(int i=0;i<(1<<chosen_num);++i)
//        nodes[i].id = i, nodes[i].size=0, nodes[i].pid = -1;
//    vector<vector<int>>node_offsets(1<<chosen_num, vector<int>());
//
//    for(int i=0;i<size;++i){
//        int new_id = SaxUtil::extendSax(FADASNode::saxes + (long)offsets[i] * (Const::segmentNum), bits_cardinality,chosenSegments);
//        nodes[new_id].size++;
//        node_offsets[new_id].push_back(offsets[i]);
//    }
//
//    if(this->layer > 1) vector<int>().swap(offsets);
//
//    int partNum = partition(nodes);
//    // build rest data node if any
//    for(auto &node:nodes)
//        if(node.size <= Const::th && node.pid == -1)
//            node.pid = ++partNum;
//
//
//    FADASNode* childrenList[partNum];
//    for(int i=0;i<partNum;++i)  childrenList[i] = new FADASNode(this, i);
//    children.resize(1 << chosen_num);
//    for(int i=0;i<(1 << chosen_num);++i){
//        if(nodes[i].size <= 0)  continue;
//        else if(nodes[i].pid == -1) {
//            children[i] = new FADASNode(this, nodes[i].size, i);
//            generateSaxAndCardinality(children[i], i);
//            children[i]->offsets.resize(nodes[i].size);
//            copy(node_offsets[i].begin(),  node_offsets[i].end(), children[i]->offsets.begin());
//            vector<int>().swap(node_offsets[i]);
//        }
//        else{
//            int _pid = nodes[i].pid;
//            children[i] = childrenList[_pid];
//            childrenList[_pid]->size += nodes[i].size;
//            generateSaxAndCardinality4LeafNode(children[i], i);
//            vector<int>().swap(node_offsets[i]);
//        }
//    }
//
//    vector<vector<int>>().swap(node_offsets);
//    auto end = chrono::system_clock::now();
//    GROW_CPU_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();
//
//    for(auto &child: children){
//        if(child!= nullptr && child->size > Const::th){
//            child->growIndex();
//        }
//    }
//
//}

void FADASNode::growIndex(){
    if(size <= Const::th)   return;
    auto start = chrono::system_clock::now();
//    int chosen_num = SaxUtil::findFirstGE(power_2, 1, Const::segmentNum + 1, size / Const::th + 1);
//    SAX_INFO* sax_info = statSAX();
//    chooseSegment(sax_info, chosen_num);
    determineSegments();
    int chosen_num = chosenSegments.size();
    // statistic children information in order to partition
    partUnit nodes[1<<chosen_num];
    for(int i=0;i<(1<<chosen_num);++i)
        nodes[i].id = i, nodes[i].size=0, nodes[i].pid = -1;
    vector<vector<int>>node_offsets(1<<chosen_num, vector<int>());

    for(int i=0;i<size;++i){
        int new_id = SaxUtil::extendSax(FADASNode::saxes + (long)offsets[i] * (Const::segmentNum), bits_cardinality,chosenSegments);
        nodes[new_id].size++;
        node_offsets[new_id].push_back(offsets[i]);
    }

    if(this->layer > 1) vector<int>().swap(offsets);

    int partNum = partitionNew(nodes, chosen_num);
    // build rest data node if any
//    for(auto &node:nodes)
//        if(node.size <= Const::th && node.pid == -1)
//            node.pid = partNum++;

    FADASNode* childrenList[partNum];
    for(int i=0;i<partNum;++i)  childrenList[i] = new FADASNode(this, i);
    children.resize(1 << chosen_num);
    for(int i=0;i<(1 << chosen_num);++i){
        if(nodes[i].size <= 0)  continue;
        else if(nodes[i].size > Const::th) {
            children[i] = new FADASNode(this, nodes[i].size, i);
            generateSaxAndCardinality(children[i], i);
            children[i]->offsets.resize(nodes[i].size);
            copy(node_offsets[i].begin(),  node_offsets[i].end(), children[i]->offsets.begin());
            vector<int>().swap(node_offsets[i]);
        }else if(nodes[i].pid == -1){
            children[i] = new FADASNode(this, nodes[i].size, i);
            generateSaxAndCardinality(children[i], i);
            vector<int>().swap(node_offsets[i]);
        }
        else{
            int _pid = nodes[i].pid;
            children[i] = childrenList[_pid];
            childrenList[_pid]->size += nodes[i].size;
            generateSaxAndCardinality4LeafNode(children[i], i);
            vector<int>().swap(node_offsets[i]);
        }
    }

    vector<vector<int>>().swap(node_offsets);
    auto end = chrono::system_clock::now();
    GROW_CPU_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();

    for(auto &child: children){
        if(child!= nullptr && child->size > Const::th){
//              cout << file_id << "+" << child->id<<endl;
            child->growIndex();
        }
    }

}

void FADASNode::growIndexLessPack() {
    if(size <= Const::th)   return;
    auto start = chrono::system_clock::now();
    int chosen_num = SaxUtil::findFirstGE(power_2, 1, Const::segmentNum + 1, size / Const::th + 1);
    SAX_INFO* sax_info = statSAX();
    chooseSegment(sax_info, chosen_num);

    // statistic children information in order to partition
    partUnit nodes[1<<chosen_num];
    for(int i=0;i<(1<<chosen_num);++i)
        nodes[i].id = i, nodes[i].size=0, nodes[i].pid = -1;
    vector<vector<int>>node_offsets(1<<chosen_num, vector<int>());

    for(int i=0;i<size;++i){
        int new_id = SaxUtil::extendSax(FADASNode::saxes + (long)offsets[i] * (Const::segmentNum), bits_cardinality,chosenSegments);
        nodes[new_id].size++;
        node_offsets[new_id].push_back(offsets[i]);
    }

    if(this->layer > 1) vector<int>().swap(offsets);

    int partNum = partitionLessPack(nodes, chosen_num);

    vector<FADASNode*>childrenList;
    for(int i=0;i<partNum;++i)
        childrenList.push_back(new FADASNode(this, i));

    children.resize(1 << chosen_num);
    for(int i=0;i<(1 << chosen_num);++i){
        if(nodes[i].size <= 0)  continue;
        else if(nodes[i].pid == -1) {
            children[i] = new FADASNode(this, nodes[i].size, i);
            generateSaxAndCardinality(children[i], i);
            if(nodes[i].size > Const::th){
                children[i]->offsets.resize(nodes[i].size);
                copy(node_offsets[i].begin(),  node_offsets[i].end(), children[i]->offsets.begin());
            }
            vector<int>().swap(node_offsets[i]);
        }
        else{
            int _pid = nodes[i].pid;
            children[i] = childrenList[_pid];
            childrenList[_pid]->size += nodes[i].size;
            generateSaxAndCardinality4LeafNode(children[i], i);
            vector<int>().swap(node_offsets[i]);
        }
    }

    vector<vector<int>>().swap(node_offsets);
    auto end = chrono::system_clock::now();
    GROW_CPU_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();

    for(auto &child: children){
        if(child!= nullptr && child->size > Const::th){
            child->growIndexLessPack();
        }
    }

}

// it has some bugs, but I cannot find them
FADASNode *FADASNode::BuildIndexWOPack(string &datafn, string &saxfn, string &paafn, vector<vector<int>> *g) {
    Const::logPrint("Start building index.");
    auto start_t = chrono::system_clock::now();
    FileUtil::checkDirClean(Const::fidxfn.c_str());
    auto end = chrono::system_clock::now();
//    long series_num = generateSaxAndPaaTbl();
    long series_num = loadSax(saxfn);
//    loadPaa(paafn);
    auto start = chrono::system_clock::now();
    Const::logPrint("Finish building sax and paa table.");
    SAX_PAA_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

    auto* root = new FADASNode();
    root->size = series_num;
    for(int &i:root->bits_cardinality)  i=0;
    partUnit nodeIn1stLayer[Const::vertexNum];
    int *navids = new int[series_num];
    for(int i=0;i<Const::vertexNum;++i)
        nodeIn1stLayer[i].id = i, nodeIn1stLayer[i].size=0, nodeIn1stLayer[i].pid = -1;

    // get 1st layer node size
    for(long i=0;i<series_num;++i){
        unsigned short *asax = saxes + i * Const::segmentNum;
        int nav_id = SaxUtil::invSaxHeadFromSax(asax, Const::bitsCardinality, Const::segmentNum);
        navids[i] = nav_id;
        nodeIn1stLayer[nav_id].size++;
    }

    Const::logPrint("Finish statistic size of nodes in the 1st layer.");

    // partition 1st layer
//    int partNum = partition1stLayer(nodeIn1stLayer, g, Const::filling_factor_1st);
    int partNum = 0;
    for(int i=0;i<Const::vertexNum;++i){
        if(nodeIn1stLayer[i].size <= Const::th && nodeIn1stLayer[i].size > 0){
            nodeIn1stLayer[i].pid = partNum++;
        }
    }
//    Const::logPrint("Finish partition");
//    int partNum = 0;
    FADASNode* childrenList[partNum];
    for(int i=0;i<partNum;++i)  childrenList[i] = new FADASNode(1, i);
    root->children.resize(Const::vertexNum);
    for(int i=0;i<Const::vertexNum;++i){
        if(nodeIn1stLayer[i].size <=0)  continue;
        if(nodeIn1stLayer[i].pid == -1) {
            assert(nodeIn1stLayer[i].size > Const::th);
            root->children[i] = new FADASNode(1, nodeIn1stLayer[i].size, i);
            root->children[i]->generateSaxAndCardIn1stLayer(i);
        }
        else{
            int pid = nodeIn1stLayer[i].pid;
            root->children[i] = childrenList[pid];
            childrenList[pid]->size += nodeIn1stLayer[i].size;
            childrenList[pid]->generateSaxAndCardIn1stLayer4LeafNode(i);
        }
    }
    Const::logPrint("Finish build index structure 1st layer.");

//    thread IO(materialize1stLayer, datafn, root, navids, Const::fidxfn);
    materialize1stLayer(datafn, root, navids, Const::fidxfn);

    // put data offsets to internal nodes in 1st layer
    for(int i=0;i<Const::vertexNum;++i)
        if(nodeIn1stLayer[i].size > Const::th)
            root->children[i]->offsets.reserve(nodeIn1stLayer[i].size);
    for(int i=0;i<series_num;++i){
        int nav_id = navids[i];
        if(nodeIn1stLayer[nav_id].size > Const::th) {
            root->children[nav_id]->offsets.push_back(i);
        }
    }
    Const::logPrint("data offsets have been put into nodes in the 1st layer.");
    end = chrono::system_clock::now();
    GROW_CPU_TIME_1st += chrono::duration_cast<chrono::microseconds>(end - start).count();

    int j = 0;
    Const::logPrint("start grow the index structure");
    for(int i=0;i<Const::vertexNum;++i){
        if(nodeIn1stLayer[i].size > Const::th) {
            cout << i <<endl;
            root->children[i]->growIndexWOPack();
        }
        if(++j%10000 == 0)  Const::logPrint(to_string(j) + " nodes in the 1st layer has been processed.");
    }
    start = chrono::system_clock::now();
    GROW_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

    Const::logPrint("build index skeleton finished.");

//    IO.join();
    Const::logPrint("Start materialize internal nodes in the 1st layer");
    for(int i=0;i<Const::vertexNum;++i)
        if(nodeIn1stLayer[i].size > Const::th)
            materializeInterNode(root->children[i], saxes);
    Const::logPrint("build index successfully!");
    auto end_t = chrono::system_clock::now();
    cout << "Total time is " << chrono::duration_cast<chrono::microseconds>(end_t - start_t).count() / 1000 << "ms."<<endl;
    cout << "Building sax and paa total time is "<<SAX_PAA_TOTAL_TIME / 1000 <<"ms, cpu time is "
         << SAX_PAA_CPU_TIME / 1000 <<"ms, I/O read time is " << SAX_PAA_READ_TIME / 1000 << "ms."<<endl;

    cout << "During the process of building 1st layer index structure, CPU time is "<< GROW_CPU_TIME_1st / 1000 <<"ms, "<<endl;

    cout << "During the process of materializing 1st layer nodes, total time is "<< MAT1_TOTAL_TIME / 1000
         <<"ms, I/O read time is "<< MAT1_READ_TIME / 1000 <<"ms, CPU statistic time  is " << MAT1_CPU_TIME_STAT / 1000
         << "ms, CPU copy Time is " << MAT1_CPU_TIME_COPY / 1000 << "ms, while I/O write time  is " << MAT1_WRITE_TIME / 1000 << "ms. "<< endl;

    cout << "During the process of growing index structure, total time is "<< GROW_TOTAL_TIME / 1000
         <<"ms, other CPU time is "<< GROW_CPU_TIME / 1000 <<"ms."<<endl;

    cout << "During the process of materializing internal nodes, total time is "<<MAT2_TOTAL_TIME / 1000
         <<"ms, I/O read time is "<<MAT2_READ_TIME / 1000 <<"ms, CPU time is " << MAT2_CPU_TIME / 1000
         << "ms, while I/O write time is " << MAT2_WRITE_TIME / 1000 << "ms." << endl;

    cout << "Random read count = " << RAND_READ_CNT << endl
         << "Random write count = " << RAND_WRITE_CNT << endl
         <<"Sequence read count = " << SEQ_READ_CNT << endl
         <<"Sequence write count = " << SEQ_WRITE_CNT << endl
         <<"Small file read bytes = "<<SMALL_FILES_BYTES_READ << endl
         <<"Small file write bytes = " << SMALL_FILES_BYTES_WRITE << endl;

    return root;
}

void FADASNode::growIndexWOPack() {
    if(size <= Const::th)   return;
    cout << file_id << ", " << id <<", " << partition_id <<", " <<size<<endl;
    auto start = chrono::system_clock::now();
    int chosen_num = SaxUtil::findFirstGE(power_2, 1, Const::segmentNum + 1, size / Const::th + 1);
    SAX_INFO* sax_info = statSAX();
    chooseSegment(sax_info, chosen_num);

    cout << "has chosen segments" << endl;
    // statistic children information in order to partition
    partUnit nodes[1<<chosen_num];
    for(int i=0;i<(1<<chosen_num);++i)
        nodes[i].id = i, nodes[i].size=0, nodes[i].pid = -1;
    vector<vector<int>>node_offsets(1<<chosen_num, vector<int>());

    cout <<"here1"<<endl;
    for(int i=0;i<size;++i){
        int new_id = SaxUtil::extendSax(FADASNode::saxes + (long)offsets[i] * (Const::segmentNum), bits_cardinality,chosenSegments);
        nodes[new_id].size++;
        node_offsets[new_id].push_back(offsets[i]);
    }

    cout <<"here2"<<endl;
    if(this->layer > 1) vector<int>().swap(offsets);

//    int partNum = partition(nodes);
    // build rest data node if any
    int partNum = 0;
    for(int i=0;i<(1<<chosen_num);++i){
        if(nodes[i].size >0 && nodes[i].size <= Const::th)
            nodes[i].pid = partNum++;
    }
    cout <<"here3"<<endl;
//    for(auto &node:nodes)
//        if(node.size <= Const::th && node.pid == -1)
//            node.pid = ++partNum;
//
//
    vector<FADASNode*>childrenList;
    if(partNum >0) {
        childrenList.resize(partNum);
        for (int i = 0; i <= partNum; ++i) childrenList[i] = new FADASNode(this, i);
    }
    children.resize(1 << chosen_num);
    cout <<"here4"<<endl;
    for(int i=0;i<(1 << chosen_num);++i){
        if(nodes[i].size <= 0)  continue;
        else if(nodes[i].pid == -1){
            children[i] = new FADASNode(this, nodes[i].size, i);
            generateSaxAndCardinality(children[i], i);
            children[i]->offsets.resize(nodes[i].size);
            copy(node_offsets[i].begin(),  node_offsets[i].end(), children[i]->offsets.begin());
            vector<int>().swap(node_offsets[i]);
        }
        else{
            int _pid = nodes[i].pid;
            children[i] = childrenList[_pid];
            childrenList[_pid]->size += nodes[i].size;
            generateSaxAndCardinality4LeafNode(children[i], i);
            vector<int>().swap(node_offsets[i]);
        }
    }

    cout <<"here5"<<endl;

    vector<vector<int>>().swap(node_offsets);
    auto end = chrono::system_clock::now();
    GROW_CPU_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();

    for(auto &child: children){
        if(child!= nullptr && child->size > Const::th){
            child->growIndexWOPack();
        }
    }

}


// TODO: this function may be optimized with SIMD
PAA_INFO* FADASNode::statPaa(){
    auto* r = new PAA_INFO();
    double split_line[Const::segmentNum],paa_max[Const::segmentNum],paa_min[Const::segmentNum],paa_mu[Const::segmentNum];
    // TODO: optimize
    double lb;  // ub is the new split line
    for(int i=0; i < Const::segmentNum; ++i)
        SaxUtil::getValueRange(sax[i] << 1, bits_cardinality[i] + 1, &lb, &split_line[i]);
    for(auto &i:paa_max) i = - numeric_limits<double>::max();
    for(auto &i:paa_min) i = numeric_limits<double>::max();
    for(auto &i:r->paa_up_size) i = 0;
    for(auto &i:r->paa_below_size) i = 0;
    for(auto &i:r->paa_variance) i = 0;
    for(auto &i:paa_mu) i=0;
    for(long offset:offsets){
        float * start = paas + offset * (Const::segmentNum);
        for(int i=0;i<Const::segmentNum;++i){
            double value = *(start + i);
            paa_mu[i] += value;
            paa_min[i] = min(paa_min[i], value);
            paa_max[i] = max(paa_max[i], value);
            if(value > split_line[i]) {
                r->paa_up_size[i]++;
            }
            else {
                r->paa_below_size[i]++;
            }
        }
    }
    for(double & i : paa_mu) {
        i /= size;
    }

    for(long offset:offsets){
        float * start = paas + offset * (Const::segmentNum);
        for(int i=0;i<Const::segmentNum;++i){
            double value = *(start + i);
            r->paa_variance[i] += (value - paa_mu[i]) * (value - paa_mu[i]);
        }
    }
    return r;
}

SAX_INFO* FADASNode::statSAX(){
    auto* r = new SAX_INFO();
    vector<vector<double>>numerical_sax(Const::segmentNum, vector<double>(size, 0));
    vector<double>numerical_sax_sum(Const::segmentNum, 0);
    // TODO: optimize
    for(auto &i:r->sax_variance) i = 0;
    for(int j = 0;j<size;++j){
        int offset = offsets[j];
        unsigned short* start = saxes + offset * (Const::segmentNum);
//        cout << endl;
        for(int i=0;i<Const::segmentNum;++i){
            unsigned short symbol = *(start + i);
//            cout << symbol << ",";
            double numer = SaxUtil::getMidLineFromSaxSymbolbc8(symbol);
            numerical_sax[i][j] = numer;
            numerical_sax_sum[i] += numer;
            bool isUp = (symbol >> (Const::bitsCardinality - 1 - bits_cardinality[i])) % 2;
            if(isUp)    r->up_size[i]++;
            else    r->below_size[i]++;
        }
    }
    for(double & i : numerical_sax_sum)
        i /= size;

    for(int i=0; i<Const::segmentNum;++i){
        for(int j=0;j<size;++j){
            r->sax_variance[i] += (numerical_sax[i][j] - numerical_sax_sum[i]) * (numerical_sax[i][j] - numerical_sax_sum[i]);
        }
    }
    return r;
}


struct tmp{
    int i{};
    double score{};
    tmp(int _i, double _score){i=_i;score = _score;}
    tmp(){;}

    static bool order(tmp a,tmp b){
        return a.score < b.score;
    }

    static bool orderdesc(tmp a,tmp b){
        return a.score > b.score;
    }
};

void FADASNode::determineFanout(int *lambda_min, int * lambda_max) const{
    if(size < 2 * Const::th)    {
        *lambda_min = 1;
        *lambda_max = 1;
        return;
    }
    *lambda_min = -1;
    *lambda_max = Const::segmentNum;
    double _min = size / (Const::th * Const::f_high);
    double _max = size / (Const::th * Const::f_low);
    if(Const::vertexNum < _min){
        *lambda_min = Const::segmentNum;
        *lambda_max = Const::segmentNum;
        return;
    }
    for(int i = 1; i <= Const::segmentNum; ++i){
        if(*lambda_min == -1){
            if((1<< i) >= _min){
                *lambda_min = i;
            }
        }else{
            if((1<<i) == _max){
                *lambda_max = i;
                break;
            }else if((1<<i) > _max){
                *lambda_max = max(i-1,*lambda_min);
                break;
            }
        }
    }
}

// determine fan-out and choose segments
void FADASNode::determineSegmentsAvgVariance(){
    int lambda_min, lambda_max;
    determineFanout(&lambda_min, &lambda_max);

    vector<unordered_map<unsigned short, int>>data_seg_symbols(Const::segmentNum);
    for(int offset:offsets){
        unsigned short* cur_sax = saxes + offset * Const::segmentNum;
        for(int i=0;i<Const::segmentNum;++i){
            data_seg_symbols[i][cur_sax[i]]++;
        }
    }

    // compute stdev of each segment
    vector<double>data_seg_mean(Const::segmentNum, 0);
    vector<double>data_seg_stdev(Const::segmentNum, 0);
    for(int i=0;i<Const::segmentNum;++i){
        auto& map = data_seg_symbols[i];
        for(auto &iter:map){
            unsigned short symbol = iter.first;
            data_seg_mean[i] += (SaxUtil::getMidLineFromSaxSymbolbc8(symbol) * iter.second);
        }
        data_seg_mean[i] /= size;
        for(auto &iter:map){
            unsigned short symbol = iter.first;
            double mid_value = SaxUtil::getMidLineFromSaxSymbolbc8(symbol);
            data_seg_stdev[i] += (iter.second * ((mid_value - data_seg_mean[i]) * (mid_value - data_seg_mean[i])));
        }
//        data_seg_stdev[i] = sqrt(data_seg_stdev[i] / size);
        data_seg_stdev[i] = data_seg_stdev[i] / size;
    }

    vector<double>().swap(data_seg_mean);
    vector<unordered_map<unsigned short, int>>().swap(data_seg_symbols);

    // start to compute the size of each node in each plan

    unordered_set<int>visited;
    double max_score = 0;
    vector<int> best_plan;
    for(int lambda = lambda_min; lambda <= lambda_max;++lambda){
        int plan_num = FADASNode::combine_num[lambda];
        for(int i=0;i<plan_num;++i){
            int *plan = FADASNode::combines[lambda][i];
            double score = 0;
            for(int j = 0; j < lambda; ++j){
                score += data_seg_stdev[plan[j]];
            }
            score /= lambda;
            if(score > max_score){
                max_score = score;
                best_plan.clear();
                for(int j = 0; j<lambda;++j)
                    best_plan.push_back(plan[j]);
            }
        }
    }

    unordered_set<int>().swap(visited);
    chosenSegments = best_plan;
}

#include <cstdlib>
void FADASNode::determineSegmentsNaive() {
    int lambda_min, lambda_max;
    determineFanout(&lambda_min, &lambda_max);
    unordered_set<int>chosen;
    srand(time(nullptr));
    for(int i=0;i<lambda_max;++i){
        int tmp;
        while(1) {
            tmp = rand() % Const::segmentNum;
            if(chosen.count(tmp) <=0){
                chosen.insert(tmp);
                break;
            }
        }
    }
    chosenSegments.assign(chosen.begin(),  chosen.end());
}
// determine fan-out and choose segments
void FADASNode::determineSegments(){
    int lambda_min, lambda_max;
    determineFanout(&lambda_min, &lambda_max);
    if(lambda_min == Const::segmentNum && lambda_max == Const::segmentNum){
        for(int i=0;i<Const::segmentNum;++i)
            chosenSegments.push_back(i);
        return;
    }
    vector<int>unit_size(Const::vertexNum, 0);

    vector<unordered_map<unsigned short, int>>data_seg_symbols(Const::segmentNum);
    for(int offset:offsets){
        unsigned short* cur_sax = saxes + offset * Const::segmentNum;
        for(int i=0;i<Const::segmentNum;++i){
            data_seg_symbols[i][cur_sax[i]]++;
        }
        int head = SaxUtil::extendSax(cur_sax, bits_cardinality);
        unit_size[head]++;
    }

    // compute stdev of each segment
    vector<double>data_seg_mean(Const::segmentNum, 0);
    vector<double>data_seg_stdev(Const::segmentNum, 0);
    for(int i=0;i<Const::segmentNum;++i){
        auto& map = data_seg_symbols[i];
        for(auto &iter:map){
            unsigned short symbol = iter.first;
            data_seg_mean[i] += (SaxUtil::getMidLineFromSaxSymbolbc8(symbol) * iter.second);
        }
        data_seg_mean[i] /= size;
        for(auto &iter:map){
            unsigned short symbol = iter.first;
            double mid_value = SaxUtil::getMidLineFromSaxSymbolbc8(symbol);
            data_seg_stdev[i] += (iter.second * ((mid_value - data_seg_mean[i]) * (mid_value - data_seg_mean[i])));
        }
//        data_seg_stdev[i] = sqrt(data_seg_stdev[i] / size);
        data_seg_stdev[i] /= size;
    }

    vector<double>().swap(data_seg_mean);
    vector<unordered_map<unsigned short, int>>().swap(data_seg_symbols);

    // start to compute the size of each node in each plan
    int plan_num;
    if(lambda_max < Const::segmentNum)
        plan_num = FADASNode::combine_num[lambda_max];
    else
        plan_num = 1;
    unordered_set<int>visited;
    double max_score = 0;
    vector<int> best_plan;
    for(int i=0;i<plan_num;++i){
        int *plan = FADASNode::combines[lambda_max][i];
        // first evaluate the whole plan
        vector<int>plan_node_sizes(1<<lambda_max, 0);
        int mask_code = MathUtil::generateMaskSettingKbits(plan, lambda_max, Const::segmentNum);
        map<int,int>max_node_size;
        for(int j=0;j<Const::vertexNum;++j){
            max_node_size[mask_code & j] += unit_size[j];
        }
//        assert(max_node_size.size() == (1 << lambda_max));
        int _ = 0;
        for(auto & iter: max_node_size){
            plan_node_sizes[_++] = iter.second;
        }
        map<int,int>().swap(max_node_size);

        double score = compute_score(plan_node_sizes, plan, lambda_max, data_seg_stdev);
        if(score > max_score){
            max_score = score;
            best_plan.clear();
            for(int j = 0; j<lambda_max;++j)
                best_plan.push_back(plan[j]);
        }

        if(lambda_min <= lambda_max - 1)
            visitPlanFromBaseTable(visited, lambda_max - 1, plan, plan_node_sizes,
                                   &max_score, best_plan, lambda_min, mask_code, data_seg_stdev, score);
        vector<int>().swap(plan_node_sizes);
    }

    unordered_set<int>().swap(visited);
    vector<int>().swap(unit_size);
    chosenSegments = best_plan;
}

void FADASNode::determineSegments(unsigned short *node_saxes) {
    int lambda_min, lambda_max;
    determineFanout(&lambda_min, &lambda_max);

    vector<int>unit_size(Const::vertexNum, 0);

    vector<unordered_map<unsigned short, int>>data_seg_symbols(Const::segmentNum);
    for(int offset: offsets){
        unsigned short* cur_sax = node_saxes + (long)offset * Const::segmentNum;
        for(int i=0;i<Const::segmentNum;++i){
            assert(cur_sax[i] >=0 && cur_sax[i] <= 255);
            data_seg_symbols[i][cur_sax[i]]++;
        }
        int head = SaxUtil::extendSax(cur_sax, bits_cardinality);
        unit_size[head]++;
    }

    // compute stdev of each segment
    vector<double>data_seg_mean(Const::segmentNum, 0);
    vector<double>data_seg_stdev(Const::segmentNum, 0);
    for(int i=0;i<Const::segmentNum;++i){
        auto& map = data_seg_symbols[i];
        for(auto &iter:map){
            unsigned short symbol = iter.first;
            data_seg_mean[i] += (SaxUtil::getMidLineFromSaxSymbolbc8(symbol) * iter.second);
        }
        data_seg_mean[i] /= size;
        for(auto &iter:map){
            unsigned short symbol = iter.first;
            double mid_value = SaxUtil::getMidLineFromSaxSymbolbc8(symbol);
            data_seg_stdev[i] += (iter.second * ((mid_value - data_seg_mean[i]) * (mid_value - data_seg_mean[i])));
        }
//        data_seg_stdev[i] = sqrt(data_seg_stdev[i] / size);
        data_seg_stdev[i] /= size;
    }

    vector<double>().swap(data_seg_mean);
    vector<unordered_map<unsigned short, int>>().swap(data_seg_symbols);

    // start to compute the size of each node in each plan
    int plan_num = FADASNode::combine_num[lambda_max];
    unordered_set<int>visited;
    double max_score = 0;
    vector<int> best_plan;
    for(int i=0;i<plan_num;++i){
        int *plan = FADASNode::combines[lambda_max][i];
        // first evaluate the whole plan
        vector<int>plan_node_sizes(1<<lambda_max, 0);
        int mask_code = MathUtil::generateMaskSettingKbits(plan, lambda_max, Const::segmentNum);
        map<int,int>max_node_size;
        for(int j=0;j<Const::vertexNum;++j){
            max_node_size[mask_code & j] += unit_size[j];
        }
//        assert(max_node_size.size() == (1 << lambda_max));
        int _ = 0;
        for(auto & iter: max_node_size){
            plan_node_sizes[_++] = iter.second;
        }
        map<int,int>().swap(max_node_size);

        double score = compute_score(plan_node_sizes, plan, lambda_max, data_seg_stdev);
        if(score > max_score){
            max_score = score;
            best_plan.clear();
            for(int j = 0; j<lambda_max;++j)
                best_plan.push_back(plan[j]);
        }

        if(lambda_min <= lambda_max - 1)
            visitPlanFromBaseTable(visited, lambda_max - 1, plan, plan_node_sizes,
                                   &max_score, best_plan, lambda_min, mask_code, data_seg_stdev, score);
        vector<int>().swap(plan_node_sizes);
    }

    unordered_set<int>().swap(visited);
    vector<int>().swap(unit_size);
    chosenSegments = best_plan;
}

void
FADASNode::visitPlanFromBaseTable(unordered_set<int> &visited, int cur_lambda, const int *plan, vector<int> &base_tbl,
                                  double *max_score, vector<int> &best_plan, int lambda_min, int mask_code,
                                  vector<double> &data_seg_stdev, double base_score) {
    // base mask is used to detect base tbl
    int base_mask = 1;
    for(int i=0;i<cur_lambda;++i)
        base_mask = (base_mask << 1) +1;

    for(int i=0;i<cur_lambda + 1 ;++i){
        int reset_pos = plan[i];
        // get the whole plan code
        int cur_whole_mask = mask_code - (1 << (Const::segmentNum - 1  - reset_pos));
        if(visited.contains(cur_whole_mask))
            continue;
        visited.insert(cur_whole_mask);
        // get the new plan
        int *new_plan = new int[cur_lambda];
        int _=0;
        for(int j=0;j<cur_lambda+1;++j){
            if(i != j)  new_plan[_++] = plan[j];
        }

        // get the current base mask
        int cur_base_mask = base_mask - (1 << (cur_lambda - i));
        map<int,int>node_size_map;
        for(int j=0;j<base_tbl.size();++j)
            node_size_map[cur_base_mask &j] += base_tbl[j];
//        assert(node_size_map.size() == (1<<cur_lambda));
        vector<int>new_tbl(1<<cur_lambda, 0);
        _ = 0;
        for(auto &iter:node_size_map)
            new_tbl[_++] = iter.second;
        map<int,int>().swap(node_size_map);
        double score  = compute_score(new_tbl, new_plan, cur_lambda, data_seg_stdev);
        if(score > *max_score){
            *max_score = score;
            best_plan.clear();
            for(_ = 0; _<cur_lambda;++_)
                best_plan.push_back(new_plan[_]);
        }
//        if(score > base_score && cur_lambda > 1)
        if(cur_lambda > lambda_min)
            visitPlanFromBaseTable(visited, cur_lambda - 1, new_plan, new_tbl,
                                   max_score, best_plan, lambda_min, cur_whole_mask, data_seg_stdev, score);
        vector<int>().swap(new_tbl);
        delete[] new_plan;
    }
}

double FADASNode::compute_score(vector<int>&node_sizes, int *plan, int lambda, vector<double>&data_seg_stdev) const{
    if(size < 2*Const::th){
        if(node_sizes[0] > Const::th || node_sizes[1] > Const::th)
            return (double)min(node_sizes[0],node_sizes[1]) / Const::th;
        return data_seg_stdev[plan[0]] * 100;
    }
    int over_th_nodes_no = 0;
    for(int _:node_sizes){
        if(_ > Const::th)
            over_th_nodes_no++;
    }
    double w = ((double)over_th_nodes_no) / ((int)node_sizes.size());
    double sum_seg = 0;
    for(int i=0;i<lambda;++i){
        sum_seg += data_seg_stdev[plan[i]];
    }
    sum_seg  = sqrt(sum_seg / lambda);
    sum_seg = exp(1+sum_seg);
//    sum_seg  = sum_seg / lambda;
    auto *tmp = new double[node_sizes.size()];
    for(int i=0;i<node_sizes.size();++i){
        tmp[i] = ((double)node_sizes[i]) / Const::th;
    }
    double stdev_fill_factor = MathUtil::deviation(tmp, node_sizes.size());

//    double balance = ((1+w) * stdev_fill_factor);
//    double balance = 1.0 / stdev_fill_factor;
    double balance = exp(-(1+w) * stdev_fill_factor);
    double ret = sum_seg + Const::alpha* balance;
    delete[] tmp;
    return ret;
}

// determine fan-out and choose segments
void FADASNode::determineSegmentsCluster(){
    int lambda_min, lambda_max;
    determineFanout(&lambda_min, &lambda_max);

    vector<int>unit_size(Const::vertexNum, 0);
    vector<vector<double>>unit_seg_sum(Const::vertexNum , vector<double>(Const::segmentNum, 0));
    vector<unordered_map<unsigned short, int>>data_seg_symbols(Const::segmentNum);

    for(int offset:offsets){
        unsigned short* cur_sax = saxes + offset * Const::segmentNum;
        int head = SaxUtil::extendSax(cur_sax, bits_cardinality);
        unit_size[head]++;
        for(int i=0;i<Const::segmentNum;++i) {
            unit_seg_sum[head][i] += SaxUtil::getMidLineFromSaxSymbolbc8(cur_sax[i]);
            data_seg_symbols[i][cur_sax[i]]++;
        }
    }
    vector<double>data_seg_mean(Const::segmentNum, 0);
    for(int i=0;i<Const::segmentNum;++i){
        auto& map = data_seg_symbols[i];
        for(auto &iter:map){
            unsigned short symbol = iter.first;
            data_seg_mean[i] += (SaxUtil::getMidLineFromSaxSymbolbc8(symbol) * iter.second);
        }
        data_seg_mean[i] /= size;
    }
    vector<unordered_map<unsigned short, int>>().swap(data_seg_symbols);

    // start to compute the size of each node in each plan
    int plan_num = FADASNode::combine_num[lambda_max];
    unordered_set<int>visited;
    double max_score = 0;
    vector<int> best_plan;
    for(int i=0;i<plan_num;++i){
        int *plan = FADASNode::combines[lambda_max][i];
        // first evaluate the whole plan
        vector<int>plan_node_sizes(1<<lambda_max, 0);
        vector<vector<double>>plan_node_seg_mean(1<<lambda_max, vector<double>(Const::segmentNum,0));
        vector<vector<double>>plan_node_seg_sum(1<<lambda_max, vector<double>(Const::segmentNum, 0));
        int mask_code = MathUtil::generateMaskSettingKbits(plan, lambda_max, Const::segmentNum);
        map<int,int>max_node_size;
        map<int,vector<double>>max_node_seg_sum;
        for(int j=0;j<Const::vertexNum;++j){
            max_node_size[mask_code & j] += unit_size[j];
            if(max_node_seg_sum[mask_code & j].empty()){
                max_node_seg_sum[mask_code & j].resize(Const::segmentNum, 0);
            }
            for(int k=0;k<Const::segmentNum;++k)
                max_node_seg_sum[mask_code & j][k] += unit_seg_sum[j][k];
        }
//        assert(max_node_size.size() == (1 << lambda_max));
        int _ = 0;
        for(auto & iter: max_node_size) plan_node_sizes[_++] = iter.second;
        map<int,int>().swap(max_node_size);
        _ = 0;
        for(auto & iter: max_node_seg_sum)  {
            vector<double>&tmp = iter.second;
            for(int k=0;k<Const::segmentNum;++k) {
                plan_node_seg_sum[_][k] = tmp[k];
                plan_node_seg_mean[_][k] = tmp[k] / plan_node_sizes[_];
            }
            ++_;
        }
        map<int,vector<double>>().swap(max_node_seg_sum);

        double score = compute_score_cluster(plan_node_sizes, plan_node_seg_mean, data_seg_mean,
                                             plan, lambda_max);
        if(score > max_score){
            max_score = score;
            best_plan.clear();
            for(int j = 0; j<lambda_max;++j)
                best_plan.push_back(plan[j]);
        }

        if(lambda_min <= lambda_max - 1)
            visitPlanFromBaseTableCluster(visited, lambda_max - 1, plan,
                                          plan_node_sizes, plan_node_seg_sum,
                                   &max_score, best_plan, lambda_min, mask_code, data_seg_mean);
        vector<int>().swap(plan_node_sizes);
        vector<vector<double>>().swap(plan_node_seg_mean);
        vector<vector<double>>().swap(plan_node_seg_sum);
    }

    unordered_set<int>().swap(visited);
    vector<int>().swap(unit_size);
    vector<vector<double>>().swap(unit_seg_sum);
    chosenSegments = best_plan;
}

void
FADASNode::visitPlanFromBaseTableCluster(unordered_set<int> &visited, int cur_lambda, const int *plan,
                                         vector<int> &base_tbl_size, vector<vector<double>>&base_tbl_seg_sum,
                                   double *max_score, vector<int> &best_plan, int lambda_min, int mask_code,
                                  vector<double> &data_seg_mean) {
    // base mask is used to detect base tbl
    int base_mask = 1;
    for(int i=0;i<cur_lambda;++i)
        base_mask = (base_mask << 1) +1;

    for(int i=0;i<cur_lambda + 1 ;++i){
        int reset_pos = plan[i];
        // get the whole plan code
        int cur_whole_mask = mask_code - (1 << (Const::segmentNum - 1  - reset_pos));
        if(visited.contains(cur_whole_mask))
            continue;
        visited.insert(cur_whole_mask);
        // get the new plan
        int *new_plan = new int[cur_lambda];
        int _=0;
        for(int j=0;j<cur_lambda+1;++j){
            if(i != j)  new_plan[_++] = plan[j];
        }

        // get the current base mask
        int cur_base_mask = base_mask - (1 << (cur_lambda - i));
        map<int,int>node_size_map;
        map<int,vector<double>>node_seg_sum_map;
        for(int j=0;j<base_tbl_size.size();++j) {
            int after_mask = cur_base_mask & j;
            node_size_map[after_mask] += base_tbl_size[j];
            if(node_seg_sum_map[after_mask].empty())
                node_seg_sum_map[after_mask].resize(Const::segmentNum, 0);
            for(int k=0;k<Const::segmentNum;++k)
                node_seg_sum_map[after_mask][k] += base_tbl_seg_sum[j][k];
        }
//        assert(node_size_map.size() == (1<<cur_lambda));
        vector<int>new_tbl(1<<cur_lambda, 0);
        vector<vector<double>>new_tbl_seg_sum(1<<cur_lambda, vector<double>(Const::segmentNum, 0));
        vector<vector<double>>new_tbl_seg_mean(1<<cur_lambda, vector<double>(Const::segmentNum, 0));

        _ = 0;
        for(auto &iter:node_size_map)
            new_tbl[_++] = iter.second;
        map<int,int>().swap(node_size_map);
        _ = 0;
        for(auto & iter: node_seg_sum_map)  {
            vector<double>&tmp = iter.second;
            for(int k=0;k<Const::segmentNum;++k) {
                new_tbl_seg_sum[_][k] = tmp[k];
                new_tbl_seg_mean[_][k] = tmp[k] / new_tbl[_];
            }
            ++_;
        }
        map<int,vector<double>>().swap(node_seg_sum_map);

        double score  = compute_score_cluster(new_tbl, new_tbl_seg_mean, data_seg_mean, new_plan, cur_lambda);
        if(score > *max_score){
            *max_score = score;
            best_plan.clear();
            for(_ = 0; _<cur_lambda;++_)
                best_plan.push_back(new_plan[_]);
        }
        if(lambda_min <= cur_lambda  -1)
            visitPlanFromBaseTableCluster(visited, cur_lambda - 1, new_plan, new_tbl, new_tbl_seg_sum,
                                   max_score, best_plan, lambda_min, cur_whole_mask, data_seg_mean);
        vector<int>().swap(new_tbl);
        vector<vector<double>>().swap(new_tbl_seg_mean);
        vector<vector<double>>().swap(new_tbl_seg_sum);
        delete[] new_plan;
    }
}

double FADASNode::compute_score_cluster(vector<int>&node_sizes, vector<vector<double>>&plan_node_seg_mean, vector<double>&seg_mean,
                                        int *plan, int lambda) const{
//    int over_th_nodes_no = 0;
//    for(int _:node_sizes){
//        if(_ > Const::th)
//            over_th_nodes_no++;
//    }
//    double w = ((double)over_th_nodes_no) / ((int)node_sizes.size());
    int cluster_num = 1 << lambda;
    double ssb = 0;
    for(int i=0;i<cluster_num;++i){
        double cluster_size_ratio = (double)node_sizes[i]/size;
        double dist =0;
        for(int k = 0;k<Const::segmentNum;++k){
            dist += ((plan_node_seg_mean[i][k] - seg_mean[k]) * (plan_node_seg_mean[i][k] - seg_mean[k]));
        }
        ssb += (cluster_size_ratio * sqrt(dist));
    }

    double *tmp = new double[node_sizes.size()];
    for(int i=0;i<node_sizes.size();++i){
        tmp[i] = ((double)node_sizes[i]) / Const::th;
    }
    double stdev_fill_factor = MathUtil::deviation(tmp, node_sizes.size());
    delete[] tmp;
//    double balance = ((1+w) * stdev_fill_factor);
//    double balance = sqrt(stdev_fill_factor);
    double balance = stdev_fill_factor;

    double ret = ssb / 1;
    return ret;
}

// determine fan-out and choose segments
void FADASNode::determineSegmentsWeakCluster(){
    int lambda_min, lambda_max;
    determineFanout(&lambda_min, &lambda_max);

    vector<int>unit_size(Const::vertexNum, 0);
    vector<vector<double>>unit_seg_sum(Const::vertexNum , vector<double>(Const::segmentNum, 0));
    vector<unordered_map<unsigned short, int>>data_seg_symbols(Const::segmentNum);

    for(int offset:offsets){
        unsigned short* cur_sax = saxes + offset * Const::segmentNum;
        int head = SaxUtil::extendSax(cur_sax, bits_cardinality);
        unit_size[head]++;
        for(int i=0;i<Const::segmentNum;++i) {
            unit_seg_sum[head][i] += SaxUtil::getMidLineFromSaxSymbolbc8(cur_sax[i]);
            data_seg_symbols[i][cur_sax[i]]++;
        }
    }
    vector<double>data_seg_mean(Const::segmentNum, 0);
    for(int i=0;i<Const::segmentNum;++i){
        auto& map = data_seg_symbols[i];
        for(auto &iter:map){
            unsigned short symbol = iter.first;
            data_seg_mean[i] += (SaxUtil::getMidLineFromSaxSymbolbc8(symbol) * iter.second);
        }
        data_seg_mean[i] /= size;
    }
    vector<unordered_map<unsigned short, int>>().swap(data_seg_symbols);

    // start to compute the size of each node in each plan
    int plan_num = FADASNode::combine_num[lambda_max];
    unordered_set<int>visited;
    double max_score = 0;
    vector<int> best_plan;
    for(int i=0;i<plan_num;++i){
        int *plan = FADASNode::combines[lambda_max][i];
        // first evaluate the whole plan
        vector<int>plan_node_sizes(1<<lambda_max, 0);
        vector<vector<double>>plan_node_seg_mean(1<<lambda_max, vector<double>(lambda_max,0));
        vector<vector<double>>plan_node_seg_sum(1<<lambda_max, vector<double>(lambda_max, 0));
        int mask_code = MathUtil::generateMaskSettingKbits(plan, lambda_max, Const::segmentNum);
        map<int,int>max_node_size;
        map<int,vector<double>>max_node_seg_sum;
        for(int j=0;j<Const::vertexNum;++j){
            max_node_size[mask_code & j] += unit_size[j];
            if(max_node_seg_sum[mask_code & j].empty()){
                max_node_seg_sum[mask_code & j].resize(lambda_max, 0);
            }
            for(int k=0;k<lambda_max;++k)
                max_node_seg_sum[mask_code & j][k] += unit_seg_sum[j][k];
        }
//        assert(max_node_size.size() == (1 << lambda_max));
        int _ = 0;
        for(auto & iter: max_node_size) plan_node_sizes[_++] = iter.second;
        map<int,int>().swap(max_node_size);
        _ = 0;
        for(auto & iter: max_node_seg_sum)  {
            vector<double>&tmp = iter.second;
            for(int k=0;k<lambda_max;++k) {
                plan_node_seg_sum[_][k] = tmp[k];
                plan_node_seg_mean[_][k] = tmp[k] / plan_node_sizes[_];
            }
            ++_;
        }
        map<int,vector<double>>().swap(max_node_seg_sum);

        double score = compute_score_cluster_weak(plan_node_sizes, plan_node_seg_mean, data_seg_mean,
                                             plan, lambda_max);
        if(score > max_score){
            max_score = score;
            best_plan.clear();
            for(int j = 0; j<lambda_max;++j)
                best_plan.push_back(plan[j]);
        }

        if(lambda_min <= lambda_max - 1)
            visitPlanFromBaseTableWeakCluster(visited, lambda_max - 1, plan,
                                          plan_node_sizes, plan_node_seg_sum,
                                          &max_score, best_plan, lambda_min, mask_code, data_seg_mean);
        vector<int>().swap(plan_node_sizes);
        vector<vector<double>>().swap(plan_node_seg_mean);
        vector<vector<double>>().swap(plan_node_seg_sum);
    }

    unordered_set<int>().swap(visited);
    vector<int>().swap(unit_size);
    vector<vector<double>>().swap(unit_seg_sum);
    chosenSegments = best_plan;
}

void
FADASNode::visitPlanFromBaseTableWeakCluster(unordered_set<int> &visited, int cur_lambda, const int *plan,
                                         vector<int> &base_tbl_size, vector<vector<double>>&base_tbl_seg_sum,
                                         double *max_score, vector<int> &best_plan, int lambda_min, int mask_code,
                                         vector<double> &data_seg_mean) {
    // base mask is used to detect base tbl
    int base_mask = 1;
    for(int i=0;i<cur_lambda;++i)
        base_mask = (base_mask << 1) +1;

    for(int i=0;i<cur_lambda + 1 ;++i){
        int reset_pos = plan[i];
        // get the whole plan code
        int cur_whole_mask = mask_code - (1 << (Const::segmentNum - 1  - reset_pos));
        if(visited.contains(cur_whole_mask))
            continue;
        visited.insert(cur_whole_mask);
        // get the new plan
        int *new_plan = new int[cur_lambda];
        int _=0;
        for(int j=0;j<cur_lambda+1;++j){
            if(i != j)  new_plan[_++] = plan[j];
        }

        // get the current base mask
        int cur_base_mask = base_mask - (1 << (cur_lambda - i));
        map<int,int>node_size_map;
        map<int,vector<double>>node_seg_sum_map;
        for(int j=0;j<base_tbl_size.size();++j) {
            int after_mask = cur_base_mask & j;
            node_size_map[after_mask] += base_tbl_size[j];
            if(node_seg_sum_map[after_mask].empty())
                node_seg_sum_map[after_mask].resize(cur_lambda, 0);
            _ = 0;
            for(int k=0;k<cur_lambda+1;++k) {
                if(i == k)  continue;
                node_seg_sum_map[after_mask][_++] += base_tbl_seg_sum[j][k];
            }
        }
//        assert(node_size_map.size() == (1<<cur_lambda));
        vector<int>new_tbl(1<<cur_lambda, 0);
        vector<vector<double>>new_tbl_seg_sum(1<<cur_lambda, vector<double>(cur_lambda, 0));
        vector<vector<double>>new_tbl_seg_mean(1<<cur_lambda, vector<double>(cur_lambda, 0));

        _ = 0;
        for(auto &iter:node_size_map)
            new_tbl[_++] = iter.second;
        map<int,int>().swap(node_size_map);
        _ = 0;
        for(auto & iter: node_seg_sum_map)  {
            vector<double>&tmp = iter.second;
            for(int k=0;k<cur_lambda;++k) {
                new_tbl_seg_sum[_][k] = tmp[k];
                new_tbl_seg_mean[_][k] = tmp[k] / new_tbl[_];
            }
            ++_;
        }
        map<int,vector<double>>().swap(node_seg_sum_map);

        double score  = compute_score_cluster_weak(new_tbl, new_tbl_seg_mean, data_seg_mean, new_plan, cur_lambda);
        if(score > *max_score){
            *max_score = score;
            best_plan.clear();
            for(_ = 0; _<cur_lambda;++_)
                best_plan.push_back(new_plan[_]);
        }
        if(lambda_min <= cur_lambda  -1)
            visitPlanFromBaseTableWeakCluster(visited, cur_lambda - 1, new_plan, new_tbl, new_tbl_seg_sum,
                                          max_score, best_plan, lambda_min, cur_whole_mask, data_seg_mean);
        vector<int>().swap(new_tbl);
        vector<vector<double>>().swap(new_tbl_seg_mean);
        vector<vector<double>>().swap(new_tbl_seg_sum);
        delete[] new_plan;
    }
}

double FADASNode::compute_score_cluster_weak(vector<int>&node_sizes, vector<vector<double>>&plan_node_seg_mean, vector<double>&seg_mean,
                                        int *plan, int lambda) const{
//    int over_th_nodes_no = 0;
//    for(int _:node_sizes){
//        if(_ > Const::th)
//            over_th_nodes_no++;
//    }
//    double w = ((double)over_th_nodes_no) / ((int)node_sizes.size());
    int cluster_num = 1 << lambda;
    double ssb = 0;
    for(int i=0;i<cluster_num;++i){
        double cluster_size_ratio = (double)node_sizes[i]/size;
        double dist =0;
        for(int k = 0;k<lambda;++k){
            dist += ((plan_node_seg_mean[i][k] - seg_mean[plan[k]]) * (plan_node_seg_mean[i][k] - seg_mean[plan[k]]));
        }
        ssb += (cluster_size_ratio * sqrt(dist));
    }

    double *tmp = new double[node_sizes.size()];
    for(int i=0;i<node_sizes.size();++i){
        tmp[i] = ((double)node_sizes[i]) / Const::th;
    }
    double stdev_fill_factor = MathUtil::deviation(tmp, node_sizes.size());
    delete[] tmp;
//    double balance = ((1+w) * stdev_fill_factor);
//    double balance = sqrt(stdev_fill_factor);
    double balance = stdev_fill_factor;

    double ret = ssb / 1;
    return ret;
}


int FADASNode::chooseOneSegment(PAA_INFO* node){
    int min = numeric_limits<int>::max(), min_index = -1;
    for(int i = 0; i < Const::segmentNum;++i){
        int big = max(node->paa_up_size[i], node->paa_below_size[i]);
        if(big < min){
            min = big;
            min_index = i;
        }
    }
    return min_index;
}

int FADASNode::chooseOneSegment(SAX_INFO* node){
    int min = numeric_limits<int>::max(), min_index = -1;
    for(int i = 0; i < Const::segmentNum;++i){
        int big = max(node->up_size[i], node->below_size[i]);
        if(big < min){
            min = big;
            min_index = i;
        }
    }
    return min_index;
}

void FADASNode::chooseSegment(PAA_INFO *paa, int chosen_num) {
    chosenSegments.resize(chosen_num);
    if(chosen_num == 1) { chosenSegments[0]=chooseOneSegment(paa) ; return;}

    tmp scores[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        if(bits_cardinality[i] >= Const::bitsCardinality)
            scores[i] = tmp(i, -1);
        else
            scores[i] = tmp(i, paa->paa_variance[i]);
    sort(scores, scores+Const::segmentNum, tmp::orderdesc);
    for(int i=0;i<chosen_num;++i)
        chosenSegments[i] = scores[i].i;
    sort(chosenSegments.begin(), chosenSegments.end());
}

void FADASNode::chooseSegment(SAX_INFO *sax_info, int chosen_num) {
    chosenSegments.resize(chosen_num);
    if(chosen_num == 1) { chosenSegments[0]=chooseOneSegment(sax_info) ; return;}

    tmp scores[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        if(bits_cardinality[i] >= Const::bitsCardinality)
            scores[i] = tmp(i, -1);
        else
            scores[i] = tmp(i, sax_info->sax_variance[i]);
    sort(scores, scores+Const::segmentNum, tmp::orderdesc);
    for(int i=0;i<chosen_num;++i)
        chosenSegments[i] = scores[i].i;
    sort(chosenSegments.begin(), chosenSegments.end());
}

void FADASNode::generateSaxAndCardIn1stLayer(int new_id){
    for(int i = Const::segmentNum - 1; i >=0 ;--i){
        unsigned short t = new_id % 2 ;
        new_id >>= 1;
        sax[i] =  t;
    }
}

void FADASNode::generateSaxAndCardinality(FADASNode* node, int new_id){
    copy(sax, sax + Const::segmentNum, node->sax);
    copy(bits_cardinality, bits_cardinality + Const::segmentNum, node->bits_cardinality);
    for(int i = chosenSegments.size() - 1; i >=0 ;--i){
        int seg = (chosenSegments)[i];
        node->bits_cardinality[seg]++;
        int t = new_id % 2 ;
        new_id >>= 1;
        node->sax[seg] = (node->sax[seg] << 1) + t;
    }
}

void FADASNode::generateSaxAndCardIn1stLayer4LeafNode(int new_id){
    if(bits_cardinality[0] == -1){
        for(int &i:bits_cardinality)    i=1;
        generateSaxAndCardIn1stLayer(new_id);
        return;
    }
    for(int i = Const::segmentNum - 1; i >=0 ;--i){
        int t = new_id % 2 ;
        new_id >>= 1;
        if(bits_cardinality[i] == 1 && sax[i] != t){
            bits_cardinality[i] = 0;
        }
    }
}

void FADASNode::generateSaxAndCardinality4LeafNode(FADASNode* node, int new_id){
    if(node->bits_cardinality[0] == -1){
        generateSaxAndCardinality(node, new_id);
        return;
    }
    for(int i = chosenSegments.size() - 1; i >=0 ;--i){
        int seg = chosenSegments[i];
        int t = new_id % 2 ;
        new_id >>= 1;
        if(node->bits_cardinality[seg] == bits_cardinality[seg] + 1 && node->sax[seg] % 2 != t){
            node->bits_cardinality[seg]--;
            node->sax[seg] >>= 1;
        }
    }
}

partUnit* find_node(vector<partUnit*>&nodes, int target_id){
    for(partUnit* node: nodes){
        if(node->id == target_id)
            return node;
    }
    return nullptr;
}

struct failUnit{
    partUnit *node{};
    int neighbor_size{};

    failUnit(partUnit* a, int b){ node = a; neighbor_size = b;}
};

static bool comp_fail_node(failUnit *x, failUnit *y){
    return x->neighbor_size > y->neighbor_size;
}

int FADASNode::partition1stLayer(partUnit *nodes_map, vector<vector<int>> *g, double filling_factor) {
    vector<partUnit*>nodes;
    int total_size = 0, node_number;
    for(int i=0;i<Const::vertexNum;++i){
        if(nodes_map[i].size <= Const::th && nodes_map[i].size > 0) {
            nodes.push_back(&nodes_map[i]); 
            total_size += nodes_map[i].size;
        }
    }

    node_number = nodes.size();
    if(node_number < 2) return 0;
    if(node_number == 2 && nodes[0]->size + nodes[1]->size > Const::th) return 0;
//    cout << "node number = " << node_number <<", total size = " << total_size << endl;
    sort(nodes.begin(),  nodes.end(), partUnit::comp_size);
    int pid = 0;
    int k = total_size / Const::th + 1;

    int _id, finish_num = 0, finish_size = 0, temp_size = 0, fail_partition_size = 0;
    partUnit *temp_node;
    vector<partUnit*>temp, candidates;
    // first loop, build 2-clique
    {
        for(int i=0;i<node_number;++i){
            if(nodes[i]->pid != -1 || nodes[i]->size > Const::th)   continue;
            temp_size = nodes[i]->size;  temp.clear(); candidates.clear();  temp.push_back(nodes[i]);
            _id = nodes[i]->id;
            for(int j=0; j<a; ++j)
            {
                temp_node = &nodes_map[(*g)[_id][j]];
                if(temp_node->size > Const::th)    continue;
                if(temp_node->pid == -1 && temp_node->size < Const::th)
                    candidates.push_back(temp_node);
            }

            sort(candidates.begin(),candidates.end(), partUnit::comp_size);
            for(partUnit* node:candidates){
                if(node->size + temp_size > Const::th) continue;
                temp.push_back(node);
                temp_size += node->size;
            }

            // fulfill the partition requirement
            if(temp_size >= filling_factor * Const::th){
                nodes[i]->pid = pid;
                for(partUnit* cur: temp)  cur->pid = pid; 
                finish_num += temp.size();
                finish_size += temp_size;
//                cout << pid << ":" << temp_size << endl;
                ++pid;
            }
        }
    }

    if(finish_num >= node_number)   return pid;
//    cout<< "After first loop, finish nodes number = " << finish_num <<", finish series number = " << finish_size <<endl
//    <<"Built " << pid <<" partitions, with avg size = " << accumulate(partition_size.begin(), partition_size.end(), 0.0) / (double)pid<<endl;

    vector<partUnit*>candidates2;
    // second loop, build 4-clique
    {
        for(int i=0;i<node_number;++i){
            if(nodes[i]->pid != -1 || nodes[i]->size > Const::th)   continue;
            _id = nodes[i]->id;
            temp_size = nodes[i]->size;  temp.clear(); temp.push_back(nodes[i]);  candidates2.clear();
            for(int j=0; j<a; ++j) {
                temp_node = &nodes_map[(*g)[_id][j]];
                if(temp_node->pid == -1 && temp_node->size + temp_size <= Const::th){
                    temp.push_back(temp_node);
                    temp_size += temp_node->size;
                }
            }
            for(int j=a;j<a+b;++j){
                temp_node = &nodes_map[(*g)[_id][j]];
                if(temp_node->pid == -1 && temp_node->size < Const::th){
                    candidates2.push_back(temp_node);
                }
            }

            sort(candidates2.begin(),candidates2.end(), partUnit::comp_size);
            for(partUnit* node:candidates2){
                if(node->size + temp_size > Const::th) continue;
                temp.push_back(node);
                temp_size += node->size;
            }

            // fulfill the partition requirement
            if(temp_size >= filling_factor * Const::th){
                nodes[i]->pid = pid;
                for(partUnit* cur: temp)  cur->pid = pid;
                finish_num += temp.size();
                finish_size += temp_size;
                ++pid;
            }
        }
    }

    if(finish_num >= node_number)   return pid;
//    cout<< "After second loop, finish nodes number = " << finish_num <<", finish series number = " << finish_size <<endl
//    <<"Built " << pid <<" partitions, with avg size = " << accumulate(partition_size.begin() + one_loop_pid, partition_size.end(), 0.0) / (double)(pid-one_loop_pid)<<endl;
    vector<failUnit*>fail_partition_node;

    {
        int no_part;
        for(int i=0;i<node_number;++i){
            if(nodes[i]->pid != -1)   continue;
            no_part = 0;
            _id = nodes[i]->id;
            for(int j=0;j<a+b+c;++j){
                temp_node = &nodes_map[(*g)[_id][j]];
                if(temp_node->size > Const::th)    continue;
                if(temp_node->pid == -1){
                    no_part+= temp_node->size;
                }
            }

            fail_partition_node.push_back(new failUnit(nodes[i], no_part + nodes[i]->size));
            fail_partition_size += nodes[i]->size;
        }
    }
//    cout<< "After filling stage, finish nodes number = " << finish_num <<", finish series number = " << finish_size <<endl
//    <<"Built " << pid <<" partitions, with avg size = " << accumulate(partition_size.begin(), partition_size.end(), 0.0) / (double)pid<<endl;
//    cout <<"Fail partition node number = " << fail_partition_node.size() << " , total size = " << fail_partition_size << endl;

    sort(fail_partition_node.begin(),  fail_partition_node.end(), comp_fail_node);
    partUnit* debug_node = &nodes_map[33];
    vector<bool>flag(Const::vertexNum, false);
    for(failUnit* node:fail_partition_node){
        if(flag[node->node->id])    continue;
        temp_size = node->node->size;   temp.clear();   temp.push_back(node->node);
        node->node->pid = pid;  flag[node->node->id] = true;
        finish_num++; finish_size+=node->node->size;
        _id = node->node->id;

        int j;
        for(j=0; j<a +b +c && temp_size < Const::th; ++j)
        {
            temp_node = &nodes_map[(*g)[_id][j]];
            if(temp_node->pid == -1&& temp_node->size + temp_size <= Const::th)
                temp_node->pid = pid, flag[temp_node->id] = true, finish_num++, finish_size += temp_node->size, temp_size+=temp_node->size, temp.push_back(temp_node);
        }

        if(temp_size >= filling_factor * Const::th) {
            ++pid; continue;
        }
        for(int i=4; i <= Const::max_radius&& temp_size < Const::th; ++i){
            int mask = (1U << i) - 1, r, last1;
            do{
                {
                    int cur = mask ^ node->node->id;
                    temp_node = &nodes_map[cur];
                    if(temp_node->pid == -1 && temp_node->size + temp_size < Const::th)
                        temp_node->pid = pid, flag[temp_node->id] = true, finish_num++, finish_size += temp_node->size, temp_size+=temp_node->size, temp.push_back(temp_node);
                    if(temp_size >= Const::th) break;
                }
                last1 = mask & -mask;
                r = mask + last1;
                mask = (((r ^ mask) >> 2) / last1) | r;
            } while (r < (1 << Const::segmentNum));
        }
        ++pid;
    }

//    cout<< "After stretch stage, finish nodes number = " << finish_num <<", finish series number = " << finish_size <<endl
//    <<"Built " << pid <<" partitions, with avg size = " << accumulate(partition_size.begin(), partition_size.end(), 0.0) / (double)pid<<endl;

    for(auto* _:fail_partition_node)    delete _;
    return pid;
}

int FADASNode::partitionLessPack(partUnit* nodes_map, int chosen_segment_number){
    unordered_set<partUnit*>nodes;
    int input_node_number = 1 << chosen_segment_number;
    int total_size = 0, node_number;
    for(int i=0;i<input_node_number;++i){
        if(nodes_map[i].size <= Const::th * Const::small_perc && nodes_map[i].size > 0) {
            nodes.insert(&nodes_map[i]);
            total_size += nodes_map[i].size;
        }
    }

    node_number = nodes.size();
    if(node_number <= 5) return 0;
//    cout << "node number = " << node_number <<", total size = " << total_size << endl;
//    sort(nodes.begin(),  nodes.end(), partUnit::comp_size);
    int pid = 0;
    for(partUnit* cur_node:nodes){
        if(cur_node->pid!=-1)   continue;   // the node has been packed
        int cur_id = cur_node->id;
        int _mask = 1;
        int cur_size = cur_node->size;
        unordered_set<partUnit*>pack;
        pack.insert(cur_node);
        for(int j=1;j<chosen_segment_number;++j){
            int cand_id = cur_id ^ _mask;
            if(nodes_map[cand_id].size <= Const::th * Const::small_perc && nodes_map[cand_id].size > 0
            && nodes_map[cand_id].pid==-1
            && cur_size + nodes_map[cand_id].size <= Const::th){
                pack.insert(&nodes_map[cand_id]);
                cur_size += nodes_map[cand_id].size;
            }
            _mask <<= 1;
        }
        if(pack.size() > 1){
            for(partUnit* cand:pack){
                cand->pid = pid;
            }
            ++pid;
        }
    }

    return pid;
}

struct pack{
    vector<bool> cur_bits;
    vector<bool> cur_mask;
    int total_size;
    int pid{};
    int masked_bits_num;
    bool disabled;

    pack(){
        total_size = 0;
        masked_bits_num = 0;
        disabled = false;
    }
    pack(partUnit* node, int chosen_segment_num, int _pid){
        total_size = node->size;
        masked_bits_num = 0;
        cur_bits.resize(chosen_segment_num, false);
        cur_mask.resize(chosen_segment_num, false);
        int _id =  node->id;
        for(int i=0;i<chosen_segment_num;++i){
            cur_bits[chosen_segment_num - 1- i] = _id % 2;
            _id >>= 1;
        }
        pid = _pid;
        node->pid = pid;
        disabled = false;
    }

    int calc_cost(int _id, int chosen_seg_num){
        int cost = 0;
        for(int i=0;i<chosen_seg_num;++i){
            if(!cur_mask[chosen_seg_num-1-i] && cur_bits[chosen_seg_num - 1 -i]!= (_id %2)){
                ++cost;
            }
            _id >>=1;
        }
        return cost;
    }

    int calc_pack_merge_cost(const pack & p, int chosen_seg_num, int *cur_cost, int *tar_cost){
        *cur_cost = 0; *tar_cost = 0;
        int cost = 0;
        for(int i=0;i<chosen_seg_num;++i){
            if(cur_mask[i] && p.cur_mask[i])    continue;
            if(cur_mask[i] && !p.cur_mask[i]){
                (*tar_cost)++;
                ++cost;
            }
            else if (!cur_mask[i] && p.cur_mask[i]) { (*cur_cost)++; ++cost;}
            else if(cur_bits[i] != p.cur_bits[i])   {
                (*cur_cost)++;
                (*tar_cost)++;
                ++cost;
            }
        }
        return cost;
    }

    void merge_pack(pack* p, int chosen_seg_num){
        pack * dis_one, * res_one;
        if(pid < p->pid){
            dis_one = p;
            res_one = this;
        }else{
            dis_one = this;
            res_one = p;
        }
        dis_one->disabled = true;
        res_one->total_size = total_size + p->total_size;
        for(int i=0;i<chosen_seg_num;++i){
            if(dis_one->cur_mask[i] && res_one->cur_mask[i])    continue;
            if(dis_one->cur_mask[i] && !res_one->cur_mask[i]) { res_one->cur_mask[i] = true; res_one->masked_bits_num++;}
            else if(!dis_one->cur_mask[i] && res_one->cur_mask[i]) continue;
            else if(dis_one->cur_bits[i] != res_one->cur_bits[i]){
                res_one->masked_bits_num++;
                res_one->cur_mask[i] = true;
            }
        }

    }

    void insert(partUnit* node, int chosen_seg_num){
        node->pid = pid;
        int _id = node->id;
        total_size += node->size;
        for(int i=0;i<chosen_seg_num;++i){
            if(!cur_mask[chosen_seg_num-1-i] && cur_bits[chosen_seg_num - 1 -i]!= (_id %2)){
                cur_mask[chosen_seg_num-1-i] = true;
                masked_bits_num++;
            }
            _id>>=1;
        }
    }

    static bool comp_size(const pack&x, const pack&y){
        return x.total_size < y.total_size;
    }
};

int FADASNode::partitionNew(vector<partUnit>& nodes_map, int chosen_segment_number){
    vector<partUnit*>nodes;
    int input_node_number = 1 << chosen_segment_number;
    int total_size = 0, node_number;
    for(int i=0;i<input_node_number;++i){
        if(nodes_map[i].size < Const::th * Const::small_perc && nodes_map[i].size > 0) {
            nodes.push_back(&nodes_map[i]);
            total_size += nodes_map[i].size;
        }
    }

    node_number = nodes.size();
    if(node_number <= 3) return 0;
//    cout << "node number = " << node_number <<", total size = " << total_size << endl;
//    sort(nodes.begin(),  nodes.end(), partUnit::comp_size);
    int pid = 0;
    int first_round_pack_num = (int)floor(total_size / (Const::th));
    if(node_number <= first_round_pack_num) return 0;
    int max_mask_num = floor(chosen_segment_number * Const::max_mask_bit_percentage);
    // 0.5 is a small number, generate more packs in the first round
    vector<pack>packs(first_round_pack_num);
    sort(nodes.begin(), nodes.end(), partUnit::comp_size);
    for(int i=0;i<first_round_pack_num;++i){
        packs[i] = pack(nodes[i], chosen_segment_number, i);
    }

    for(partUnit* cur_node:nodes) {
        if (cur_node->pid != -1) continue;   // the node has been packed
        int cur_id = cur_node->id;

        int min_cost = chosen_segment_number;
        int pack_id = -1;
        for (pack &p: packs) {
            if (p.total_size + cur_node->size > Const::th || p.masked_bits_num >= max_mask_num) continue;
            int cost = p.calc_cost(cur_id, chosen_segment_number);
            if(cost + p.masked_bits_num >= max_mask_num)    continue;
            if (cost < min_cost) {
                min_cost = cost;
                pack_id = p.pid;
            }
        }
        if (pack_id == -1) {
            packs.emplace_back(cur_node, chosen_segment_number, packs.size());
        }else {
            packs[pack_id].insert(cur_node, chosen_segment_number);
        }
    }

    // merge packs
    unordered_map<int,int>pid_map;
    for(int i = 0;i<packs.size();++i)
        pid_map[i] = i;
    sort(packs.begin(),packs.end(), pack::comp_size);
    for(int i=0;i<packs.size();++i){
        pack& cur_pack = packs[i];
        if(cur_pack.pid != pid_map[cur_pack.pid])    continue;

        int min_cost = chosen_segment_number;
        int min_size = numeric_limits<int>::max();
        int min_pack_id = -1;
        int cur_cost, tar_cost;
        for(int j=0;j<packs.size();++j){
            pack&target_pack = packs[j];
            if(target_pack.disabled || cur_pack.pid == target_pack.pid || cur_pack.total_size + target_pack.total_size > Const::th
            || cur_pack.masked_bits_num >= max_mask_num || target_pack.masked_bits_num >= max_mask_num) continue;
            cur_pack.calc_pack_merge_cost(target_pack, chosen_segment_number, &cur_cost, &tar_cost);
            if(cur_cost + cur_pack.masked_bits_num >= max_mask_num ||
               tar_cost + target_pack.masked_bits_num >= max_mask_num) continue;
            int cost = cur_cost + tar_cost;
            if(cost < min_cost || (cost == min_cost && cur_pack.total_size < min_size)){
                min_cost = cost;
                min_size = target_pack.total_size;
                min_pack_id = j;
            }
        }

        if(min_size < numeric_limits<int>::max()){
            cur_pack.merge_pack(&packs[min_pack_id], chosen_segment_number);
        }
    }

    // re-assign the pids to the nodes
    int max_pid = 0;
    for(partUnit* node:nodes) {
        node->pid = pid_map[node->pid];
        max_pid = max(max_pid, node->pid);
    }

    return max_pid + 1;
}

int FADASNode::partitionNew(partUnit* nodes_map, int chosen_segment_number){
    vector<partUnit*>nodes;
    int input_node_number = 1 << chosen_segment_number;
    int total_size = 0, node_number;
    for(int i=0;i<input_node_number;++i){
        if(nodes_map[i].size < Const::th * Const::small_perc && nodes_map[i].size > 0) {
            nodes.push_back(&nodes_map[i]);
            total_size += nodes_map[i].size;
        }
    }

    node_number = nodes.size();
    if(node_number <= 3) return 0;
//    cout << "node number = " << node_number <<", total size = " << total_size << endl;
//    sort(nodes.begin(),  nodes.end(), partUnit::comp_size);
    int pid = 0;
    int first_round_pack_num = (int)floor(total_size / (Const::th));
    if(node_number <= first_round_pack_num) return 0;
    int max_mask_num = floor(chosen_segment_number * Const::max_mask_bit_percentage);
    // 0.5 is a small number, generate more packs in the first round
    vector<pack>packs(first_round_pack_num);
    sort(nodes.begin(), nodes.end(), partUnit::comp_size);
    for(int i=0;i<first_round_pack_num;++i){
        packs[i] = pack(nodes[i], chosen_segment_number, i);
    }

    for(partUnit* cur_node:nodes) {
        if (cur_node->pid != -1) continue;   // the node has been packed
        int cur_id = cur_node->id;

        int min_cost = chosen_segment_number;
        int pack_id = -1;
        for (pack &p: packs) {
            if (p.total_size + cur_node->size > Const::th || p.masked_bits_num >= max_mask_num) continue;
            int cost = p.calc_cost(cur_id, chosen_segment_number);
            if(cost + p.masked_bits_num >= max_mask_num)    continue;
            if (cost < min_cost) {
                min_cost = cost;
                pack_id = p.pid;
            }
        }
        if (pack_id == -1) {
            packs.emplace_back(cur_node, chosen_segment_number, packs.size());
        }else {
            packs[pack_id].insert(cur_node, chosen_segment_number);
        }
    }

    // merge packs
    unordered_map<int,int>pid_map;
    for(int i = 0;i<packs.size();++i)
        pid_map[i] = i;
    sort(packs.begin(),packs.end(), pack::comp_size);
    for(int i=0;i<packs.size();++i){
        pack& cur_pack = packs[i];
        if(cur_pack.pid != pid_map[cur_pack.pid])    continue;

        int min_cost = chosen_segment_number;
        int min_size = numeric_limits<int>::max();
        int min_pack_id = -1;
        int cur_cost, tar_cost;
        for(int j=0;j<packs.size();++j){
            pack&target_pack = packs[j];
            if(target_pack.disabled || cur_pack.pid == target_pack.pid || cur_pack.total_size + target_pack.total_size > Const::th
               || cur_pack.masked_bits_num >= max_mask_num || target_pack.masked_bits_num >= max_mask_num) continue;
            cur_pack.calc_pack_merge_cost(target_pack, chosen_segment_number, &cur_cost, &tar_cost);
            if(cur_cost + cur_pack.masked_bits_num >= max_mask_num ||
               tar_cost + target_pack.masked_bits_num >= max_mask_num) continue;
            int cost = cur_cost + tar_cost;
            if(cost < min_cost || (cost == min_cost && cur_pack.total_size < min_size)){
                min_cost = cost;
                min_size = target_pack.total_size;
                min_pack_id = j;
            }
        }

        if(min_size < numeric_limits<int>::max()){
            cur_pack.merge_pack(&packs[min_pack_id], chosen_segment_number);
        }
    }

    // re-assign the pids to the nodes
    int max_pid = 0;
    for(partUnit* node:nodes) {
        node->pid = pid_map[node->pid];
        max_pid = max(max_pid, node->pid);
    }

    return max_pid + 1;
}


int findFirstLE(vector<partUnit*>&nodes, int start, int end, int target){
    int mid, nn = end + 1;
    while (start < end){
        mid = (start + end) / 2;
        if(nodes[mid]->size > target){
            start = mid + 1;
        }else{
            end = mid;
        }
    }
    if(start == end && nodes[start]->size <= target)    return  start;
    return  nn;
}

int FADASNode::partition(partUnit* nodes_map, int chosen_num){
    vector<partUnit*>nodes;
    int input_node_num = 1 << chosen_num;
    int total_size = 0, node_number;
    for(int i=0;i<input_node_num;++i){
        if(nodes_map[i].size <= Const::th) {
            nodes.push_back(&nodes_map[i]);
            total_size += nodes_map[i].size;
        }
    }

    node_number = nodes.size();
    if(node_number < 2) return 0;
    if(node_number == 2 && nodes[0]->size + nodes[1]->size > Const::th) return 0;
//    cout << "node number = " << node_number <<", total size = " << total_size << endl;
    sort(nodes.begin(),  nodes.end(), partUnit::comp_size);
    int pid = 0;
    for(int i=0;i<node_number;++i){
        if(nodes[i]->pid != -1) continue;
        int target = Const::th - nodes[i]->size;
        int start = findFirstLE(nodes, i + 1, node_number - 1, target);
        nodes[i]->pid = pid;
        for(int j=start; j < node_number; ++j){
            if(nodes[j]->pid != -1 || nodes[j]->size > target) continue;
            nodes[j]->pid = pid;
            target -= nodes[j]->size;
        }
        ++pid;
    }
    return pid;
}

void FADASNode::getIndexStats(){
    int total_leaf_node_num = getLeafNodeNum();
    int total_size = getTotalSize();
    cout << "Total size = " << total_size << endl;
    cout <<"Total nodes number = " << getNodeNum() << endl;
    cout << "Leaf node number = " << total_leaf_node_num << endl;
    cout << "1st layer node number = " << get1stLayerNodesNo() <<endl;
    cout << "1st layer internal node number = " << get1stLayerInterNodesNo() << endl;
    cout << "1st layer internal series number = " << get1stLayerInterNodeSeriesNo() << endl;
    cout << "Max. height = " << getMaxHeight() - 1 <<endl;
    cout << "Avg. Height = " << getSumHeight() / (double) total_leaf_node_num << endl;
    cout <<"Avg. Filling Factor = "<< total_size / (double)total_leaf_node_num / Const::th << endl;
    cout << "Bias leaf node ratio = " << (double)getBiasLeafNodeNum() / total_leaf_node_num << endl;
    double sum = 0, sum_square  = 0, sum_dist = 0;
    int deep_leaf_num = 0;
    getBoundRange(&sum, &sum_square, &deep_leaf_num, &sum_dist);
    cout << "deep layer leaf number = " << deep_leaf_num <<endl;
    cout << "total sum = " << sum <<", avg sum = " << sum / deep_leaf_num << endl;
    cout << "total sum_square = " << sum_square << ", avg sum square = " << sum_square / deep_leaf_num << endl;
    cout << "total sum_distance = " << sum_dist << ", avg sum dist = " << sum_dist / deep_leaf_num << endl;
}

int FADASNode::get1stLayerInterNodesNo(){
    unordered_set<FADASNode*>node;
    for(FADASNode* child:children){
        if(child== nullptr || child->size <= Const::th || node.find(child) != node.end())
            continue;
        node.insert(child);
    }
    return node.size();
}

int FADASNode::get1stLayerNodesNo(){
    unordered_set<FADASNode*>node;
    for(FADASNode* child:children){
        if(child== nullptr || node.find(child) != node.end())
            continue;
        node.insert(child);
    }
    return node.size();
}

int FADASNode::get1stLayerInterNodeSeriesNo(){
    unordered_set<FADASNode*>node;
    int ret = 0;
    for(FADASNode* child:children){
        if(child== nullptr || child->size <= Const::th || node.find(child) != node.end())
            continue;
        node.insert(child);
        ret += child->size;
    }
    return ret;
}

int FADASNode::getMaxHeight(){
    if(!isInternalNode())    return 1;
    int max_height = 0;
    unordered_set<FADASNode*>hash_map;
    for(FADASNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        max_height = max(child->getMaxHeight(), max_height);
        hash_map.insert(child);
    }
    return max_height + 1;
}

int FADASNode::getLeafNodeNum(){
    if(!isInternalNode())    return 1;
    int sum = 0;
    unordered_set<FADASNode*>hash_map;
    for(FADASNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum += child->getLeafNodeNum();
        hash_map.insert(child);
    }
    return sum;
}

void FADASNode::getBoundRange(double *sum, double *sum_square, int *leafNum, double *sum_dist) {
    int len = 1 << chosenSegments.size();
    unordered_set<FADASNode*>visited;

    for(int i=0;i<len;++i){
        if(children[i] == nullptr)    continue;
        if(visited.count(children[i]) > 0)  continue;
        if(children[i]->isInternalNode()){
            children[i]->getBoundRange(sum, sum_square, leafNum, sum_dist);
        }else{
            auto child = children[i];
            visited.insert(child);
            *leafNum += 1;

            double dist = 0;
            int bc[Const::segmentNum];
            unsigned short s[Const::segmentNum];
            SaxUtil::extendSax(sax, bits_cardinality, chosenSegments, i, s,bc);
            for(int j=0; j < Const::segmentNum; ++j){
                double lb, ub;
                SaxUtil::getValueRange(s[j], bc[j], &lb, &ub);
                if(lb == -numeric_limits<double>::max())    lb = -4;
                if(ub == numeric_limits<double>::max()) ub = 4;
                *sum += (ub - lb);
                *sum_square += ((ub-lb) * (ub -lb));
                dist += Const::tsLengthPerSegment * ((ub-lb) * (ub -lb));
            }
            *sum_dist += sqrt(dist);
            cout << sqrt(dist) << ",";
        }
    }
}

// isax bias leaf nodes number
int FADASNode::getBiasLeafNodeNum(){
    if(!isInternalNode()){
        int max_b = 0, min_b=8;
        for(int &bc:bits_cardinality){
            max_b = max(max_b, bc);
            min_b = min(min_b, bc);
        }
        return (max_b - min_b >= 4);
    }
    int sum = 0;
    unordered_set<FADASNode*>hash_map;
    for(FADASNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum += child->getBiasLeafNodeNum();
        hash_map.insert(child);
    }
    return sum;
}

int FADASNode::getTotalSize(){
    if(!isInternalNode())    return size;
    int sum = 0;
    unordered_set<FADASNode*>hash_map;
    for(FADASNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum += child->getTotalSize();
        hash_map.insert(child);
    }
    return sum;
}

int FADASNode::getNodeNum(){
    if(!isInternalNode())    return 1;
    int sum = 0;
    unordered_set<FADASNode*>hash_map;
    for(FADASNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum += child->getNodeNum();
        hash_map.insert(child);
    }
    return sum + 1;
}

int FADASNode::getSumHeight(){
    if(!isInternalNode())    return layer;
    int sum_height = 0;
    unordered_set<FADASNode*>hash_map;
    for(FADASNode*child:children){
        if(child == nullptr || hash_map.find(child) != hash_map.end())    continue;
        sum_height += child->getSumHeight();
        hash_map.insert(child);
    }
    return sum_height;
}

int FADASNode::loadSax(const string & saxfn){
    long f_size = FileUtil::getFileSize(saxfn.c_str()), series_num = f_size / (sizeof(unsigned short) * Const::segmentNum);
    saxes = new unsigned short [f_size / sizeof(unsigned short )];
    FILE *f = fopen(saxfn.c_str(), "rb");
    fread(saxes, sizeof(unsigned short ), f_size / sizeof(unsigned short), f);
    fclose(f);
    Const::logPrint("Finish loading sax");
    return series_num;
}

void FADASNode::loadPaa(const string & paafn){
    long f_size = FileUtil::getFileSize(paafn.c_str());
    paas = new float [f_size / sizeof(float )];
    FILE *f = fopen(paafn.c_str(), "rb");
    fread(paas, sizeof(float ), f_size / sizeof(float ), f);
    fclose(f);
    Const::logPrint( "Finish loading paa");
}

long FADASNode::generateSaxAndPaaTbl(){
    string fn = Const::datafn;
    long fs = FileUtil::getFileSize(fn.c_str());
    long series_num = fs / Const::tsLengthBytes;
    cout << "Total Series Number is "<<series_num <<endl;
    float ts[Const::tsLength];
    saxes = new unsigned short[series_num * Const::segmentNum];
    paas = new float[series_num * Const::segmentNum];
    long rest = series_num, cur = 0;
    FILE *f = fopen(fn.c_str(), "rb");

    RAND_READ_CNT++;
    SEQ_READ_CNT += series_num;
    RAND_WRITE_CNT+=2;
    SEQ_WRITE_CNT += series_num;

    while(rest > 0){
        int num;
        if(rest > 400000)    num = 400000;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto start = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto end = chrono::system_clock::now();
        SAX_PAA_READ_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();

        for(int i=0;i<num;++i){
            if(isnan(tss[i * Const::tsLength])){
                for(int j = 0; j < Const::segmentNum; ++j)
                    saxes[i * Const::segmentNum + j] = 0;
                cout << "Dirty data: "<<i << "," <<endl;
            }
            else{
                SaxUtil::paaAndSaxFromTs(tss + i * Const::tsLength,
                                         paas + (cur+ i) * Const::segmentNum, saxes + (cur+ i) * Const::segmentNum,
                                         Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
            }
        }
        delete[] tss;
        rest -=num;
        cur+=num;
        start = chrono::system_clock::now();
        SAX_PAA_CPU_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();
    }

    fclose(f);
    return series_num;
}

void FADASNode::getFileNameInsert(const string &index_dir, string &sax_file, string &data_file) const {
    sax_file = index_dir;
    data_file = index_dir;
    if(layer == 1){
        sax_file += "1_";
        data_file += "1_";
        if(partition_id == -1){
            sax_file += to_string(id) + "_sax_L";
            data_file += to_string(id) + "_L";
        }else{
            sax_file += to_string(partition_id) + "_sax";
            data_file += to_string(partition_id);
        }
    }else{
        sax_file += getFileName() + "_sax";
        if(partition_id == -1) sax_file += "_L";
        data_file += getFileName();
        if(partition_id == -1) data_file += "_L";
    }
}

long FADASNode::generateSaxTbl(){
    string fn = Const::datafn;
    long series_num;
    if(Const::series_num == -1) {
        long fs = FileUtil::getFileSize(fn.c_str());
        series_num = fs / Const::tsLengthBytes;
    }else{
        series_num = Const::series_num;
    }
    cout << "Total Series Number is "<<series_num <<endl;
    float ts[Const::tsLength];
    saxes = new unsigned short[series_num * Const::segmentNum];
    long rest = series_num, cur = 0;
    FILE *f = fopen(fn.c_str(), "rb");

    RAND_READ_CNT++;
    SEQ_READ_CNT += series_num;
    RAND_WRITE_CNT+=2;
    SEQ_WRITE_CNT += series_num;

    while(rest > 0){
        unsigned num;
        if(rest > 2000000)    num = 2000000;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto start = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto end = chrono::system_clock::now();
        SAX_PAA_READ_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();

        for(int i=0;i<num;++i){
            if(isnan(tss[i * Const::tsLength])){
                for(int j = 0; j < Const::segmentNum; ++j)
                    saxes[i * Const::segmentNum + j] = 0;
                cout << "Dirty data: "<<i << "," <<endl;
            }
            else{
                SaxUtil::saxFromTs(tss + i * Const::tsLength, saxes + (cur+ i) * Const::segmentNum,
                                         Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
            }
        }
        delete[] tss;
        rest -=num;
        cur+=num;
        start = chrono::system_clock::now();
        SAX_PAA_CPU_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();
    }

    fclose(f);
    return series_num;
}

void FADASNode::generateSaxTbl(const float*tss, int series_num){
    saxes = new unsigned short[series_num * Const::segmentNum];

    for(long i = 0; i< series_num;++i)
        SaxUtil::saxFromTs(tss + i * Const::tsLength, saxes + i * Const::segmentNum,
                           Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);

}

FADASNode *FADASNode::loadFromDisk(const string &saxfn, const string &idxfn, bool need_sax) {
    if(need_sax)
        loadSax(saxfn);
    ifstream ifs(idxfn, ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    auto *g = new FADASNode();
    ia >> (*g);
    ifs.close();
    return g;
}

int FADASNode::assignLeafNum() {
    if(!isInternalNode()) {
        leaf_num = 1;
        return 1;
    }

    unordered_set<FADASNode*>visited;
    for(FADASNode* child: children){
        if(child == nullptr || visited.count(child) > 0)    continue;
        visited.insert(child);
        leaf_num += child->assignLeafNum();
    }

    return leaf_num;
}
