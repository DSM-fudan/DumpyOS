//
// Created by wzy on 2021/9/24.
//

#include "../../include/Tardis/TardisTreeNode.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Utils/TimeSeriesUtil.h"
#include "../../include/Const.h"
#include "../../include/DataStructures/FADASNode.h"
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <chrono>
#include <cmath>
#include <iostream>

unsigned short * TardisTreeNode::saxes = nullptr;
float *TardisTreeNode::dataset = nullptr;
int TardisTreeNode::a = MathUtil::nChooseK(Const::segmentNum, 1), TardisTreeNode::b = MathUtil::nChooseK(Const::segmentNum, 2), TardisTreeNode::c = MathUtil::nChooseK(Const::segmentNum, 3);

static long MAT1_CPU_TIME_STAT = 0, MAT1_CPU_TIME_COPY = 0, MAT1_WRITE_TIME = 0, MAT1_TOTAL_TIME = 0, MAT1_READ_TIME = 0,
        MAT2_CPU_TIME = 0, MAT2_WRITE_TIME = 0, MAT2_TOTAL_TIME = 0, SAX_PAA_TOTAL_TIME = 0, MAT2_READ_TIME = 0,
        GROW_CPU_TIME = 0,  GROW_CPU_TIME_1st = 0, GROW_TOTAL_TIME = 0,
        RAND_READ_CNT = 0, RAND_WRITE_CNT = 0, SEQ_READ_CNT = 0, SEQ_WRITE_CNT = 0, SAX_PAA_READ_TIME, SAX_PAA_CPU_TIME;


void TardisTreeNode::materialize1stLayer(const string & index_dir) const{
    auto start_t = chrono::system_clock::now();
    Const::logPrint("Start move data to disk file in 1st layer.");
    FILE *f = fopen(Const::datafn.c_str(), "r");
    long rest = size, total = size, cur = 0;
    unordered_map<int, LBL_UNIT>fbl_inter, fbl_leaves;

    RAND_READ_CNT++;

    // There is another implementation method that fbl in each node stores a pointer vector where each pointer points to a series.
    // This may incur many write calls.
    while(rest > 0){
        fbl_inter.clear();
        fbl_leaves.clear();
        long num;
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto end = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto start = chrono::system_clock::now();
        MAT1_READ_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();

        // statistic the size of each node fbl size, and allocate memory for them
        for(long i=cur;i<cur+num;++i){
            int _id = SaxUtil::invSaxHeadFromSax(saxes + i * Const::segmentNum, Const::bitsCardinality, Const::segmentNum);
            TardisTreeNode* target = (*children)[_id];
            if(target->isLeafNode()){
                fbl_leaves[target->partitionId].buffer.push_back(tss + (i-cur) * Const::tsLength);
            }else{
                fbl_inter[_id].buffer.push_back(tss + (i-cur) * Const::tsLength);
            }
        }

        end = chrono::system_clock::now();
        MAT1_CPU_TIME_STAT += chrono::duration_cast<chrono::microseconds>(end - start).count();

        // write series in order to node file from node fbl
        for(auto &iter:fbl_inter){
            string outfile = Const::tardisfn + "U_" + to_string(iter.first);
            FILE *outf = fopen(outfile.c_str(), "a");

            RAND_WRITE_CNT++;

            for(float *dat:iter.second.buffer)
                fwrite(dat, sizeof(float), Const::tsLength, outf);
            fclose(outf);
        }
        for(auto &iter:fbl_leaves){
            string outfile = Const::tardisfn + "1_" + to_string(iter.first);
            FILE *outf = fopen(outfile.c_str(), "a");
            RAND_WRITE_CNT++;

            for(float *dat:iter.second.buffer)
                fwrite(dat, sizeof(float), Const::tsLength, outf);
            fclose(outf);
        }

        start = chrono::system_clock::now();
        MAT1_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();
        delete[] tss;

        rest-=num;
        cur += num;
        Const::logPrint("Now in 1st layer " + to_string((double)cur / (double)total * 100) + "% series have been written to disk.");
    }

    fclose(f);
    auto end_t = chrono::system_clock::now();
    MAT1_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(end_t - start_t).count();
}

void TardisTreeNode::materializeInterNode(const string & index_dir){
    auto start_t = chrono::system_clock::now();

    FILE *f = fopen((Const::tardisfn + "U_" + to_string(id)).c_str(), "r");
    long rest = size, cur = 0, num;
    unordered_map<string, LBL_UNIT>lbl;

    RAND_READ_CNT++;

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
            unsigned short *_sax = saxes + (long)offsets[i] * Const::segmentNum;
            TardisTreeNode* target = route(_sax);
            assert(target->size > 0);
            while(!target->isLeafNode()){
                assert(target->size > 0);
                target= target->route(_sax);
            }
            lbl[target->getFileName()].buffer.push_back(tss + (i - cur) * Const::tsLength);
        }

        end = chrono::system_clock::now();
        MAT2_CPU_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();

        for(auto &iter:lbl){
            string outfile = Const::tardisfn + iter.first;
            FILE *outf = fopen(outfile.c_str(), "a");
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
    FileUtil::FileRemove((Const::tardisfn + "U_" + to_string(id)).c_str());
    auto end_t = chrono::system_clock::now();
    MAT2_TOTAL_TIME += chrono::duration_cast<chrono::microseconds>(end_t - start_t).count();
}

void TardisTreeNode::materialize(const string &index_dir) const {
    materialize1stLayer(Const::tardisfn);
    for(auto & iter:*children){
        if(iter.second != nullptr && iter.second->size > Const::th && !iter.second->isLeafNode())
            iter.second->materializeInterNode(Const::tardisfn);
    }

    cout << "During the process of materializing 1st layer nodes, total time is "<< MAT1_TOTAL_TIME / 1000
         <<"ms, I/O read time is "<< MAT1_READ_TIME / 1000 <<"ms, CPU statistic time  is " << MAT1_CPU_TIME_STAT / 1000
         << "ms, CPU copy Time is " << MAT1_CPU_TIME_COPY / 1000 << "ms, while I/O write time  is " << MAT1_WRITE_TIME / 1000 << "ms. "<< endl;

    cout << "During the process of materializing internal nodes, total time is "<<MAT2_TOTAL_TIME / 1000
         <<"ms, I/O read time is "<<MAT2_READ_TIME / 1000 <<"ms, CPU time is " << MAT2_CPU_TIME / 1000
         << "ms, while I/O write time is " << MAT2_WRITE_TIME / 1000 << "ms." << endl;

    cout << "Random read count = " << RAND_READ_CNT << endl
         << "Random write count = " << RAND_WRITE_CNT << endl;
}

TardisTreeNode * TardisTreeNode::BuildTardisTree(string &saxfn, bool on_disk) {
    long series_num;
    if(on_disk){
//        series_num = generateSaxTbl();;
        long f_size = FileUtil::getFileSize(Const::datafn.c_str());
        series_num = f_size / Const::tsLengthBytes;
        saxes = new unsigned short[series_num * Const::segmentNum];
        FILE *f = fopen(saxfn.c_str(), "rb");
        fread(saxes, sizeof(unsigned short), series_num * Const::segmentNum, f);
        fclose(f);
        cout << "Finish loading sax"<<endl;
    }else{
        series_num = Const::series_num;
        saxes = new unsigned short[series_num * Const::segmentNum];
        FILE *f = fopen(saxfn.c_str(), "rb");
        fread(saxes, sizeof(unsigned short), series_num * Const::segmentNum, f);
        fclose(f);
        cout << "Finish loading sax"<<endl;
    }

//    saxes = new unsigned short[series_num * Const::segmentNum];
//    FILE *f = fopen(saxfn.c_str(), "rb");
//    fread(saxes, sizeof(unsigned short), series_num * Const::segmentNum, f);
//    fclose(f);
//    cout << "Finish loading sax"<<endl;

    auto* root = new TardisTreeNode();
    root->children = new unordered_map<int, TardisTreeNode*>();

//    for(long i=0;i<series_num;++i){
//        int _id = SaxUtil::invSaxHeadFromSax(saxes + i * Const::segmentNum, Const::bitsCardinality, Const::segmentNum);
//        TardisTreeNode* target = ((*root->children)[_id]);
//        if(target == nullptr){
//            (*root->children)[_id] = new TardisTreeNode(_id, 1, saxes + i * Const::segmentNum, nullptr);
//            target = (*root->children)[_id];
//        }
//        target->insert(i);
//    }

    for(int i=0;i<series_num;++i) {
        if(i%10000000 == 0)  cout << i << endl;
        root->insert(i);
    }

    return root;
}

void TardisTreeNode::insert(int offset) {
    unsigned short *asax = saxes + (long)offset * Const::segmentNum;
    ++size;
    if(!isLeafNode()) {
        TardisTreeNode *target = routeToTarget(asax);
        if(target->layer == 1 && !target->isLeafNode())
            target->offsets.push_back(offset);
        target->insert(offset);
    } else if(size <= Const::th || layer >= Const::bitsCardinality)   offsets.push_back(offset);
    else{   //split
        offsets.push_back(offset);
        children = new unordered_map<int, TardisTreeNode*>();
        for(int off:offsets){
            unsigned short* _sax = saxes + (long)off * Const::segmentNum;;
            TardisTreeNode* target = routeToTarget(_sax);
            target->insert(off);
        }
        if(layer > 1)
            vector<int>().swap(offsets);
    }
}

TardisTreeNode* TardisTreeNode::routeToTarget(unsigned short *asax){
    assert(!isLeafNode());
    int nav_id = SaxUtil::invSaxHeadkFromSax(asax, Const::bitsCardinality, Const::segmentNum, layer + 1);
    if(children->find(nav_id) == children->end()) {
        auto *target =  new TardisTreeNode(nav_id, layer + 1, sax, this);
        (*children)[nav_id] = target;
    }
    return (*children)[nav_id];
}

TardisTreeNode* TardisTreeNode::route(unsigned short *asax) const{
    assert(!isLeafNode());
    int nav_id = SaxUtil::invSaxHeadkFromSax(asax, Const::bitsCardinality, Const::segmentNum, layer + 1);
    return (*children)[nav_id];
}

vector<OffsetDist *> *TardisTreeNode::preparePqUsingSaxLbNew(double bsf, float *queryPaa) const {
    auto*pq = new vector<OffsetDist*>();
    auto v_off = offsets;
    for(int i=0; i<size; ++i){
        double lbDist = SaxUtil::LowerBound_Paa_iSax(queryPaa, saxes + v_off[i] * Const::segmentNum,
                                                     Const::bitsCardinality);
        if(lbDist >= bsf)   continue;
        pq->push_back(new OffsetDist(v_off[i], lbDist));
    }
    sort(pq->begin(),  pq->end(), OffsetDistMinHeap());
    return pq;
}

void TardisTreeNode::exactSearchKnnInMemory(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap) const {
    if(isLeafNode()){
        double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
        float *ts;
        for(int offset: offsets){
            ts = dataset + (long)offset * Const::tsLength;
            double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts, Const::tsLength);
            if(heap.size() < k){
                heap.push_back(new PqItemSeries(ts, dist));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            }else if(dist < bsf){
                pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
                delete heap.back();
                heap.pop_back();
                heap.push_back(new PqItemSeries(ts, dist));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            }
            if(heap.size() >= k)    bsf = heap[0]->dist;
        }
    }
    else{
        for(auto & iter : *children)
            if(iter.second != nullptr)
                iter.second->exactSearchKnnInMemory(k, queryTs, heap);
    }
}

void TardisTreeNode::exactSearchKnnInMemory(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, int &threshold) const {
    if(isLeafNode()){
        double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
        float *ts;
        int bound = min(threshold, size);
        for(int i=0;i<bound;++i){
            ts = dataset + (long)offsets[i] * Const::tsLength;
            double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts, Const::tsLength);
            if(heap.size() < k){
                heap.push_back(new PqItemSeries(ts, dist));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            }else if(dist < bsf){
                pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
                delete heap.back();
                heap.pop_back();
                heap.push_back(new PqItemSeries(ts, dist));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            }
            if(heap.size() >= k)    bsf = heap[0]->dist;
        }
        threshold -=bound;
    }else{
        for(auto & iter : *children){
            if(threshold <=0) break;
            if(iter.second != nullptr)
                iter.second->exactSearchKnnInMemory(k, queryTs, heap, threshold);
        }
    }

}

void TardisTreeNode::exactSearchKnn(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap) const {
    if(isLeafNode()){
        double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
        FILE *f = fopen(Const::datafn.c_str(), "rb");
        float *ts;
        for(int offset: offsets){
            fseek(f, (long)offset * Const::tsLengthBytes, SEEK_SET);
            ts = new float[Const::tsLength];
            fread(ts, sizeof(float), Const::tsLength, f);
            double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts, Const::tsLength);
            if(heap.size() < k){
                heap.push_back(new PqItemSeries(ts, dist, true));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            }else if(dist < bsf){
                pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
                delete heap.back();
                heap.pop_back();
                heap.push_back(new PqItemSeries(ts, dist, true));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            } else delete[] ts;
            if(heap.size() >= k)    bsf = heap[0]->dist;
        }
        fclose(f);
    }
    else{
        for(auto & iter : *children)
            if(iter.second != nullptr)
                iter.second->exactSearchKnn(k, queryTs, heap);
    }
}

// only support nodes under this are less than node number
void TardisTreeNode::exactSearchKnn(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, int& node_number) const {
    if(node_number <= 0)    return;
    if(isLeafNode()){
        double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
        FILE *f = fopen(Const::datafn.c_str(), "rb");
        float *ts;
        for(int offset: offsets){
            fseek(f, (long)offset * Const::tsLengthBytes, SEEK_SET);
            ts = new float[Const::tsLength];
            fread(ts, sizeof(float), Const::tsLength, f);
            double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts, Const::tsLength);
            if(heap.size() < k){
                heap.push_back(new PqItemSeries(ts, dist, true));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            }else if(dist < bsf){
                pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
                delete heap.back();
                heap.pop_back();
                heap.push_back(new PqItemSeries(ts, dist, true));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            } else delete[] ts;
            if(heap.size() >= k)    bsf = heap[0]->dist;
        }
        fclose(f);
    }
    else{
        int max_pid = -1;
        for(auto & iter : *children)
            if(iter.second != nullptr)
                max_pid = max(max_pid, iter.second->partitionId);
//        assert(max_pid < node_number);
        // subtract all leaf node number
        node_number -= (max_pid + 1);
        for(auto & iter : *children)
            if(iter.second != nullptr)
                iter.second->exactSearchKnn(k, queryTs, heap);
    }
}


void TardisTreeNode::exactSearchKnnLeaf(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap) const {
    assert(isLeafNode());

    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    FILE *f = fopen(Const::datafn.c_str(), "rb");
    float *ts;
    for(int offset:offsets){
        fseek(f, (long)offset * Const::tsLengthBytes, SEEK_SET);
        ts = new float[Const::tsLength];
        fread(ts, sizeof(float), Const::tsLength, f);
        double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts, Const::tsLength);
        if(heap.size() < k){
            heap.push_back(new PqItemSeries(ts, dist, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(ts, dist, true));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        } else delete[] ts;

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }
    fclose(f);

}

TardisTreeNode *TardisTreeNode::loadFromDisk(const string &saxfn, const string &idxfn, bool on_disk) {
    long series_num;
    if(!on_disk){
        series_num  = Const::series_num;
        dataset = new float[series_num * Const::tsLength];
        FILE * df = fopen(Const::datafn.c_str(), "rb");
        long num = fread(dataset, sizeof(float), series_num * Const::tsLength, df);
        if(num != series_num * Const::tsLength){
            cout << "Cannot read data into memory "<<endl;
            exit(-1);
        }
        fclose(df);
    }else{
        long f_size = FileUtil::getFileSize(saxfn.c_str());
        series_num = f_size / (4 * Const::segmentNum);
    }

    saxes = new unsigned short[series_num * Const::segmentNum];
    FILE *f = fopen(saxfn.c_str(), "rb");
    fread(saxes, sizeof(unsigned short), series_num * Const::segmentNum, f);
    fclose(f);
    cout << "sax load finish. now load tree index." << endl;

    ifstream ifs(idxfn.c_str(), ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    auto *g = new TardisTreeNode();
    ia >> (*g);
    ifs.close();
    return g;

}

TardisTreeNode *TardisTreeNode::loadFromDisk(const string &idxfn) {
    ifstream ifs(idxfn.c_str(), ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    auto *g = new TardisTreeNode();
    ia >> (*g);
    ifs.close();
    return g;

}

//double TardisTreeNode::scoreGraph(unordered_map<int, TardisTreeNode*>*nodes_map){
//    int node_number = 0, total_size = 0;
//    for(auto& iter:*nodes_map){
//        if(iter.second != nullptr && iter.second->isLeafNode())
//            node_number++, total_size += iter.second->size;
//    }
//    return Const::node_number_weight * node_number / Const::vertexNum + Const::node_size_weight * total_size / (Const::th * node_number);
//}
//
//static int gid = 0;
//
//void TardisTreeNode::partition(TardisTreeNode* root, vector<vector<int>> *g){
//    if(root == nullptr || root->children == nullptr || root->size <= Const::th)
//        return;
//    if(scoreGraph(root->children) < Const::dense_threshold)
//        Partition::bubble(root->children);
//    else
//        Partition::growingGraph(root->children, g);
//    cout << gid++ << endl;
//    for(auto& iter:*root->children)
//        partition(iter.second, g);
//}

void TardisTreeNode::mergingTardis(TardisTreeNode* root){
    if(root == nullptr || root->children == nullptr || root->size <= Const::th)   return;

    // do merge
    vector<TardisTreeNode*>candidates;
    for(auto &iter:*root->children){
        TardisTreeNode *tmp = iter.second;
        if(tmp != nullptr && tmp->isLeafNode() && tmp->size <= Const::th)
            candidates.push_back(tmp);
    }
    sort(candidates.begin(), candidates.end(), TardisTreeNode::order);
    int pid = 0, cur_size = 0;
    for(TardisTreeNode*_:candidates){
        if(cur_size + _->size <= Const::th) {
            _->partitionId = pid;
            cur_size += _->size;
        }else{
            ++pid;
            cur_size = _->size;
            _->partitionId = pid;
        }
    }

    // recursive
    for(auto &iter:*root->children){
        TardisTreeNode *tmp = iter.second;
        if(tmp!= nullptr && !tmp->isLeafNode())
            mergingTardis(tmp);
    }
}

int TardisTreeNode::getLeafNodeNumber(TardisTreeNode * root){
    if(root == nullptr) return 0;
    if(root->size <= Const::th) return 1;
    int sum = 0;
    unordered_set<int>seen;
    for(auto& iter:*root->children){
        if(iter.second == nullptr) continue;
        if(iter.second->size <= Const::th){
            if(seen.find(iter.second->partitionId) != seen.end())   continue;
            sum++;
            seen.insert(iter.second->partitionId);
        }else{
            sum += getLeafNodeNumber(iter.second);
        }
    }
    return sum;
}

void TardisTreeNode::getLeafNodeSize(TardisTreeNode * root, ofstream &f){
    if(root == nullptr) return;
    if(root->size <= Const::th){
        int s = root->size;
        f << to_string(s) + ",";
        return;
    }
    int sum = 0;
    for(auto& iter:*root->children){
        if(iter.second == nullptr) continue;
        if(iter.second->size <= Const::th){
            f << to_string(iter.second->size) + ",";
        }else{
            getLeafNodeSize(iter.second, f);
        }
    }
    return;
}

int TardisTreeNode::getNodeNumber(TardisTreeNode * root){
    if(root == nullptr) return 0;
    if(root->children == nullptr) return 1;
    int sum = 0;
    unordered_set<int>seen;
    for(auto& iter:*root->children){
        if(iter.second == nullptr) continue;
        if(iter.second->size <= Const::th){
            if(seen.find(iter.second->partitionId) != seen.end())   continue;
            sum++;
            seen.insert(iter.second->partitionId);
        }else{
            sum += getNodeNumber(iter.second);
        }
    }
    return sum + 1;
}

int TardisTreeNode::getTotalSize(TardisTreeNode * root){
    if(root == nullptr) return 0;
    if(root->children == nullptr) return root->size;
    int sum = 0;
    unordered_set<int>seen;
    for(auto& iter:*root->children){
        if(iter.second == nullptr) continue;   // may not happen
        sum += getTotalSize(iter.second);
    }
    return sum;
}

int TardisTreeNode::getHeight(TardisTreeNode *root){
    if(root == nullptr) return 0;
    if(root->children == nullptr) return 1;
    int max = 0;
    for(auto& iter:*root->children){
        if(iter.second == nullptr) continue;
        int h = getHeight(iter.second);
        max = max >= h? max : h;
    }
    return max + 1;
}

int TardisTreeNode::getSumHeight(TardisTreeNode * root) const{
    if(root == nullptr) return 0;
    if(root->children == nullptr) return root->layer;
    int sum = 0;
    unordered_set<int>seen;
    for(auto& iter:*root->children){
        if(iter.second == nullptr || seen.find(iter.second->partitionId) != seen.end()) continue;   // may not happen
        sum += getSumHeight(iter.second);
        seen.insert(iter.second->partitionId);
    }
    return sum;
}

int TardisTreeNode::getGraphNum(TardisTreeNode* root){
    if(root == nullptr || root->children == nullptr || root->size <= Const::th)   return 0;
    int sum = 1;
    for(auto iter:*root->children){
        sum += getGraphNum(iter.second);
    }
    return sum;
}

void TardisTreeNode::getIndexStats(){
    ofstream outfile;
    outfile.open("/mnt/c/codes/rand-leaf.out", ios::out | ios::app);//输入文件的路径
    getLeafNodeSize(this, outfile);
    int total_leaf_node_num = getLeafNodeNumber(this);
    int total_size = getTotalSize(this);
    cout << "Total size = " << total_size << endl;
    cout <<"Total nodes number = " << getNodeNumber(this) << endl;
    cout << "Leaf node number = " << total_leaf_node_num << endl;
    cout << "Max. height = " << getHeight(this) - 1 <<endl;
    cout << "Avg. Height = " << (double)getSumHeight(this) / (double) total_leaf_node_num << endl;
    cout <<"Avg. Filling Factor = "<< total_size / (double)total_leaf_node_num / Const::th << endl;
    outfile.close();
}

long TardisTreeNode::generateSaxTbl(){
    string fn = Const::datafn;
    long fs = FileUtil::getFileSize(fn.c_str());
    long series_num = fs / Const::tsLengthBytes;
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
        long num;
        if(rest > 4000000)    num = 4000000;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto start = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto end = chrono::system_clock::now();
        SAX_PAA_READ_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();

        for(long i=0;i<num;++i){
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