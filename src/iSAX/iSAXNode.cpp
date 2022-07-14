//
// Created by wzy on 2021/11/20.
//

#include "../../include/DataStructures/iSAXNode.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Utils/SaxUtil.h"
#include "../../include/Utils/TimeSeriesUtil.h"
#include "../../include/DataStructures/TimeSeries.h"
#include "../../include/Const.h"
#include <iostream>
#include <cmath>

unsigned short * iSAXRoot::saxes = nullptr;
float *iSAXRoot::paas = nullptr;
float *iSAXRoot::dataset = nullptr;

iSAXRoot * iSAXRoot::buildIndex(const string &saxfn, long series_num) {
    loadSax(saxfn, series_num);
    cout << "Finish loading sax"<<endl;

    auto* root = new iSAXRoot();

    for(int i=0;i<series_num;++i)
        root->insertPrev(i);

    root->bulkInsert(false);

    delete[] saxes;
    return root;
}

iSAXRoot * iSAXRoot::buildIndexFuzzy(const string &saxfn, const string &paafn) {
    long series_num = loadSax(saxfn, -1);
    cout << "Finish loading sax"<<endl;

    loadPaa(paafn);


    auto* root = new iSAXRoot();

    for(int i=0;i<series_num;++i)
        root->insertPrev(i);

    root->bulkInsert(true);

    delete[] saxes;
    return root;
}

void iSAXRoot::bulkInsert(bool fuzzy) {
    for(auto node : children){
        if(node == nullptr) continue;
        for(int offset:node->fbl){
            node->insert(offset, true, fuzzy);
        }
        vector<int>().swap(node->fbl);
    }
}

void iSAXRoot::insertPrev(int offset){
    iSAXNode * child  = route2firstLayer(saxes + offset * Const::segmentNum);
    child->fbl.push_back(offset);
}

void iSAXNode::insert(int offset, bool isReal, bool fuzzy) {
    iSAXNode* node = route2Target(iSAXRoot::saxes + offset * Const::segmentNum, isReal);
    node->append(offset, isReal);
    if(node->size > Const::th)
        node->split(fuzzy);
}

void iSAXNode::split(bool fuzzy) {
    segment_id = chooseSegment(false);
    left = new iSAXNode(this, true);
    right = new iSAXNode(this, false);
    int bc = bits_cardinality[segment_id];
    for(int offset:offsets){
        unsigned short *_sax = iSAXRoot::saxes + offset * Const::segmentNum;
        int res = _sax[segment_id] >> (Const::bitsCardinality - bc - 1);
        if(res % 2) right->append(offset, true);
        else left->append(offset, true);
    }
    vector<int>().swap(offsets);
    if(fuzzy)   fuzzyDuplicate(this);
}

int iSAXNode::chooseSegment(bool is1st) {
    int ret = -1, sum = 0; double dist = numeric_limits<double>::max(), lb, ub, mean;
    for(int i=0;i<Const::segmentNum;++i) {
        if(bits_cardinality[i] >= Const::bitsCardinality)  continue;
        sum = 0;
        if(is1st){
            for(int offset:fbl)
                sum += iSAXRoot::saxes[offset * Const::segmentNum];
            sum /= fbl.size();
        }else{
            for(int offset:offsets)
                sum += iSAXRoot::saxes[offset * Const::segmentNum];
            sum /= size;
        }

        SaxUtil::getValueRange(sum, Const::bitsCardinality, &lb, &ub);
        if(sum == 0)    mean = ub;
        else if(sum == (1 << Const::bitsCardinality) - 1) mean = lb;
        else mean = (ub + lb) / 2.0;
        SaxUtil::getValueRange(sax[i] << 1, bits_cardinality[i] + 1, &lb, &ub);
        mean = abs(ub - mean);
        if(mean < dist){
            ret = i;
            dist = mean;
        }
    }
    assert(ret!=-1);
    return ret;
}

void iSAXNode::append(int offset, bool isReal) {
    size++;
    if(isReal)  actual_size++;
    offsets.push_back(offset);
}

void iSAXNode::fuzzyDuplicate(iSAXNode* parent){
    double range, lb, ub, split_line;
    float *paa;
    iSAXNode*node;
    int new_id, id, seg = parent->segment_id;
    bool flag = false;
    SaxUtil::getValueRange(parent->sax[seg] << 1, parent->bits_cardinality[seg] + 1, &lb, &split_line);

    {
        node = parent->left;
        if(node->sax[seg] == (1 << node->bits_cardinality[seg]) - 1)
            SaxUtil::getValueRange(node->sax[seg] - 1, node->bits_cardinality[seg], &lb, &ub);
        else if(node->sax[seg] == 0)
            SaxUtil::getValueRange(1, node->bits_cardinality[seg], &lb, &ub);
        else
            SaxUtil::getValueRange(node->sax[seg] , node->bits_cardinality[seg], &lb, &ub);
        range = (ub -lb) * Const::boundary;

        if(node->size > Const::filling_factor_1st * Const::th)  {
            for(int j=0;j<node->actual_size;++j){
                // iterate all actual series
                flag = false;
                paa = iSAXRoot::paas + node->offsets[j] * Const::segmentNum;
                if(node->sax[seg] == (1 << node->bits_cardinality[seg]) - 1){
                    if(paa[seg] - ub <= range) flag = true;
                }else if(node->sax[seg] == 0){
                    if(lb - paa[seg] <= range)  flag = true;
                }else if(node->sax[seg] % 2 == 0){
                    if(abs(ub - paa[seg]) <= range) flag = true;
                }else{
                    if(abs(paa[seg] - lb) <= range)  flag = true;
                }
                if(!flag)   continue;
                if(parent->right->actual_size <= Const::th && parent->right->size >= Const::th);
                else parent->right->insert(node->offsets[j], false, false);
            }
        }
    }

    {
        node = parent->right;
        if(node->sax[seg] == (1 << node->bits_cardinality[seg]) - 1)
            SaxUtil::getValueRange(node->sax[seg] - 1, node->bits_cardinality[seg], &lb, &ub);
        else if(node->sax[seg] == 0)
            SaxUtil::getValueRange(1, node->bits_cardinality[seg], &lb, &ub);
        else
            SaxUtil::getValueRange(node->sax[seg] , node->bits_cardinality[seg], &lb, &ub);
        range = (ub -lb) * Const::boundary;

        if(node->size > Const::filling_factor_1st * Const::th)  {
            for(int j=0;j<node->actual_size;++j){
                // iterate all actual series
                flag = false;
                paa = iSAXRoot::paas + node->offsets[j] * Const::segmentNum;
                if(node->sax[seg] == (1 << node->bits_cardinality[seg]) - 1){
                    if(paa[seg] - ub <= range) flag = true;
                }else if(node->sax[seg] == 0){
                    if(lb - paa[seg] <= range)  flag = true;
                }else if(node->sax[seg] % 2 == 0){
                    if(abs(ub - paa[seg]) <= range) flag = true;
                }else{
                    if(abs(paa[seg] - lb) <= range)  flag = true;
                }
                if(!flag) continue;
                if(parent->left->actual_size <= Const::th && parent->left->size >= Const::th);
                else parent->left->insert(node->offsets[j], false, false);
            }
        }
    }
}

iSAXNode * iSAXNode::route2Target(const unsigned short *_sax, bool update) {
    if(isLeafNode())    return this;
    iSAXNode*cur = this;
    while (!cur->isLeafNode()){
        if(update)  cur->size++;
        cur = (_sax[cur->segment_id] >> (Const::bitsCardinality - cur->bits_cardinality[cur->segment_id] - 1)) % 2 == 0 ? cur->left : cur->right;
    }
    return cur;
}

bool iSAXNode::route2Left(const unsigned short *_sax){
    return (_sax[segment_id] >> (Const::bitsCardinality - bits_cardinality[segment_id] - 1)) % 2 == 0;
}

iSAXNode * iSAXNode::route2Target(const unsigned short *_sax, int threshold) {
    if(isLeafNode())    return this;
    iSAXNode*cur = this;
    while (!cur->isLeafNode()){
        if(cur->size <= threshold)    return cur;
        cur = (_sax[cur->segment_id] >> (Const::bitsCardinality - cur->bits_cardinality[cur->segment_id] - 1)) % 2 == 0 ? cur->left : cur->right;
    }
    return cur;
}

iSAXNode * iSAXNode::route1Step(const unsigned short *_sax) {
    assert(!isLeafNode());
    return (_sax[segment_id] >> (Const::bitsCardinality - bits_cardinality[segment_id] - 1)) % 2 == 0 ? left : right;

}

iSAXNode * iSAXRoot::route2Target(unsigned short *_sax, bool update) {
    int id1 = SaxUtil::invSaxHeadFromSax(_sax, Const::bitsCardinality, Const::segmentNum);
    if(children[id1] == nullptr) {
        children[id1] = new iSAXNode(1, _sax);
        for (auto &_:children[id1]->bits_cardinality) _ = 1;
        return children[id1];
    }
    return children[id1]->route2Target(_sax, update);
}

iSAXNode * iSAXRoot::route2Target(unsigned short *_sax, int threshold) {
    int id1 = SaxUtil::invSaxHeadFromSax(_sax, Const::bitsCardinality, Const::segmentNum);
    if(children[id1] == nullptr) {
        return nullptr;
    }
    return children[id1]->route2Target(_sax, threshold);
}

iSAXNode * iSAXRoot::route2firstLayer(unsigned short *_sax) {
    int id1 = SaxUtil::invSaxHeadFromSax(_sax, Const::bitsCardinality, Const::segmentNum);
    if(children[id1] == nullptr) {
        children[id1] = new iSAXNode(1, _sax);
        for (auto &_:children[id1]->bits_cardinality) _ = 1;
    }
    return children[id1];
}

void iSAXNode::exactSearchKnn(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap) {
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

void iSAXNode::exactSearchKnnInMemory(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, int &threshold) {
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    float *ts;
    int bound = min(threshold, (int)offsets.size());
    for(int i=0;i<bound;++i){
        ts = iSAXRoot::dataset + (long)offsets[i] * Const::tsLength;
        double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts, Const::tsLength);
        if(heap.size() < k){
            heap.push_back(new PqItemSeries(ts, dist, false, false));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(ts, dist, false, false));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }
    threshold -= bound;

}

void iSAXNode::exactSearchKnnInMemory(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap) {
    int x=numeric_limits<int>::max();
    if(isLeafNode()) exactSearchKnnInMemory(k, queryTs, heap, x);
    if(left != nullptr) left->exactSearchKnnInMemory(k, queryTs, heap);
    if(right != nullptr)    right->exactSearchKnnInMemory(k, queryTs, heap);
}

void iSAXNode::exactSearchKnnInMemoryEntry(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, int &rest) {
    if(rest <= 0)   return;
    if(isLeafNode()) { exactSearchKnnInMemory(k, queryTs, heap, rest);  }
    if(rest <= 0)   return;
    if(left != nullptr) left->exactSearchKnnInMemoryEntry(k, queryTs, heap, rest);
    if(rest <= 0)   return;
    if(right != nullptr)    right->exactSearchKnnInMemoryEntry(k, queryTs, heap,rest);
}


long iSAXRoot::loadSax(const string &saxfn, long series_num) {
    if(series_num == -1){
        long f_size = FileUtil::getFileSize(saxfn.c_str());
        series_num = f_size / (sizeof(unsigned short ) * Const::segmentNum);
    }
    saxes = new unsigned short[series_num * Const::segmentNum];
    FILE *f = fopen(saxfn.c_str(), "rb");
    fread(saxes, sizeof(unsigned short), series_num * Const::segmentNum, f);
    fclose(f);
    cout << "Finish loading sax"<<endl;
    return series_num;
}

void iSAXRoot::loadPaa(const string & paafn){
    long f_size = FileUtil::getFileSize(paafn.c_str());
    paas = new float [f_size / sizeof(float)];
    FILE *f = fopen(paafn.c_str(), "rb");
    fread(paas, sizeof(float), f_size / sizeof(float), f);
    fclose(f);
    cout << "Finish loading paa"<<endl;
}

iSAXRoot *iSAXRoot::loadFromDisk(const string &saxfn, const string &idxfn, long series_num) {
    loadSax(saxfn, series_num);
    if(series_num != -1){
        dataset = new float[series_num * Const::tsLength];
        FILE *df = fopen(Const::datafn.c_str(), "rb");
        fread(dataset, sizeof(float), series_num * Const::tsLength, df);
        fclose(df);
    }
    ifstream ifs(idxfn.c_str(), ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    auto *g = new iSAXRoot();
    ia >> (*g);
    ifs.close();
    return g;
}

iSAXRoot *iSAXRoot::loadFromDisk(const string &saxfn, const string &paafn,const string &idxfn) {
    loadSax(saxfn, 0);
    loadPaa(paafn);
    ifstream ifs(idxfn.c_str(), ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    auto *g = new iSAXRoot();
    ia >> (*g);
    ifs.close();
    return g;
}
