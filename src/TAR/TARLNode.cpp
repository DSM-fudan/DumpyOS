//
// Created by wzy on 2022/7/1.
//
#include <cstdio>
#include <chrono>
#include "../../include/TAR/TARLNode.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Utils/SaxUtil.h"
#include "../../include/Const.h"
#include "../../include/Utils/TimeSeriesUtil.h"

extern long TAR_WRITE_TIME, TAR_READ_TIME;

TARLNode * TARLNode::buildLocalIndex(int pid){
    long unitSize = Const::tsLengthBytes + 8 + (Const::segmentNum / 4) * Const::bitsCardinality * 2;
    long series_num = FileUtil::getFileSize((Const::tardisfn + to_string(pid)).c_str()) / unitSize;
    if(series_num <= 0) return nullptr;
    vector<ser_item>datas(series_num);
    unsigned invSax[Const::bitsCardinality];
    auto start = chrono::system_clock::now();
    FILE *f = fopen((Const::tardisfn + to_string(pid)).c_str(),"rb");
    for(long i=0; i<series_num;++i) {
        datas[i].ts.resize(Const::tsLength);
        fread(&datas[i].id, sizeof(long), 1, f);
        fread(&datas[i].ts[0], sizeof(float), Const::tsLength, f);
        fread(invSax, sizeof(unsigned), Const::bitsCardinality, f);
        datas[i].invsax = SaxUtil::invSax2String(invSax);
    }
    fclose(f);
    auto end = chrono::system_clock::now();
    TAR_READ_TIME +=  chrono::duration_cast<chrono::microseconds>(end - start).count();

    auto root = new TARLNode();
    for(ser_item&data: datas){
        root->addRecordORIGIN(data, pid);
    }

    root = prueSingleElementLayer(root);
    if(root->rcdNbr <= Const::th){
        root->fetchAllRecords(root->buffer);
        root->deleteDescendants();
        root->descendants.clear();
    }else
        root->shrinkNodeLayer();
    start = chrono::system_clock::now();
    root->save2Disk((Const::tardisfn + to_string(pid) + "_tmp"));
    end = chrono::system_clock::now();
    TAR_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();

    // traverse through all leaf nodes and write them
//    long cur_offset = 0;
//    vector<TARLNode*>leaves;
//    collectAllLeaves(leaves);
//    FILE *outf = fopen((Const::tardisfn + to_string(pid) + "_tmp").c_str(), "wb");
//    for(TARLNode* leaf: leaves){
//        leaf->file_offset = cur_offset;
//        leaf->file_length = leaf->buffer.size();
//        auto local_start = chrono::system_clock::now();
//        for(item& unit: leaf->buffer){
//            fwrite(&unit.id, sizeof(long ), 1, outf);
//            fwrite(unit.ts, sizeof(float), Const::tsLength, outf);
//            SaxUtil::string2invSAX(unit.invsax, invSax);
//            delete[] unit.ts;
//            fwrite(invSax, sizeof(unsigned ), Const::bitsCardinality, outf);
//        }
//        auto local_end = chrono::system_clock::now();
//        TAR_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(local_end - local_start).count();
//        cur_offset += leaf->buffer.size();
//        leaf->buffer.clear();
//    }
//
//    fclose(outf);

    FileUtil::FileRemove((Const::tardisfn + to_string(pid)).c_str());
    if(rename((Const::tardisfn + to_string(pid) + "_tmp").c_str(), (Const::tardisfn + to_string(pid)).c_str())){
        Const::logPrint("rename error!");
        exit(-6);
    }
    return root;
}

//void TARLNode::insertBatch(vector<item*>&datas, int pid){
//    for(item*data: datas){
//        addRecord(data, pid);
//    }
//
//    // traverse through all leaf nodes and write them
//    long cur_offset = 0;
//    vector<TARLNode*>leaves;
//    collectAllLeaves(leaves);
//    FILE *outf = fopen((Const::tardisfn + to_string(pid) + "_tmp").c_str(), "wb");
//    unsigned invSax[Const::bitsCardinality];
//    for(TARLNode* leaf: leaves){
//        if(leaf->file_length > 0){
//            long unit_size= sizeof(unsigned ) * Const::bitsCardinality + Const::tsLengthBytes + sizeof(long );
//            long id;
//            string pack_file = Const::tardisfn + to_string(pid);
//            auto start = chrono::system_clock::now();
//            FILE *f = fopen(pack_file.c_str(),"rb");
//            fseek(f, file_offset*unit_size, SEEK_SET);
//            for(int i=0; i<leaf->file_length;++i) {
//                item cur;
//                cur.ts = new float[Const::tsLength];
//                fread(&id, sizeof(long), 1, f);
//                fread(cur.ts, sizeof(float), Const::tsLength, f);
//                fread(invSax, sizeof(unsigned), Const::bitsCardinality, f);
//                cur.invsax = SaxUtil::invSax2String(invSax);
//                cur.id = id;
//                leaf->buffer.push_back(cur);
//            }
//            fclose(f);
//            auto end = chrono::system_clock::now();
//            TAR_READ_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();
//        }
//        leaf->file_offset = cur_offset;
//        leaf->file_length = leaf->buffer.size();
//        auto start = chrono::system_clock::now();
//        for(item& unit: leaf->buffer){
//            fwrite(&unit.id, sizeof(long ), 1, outf);
//            fwrite(unit.ts, sizeof(float), Const::tsLength, outf);
//            SaxUtil::string2invSAX(unit.invsax, invSax);
//            delete[] unit.ts;
//            fwrite(invSax, sizeof(unsigned ), Const::bitsCardinality, outf);
//        }
//        auto end = chrono::system_clock::now();
//        TAR_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();
//        cur_offset += leaf->buffer.size();
//        leaf->buffer.clear();
//    }
//
//    fclose(outf);
//    FileUtil::FileRemove((Const::tardisfn + to_string(pid)).c_str());
//    if(rename((Const::tardisfn + to_string(pid) + "_tmp").c_str(), (Const::tardisfn + to_string(pid)).c_str())){
//        Const::logPrint("rename error!");
//        exit(-6);
//    }
//}

void TARLNode::collectAllLeaves(vector<TARLNode*>&leaves){
    if(descendants.empty()) leaves.push_back(this);
    for(auto &iter:descendants){
        iter.second->collectAllLeaves(leaves);
    }
}

//void TARLNode::addRecord(item *data, int pid) {
//    rcdNbr++;
//    if(rcdNbr <= Const::th){
//        item newitem = *data;
//        newitem.ts = new float[Const::tsLength];
//        copy(data->ts, data->ts+Const::tsLength, newitem.ts);
//        buffer.push_back(newitem);
//    }else{
//        if(descendants.empty()){
//            split(pid);
//        }
//        string new_key = SaxUtil::invSaxHeadKFromInvSax(data->invsax, layer + 1);
//        if(!descendants.contains(new_key)){
//            auto node = new TARLNode();
//            node->layer = layer + 1;
//            node->invSAX = new_key;
//            descendants[new_key] = node;
//        }
//        descendants[new_key]->addRecord(data, pid);
//    }
//
//}

void TARLNode::addRecordORIGIN(const ser_item &data, int pid) {
    rcdNbr++;
    if(this->layer == Const::bitsCardinality){
        ser_item newitem(data);
        buffer.push_back(newitem);
    }else{
        string new_key = SaxUtil::invSaxHeadKFromInvSax(data.invsax, layer + 1);
        if(!descendants.contains(new_key)){
            auto node = new TARLNode();
            node->layer = layer + 1;
            node->invSAX = new_key;
            descendants[new_key] = node;
        }
        descendants[new_key]->addRecordORIGIN(data, pid);
    }
}

void TARLNode::shrinkNodeLayer(){
    for(auto &iter: descendants){
        if(!iter.second->descendants.empty()){
            if(iter.second->rcdNbr <= Const::th){
                iter.second->fetchAllRecords(iter.second->buffer);
                iter.second->deleteDescendants();
                iter.second->descendants.clear();
            }else{
                iter.second->shrinkNodeLayer();
            }
        }
    }
}

void TARLNode::deleteDescendants(){
    for(auto &iter: descendants){
        if(iter.second->descendants.empty()){
            delete iter.second;
        }else{
            iter.second->deleteDescendants();
            delete iter.second;
        }
    }
}

void TARLNode::fetchAllRecords(vector<ser_item>&buf){
    if(descendants.empty()){
        for(auto &_:buffer){
            buf.push_back(_);
        }
    }else{
        for(auto &iter:descendants){
            iter.second->fetchAllRecords(buf);
        }
    }
}

TARLNode* TARLNode::prueSingleElementLayer(TARLNode* root){
    TARLNode* cur = root, *before;
    while (cur->descendants.size() == 1){
        before = cur;
        cur = cur->descendants.begin()->second;
        delete before;
    }
    return cur;
}

//void TARLNode::split(int pid) {
//    vector<item>file_data;
//    unsigned invSax[Const::bitsCardinality];
//    long unit_size= sizeof(unsigned ) * Const::bitsCardinality + Const::tsLengthBytes + sizeof(long );
//    long id;
//    auto start = chrono::system_clock::now();
//    if(file_length > 0){
//        FILE *f = fopen((Const::tardisfn + to_string(pid)).c_str(),"rb");
//        fseek(f, file_offset*unit_size, SEEK_SET);
//        for(int i=0;i<file_length;++i) {
//            item cur;
//            cur.ts = new float[Const::tsLength];
//            fread(&id, sizeof(long), 1, f);
//            fread(cur.ts, sizeof(float), Const::tsLength, f);
//            fread(invSax, sizeof(unsigned), Const::bitsCardinality, f);
//            cur.invsax = SaxUtil::invSax2String(invSax);
//            cur.id = id;
//            file_data.push_back(cur);
//        }
//        fclose(f);
//    }
//    auto end = chrono::system_clock::now();
//    TAR_READ_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();
//
//    //split file data to children
//    for(auto&cur:file_data){
//        string new_key = SaxUtil::invSaxHeadKFromInvSax(cur.invsax, layer + 1);
//        if(!descendants.contains(new_key)){
//            auto node = new TARLNode();
//            node->layer = layer + 1;
//            node->invSAX = new_key;
//            descendants[new_key] = node;
//        }
//        descendants[new_key]->buffer.push_back(cur);
//    }
//
//    // split buffer data
//    for(auto &cur:buffer){
//        string new_key = SaxUtil::invSaxHeadKFromInvSax(cur.invsax, layer + 1);
//        if(!descendants.contains(new_key)){
//            auto node = new TARLNode();
//            node->layer = layer + 1;
//            node->invSAX = new_key;
//            descendants[new_key] = node;
//        }
//        descendants[new_key]->buffer.push_back(cur);
//    }
//
//    file_offset = -1;
//    file_length = 0;
//    buffer.clear();
//}

TARLNode* TARLNode::route2Leaf(string& invsax){
    if(descendants.empty()) return this;
    string key = SaxUtil::invSaxHeadKFromInvSax(invsax, layer + 1);
    if(descendants.contains(key)){
        return descendants[key]->route2Leaf(invsax);
    }else{
        // if no target leaf node, then search its sibling randomly
        return descendants.begin()->second->route2Leaf(invsax);
    }
}

extern int _search_num, layer;
extern long READ_TIME;
void TARLNode::search(int k, const float *query, vector<PqItemSeries *> &heap) const{
    assert(descendants.empty());
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;

    _search_num += rcdNbr; ::layer = layer;
    for(int i=0;i<rcdNbr;++i){
//        struct timeval start{};
//        Const::timer_start(&start);
        double dist = TimeSeriesUtil::euclideanDist(query, buffer[i].ts, Const::tsLength, bsf);
//        DIST_CALC_TIME += Const::timer_end(&start);

        if(heap.size() < k){
            heap.push_back(new PqItemSeries(buffer[i].ts, dist, true, false));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(buffer[i].ts, dist, true, false));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
}

void TARLNode::search_dtw(int k, float *query, vector<PqItemSeries *> &heap, int window_size) const{
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
//
    _search_num += rcdNbr; ::layer = layer;
    for(int i=0;i<rcdNbr;++i){
//        struct timeval start{};
//        Const::timer_start(&start);
        double dist = TimeSeriesUtil::dtw(query, buffer[i].ts, Const::tsLength, window_size, bsf);
//        DIST_CALC_TIME += Const::timer_end(&start);

        if(heap.size() < k){
            heap.push_back(new PqItemSeries(buffer[i].ts, dist, true, false));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }else if(dist < bsf){
            pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            delete heap.back();
            heap.pop_back();
            heap.push_back(new PqItemSeries(buffer[i].ts, dist, true, false));
            push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
        }

        if(heap.size() >= k)    bsf = heap[0]->dist;
    }

    for(PqItemSeries*s: heap){
        if(s->needDeepCopy) s->copyData();
    }
}

int TARLNode::getLeafNodeNbr(){
    if(descendants.empty()) return 1;
    int ret = 0;
    for(auto &iter:descendants){
        ret += iter.second->getLeafNodeNbr();
    }
    return ret;
}

int TARLNode::getNodeNbr(){
    if(descendants.empty()) return 1;
    int ret = 1;
    for(auto &iter:descendants){
        ret += iter.second->getNodeNbr();
    }
    return ret;
}

int TARLNode::getMaxHeight(){
    if(descendants.empty()) return 1;
    int ret = 0;
    for(auto &iter:descendants){
        ret = max(ret, iter.second->getMaxHeight());
    }
    return ret + 1;
}

long TARLNode::getSumHeight(){
    if(descendants.empty()) return 1;
    long ret = 0;
    for(auto &iter:descendants){
        ret += iter.second->getSumHeight() + 1;
    }
    return ret;
}

void TARLNode::save2Disk(const string &output) {
    ofstream ofs(output.c_str(), ios::binary);
    boost::archive::binary_oarchive oa(ofs);
    oa << (*this);
    ofs.close();
}

TARLNode *TARLNode::loadFromDisk(const string &idxfn) {
    ifstream ifs(idxfn, ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    auto *g = new TARLNode();
    ia >> (*g);
    ifs.close();
    return g;
}
