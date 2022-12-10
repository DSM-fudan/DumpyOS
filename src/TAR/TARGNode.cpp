//
// Created by wzy on 2022/7/1.
//
#include <unordered_set>
#include <chrono>
#include "../../include/TAR/TARGNode.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Utils/SaxUtil.h"

int TARGNode::pid_factory = 0;
int TARGNode::block_threshold = 0;

extern long TAR_READ_TIME = 0, TAR_CPU_TIME = 0, TAR_WRITE_TIME = 0;

TARGNode* TARGNode::buildIndex(){
    FileUtil::checkDirClean(Const::tardisfn.c_str());

    auto start = chrono::system_clock::now();
    long series_num = FileUtil::getFileSize(Const::datafn.c_str()) / Const::tsLengthBytes;
    Const::logPrint("Start building. Total series number = " + to_string(series_num));

    auto root = new TARGNode();
    auto sample_invSAX_tbl = loadSampling();
    Const::logPrint("loaded sampling finished!");
    auto*invsax_freq = new unordered_map<string, int>();
    computeFreqTbl(invsax_freq, sample_invSAX_tbl, series_num * Const::tardis_sample_percent);
    delete[] sample_invSAX_tbl;
    vector<unordered_map<string, int>*>* node_list  = buildTARGbySampling(invsax_freq);
    Const::logPrint("get node list finished!");
    loadCollectionSampling(node_list, root);
    for(auto*_:*node_list)  delete _;
    delete node_list;
    Const::logPrint("build TARDIS-G finished!");

    long blockSpace = 0.95 * 128 * 1024 * 1024;
    long unitSize = Const::tsLengthBytes + 8 + (Const::segmentNum / 4) * Const::bitsCardinality * 2;
    block_threshold = blockSpace / unitSize;
    assignPIdSampling(root);
    Const::logPrint("finish leaf node packing for TARDIS-G");

    pid_factory++;
    vector<TARLNode*>local_roots(pid_factory);
    Const::logPrint("Totally "+ to_string(pid_factory) + " packs.");
    partition(root, series_num);

    for(int i=0; i< pid_factory; ++i)
        TARLNode::buildLocalIndex(i);
    Const::logPrint("Finish building local index");


//    root->assignLocalRootForAllLeaves(local_roots);
    auto end = chrono::system_clock::now();
    TAR_CPU_TIME = chrono::duration_cast<chrono::microseconds>(end - start).count() - TAR_READ_TIME - TAR_WRITE_TIME;

    Const::logPrint("build index finished");
    cout << "WRITE TIME = " << TAR_WRITE_TIME / 1000.0 <<" ms."<<endl
    <<"READ TIME = " << TAR_READ_TIME / 1000.0 <<" ms." <<endl
    <<"CPU Time = " << TAR_CPU_TIME / 1000.0 <<"ms" <<endl;
    return root;
}

void TARGNode::partition(TARGNode*root, long series_num){
    long rest = series_num, cur = 0;
    long unitSize = Const::tsLengthBytes + 8 + (Const::segmentNum / 4) * Const::bitsCardinality * 2;
    long buffer_num = Const::fbl_size * 1024L * 1024L / unitSize;
    FILE *f  = fopen(Const::datafn.c_str(), "rb");
    unsigned invSax[Const::bitsCardinality];

    while(rest > 0) {
        long num;
        if (rest > buffer_num) num = buffer_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];
        rest -= num;

        auto local_start = chrono::system_clock::now();
        fread(tss, sizeof(float ), num * Const::tsLength, f);
        auto local_end = chrono::system_clock::now();
        TAR_READ_TIME += chrono::duration_cast<chrono::microseconds>(local_end - local_start).count();
        vector<item>saxTsPairRdd(num);
        unsigned invsax[Const::cardinality];
        for(long i =0; i< num;++i){
            saxTsPairRdd[i].id = cur + i;
            saxTsPairRdd[i].ts = tss + i * Const::tsLength;
            auto _ = SaxUtil::saxFromTs(tss + i *Const::tsLength, Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
            SaxUtil::invSaxFromSax(_, invsax, Const::bitsCardinality, Const::segmentNum);
            delete _;
            saxTsPairRdd[i].invsax = SaxUtil::invSax2String(invsax);
        }
        unordered_map<int, vector<item*> > saxTsIdRdd;

        TARGNode* parent;
        int cur_pid;
        for(long i=0;i<num;++i){
            parent = root->route(saxTsPairRdd[i].invsax, &cur_pid);
            if(cur_pid == -1){
                auto node  = new TARGNode();
                node->layer = parent->layer + 1;
                node->pid = pid_factory - 1;
                node->invSAX = saxTsPairRdd[i].invsax;
                node->ancestor = parent;
                parent->children[SaxUtil::invSaxHeadKFromInvSax(saxTsPairRdd[i].invsax, node->layer)] = node;
                cur_pid = node->pid;
            }
            saxTsIdRdd[cur_pid].push_back(&saxTsPairRdd[i]);
        }

        auto start = chrono::system_clock::now();
        for(auto &partition: saxTsIdRdd){
            string out_file_name = Const::tardisfn + to_string(partition.first);
            FILE *outf = fopen(out_file_name.c_str(), "ab");
            for(auto &unit:partition.second){
                fwrite(&unit->id, sizeof(long ), 1, outf);
                fwrite(unit->ts, sizeof(float), Const::tsLength, outf);
                SaxUtil::string2invSAX(unit->invsax, invSax);
                fwrite(invSax, sizeof(unsigned ), Const::bitsCardinality, outf);
            }
            fclose(outf);
        }
        auto end = chrono::system_clock::now();
        TAR_WRITE_TIME += chrono::duration_cast<chrono::microseconds>(end- start).count();

        delete[] tss;
        cur += num;
        Const::logPrint("Pack Partition, progress" + to_string((double)cur / series_num));

    }
}

void TARGNode::getIndexStats(){
    int total_leaf_node_num = getLeafNodeNum();
    int total_size = getTotalSize();
    cout << "Total size = " << total_size << endl;
    cout <<"Total nodes number = " << getNodeNum() << endl;
    cout << "Leaf node number = " << total_leaf_node_num << endl;
//    cout << "1st layer node number = " << get1stLayerNodesNo() <<endl;
//    cout << "1st layer internal node number = " << get1stLayerInterNodesNo() << endl;
//    cout << "1st layer internal series number = " << get1stLayerInterNodeSeriesNo() << endl;
    cout << "Max. height = " << getMaxHeight()  <<endl;
    cout << "Avg. Height = " << getAvgHeight() / (double) total_leaf_node_num << endl;
    cout <<"Avg. Filling Factor = "<< total_size / (double)total_leaf_node_num / Const::th << endl;
//    cout << "Bias leaf node ratio = " << (double)getBiasLeafNodeNum() / total_leaf_node_num << endl;
}

int TARGNode::getLeafNodeNum(){
    int sum = 0;
    vector<string>file_names;
    FileUtil::getFiles(Const::tardisfn, file_names);
    for(auto &file_name: file_names){
        if(file_name.substr(file_name.size() - 8) == "root.idx") continue;
        TARLNode* l_root = TARLNode::loadFromDisk(file_name);
        sum += l_root->getLeafNodeNbr();

        l_root->deleteDescendants();
        delete l_root;
    }
    return sum;
}

int TARGNode::getTotalSize(){
    int sum = 0;
    vector<string>file_names;
    FileUtil::getFiles(Const::tardisfn, file_names);
    for(auto &file_name: file_names){
        if(file_name.substr(file_name.size() - 8) == "root.idx" ) continue;
        TARLNode* l_root = TARLNode::loadFromDisk(file_name);
        sum += l_root->rcdNbr;

        l_root->deleteDescendants();
        delete l_root;
    }
    return sum;
}

int TARGNode::getNodeNum(){
    int sum = 0;
    vector<string>file_names;
    FileUtil::getFiles(Const::tardisfn, file_names);
    for(auto &file_name: file_names){
        if(file_name.substr(file_name.size() - 8) == "root.idx") continue;
        TARLNode* l_root = TARLNode::loadFromDisk(file_name);
        sum += l_root->getNodeNbr();

        l_root->deleteDescendants();
        delete l_root;
    }

    return sum + getGlobalNodeNum();
}

int TARGNode::getMaxHeight(){
    vector<int>local_heights;
    vector<string>file_names;
    FileUtil::getFiles(Const::tardisfn, file_names);
    local_heights.resize(file_names.size() - 1, 0);
    for(auto &file_name: file_names){
        if(file_name.substr(file_name.size() - 8) == "root.idx") continue;
        TARLNode* l_root = TARLNode::loadFromDisk(file_name);
        int pos = file_name.rfind('/');
        int _pid = atoi(file_name.substr(pos).c_str());

        local_heights[_pid] = l_root->getMaxHeight();
        l_root->deleteDescendants();
        delete l_root;
    }

    return getMaxHeightGlobal(local_heights);
}

int TARGNode::getMaxHeightGlobal(vector<int>&local_heights) {
    if(children.empty())    return local_heights[pid];
    int ret = 0;
    for(auto &iter: children){
        ret = max(ret, iter.second->getMaxHeightGlobal(local_heights));
    }
    return  ret + 1;
}

long TARGNode::getAvgHeight(){
    vector<int>local_heights;
    vector<int>leaf_nbr;
    vector<bool>visited;
    vector<string>file_names;
    FileUtil::getFiles(Const::tardisfn, file_names);
    local_heights.resize(file_names.size() - 1, 0);
    leaf_nbr.resize(file_names.size() - 1, 0);
    visited.resize(file_names.size() - 1, false);
    for(auto &file_name: file_names){
        if(file_name.substr(file_name.size() - 8) == "root.idx") continue;
        TARLNode* l_root = TARLNode::loadFromDisk(file_name);
        int pos = file_name.rfind('/');
        int _pid = atoi(file_name.substr(pos).c_str());
        local_heights[_pid] = l_root->getSumHeight();
        leaf_nbr[_pid] = l_root->getLeafNodeNbr();

        l_root->deleteDescendants();
        delete l_root;
    }

    return getAvgHeightGlobal(local_heights,leaf_nbr);
}

long TARGNode::getAvgHeightGlobal(vector<int>&local_heights, vector<int>&leaf_nbr) {
    if(children.empty()){
        if(leaf_nbr[pid] > 0) {
            long to_ret = local_heights[pid] + leaf_nbr[pid] * layer;
            leaf_nbr[pid] = -1;
            return to_ret;
        }
        else
            return 0 ;
    }
    long ret = 0;
    for(auto &iter: children){
        ret += iter.second->getAvgHeightGlobal(local_heights, leaf_nbr);
    }
    return  ret;
}

int TARGNode::getGlobalNodeNum(){
    if(children.empty())    return 0;
    int sum = 1;
    for(auto &iter: children){
        sum += iter.second->getGlobalNodeNum();
    }
    return sum;
}

//void TARGNode::assignLocalRootForAllLeaves(vector<TARLNode*>&local_roots){
//    if(children.empty()){
//        local_root = local_roots[pid];
//    }else{
//        for(auto &iter: children){
//            iter.second->assignLocalRootForAllLeaves(local_roots);
//        }
//    }
//}

// return pid if a leaf pack found, otherwise the deepest internal node
TARGNode * TARGNode::route(string& invsax, int *ret_pid){
    if(children.empty()) {
        *ret_pid = pid;
        return this;
    }
    string key = SaxUtil::invSaxHeadKFromInvSax(invsax, layer + 1);
    if(children.contains(key)){
        return children[key]->route(invsax, ret_pid);
    }else{
        *ret_pid = -1;
        return this;
    }
}

void TARGNode::assignPIdSampling(TARGNode* node){
    if(node->layer > 0 && node->size <= Const::th){
        node->pid = pid_factory++;
        node->children.clear();
    }else{
        vector<TARGNode*>tless;
        for(auto &iter: node->children){
            if(iter.second->size <= Const::th)
                tless.push_back(iter.second);
        }
        sort(tless.begin(), tless.end(), TARGNode::order);
        int cur_size = 0;
        for(TARGNode*_:tless){
            if(cur_size + _->size <= block_threshold){
                _->pid = pid_factory;
                cur_size += _->size;
                _->children.clear();
            }else{
                ++pid_factory;
                cur_size = _->size;
                _->pid = pid_factory;
                _->children.clear();
            }
        }

        ++pid_factory;

        for(auto &iter: node->children){
            if(iter.second->size > Const::th)
                assignPIdSampling(iter.second);
        }
    }
}

void TARGNode::loadCollectionSampling(vector<unordered_map<string, int>*>* node_list, TARGNode* root){
    for(int tree_layer = 0; tree_layer< node_list->size();++tree_layer){
        for(auto &iter: *(*node_list)[tree_layer]){
            auto* node = new TARGNode();
            node->layer = tree_layer + 1;
            node->size = iter.second / Const::tardis_sample_percent;
            node->invSAX = iter.first;
            node->addTreeNodeHexRobust(root);
        }
    }
}

void TARGNode::addTreeNodeHexRobust(TARGNode* start){
    if(layer == start->layer + 1){
        ancestor = start;
        start->children[invSAX] = this;
    }else{
        if(layer > start->layer){
            if(!start->children.contains(SaxUtil::invSaxHeadKFromInvSax(invSAX, start->layer + 1))){
                string newkey = SaxUtil::invSaxHeadKFromInvSax(invSAX, start->layer + 1);
                auto newNode = new TARGNode();
                newNode->layer = start->layer + 1;
                newNode->invSAX = newkey;
                newNode->size = size;
                newNode->ancestor = start;
                start->children[newkey] = newNode;
            }
            string key = SaxUtil::invSaxHeadKFromInvSax(invSAX, start->layer + 1);
            addTreeNodeHexRobust(start->children[key]);
        }else{
            Const::logPrint("add node error");
            exit(-2);
        }
    }
}

unsigned * TARGNode::loadSampling(){
    FILE *f  = fopen(Const::datafn.c_str(), "rb");
    long series_num = FileUtil::getFileSize(Const::datafn.c_str()) / Const::tsLengthBytes;
    long load_num = series_num * Const::tardis_sample_percent;
    auto ret = new unsigned[load_num * Const::bitsCardinality];
    long rest = series_num;
    long cur = 0;
    while (rest > 0){
        long num, take_num;
        if(rest > Const::fbl_series_num) {
            num = Const::fbl_series_num;
            take_num = num * Const::tardis_sample_percent;
        }
        else {
            num = rest;
            take_num = load_num - cur;
        }

        auto *tss = new float[take_num * Const::tsLength];
        auto start = chrono::system_clock::now();
        fseek(f,(series_num - rest) * Const::tsLengthBytes,SEEK_SET);
        fread(tss, sizeof(float),take_num * Const::tsLength,  f);
        auto end = chrono::system_clock::now();
        TAR_READ_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();

        for(long i=0;i<take_num;++i){
            auto _ = SaxUtil::saxFromTs(tss + i *Const::tsLength, Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
            SaxUtil::invSaxFromSax(_, ret + (cur + i) * Const::bitsCardinality, Const::bitsCardinality, Const::segmentNum);
            delete _;
        }

        cur+=take_num;
        Const::logPrint(to_string(cur) + " series have been read.");
        rest-=num;
        delete[] tss;
    }
    fclose(f);
    assert(cur <= load_num);
    return ret;
}

void TARGNode::computeFreqTbl(unordered_map<string, int> *invsax_freq,
                              unsigned *sample_invSAX_tbl, long tbl_size){
    for(long i =0;i<tbl_size;++i){
        string invsax = SaxUtil::invSax2String(sample_invSAX_tbl + i * Const::bitsCardinality);
        (*invsax_freq)[invsax]++;
    }
}


void TARGNode::computeFreqTblSpecLayer(unordered_map<string, int>*pairRDD, unordered_map<string,dataRDDItem>*dataRDD,
                             unordered_map<string,int>*saxNbrPairRDD, int layer){
    for(auto &iter:*pairRDD){
        string key = iter.first;
        string new_key = SaxUtil::invSaxHeadKFromInvSax(key, layer);
        dataRDD->emplace(key, dataRDDItem(new_key, iter.second));
        (*saxNbrPairRDD)[new_key] += iter.second;
    }
}

vector<unordered_map<string, int>*>* TARGNode::buildTARGbySampling(unordered_map<string, int>*invsax_freq){
    // invsax_freq : invsax(full), freq
    int tree_layer = 1;
    bool reachBottom = false;
    auto *output = new vector<unordered_map<string, int>*>();
    unordered_map<string, int>*pairRDD = invsax_freq;

    while(!reachBottom){
        auto dataRDD = new unordered_map<string,dataRDDItem>();
        auto saxNbrPairRDD = new unordered_map<string,int>();
        // pairRDD: invsax(full), freq
        // dataRDD: invSAX(full), invSAX(b), freq
        // saxNbrPairRDD: invsax(b), freq
        computeFreqTblSpecLayer(pairRDD, dataRDD, saxNbrPairRDD, tree_layer);
        output->push_back(saxNbrPairRDD);

        int max_freq  =0;
        for(auto &iter:*saxNbrPairRDD)
            max_freq = max(max_freq, iter.second);
        if(max_freq / Const::tardis_sample_percent <= Const::th) {
            reachBottom = true;
            delete pairRDD;
            delete dataRDD;
        }
        else{
            // smallHashSet: invSAX(full)
            unordered_set<string>smallHashSet;
            for(auto &iter: *saxNbrPairRDD){
                if(iter.second / Const::tardis_sample_percent <= Const::th)
                    smallHashSet.insert(iter.first);
            }
            delete pairRDD;
            pairRDD = new unordered_map<string, int>();
            for(auto &iter: *dataRDD){
                if(!smallHashSet.contains(iter.second.invsax_full)){
                    (*pairRDD)[iter.first] = iter.second.freq;
                }
            }
            delete dataRDD;
            tree_layer++;
        }
    }
    return output;
}
