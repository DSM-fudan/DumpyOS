//
// Created by wzy on 2021/11/8.
//

#include "../../include/DataStructures/IPGNode.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Utils/SaxUtil.h"
#include "../../include/Utils/TimeSeriesUtil.h"
#include "../../include/DataStructures/TimeSeries.h"
#include "../../include/DataStructures/IPGPartition.h"
#include "../../include/Const.h"
#include <iostream>
#include <set>
#include <algorithm>
#include <fstream>
#include <cmath>
const int IPGNode::power_2[]{1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536};
const int IPGNode::mask_16[]{
    0b1000000000000000,
    0b0100000000000000,
    0b0010000000000000,
    0b0001000000000000,
    0b0000100000000000,
    0b0000010000000000,
    0b0000001000000000,
    0b0000000100000000,
    0b0000000010000000,
    0b0000000001000000,
    0b0000000000100000,
    0b0000000000010000,
    0b0000000000001000,
    0b0000000000000100,
    0b0000000000000010,
    0b0000000000000001,
};

int*** IPGNode::combines = readCombines();

unsigned short*IPGNode::saxes = nullptr;
float *IPGNode::paas = nullptr;
string IPGNode::rowDataFileName = "";
const int combine_num[] = {0, 16, 120, 560, 1820, 4368, 8008, 11440, 12870, 11440, 8008, 4368, 1820, 560, 12, 16};
static int fuzzy_num = 0;

int *** IPGNode::readCombines(){
    string base = "../combines/" + to_string(Const::segmentNum) + "-";
    auto ret = new int**[15];
    for(int i=1;i<=14;++i){
        ret[i] = new int*[combine_num[i]];
        ifstream f(base + to_string(i) + ".txt", ios::in);
        for(int j=0;j<combine_num[i];++j){
            ret[i][j] = new int[i];
            for(int k=0;k<i;++k)
                f >> ret[i][j][k];
        }
        f.close();
    }
    return ret;
}

IPGNode * IPGNode::BuildIPG(string &saxfn, string &paafn, vector<vector<int>> *g, const string &method) {

    long series_num = loadSax(saxfn);
    loadPaa(paafn);
//    long series_num = generateSaxAndPaaTbl();
    auto* root = new IPGNode();
    root->children = new unordered_map<int, IPGNode*>();
    root->size = series_num;
    for(int &i:root->bits_cardinality)  i=0;
    // the same 1st layer
    cout << "start build first layer"<<endl;
    root->build1stLayer(series_num);

    cout << "start partition first layer." << endl;
    IPGPartition::growingGraph(root->children, g, Const::filling_factor_1st, Const::segmentNum);
//    GraphPartitioner::doPartitionStatic(root->children);

    cout << "start stat paa in first layer:" << root->children->size()<< endl;
    int j=0;
    for(auto &iter:*root->children){
        if(iter.second->size > Const::th)
            iter.second->statPaa1stLayer();
        if(++j%10000==0) cout << j << endl;
    }

    j = 0;
    cout << "start grow the index structure" << endl;
    for(auto &iter:*root->children){
        if(iter.second->size > Const::th){
            iter.second->growIPGIndex(method);
        }
        if(++j%10000 == 0)  cout << j << endl;
    }

    cout << "build index finished." << endl;

    return root;
}

void IPGNode::growIPGIndex(const string &method) {
    if(size <= Const::th)   return;
    int chosen_num = SaxUtil::findFirstGE(power_2, 1, Const::segmentNum + 1, size / Const::th + 1);
    vector<int>chosen_segments = chooseSegment(this, chosen_num, method);
    if(children == nullptr) children = new unordered_map<int, IPGNode*>();
    for(int i=0;i<actual_size;++i){
        int new_id = SaxUtil::extendSax(IPGNode::paas + (long)offsets[i] * (Const::segmentNum), bits_cardinality,chosen_segments);
        IPGNode* node = (*children)[new_id];
        if(node == nullptr) {
            node = new IPGNode(new_id, layer + 1);
            (*children)[new_id] = node;
        }
        node->insert(offsets[i], true);
    }

    // resolve boundary nodes
    int sw, cardinality,offset, res = 0;
    float *paa;
    for(int i=actual_size;i<size;++i){
        paa = IPGNode::paas + (long)offsets[i] * Const::segmentNum;
        res = 0;
        for(int segment:chosen_segments){
            cardinality = 1 << (bits_cardinality[segment] + 1);
            offset = ((cardinality - 1) * (cardinality - 2)) / 2;
            int index = SaxUtil::findFirstGE(SaxUtil::breakpoints, offset, cardinality -1, paa[segment]);
            if(index >= 0)  sw = (index - offset);
            else    cout<<"ERROR!!!!!!!";

            if(sw >> 1 == sax[segment]) res = (res << 1) + (sw % 2);
            else{
                if(sw >> 1 > sax[segment])  res = (res << 1) + 1;
                else    res <<= 1;
            }
        }
        IPGNode* node = (*children)[res];
        if(node == nullptr) {
            node = new IPGNode(res, layer + 1);
            (*children)[res] = node;
        }
        (*children)[res]->insert(offsets[i], false);
    }

    // generate sax and bits cardinality for each segment
    for(auto &iter: *children){
        if(iter.second != nullptr) {
            generateSaxAndCardinality(iter.second, iter.first, chosen_segments);
        }
    }

    vector<int>().swap(offsets);

    //partition
    IPGPartition::balanceMatch(children);
    //GraphPartitioner::doPartitionStatic(children);

    // build rest data node if any
    int pid = children->size();
    for(auto &iter:*children){
        if(iter.second->actual_size <= Const::th && iter.second->partitionId == -1){
            assert(iter.second->dataNode == nullptr);
            buildDataNode(iter.second);
            iter.second->partitionId = ++pid;
        }
    }

    // stat paa for internal node
    for(auto &iter:*children){
        if(iter.second->partitionId == -1)
            iter.second->statPaa();
    }

    for(auto &iter: *children){
        if(iter.second->size > Const::th){
//            assert(iter.second->partitionId == -1);
//            assert(iter.second->dataNode == nullptr);
//            assert(iter.second->actual_size > Const::tardis_threshold);
            iter.second->growIPGIndex(method);
        }
    }

}

IPGNode * IPGNode::BuildIPGFuzzy(string &saxfn, string &paafn, vector<vector<int>> *g, const string &method) {
    int series_num = loadSax(saxfn);
    loadPaa(paafn);
    auto* root = new IPGNode();
    root->children = new unordered_map<int, IPGNode*>();
    root->size = series_num;
    for(int &i:root->bits_cardinality)  i=0;
    // the same 1st layer
    cout << "start build first layer"<<endl;
    root->build1stLayer(series_num);

    cout << "start partition first layer." << endl;
    IPGPartition::growingGraph(root->children, g, Const::filling_factor_1st, Const::segmentNum);
//    GraphPartitioner::doPartitionStatic(root->children);

    cout << "start stat paa in first layer:" << root->children->size()<< endl;
    int j=0;
    for(auto &iter:*root->children){
        if(iter.second->size > Const::th)
            iter.second->statPaa1stLayer();
        else
            iter.second->statPaaMean();
        if(++j%10000==0) cout << j << endl;
    }

    cout << "start replica first layer" << endl;
    fuzzyFirstLayer(root);
    cout << "replica number = " << fuzzy_num << endl;

    j = 0;
    cout << "start grow the index structure" << endl;
    for(auto &iter:*root->children){
        if(iter.second->size > Const::th){
            iter.second->growIPGIndexFuzzy(g, method);
        }
        if(++j%10000 == 0)  cout << j << endl;
    }

    cout << "total replica number = " << fuzzy_num << endl;
    cout << "build index finished." << endl;

    return root;
}

void IPGNode::growIPGIndexFuzzy(vector<vector<int>> *g, const string &method) {
    if(size <= Const::th)   return;
    int chosen_num = SaxUtil::findFirstGE(power_2, 1, Const::segmentNum + 1, size / Const::th + 1);
    vector<int>chosen_segments = chooseSegment(this, chosen_num, method);
    if(children == nullptr) children = new unordered_map<int, IPGNode*>();
    for(int i=0;i<actual_size;++i){
        int new_id = SaxUtil::extendSax(IPGNode::paas + offsets[i] * (Const::segmentNum), bits_cardinality,chosen_segments);
        IPGNode* node = (*children)[new_id];
        if(node == nullptr) {
            node = new IPGNode(new_id, layer + 1);
            (*children)[new_id] = node;
        }
        node->insert(offsets[i], true);
    }

    // resolve boundary series
    int sw, cardinality,offset, res = 0;
    float *paa;
    for(int i=actual_size;i<size;++i){
        paa = IPGNode::paas + offsets[i] * Const::segmentNum;
        res = 0;
        for(int segment:chosen_segments){
            cardinality = 1 << (bits_cardinality[segment] + 1);
            offset = ((cardinality - 1) * (cardinality - 2)) / 2;
            int index = SaxUtil::findFirstGE(SaxUtil::breakpoints, offset, cardinality -1, paa[segment]);
            if(index >= 0)  sw = (index - offset);
            else    cout<<"ERROR!!!!!!!";

            if(sw >> 1 == sax[segment]) res = (res << 1) + (sw % 2);
            else{
                if(sw >> 1 > sax[segment])  res = (res << 1) + 1;
                else    res <<= 1;
            }
        }
        IPGNode* node = (*children)[res];
        if(node == nullptr) {
            node = new IPGNode(res, layer + 1);
            (*children)[res] = node;
        }
        (*children)[res]->insert(offsets[i], false);
    }

    // generate sax and bits cardinality for each child node
    for(auto &iter: *children){
        if(iter.second != nullptr) {
            generateSaxAndCardinality(iter.second, iter.first, chosen_segments);
        }
    }

    vector<int>().swap(offsets);

    //partition
    IPGPartition::growingGraph(children, g, Const::filling_factor, chosen_segments.size());
    //GraphPartitioner::doPartitionStatic(children);

    int pid = children->size();
    for(auto &iter:*children){
        if(iter.second->actual_size <= Const::th && iter.second->partitionId == -1){
            buildDataNode(iter.second);
            iter.second->partitionId = ++pid;
        }
    }

    for(auto &iter:*children){
        if(iter.second->partitionId == -1)
            iter.second->statPaa();
    }

    //replica
    fuzzy(this, chosen_segments);

    for(auto &iter: *children){
        if(iter.second->size > Const::th){
//            assert(iter.second->partitionId == -1);
//            assert(iter.second->dataNode == nullptr);
//            assert(iter.second->actual_size > Const::th);
            iter.second->growIPGIndexFuzzy(g, method);
        }
    }

}

IPGNode * IPGNode::BuildIPGGrid(string &saxfn, string &paafn, vector<vector<int>> *g) {
    int series_num = loadSax(saxfn);
    loadPaa(paafn);
    auto* root = new IPGNode();
    root->children = new unordered_map<int, IPGNode*>();
    root->size = series_num;
    for(int &i:root->bits_cardinality)  i=0;
    // the same 1st layer
    cout << "start build first layer"<<endl;
    root->build1stLayer(series_num);

    cout << "start partition first layer." << endl;
    IPGPartition::growingGraph(root->children, g, Const::filling_factor_1st, Const::segmentNum);
//    GraphPartitioner::doPartitionStatic(root->children);

    int j = 0;
    cout << "start grow the index structure" << endl;
    for(auto &iter:*root->children){
        if(iter.second->size > Const::th){
            iter.second->growIPGIndexGrid(g);
        }
        if(++j%100 == 0)  cout << j << endl;
    }

    cout << "build index finished." << endl;

    return root;
}

void IPGNode::growIPGIndexGrid(vector<vector<int>> *g) {
    if(size <= Const::th)   return;
    int chosen_num = SaxUtil::findFirstGE(power_2, 1, Const::segmentNum + 1, size / Const::th + 1);
    if(children == nullptr) children = new unordered_map<int, IPGNode*>();
    if(sax_cnt == -1)   countDistinctSaxWords();

    int plan_num = combine_num[chosen_num];

    double min_score = numeric_limits<double>::max();
    int min_plan_index = -1;
    for(int i=0;i<plan_num;++i){
        int* plan = combines[chosen_num][i];
        for(int j=0; j < (1<<chosen_num); ++j){
            auto* node = new IPGNode(j, -1);
            (*children)[j] = node;
        }

        trySplitting2Children(plan, chosen_num);
        int totalPartition = IPGPartition::growingGraph(children, g, Const::filling_factor,  chosen_num);
        int num = 0; double sum_invert_avg = 0;
        for(auto &iter:*children){
            if(iter.second->size > Const::th || iter.second->partitionId == -1)
            {
                iter.second->countDistinctSaxWords();
                ++num;
                sum_invert_avg += (1.0 / iter.second->sax_cnt);
            }
        }
        for(int j=0; j<totalPartition;++j){
            int res = countDistinctSaxWords(j);
            if(res > 0){
                sum_invert_avg += (1.0 / res);
                ++num;
            }
        }
        double score = num / sum_invert_avg;
        if(score < min_score){
            min_score = score;
            min_plan_index = i;
        }

        for(auto &child: *children)
            delete child.second;
        children->clear();

    }

    chosenSegments.resize(chosen_num);
    copy(combines[chosen_num][min_plan_index], combines[chosen_num][min_plan_index] + chosen_num, chosenSegments.begin());

    for(int i=0;i<size;++i){
        int new_id = SaxUtil::extendSax(IPGNode::paas + offsets[i] * (Const::segmentNum), bits_cardinality,chosenSegments);
        IPGNode* node = (*children)[new_id];
        if(node == nullptr) {
            node = new IPGNode(new_id, layer + 1);
            (*children)[new_id] = node;
        }
        node->insert(offsets[i], true);
    }

    {
        //    // resolve boundary nodes

//    int sw, cardinality,offset, res = 0;
//    double *paa;
//    for(int i=actual_size;i<size;++i){
//        paa = IPGNode::paas + offsets[i] * Const::segmentNum;
//        res = 0;
//        for(int segment:chosen_segments){
//            cardinality = 1 << (bits_cardinality[segment] + 1);
//            offset = ((cardinality - 1) * (cardinality - 2)) / 2;
//            int index = SaxUtil::findFirstGE(SaxUtil::breakpoints, offset, cardinality -1, paa[segment]);
//            if(index >= 0)  sw = (index - offset);
//            else    cout<<"ERROR!!!!!!!";
//
//            if(sw >> 1 == sax[segment]) res = (res << 1) + (sw % 2);
//            else{
//                if(sw >> 1 > sax[segment])  res = (res << 1) + 1;
//                else    res <<= 1;
//            }
//        }
//        IPGNode* node = (*children)[res];
//        if(node == nullptr) {
//            node = new IPGNode(res, layer + 1);
//            (*children)[res] = node;
//        }
//        (*children)[res]->insert(offsets[i], false);
//    }
//
//    //GraphPartitioner::doPartitionStatic(children);
//
}

    // generate sax and bits cardinality for each segment
    for(auto &iter: *children){
        if(iter.second != nullptr) {
            generateSaxAndCardinality(iter.second, iter.first, chosenSegments);
        }
    }

    vector<int>().swap(offsets);

    //partition
    IPGPartition::growingGraph(children, g, Const::filling_factor, chosen_num);

    // build rest data node if any
    int pid = children->size();
    for(auto &iter:*children){
        if(iter.second->actual_size <= Const::th && iter.second->partitionId == -1){
            assert(iter.second->dataNode == nullptr);
            buildDataNode(iter.second);
            iter.second->partitionId = ++pid;
        }
    }
    for(auto &iter: *children){
        if(iter.second->size > Const::th){
            iter.second->growIPGIndexGrid(g);
        }
    }

}

IPGNode * IPGNode::BuildIPGCluster(string &saxfn, string &paafn, vector<vector<int>> *g) {
    int series_num = loadSax(saxfn);
    loadPaa(paafn);
    auto* root = new IPGNode();
    root->children = new unordered_map<int, IPGNode*>();
    root->size = series_num;
    for(int &i:root->bits_cardinality)  i=0;
    // the same 1st layer
    cout << "start build first layer"<<endl;
    root->build1stLayer(series_num);

    cout << "start partition first layer." << endl;
    IPGPartition::growingGraph(root->children, g, Const::filling_factor_1st, Const::segmentNum);
//    GraphPartitioner::doPartitionStatic(root->children);

    int j = 0;
    Const::logPrint("start grow the index structure");
    for(auto &iter:*root->children){
        if(iter.second->size > Const::th ){
            iter.second->growIPGIndexCluster();
        }
        if(++j%10000 == 0)  cout << j << endl;
    }
    Const::logPrint("build index finished.");

    return root;
}

struct HYP_NODE{
    int size;
    double sum_paa[Const::segmentNum];

    HYP_NODE(){size = 0; for(auto &i:sum_paa)   i=0;}
};

void IPGNode::growIPGIndexCluster() {
    if(size <= Const::th)   return;
    int chosen_num = SaxUtil::findFirstGE(power_2, 1, Const::segmentNum + 1, size / Const::th + 1);
    if(children == nullptr) children = new unordered_map<int, IPGNode*>();
    vector<int>chosen_segments;
    if(chosen_num < 1) {
        statPaaOne();
        chosen_segments.push_back(chooseOneSegment(this));
    }
    else{
        statPaaMaxMinMean();
        double lb, ub;  // ub is the new split line
        for(int i=0;i<Const::segmentNum;++i){
            if(bits_cardinality[i] >= Const::bitsCardinality)   continue;
            SaxUtil::getValueRange(sax[i] << 1, bits_cardinality[i] + 1, &lb, &ub);
            if(paa_max[i] <= ub || paa_min[i] >= ub){
                bits_cardinality[i]++;
                sax[i] <<= 1;
            }
        }

        int bc[Const::segmentNum];
        copy(bits_cardinality, bits_cardinality + Const::segmentNum, bc);
        sort(bc, bc + Const::segmentNum);
        int sum = 0;
        for(int & i:bits_cardinality)   sum += i;
        sum /= Const::segmentNum;
        int min_bc = bc[0], max_bc = bc[Const::segmentNum - 1], range =  max_bc - min_bc;
//        double bc_baseline = 1.0;   // lower bits preference
//        for(int i=Const::segmentNum - 1; i>=Const::segmentNum - chosen_num;--i)
//            bc_baseline += bc[i];

        // this flag controls whether we omit plans containing certain segment
        bool  flag = (range >= Const::max_diff);

        vector<HYP_NODE>hyp_nodes(Const::vertexNum);
        // split to all possible hypothesis nodes
        for(int offset:offsets){
            int nav_id = SaxUtil::extendSax(saxes + offset * Const::segmentNum, bits_cardinality);
            hyp_nodes[nav_id].size++;
            float *paa_start = paas + offset * Const::segmentNum;
            for(int i=0;i<Const::segmentNum;++i)
                hyp_nodes[nav_id].sum_paa[i]+=paa_start[i];
        }

        int plan_num = combine_num[chosen_num];
        double max_score = -numeric_limits<double>::max();
        int max_plan_index = -1, new_nodes_num = power_2[chosen_num];
        for(int i=0;i<plan_num;++i){
            int* plan = combines[chosen_num][i];
            if(flag) {
                int j = 0;
                for (; j < chosen_num; ++j) {
                    if(bits_cardinality[plan[j]] >= sum + Const::max_diff)
                        break;
                }
                if(j < chosen_num)
                    continue;
            }
            vector<int>sizes(new_nodes_num, 0);
//            vector<vector<double> >centroids(new_nodes_num, vector<double>(chosen_num, 0));
            vector<vector<double> >centroids(new_nodes_num, vector<double>(Const::segmentNum, 0));

            for(int j=0;j<Const::vertexNum;++j){
                int belong_id = 0;
                for(int k=0;k<chosen_num;++k){
                    int val = j & power_2[Const::segmentNum - 1 - plan[k]];
                    if(val) belong_id |= power_2[chosen_num - 1 - k];
                }
                if(hyp_nodes[j].size > 0){
                    sizes[belong_id] += hyp_nodes[j].size;
//                    for(int k=0;k<chosen_num;++k)
//                        centroids[belong_id][k] += hyp_nodes[j].sum_paa[plan[k]];
                    for(int k=0;k<Const::segmentNum;++k)
                        centroids[belong_id][k] += hyp_nodes[j].sum_paa[k];

                }
            }
            for(int j=0;j<new_nodes_num;++j)
//                for(int k=0;k<chosen_num;++k)
//                    if(sizes[j] > 0)
//                        centroids[j][k] /= sizes[j];
                for(int k=0;k<Const::segmentNum;++k)
                    if(sizes[j] > 0)
                        centroids[j][k] /= sizes[j];

            double score = 0;
            for(int j=0;j<new_nodes_num;++j){
                if(sizes[j] <=0)    continue;
                double tmp_sum = 0;
//                max_size = sizes[j] > max_size ? sizes[j] : max_size;
//                for(int k=0;k<chosen_num;++k){
//                    double _ = centroids[j][k] - paa_mu[plan[k]];
//                    tmp_sum += (_ * _);
//                }
                for(int k=0;k<Const::segmentNum;++k){
                    double _ = centroids[j][k] - paa_mu[k];
                    tmp_sum += (_ * _);
                }

                score += (sizes[j] * tmp_sum);
            }

//            int sum_bits = 0;
//            for(int j=0;j<chosen_num;++j)
//                sum_bits += bits_cardinality[plan[j]];
//            score *= (1 - (sum_bits / bc_baseline));
//            // balance constraints has been considered in SSB equation
//            score *= (1- (max_size / (double)size));


            if(score > max_score){
                max_score = score;
                max_plan_index = i;
            }

        }

        int* plan = combines[chosen_num][max_plan_index];
        for(int j=0;j<chosen_num;++j)
            chosen_segments.push_back(plan[j]);
        vector<HYP_NODE>().swap(hyp_nodes);
    }


    for(int i=0;i<size;++i){
        int offset = offsets[i];
        int new_id = SaxUtil::extendSax(saxes + offset * Const::segmentNum, bits_cardinality,chosen_segments);
        IPGNode* node = (*children)[new_id];
        if(node == nullptr) {
            node = new IPGNode(new_id, layer + 1);
            (*children)[new_id] = node;
        }
        node->insert(offset, true);
    }

    // generate sax and bits cardinality for each segment
    for(auto &iter: *children){
        if(iter.second != nullptr) {
            generateSaxAndCardinality(iter.second, iter.first, chosen_segments);
        }
    }

    vector<int>().swap(offsets);

    //partition
    IPGPartition::balanceMatch(children);

    // build rest data node if any
    int pid = children->size();
    for(auto &iter:*children){
        if(iter.second->actual_size <= Const::th && iter.second->partitionId == -1){
            assert(iter.second->dataNode == nullptr);
            buildDataNode(iter.second);
            iter.second->partitionId = ++pid;
        }
    }
    for(auto &iter: *children){
        if(iter.second->size > Const::th){
            iter.second->growIPGIndexCluster();
        }
    }

}

IPGNode * IPGNode::BuildIPGClusterDynamic(string &saxfn, string &paafn, vector<vector<int>> *g) {
    int series_num = loadSax(saxfn);
    loadPaa(paafn);
    auto* root = new IPGNode();
    root->children = new unordered_map<int, IPGNode*>();
    root->size = series_num;
    for(int &i:root->bits_cardinality)  i=0;
    // the same 1st layer
    cout << "start build first layer"<<endl;
    root->build1stLayer(series_num);

    cout << "start partition first layer." << endl;
    IPGPartition::growingGraph(root->children, g, Const::filling_factor_1st, Const::segmentNum);
//    GraphPartitioner::doPartitionStatic(root->children);

    int j = 0;
    Const::logPrint("start grow the index structure");
    for(auto &iter:*root->children){
        if(iter.second->size > Const::th){
            iter.second->growIPGIndexClusterDynamic();
        }
        if(++j%10000 == 0)  cout << j << endl;
    }
    Const::logPrint("build index finished.");

    return root;
}

void IPGNode::growIPGIndexClusterDynamic() {
    if(size <= Const::th)   return;
    int chosen_num = SaxUtil::findFirstGE(power_2, 1, Const::segmentNum + 1, size / Const::th + 1);
    if(children == nullptr) children = new unordered_map<int, IPGNode*>();
    vector<int>chosen_segments;
    if(chosen_num < 1) {
        statPaaOne();
        chosen_segments.push_back(chooseOneSegment(this));
    }
    else{
        int bc[Const::segmentNum];
        copy(bits_cardinality, bits_cardinality + Const::segmentNum, bc);
        sort(bc, bc + Const::segmentNum);
        int min_bc = bc[0], max_bc = bc[Const::segmentNum - 1], range =  max_bc - min_bc;
        double bc_baseline = 1.0;   // lower bits preference
        for(int i=Const::segmentNum - 1; i>=Const::segmentNum - chosen_num;--i)
            bc_baseline += bc[i];

        // this flag controls whether we omit plans containing certain segment
        bool  flag = (range >= Const::max_diff);
        statPaaMean();
        vector<HYP_NODE>hyp_nodes(Const::vertexNum);
        // split to all possible hypothesis nodes
        for(int offset:offsets){
            int nav_id = SaxUtil::getNewId(paas + offset * Const::segmentNum, paa_mu);
            hyp_nodes[nav_id].size++;
            float *paa_start = paas + offset * Const::segmentNum;
            for(int i=0;i<Const::segmentNum;++i)
                hyp_nodes[nav_id].sum_paa[i]+=paa_start[i];
        }

        int plan_num = combine_num[chosen_num];
        double max_score = -numeric_limits<double>::max();
        int max_plan_index = -1, new_nodes_num = power_2[chosen_num];
        for(int i=0;i<plan_num;++i){
            int* plan = combines[chosen_num][i];
//            if(flag) {
//                int j = 0;
//                for (; j < chosen_num; ++j) {
//                    if(bits_cardinality[plan[j]] >= min_bc + Const::max_diff)
//                        break;
//                }
//                if(j < chosen_num)
//                    continue;
//            }
            vector<int>sizes(new_nodes_num, 0);
//            vector<vector<double> >centroids(new_nodes_num, vector<double>(chosen_num, 0));
            vector<vector<double> >centroids(new_nodes_num, vector<double>(Const::segmentNum, 0));

            for(int j=0;j<Const::vertexNum;++j){
                int belong_id = 0;
                for(int k=0;k<chosen_num;++k){
                    int val = j & power_2[Const::segmentNum - 1 - plan[k]];
                    if(val) belong_id |= power_2[chosen_num - 1 - k];
                }
                if(hyp_nodes[j].size > 0){
                    sizes[belong_id] += hyp_nodes[j].size;
//                    for(int k=0;k<chosen_num;++k)
//                        centroids[belong_id][k] += hyp_nodes[j].sum_paa[plan[k]];
                    for(int k=0;k<Const::segmentNum;++k)
                        centroids[belong_id][k] += hyp_nodes[j].sum_paa[k];

                }
            }
            for(int j=0;j<new_nodes_num;++j)
//                for(int k=0;k<chosen_num;++k)
//                    centroids[j][k] /= sizes[j];
                for(int k=0;k<Const::segmentNum;++k)
                    centroids[j][k] /= sizes[j];

            double score = 0;
            for(int j=0;j<new_nodes_num;++j){
                double tmp_sum = 0;
//                max_size = sizes[j] > max_size ? sizes[j] : max_size;
//                for(int k=0;k<chosen_num;++k){
//                    double _ = centroids[j][k] - paa_mu[plan[k]];
//                    tmp_sum += (_ * _);
//                }
                for(int k=0;k<Const::segmentNum;++k){
                    double _ = centroids[j][k] - paa_mu[k];
                    tmp_sum += (_ * _);
                }

                score += (sizes[j] * tmp_sum);
            }

            int sum_bits = 0;
            for(int j=0;j<chosen_num;++j)
                sum_bits += bits_cardinality[plan[j]];
            score *= (1 - (sum_bits / bc_baseline));
            // balance constraints has been considered in SSB equation
//            score *= (1- (max_size / (double)size));


            if(score > max_score){
                max_score = score;
                max_plan_index = i;
            }

        }

        int* plan = combines[chosen_num][max_plan_index];
        for(int j=0;j<chosen_num;++j)
            chosen_segments.push_back(plan[j]);
        vector<HYP_NODE>().swap(hyp_nodes);
    }

    int bc[Const::segmentNum];
    copy(bits_cardinality, bits_cardinality + Const::segmentNum, bc);
    for(int seg:chosen_segments)    bc[seg]++;
    for(int i=0;i<size;++i){
        int offset = offsets[i];
        int new_id = SaxUtil::getNewId(paas + offset * Const::segmentNum, paa_mu,chosen_segments);
        IPGNode* node = (*children)[new_id];
        if(node == nullptr) {
            node = new IPGNode(new_id, layer + 1);
            (*children)[new_id] = node;
            copy(bc, bc+Const::segmentNum, (*children)[new_id]->bits_cardinality);
        }
        node->insert(offset, true);
    }

    vector<int>().swap(offsets);

    //partition
    IPGPartition::balanceMatch(children);

    // build rest data node if any
    int pid = children->size();
    for(auto &iter:*children){
        if(iter.second->actual_size <= Const::th && iter.second->partitionId == -1){
            assert(iter.second->dataNode == nullptr);
            buildDataNode(iter.second);
            iter.second->partitionId = ++pid;
        }
    }
    for(auto &iter: *children){
        if(iter.second->size > Const::th){
            iter.second->growIPGIndexCluster();
        }
    }

}


void IPGNode::trySplitting2Children(const int *chosen_segments, int chosen_num){
    int split[chosen_num];
    for(int i=0;i<chosen_num;++i)
    {
        assert(bits_cardinality[chosen_segments[i]] < Const::bitsCardinality);
        split[i] = (sax[chosen_segments[i]] << 1) + 1;
        split[i] <<= (Const::bitsCardinality - bits_cardinality[chosen_segments[i]] - 1);
    }
    int new_id = 0;
    unsigned short*_sax;
    for(int offset:offsets){
        new_id = 0;
        _sax = saxes + offset * Const::segmentNum;
        for(int i=0;i<chosen_num;++i){
            new_id <<= 1;
            if(_sax[chosen_segments[i]] >= split[i])    ++new_id;
        }
        (*children)[new_id]->fakeInsert(offset);
    }
}

struct SaxHasher
{
    size_t operator()(const unsigned short *sax) const noexcept
    {
        size_t ret = 0;
        for(int i=0;i<Const::segmentNum;++i)
            ret += sax[i];
        return ret;
    }
};

struct SaxComparator
{
    bool operator()(const unsigned short *sax1, const unsigned short *sax2) const noexcept
    {
        for(int i=0;i<Const::segmentNum;++i)
            if(sax1[i] != sax2[i])
                return false;
        return true;
    }
};

void IPGNode::countDistinctSaxWords(){
    unordered_set<unsigned short *, SaxHasher, SaxComparator>sax_set;
    for(int offset:offsets){
        sax_set.insert(saxes + offset * Const::segmentNum);
    }
    sax_cnt = sax_set.size();
}

int IPGNode::countDistinctSaxWords(int pid) const{
    unordered_set<unsigned short *, SaxHasher, SaxComparator>sax_set;
    for(auto &iter:*children){
        if(iter.second->partitionId == pid){
            for(int offset:iter.second->offsets){
                sax_set.insert(saxes + offset * Const::segmentNum);
            }
        }
    }
    return sax_set.size();
}

IPGNode * IPGNode::BuildDynamicIPG(string &saxfn, string &paafn, vector<vector<int>> *g){
    int series_num = loadSax(saxfn);
    loadPaa(paafn);
    auto* root = new IPGNode();
    root->children = new unordered_map<int, IPGNode*>();
    root->size = series_num;
    for(int &i:root->bits_cardinality)  i=0;
    // the same 1st layer
    cout << "start build first layer"<<endl;
    root->build1stLayerDynamicMean(series_num);

    for(auto iter = (*root->children).begin(); iter != (*root->children).end();){
        if(iter->second->size <= 0)
            iter = (*root->children).erase(iter);
        else    ++iter;
    }

    cout << "nodes number in 1st layer: "<< (*root->children).size() << endl;

    cout << "start partition first layer." << endl;
    IPGPartition::growingGraph(root->children, g, Const::filling_factor_1st, Const::segmentNum);

    // stat paa for internal node
    for(auto &iter:*root->children){
        if(iter.second->partitionId == -1)
            iter.second->statPaaRange();
    }

    int j = 0;
    cout << "start grow the index structure" << endl;
    for(auto &iter:*root->children){
        if(iter.second->size > Const::th){
            iter.second->growDynamicIPGIndexRange(g);
        }
        if(++j%10000 == 0)  cout << j << endl;
    }

    cout << "build index finished." << endl;

    return root;
}

void IPGNode::growDynamicIPGIndex(vector<vector<int>> *g){
    if(size <= Const::th)   return;
    int chosen_num = SaxUtil::findFirstGE(power_2, 1, Const::segmentNum + 1, size / Const::th);
    chooseSegmentDynamic(this, chosen_num);

    if(children == nullptr) children = new unordered_map<int, IPGNode*>();
    for(int i=0;i<actual_size;++i){
        int new_id = getDynamicId(IPGNode::paas + offsets[i] * Const::segmentNum);
        IPGNode* node = (*children)[new_id];
        if(node == nullptr) {
            node = new IPGNode(new_id, layer + 1);
            (*children)[new_id] = node;
        }
        node->insert(offsets[i], true);
    }

    vector<int>().swap(offsets);

    //partition
    IPGPartition::growingGraph(children, g, Const::filling_factor, chosen_num);
    //GraphPartitioner::doPartitionStatic(children);

    // build rest data node if any
    int pid = children->size();
    for(auto &iter:*children){
        if(iter.second->actual_size <= Const::th && iter.second->partitionId == -1){
            assert(iter.second->dataNode == nullptr);
            buildDataNode(iter.second);
            iter.second->partitionId = ++pid;
        }
    }

    for(auto &iter: *children){
        if(iter.second->size > Const::th){
//            assert(iter.second->partitionId == -1);
//            assert(iter.second->dataNode == nullptr);
//            assert(iter.second->actual_size > Const::tardis_threshold);
            iter.second->growDynamicIPGIndex(g);
        }
    }
}

void IPGNode::growDynamicIPGIndexRange(vector<vector<int>> *g){
    if(size <= Const::th)   return;
    int chosen_num = SaxUtil::findFirstGE(power_2, 1, Const::segmentNum + 1, size / Const::th + 1);
    chooseSegmentDynamicRange(this, chosen_num);

    if(children == nullptr) children = new unordered_map<int, IPGNode*>();
    for(int i=0;i< (1 << chosen_num) ;++i){
        auto* node = new IPGNode(i, layer + 1);
//        copy(paa_max, paa_max + Const::segmentNum, node->paa_max);
//        copy(paa_min, paa_min + Const::segmentNum, node->paa_min);
//        int tmp = i;
//        for(int j=chosen_num - 1;j >= 0;--j, tmp >>= 1){
//            if(tmp % 2)
//                node->paa_min[chosenSegments[j]] = (paa_max[chosenSegments[j]] + paa_min[chosenSegments[j]]) / 2.0;
//            else
//                node->paa_max[chosenSegments[j]] = (paa_max[chosenSegments[j]] + paa_min[chosenSegments[j]]) / 2.0;
//        }
        (*children)[i] = node;
    }

//    vector<double>means(Const::segmentNum);
//    for(int i=0;i<Const::segmentNum;++i)
//        means[i] = (paa_min[i] + paa_max[i]) / 2.0;

    float *paa;
    for(int i=0;i<size;++i){
        int new_id = 0;
        paa = paas + offsets[i] * Const::segmentNum;
        for(int j=0; j< chosen_num;++j){
            new_id <<=1;
            if(paa[chosenSegments[j]] > paa_mu[chosenSegments[j]])
                new_id++;
        }
        (*children)[new_id]->insert(offsets[i], true);
    }

    vector<int>().swap(offsets);

    //partition
    IPGPartition::growingGraph(children, g, Const::filling_factor, chosen_num);
    //GraphPartitioner::doPartitionStatic(children);

    // build rest data node if any
    int pid = children->size();
    for(auto &iter:*children){
        if(iter.second->actual_size <= Const::th && iter.second->partitionId == -1){
            assert(iter.second->dataNode == nullptr);
            buildDataNode(iter.second);
            iter.second->partitionId = ++pid;
        }
    }

    // stat paa for internal node
    for(auto &iter:*children){
        if(iter.second->partitionId == -1)
            iter.second->statPaaRange();
    }

    for(auto &iter: *children){
        if(iter.second->size > Const::th){
//            assert(iter.second->partitionId == -1);
//            assert(iter.second->dataNode == nullptr);
//            assert(iter.second->actual_size > Const::tardis_threshold);
            iter.second->growDynamicIPGIndexRange(g);
        }
    }
}

int IPGNode::getDynamicId(const float *paa) const{
    int new_id = 0;

    for(int i=0; i< chosenSegments.size();++i){
        if(paa[chosenSegments[i]] < split_line[i])
            new_id <<=1;
        else
            new_id = (new_id << 1) + 1;
    }
    return new_id;
}

void IPGNode::generateSaxAndCardinality(IPGNode* node, int new_id, vector<int>& chosen_segments){
    copy(sax, sax + Const::segmentNum, node->sax);
    copy(bits_cardinality, bits_cardinality + Const::segmentNum, node->bits_cardinality);
    for(int i = chosen_segments.size() - 1; i >=0 ;--i){
        int seg = chosen_segments[i];
        node->bits_cardinality[seg]++;
        int t = new_id % 2 ;
        new_id >>= 1;
        node->sax[seg] = (node->sax[seg] << 1) + t;
    }
}

void IPGNode::statPaa1stLayer(){
    vector<double>split_line(Const::segmentNum, -0.67448975019608193);
    for(int i=0; i < Const::segmentNum; ++i)
        if(sax[i] == 1)
            split_line[i] = 0.67448975019608193;
    for(auto &i:paa_up_size) i = 0;
    for(auto &i:paa_below_size) i = 0;
    for(auto &i:paa_max) i = - numeric_limits<float>::max();
    for(auto &i:paa_min) i = numeric_limits<float>::max();
    for(auto &i:paa_up_min) i = numeric_limits<float>::max();
    for(auto &i:paa_below_max) i = - numeric_limits<float>::max();
    for(auto &i:paa_dist) i = 0;
    for(auto &i:paa_mu) i=0;
    for(auto &i:paa_up_avg) i=0;
    for(auto &i:paa_below_avg) i=0;
    for(int offset:offsets){
        float* start = paas + (long)offset * (Const::segmentNum);
        for(int i=0;i<Const::segmentNum;++i){
            float value = *(start + i);
            paa_mu[i] += value;
            paa_min[i] = min(paa_min[i], value);
            paa_max[i] = max(paa_max[i], value);
            paa_dist[i]  += abs(value - split_line[i]);
            if(value > split_line[i]) {
                paa_up_min[i] = min(paa_up_min[i], value);
                paa_up_avg[i] += value;
                paa_up_size[i]++;
            }
            else {
                paa_below_max[i] = max(paa_below_max[i], value);
                paa_below_avg[i]+=value;
                paa_below_size[i]++;
            }
        }
    }
    for(int i=0;i<Const::segmentNum;++i) {
        paa_mu[i] /= size;
        if(paa_below_size[i]>0)
            paa_below_avg[i] /= paa_below_size[i];
        if(paa_up_size[i] > 0)
            paa_up_avg[i] /= paa_up_size[i];
    }

    vector<double>paa_sum_square(Const::segmentNum, 0);
    vector<double>paa_sum_square_up(Const::segmentNum, 0);
    vector<double>paa_sum_square_below(Const::segmentNum, 0);
    for(int offset:offsets){
        float* start = paas + (long)offset * (Const::segmentNum);
        for(int i=0;i<Const::segmentNum;++i){
            float value = *(start + i);
            paa_sum_square[i] += (value - paa_mu[i]) * (value - paa_mu[i]);
            if(value > split_line[i])
                paa_sum_square_up[i] += (value - paa_up_avg[i]) * (value - paa_up_avg[i]);
            else
                paa_sum_square_below[i] += (value - paa_below_avg[i]) * (value - paa_below_avg[i]);
        }
    }
    // TODO: remove
    for(int i=0;i<Const::segmentNum;++i) {
        paa_sigma[i] = sqrt(paa_sum_square[i] / size);
        paa_up_sigma[i] = sqrt(paa_sum_square_up[i] / paa_up_size[i]);
        paa_below_sigma[i] = sqrt(paa_sum_square_below[i] / paa_below_size[i]);
    }
}

void IPGNode::statPaa(){
    vector<double>split_line(Const::segmentNum, 0);
    // TODO: optimize
    double lb;  // ub is the new split line
    for(int i=0; i < Const::segmentNum; ++i)
        SaxUtil::getValueRange(sax[i] << 1, bits_cardinality[i] + 1, &lb, &split_line[i]);
    for(auto &i:paa_max) i = - numeric_limits<float>::max();
    for(auto &i:paa_min) i = numeric_limits<float>::max();
    for(auto &i:paa_up_size) i = 0;
    for(auto &i:paa_below_size) i = 0;
    for(auto &i:paa_up_min) i = numeric_limits<float>::max();
    for(auto &i:paa_below_max) i = - numeric_limits<float>::max();
    for(auto &i:paa_dist) i = 0;
    for(auto &i:paa_mu) i=0;
    for(auto &i:paa_up_avg) i=0;
    for(auto &i:paa_below_avg) i=0;
    for(int offset:offsets){
        float* start = paas + (long)offset * (Const::segmentNum);
        for(int i=0;i<Const::segmentNum;++i){
            float value = *(start + i);
            paa_mu[i] += value;
            paa_min[i] = min(paa_min[i], value);
            paa_max[i] = max(paa_max[i], value);
            paa_dist[i]  += abs(value - split_line[i]);
            if(value > split_line[i]) {
//                paa_up_min[i] = min(paa_up_min[i], value);
                paa_up_avg[i] += value;
                paa_up_size[i]++;
            }
            else {
//                paa_below_max[i] = max(paa_below_max[i], value);
                paa_below_avg[i]+=value;
                paa_below_size[i]++;
            }
        }
    }
    for(int i=0;i<Const::segmentNum;++i) {
        paa_mu[i] /= size;
        if(paa_below_size[i]>0)
            paa_below_avg[i] /= paa_below_size[i];
        if(paa_up_size[i] > 0)
            paa_up_avg[i] /= paa_up_size[i];
    }

    vector<double>paa_sum_square(Const::segmentNum, 0);
    vector<double>paa_sum_square_up(Const::segmentNum, 0);
    vector<double>paa_sum_square_below(Const::segmentNum, 0);
    for(int offset:offsets){
        float* start = paas + (long)offset * (Const::segmentNum);
        for(int i=0;i<Const::segmentNum;++i){
            float value = *(start + i);
            paa_sum_square[i] += (value - paa_mu[i]) * (value - paa_mu[i]);
            if(value > split_line[i])
                paa_sum_square_up[i] += (value - paa_up_avg[i]) * (value - paa_up_avg[i]);
            else
                paa_sum_square_below[i] += (value - paa_below_avg[i]) * (value - paa_below_avg[i]);
        }
    }
    for(int i=0;i<Const::segmentNum;++i) {
        paa_sigma[i] = sqrt(paa_sum_square[i] / size);
        paa_up_sigma[i] = sqrt(paa_sum_square_up[i] / paa_up_size[i]);
        paa_below_sigma[i] = sqrt(paa_sum_square_below[i] / paa_below_size[i]);
    }
}

void IPGNode::statPaaOne(){
    vector<double>split_line(Const::segmentNum, 0);
    // TODO: optimize
    double lb;  // ub is the new split line
    for(int i=0; i < Const::segmentNum; ++i)
        SaxUtil::getValueRange(sax[i] << 1, bits_cardinality[i] + 1, &lb, &split_line[i]);
    for(auto &i:paa_up_size) i = 0;
    for(auto &i:paa_below_size) i = 0;
    for(int offset:offsets){
        float* start = paas + offset * (Const::segmentNum);
        for(int i=0;i<Const::segmentNum;++i){
            float value = *(start + i);
            if(value > split_line[i])
                paa_up_size[i]++;
            else
                paa_below_size[i]++;
        }
    }
}

void IPGNode::statPaaMean(){
    for(auto &i:paa_mu) i=0;
    for(int offset:offsets){
        float* start = paas + offset * (Const::segmentNum);
        for(int i=0;i<Const::segmentNum;++i){
            float value = *(start + i);
            paa_mu[i] += value;
        }
    }
    for(float & i : paa_mu)
        i /= size;
}

void IPGNode::statPaaMaxMinMean(){
    for(auto &i:paa_mu) i=0;
    for(auto &i:paa_max) i= - numeric_limits<float>::max();
    for(auto &i:paa_min) i= numeric_limits<float>::max();
    for(int offset:offsets){
        float* start = paas + offset * (Const::segmentNum);
        for(int i=0;i<Const::segmentNum;++i){
            float value = *(start + i);
            paa_mu[i] += value;
            paa_min[i] = min(paa_min[i], value);
            paa_max[i] = max(paa_max[i], value);
        }
    }
    for(float & i : paa_mu)
        i /= size;
}

void IPGNode::statPaaRange(){
    for(auto &i:paa_max) i = - numeric_limits<float>::max();
    for(auto &i:paa_min) i = numeric_limits<float>::max();
    for(auto &i:paa_mu) i = 0;
    for(int offset:offsets){
        float* start = paas + offset * (Const::segmentNum);
        for(int i=0;i<Const::segmentNum;++i){
            float value = *(start + i);
            paa_max[i] = max(paa_max[i], value);
            paa_min[i] = min(paa_min[i], value);
            paa_mu[i] += value;
        }
    }

    for(float & i : paa_mu)
        i /= size;
}

void IPGNode::insert(int offset, bool isReal) {
    offsets.push_back(offset);
    ++size;
    if(dataNode != nullptr) dataNode->size++;
    if(isReal) { ++actual_size;}
}

//vector<int> IPGNode::chooseSegmentsBeamSearch(IPGNode *node, int chosen_num){
//    if(chosen_num == 1)
//        return vector<int>{chooseOneSegment(node)};
//
//}

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

int IPGNode::chooseOneSegment(IPGNode* node){
    int min = node->size, min_index = -1;
    for(int i = 0; i < Const::segmentNum;++i){
        int big = max(node->paa_up_size[i], node->paa_below_size[i]);
        if(big < min){
            min = big;
            min_index = i;
        }
    }
    return min_index;
}

vector<int> (*IPGNode::method_list[6]) (IPGNode *node, int chosen_num) = {
        chooseSegmentRange, chooseSegmentStdev, chooseSegmentdistance, chooseSegmentdunn, chooseSegment3sigma, chooseSegment3sigmaNorm
};

vector<int> IPGNode::chooseSegment(IPGNode *node, int chosen_num, const string &method) {
    if(chosen_num == 1)
        return vector<int>{chooseOneSegment(node)};
    unordered_map<string, int>method_map;
    method_map["range"] = 0;
    method_map["stdev"] = 1;
    method_map["dist"] = 2;
    method_map["dunn"] = 3;
    method_map["3sigma"] = 4;
    method_map["3sigmaNorm"] = 5;
    return method_list[method_map[method]](node, chosen_num);
}

void IPGNode::chooseSegmentDynamic(IPGNode *node, int chosen_num){
    vector<float*>paa_list(node->size, nullptr);
    vector<double>lines(Const::segmentNum, 0);
    int j =0;
    for(int offset:node->offsets){
        paa_list[j++] = paas + offset * Const::segmentNum;
    }

    tmp scores[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i){
        scores[i] = tmp(i, scoreSegmentDynamic(node, i, paa_list, &lines[i]));
    }
    sort(scores, scores+Const::segmentNum, tmp::orderdesc);

    node->chosenSegments.resize(chosen_num);
    node->split_line.resize(chosen_num);
    for(int i=0;i<chosen_num;++i) {
        node->chosenSegments[i] = scores[i].i;
        node->split_line[i] = lines[scores[i].i];
    }
}

void IPGNode::chooseSegmentDynamicRange(IPGNode *node, int chosen_num){
    tmp scores[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i){
        scores[i] = tmp(i, (node->paa_max[i] - node->paa_min[i]));
    }
    sort(scores, scores+Const::segmentNum, tmp::orderdesc);

    node->chosenSegments.resize(chosen_num);
    for(int i=0;i<chosen_num;++i)
        node->chosenSegments[i] = scores[i].i;

}

double IPGNode::scoreSegmentDynamic(IPGNode* node, int i, vector<float *> &paa_list, double *line){
    int node_size = node->size;
    vector<double>paa_segment(node_size);
    for(int j=0;j<node_size;++j){
        paa_segment[j] = paa_list[j][i];
    }
    sort(paa_segment.begin(), paa_segment.end());
    int start = node_size * Const::imbalance, end = node_size - start;
    double max_gap = 0, gap;
    int index;
    for(int j = start;j<end;++j){
        gap = paa_segment[j + 1] - paa_segment[j];
        if(gap > max_gap){
            max_gap = gap;
            index = j;
        }
    }
    *line = (paa_segment[index + 1] + paa_segment[index]) / 2;
    return max_gap;
}

vector<int> IPGNode::chooseSegmentStdev(IPGNode *node, int chosen_num){
    int min_bc = 0;
    for(int &i: node->bits_cardinality)
        min_bc = min(min_bc, i);
    tmp scores[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        if(node->bits_cardinality[i] >= min_bc + Const::max_diff)
            scores[i] = tmp(i, numeric_limits<double>::max());
        else
            scores[i] = tmp(i, scoreSegmentstdev(node, i));
    sort(scores, scores+Const::segmentNum, tmp::order);

    vector<int>ret;
    int i = 0;
    for(;scores[i].score == -1 && i < Const::segmentNum;++i)
        ret.push_back(scores[i].i);
    for(;i<chosen_num;++i)
        ret.push_back(scores[i].i) ;
    sort(ret.begin(),  ret.end());
    return ret;
}

vector<int> IPGNode::chooseSegmentRange(IPGNode *node, int chosen_num){
    tmp scores[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        scores[i] = tmp(i, scoreSegmentrange(node, i));
    sort(scores, scores+Const::segmentNum, tmp::order);

    vector<int>ret;
    int _num = 0;
    for(int i=0;i<Const::segmentNum && _num < chosen_num;++i)
    {
        if(scores[i].score == -1)
            ret.push_back(scores[i].i);
        else if(scores[i].score < numeric_limits<double>::max()){
            ret.push_back(scores[i].i);
            ++_num;
        }
    }
    sort(ret.begin(),  ret.end());
    return ret;
}

vector<int> IPGNode::chooseSegmentdistance(IPGNode *node, int chosen_num){

    tmp scores[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        scores[i] = tmp(i, scoreSegmentdistance(node, i));
    sort(scores, scores+Const::segmentNum, tmp::orderdesc);

    vector<int>ret;
    int _num = 0;
    for(int i=0;i<Const::segmentNum && _num < chosen_num;++i)
    {
        if(scores[i].score == numeric_limits<double>::max())
            continue;
        if(scores[i].score == -1)
            break;
        ret.push_back(scores[i].i);
        ++_num;
    }

    for(int i = Const::segmentNum - 1; i >=0 && scores[i].score == -1; --i)
        ret.push_back(scores[i].i);
    sort(ret.begin(),  ret.end());
    return ret;
}

vector<int> IPGNode::chooseSegmentsigmarange(IPGNode *node, int chosen_num){

    tmp scores_range[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        scores_range[i] = tmp(i, scoreSegmentrange(node, i));
    sort(scores_range, scores_range + Const::segmentNum, tmp::order);

    tmp scores_sigma[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        scores_sigma[i] = tmp(i, scoreSegmentstdev(node, i));
    sort(scores_sigma, scores_sigma + Const::segmentNum, tmp::order);

    vector<int>ret;
    int num = 0;
    for(int i=0;i<Const::segmentNum && scores_range[i].score == -1;++i){
        ret.push_back(scores_range[i].i);
        ++num;
    }

    unordered_set<int>range,sigma;
    int start = num;
    int j = start;
    for(; j< start + chosen_num && j< Const::segmentNum &&  scores_range[j].score < numeric_limits<double>::max();++j){
        range.insert(scores_range[j].i);
        sigma.insert(scores_sigma[j].i);
    }

    for(auto & iter:sigma){
        if(range.find(iter) != range.end()){
            ++num;
            ret.push_back(iter);
        }
    }

    int flag = 0, cur;
    unordered_set<int> lists[2]{range, sigma};
    while (num < chosen_num && j< Const::segmentNum &&  scores_range[j].score < numeric_limits<double>::max()){
        cur = flag ? scores_sigma[j].i : scores_range[j].i;
        if(lists[1-flag].find(cur) != lists[1-flag].end()){
            ++num;
            ret.push_back(cur);
        }else{
            lists[flag].insert(cur);
        }

        if(flag)    ++j;
        flag = 1 - flag;
    }

    sort(ret.begin(),  ret.end());
    return ret;
}

vector<int> IPGNode::chooseSegmentdunn(IPGNode *node, int chosen_num){
    tmp scores[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        scores[i] = tmp(i, scoreSegmentdunn(node, i));
    sort(scores, scores+Const::segmentNum, tmp::orderdesc);

    vector<int>ret;
    int _num = 0;
    for(int i=0;i<Const::segmentNum && _num < chosen_num;++i)
    {
        if(scores[i].score == numeric_limits<double>::max())
            continue;
        if(scores[i].score == -1)
            break;
        ret.push_back(scores[i].i);
        ++_num;
    }

    for(int i = Const::segmentNum - 1; i >=0 && scores[i].score == -1; --i)
        ret.push_back(scores[i].i);
    sort(ret.begin(),  ret.end());
    return ret;
}

vector<int> IPGNode::chooseSegment3sigma(IPGNode *node, int chosen_num){
    tmp scores[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        scores[i] = tmp(i, scoreSegment3sigma(node, i));
    sort(scores, scores+Const::segmentNum, tmp::orderdesc);

    vector<int>ret;
    int _num = 0;
    for(int i=0;i<Const::segmentNum && _num < chosen_num;++i)
    {
        if(scores[i].score == numeric_limits<double>::max())
            continue;
        if(scores[i].score == -1)
            break;
        ret.push_back(scores[i].i);
        ++_num;
    }

    for(int i = Const::segmentNum - 1; i >=0 && scores[i].score == -1; --i)
        ret.push_back(scores[i].i);
    sort(ret.begin(),  ret.end());
    return ret;
}

vector<int> IPGNode::chooseSegment3sigmaNorm(IPGNode *node, int chosen_num){
    tmp scores[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        scores[i] = tmp(i, scoreSegment3sigmaNorm(node, i));
    sort(scores, scores+Const::segmentNum, tmp::orderdesc);

    vector<int>ret;
    int _num = 0;
    for(int i=0;i<Const::segmentNum && _num < chosen_num;++i)
    {
        if(scores[i].score == numeric_limits<double>::max())
            continue;
        if(scores[i].score == -1)
            break;
        ret.push_back(scores[i].i);
        ++_num;
    }

    for(int i = Const::segmentNum - 1; i >=0 && scores[i].score == -1; --i)
        ret.push_back(scores[i].i);
    sort(ret.begin(),  ret.end());
    return ret;
}

//double IPGNode::scoreSegmentdb(IPGNode* node,int id){
//    double a_mu = node->paa_mu[id], a_sigma = node->paa_sigma[id], a_min = node->paa_min[id], a_max = node->paa_max[id];
//    double lb, ub;  // ub is the new split line
//    SaxUtil::getValueRange(node->sax[id] << 1, node->bits_cardinality[id] + 1, &lb, &ub);
//    if(a_min >= ub || a_max <= ub)  return -1;
//    if(a_mu + 3* a_sigma <=ub || a_mu - 3*a_sigma >=ub) return numeric_limits<double>::max();
//    double dcen = (node->paa_up_avg[id] - node->paa_below_avg[id]);
//}

double IPGNode::scoreSegment3sigma(IPGNode* node,int i){
    if(node->bits_cardinality[i] >= Const::bitsCardinality)    return -1;
    double a_mu = node->paa_mu[i], a_sigma = node->paa_sigma[i], a_min = node->paa_min[i], a_max = node->paa_max[i];
    double lb, ub;  // ub is the new split line
    SaxUtil::getValueRange(node->sax[i] << 1, node->bits_cardinality[i] + 1, &lb, &ub);
    if(a_min >= ub || a_max <= ub)  return -1;
    if(a_mu + 3* a_sigma <=ub || a_mu - 3*a_sigma >=ub) return numeric_limits<double>::max();
//    double dmin = node->paa_up_min[id] - node->paa_below_max[id];
    return node->paa_sigma[i] / (node->paa_up_sigma[i] * node->paa_below_sigma[i]);
}

double IPGNode::scoreSegment3sigmaNorm(IPGNode* node,int i){
    if(node->bits_cardinality[i] >= Const::bitsCardinality)    return -1;
    double a_mu = node->paa_mu[i], a_sigma = node->paa_sigma[i], a_min = node->paa_min[i], a_max = node->paa_max[i];
    double lb, ub;  // ub is the new split line
    SaxUtil::getValueRange(node->sax[i] << 1, node->bits_cardinality[i] + 1, &lb, &ub);
    if(a_min >= ub || a_max <= ub)  return -1;
    if(a_mu + 3* a_sigma <=ub || a_mu - 3*a_sigma >=ub) return numeric_limits<double>::max();
//    double dmin = node->paa_up_min[id] - node->paa_below_max[id];
    return (node->paa_sigma[i] * node->paa_sigma[i]) / (node->paa_up_sigma[i] * node->paa_below_sigma[i] * (a_max - a_min)*(a_max - a_min));
}

double IPGNode::scoreSegmentdunn(IPGNode* node,int i){
    double a_mu = node->paa_mu[i], a_sigma = node->paa_sigma[i], a_min = node->paa_min[i], a_max = node->paa_max[i];
    double lb, ub;  // ub is the new split line
    SaxUtil::getValueRange(node->sax[i] << 1, node->bits_cardinality[i] + 1, &lb, &ub);
    if(a_min >= ub || a_max <= ub)  return -1;
    if(a_mu + 3* a_sigma <=ub || a_mu - 3*a_sigma >=ub) return numeric_limits<double>::max();
//    double dmin = node->paa_up_min[id] - node->paa_below_max[id];
    double dmin = node->paa_up_avg[i] - node->paa_below_avg[i];
    double diam = max((a_max - node->paa_up_min[i]), (node->paa_below_max[i]- a_min));
    double res = dmin / diam;
    return res >= 0? res : 0;
}

double IPGNode::scoreSegmentmean(IPGNode* node, int i){
    double a_mu = node->paa_mu[i], a_sigma = node->paa_sigma[i], a_min = node->paa_min[i], a_max = node->paa_max[i];
    double lb, ub;  // ub is the new split line
    SaxUtil::getValueRange(node->sax[i] << 1, node->bits_cardinality[i] + 1, &lb, &ub);
    if(a_min >= ub || a_max <= ub)  return -1;
    if(a_mu + 3* a_sigma <=ub || a_mu - 3*a_sigma >=ub) return numeric_limits<double>::max();
    return abs(ub - a_mu);

}

double IPGNode::scoreSegmentdistance(IPGNode* node, int i){
    if(node->bits_cardinality[i] >= Const::bitsCardinality)    return -1;
    double a_mu = node->paa_mu[i], a_sigma = node->paa_sigma[i], a_min = node->paa_min[i], a_max = node->paa_max[i];
    double lb, ub;  // ub is the new split line
    SaxUtil::getValueRange(node->sax[i] << 1, node->bits_cardinality[i] + 1, &lb, &ub);
    if(a_min >= ub || a_max <= ub)  return -1;
    if(a_mu + 3* a_sigma <=ub || a_mu - 3*a_sigma >=ub) return numeric_limits<double>::max();
    return node->paa_dist[i];
}

double IPGNode::scoreSegmentstdev(IPGNode* node, int i){
    if(node->bits_cardinality[i] >= Const::bitsCardinality)    return -1;
    double a_mu = node->paa_mu[i], a_sigma = node->paa_sigma[i], a_min = node->paa_min[i], a_max = node->paa_max[i];
    double lb, ub;  // ub is the new split line
    SaxUtil::getValueRange(node->sax[i] << 1, node->bits_cardinality[i] + 1, &lb, &ub);
    if(a_min >= ub || a_max <= ub)  return -1;
    if(a_mu + 3* a_sigma <=ub || a_mu - 3*a_sigma >=ub) return numeric_limits<double>::max();
    return 999 - a_sigma;

}

double IPGNode::scoreSegmentrange(IPGNode* node, int i){
    if(node->bits_cardinality[i] >= Const::bitsCardinality)    return -1;
    double a_mu = node->paa_mu[i], a_sigma = node->paa_sigma[i], a_min = node->paa_min[i], a_max = node->paa_max[i];
    double lb, ub;  // ub is the new split line
    SaxUtil::getValueRange(node->sax[i] << 1, node->bits_cardinality[i] + 1, &lb, &ub);
    if(a_min >= ub || a_max <= ub)  return -1;
    if(a_mu + 3* a_sigma <=ub || a_mu - 3*a_sigma >=ub) return numeric_limits<double>::max();
    return 999 - (a_max - a_min);

}

void IPGNode::build1stLayer(long series_number) const{
    for(long i=0;i<series_number;++i){
        if(i%10000000 == 0) cout << i << endl;
        unsigned short *asax = saxes + i * Const::segmentNum;
        int nav_id = SaxUtil::invSaxHeadFromSax(asax, Const::bitsCardinality, Const::segmentNum);
        IPGNode* node = (*children)[nav_id];
        if(node== nullptr){
            node = new IPGNode(nav_id, layer + 1);
            for(int &_:node->bits_cardinality)    _ = 1;
            int tmp = nav_id;
            for(int j=Const::segmentNum - 1;j>=0;--j) {
                node->sax[j] = tmp % 2;
                tmp >>= 1;
            }
            (*children)[nav_id] = node;
        }
        node->insert(i, true);
    }
}

void IPGNode::build1stLayerDynamic(int series_number) {
    for(auto &_:paa_max) _  = -numeric_limits<float>::max();
    for(auto &_:paa_min) _  = numeric_limits<float>::max();
    float * paa = paas;
    for(int i=0;i<series_number;++i){
        if(i%10000000 == 0) cout << i << endl;
        for(int j = 0; j<Const::segmentNum;++j, ++paa){
            paa_max[j] = max(paa_max[j], *paa);
            paa_min[j] = min(paa_min[j], *paa);
        }
    }

    vector<double>means(Const::segmentNum);
    for(int i=0;i<Const::segmentNum;++i)
        means[i] = (paa_max[i] + paa_min[i]) / 2.0;

    for(int i=0;i<Const::vertexNum;++i){
        auto* node = new IPGNode(i, layer + 1);
//        int tmp = i;
//        for(int j=Const::segmentNum - 1;j >= 0;--j, tmp >>= 1){
//            if(tmp % 2) {
//                node->paa_max[j] = paa_max[j];
//                node->paa_min[j] = means[j];
//            }else{
//                node->paa_max[j] = means[j];
//                node->paa_min[j] = paa_min[j];
//            }
//        }
        (*children)[i] = node;
    }

    paa = paas;
    for(int i=0;i<series_number;++i){
        if(i%10000000 == 0) cout << i << endl;
        int nav_id = 0;
        for(int j=0;j<Const::segmentNum;++j, ++paa){
            nav_id <<= 1;
            if(*paa > means[j])  ++nav_id;
        }
        (*children)[nav_id]->insert(i, true);
    }
}

void IPGNode::build1stLayerDynamicMean(int series_number) {
    for(auto &_:paa_mu) _  = 0;
    float * paa = paas;
    for(int i=0;i<series_number;++i){
        if(i%10000000 == 0) cout << i << endl;
        for(int j = 0; j<Const::segmentNum;++j, ++paa){
            paa_mu[j] += *paa;
        }
    }

    for(float & i : paa_mu)
        i /= series_number;

    for(int i=0;i<Const::vertexNum;++i){
        auto* node = new IPGNode(i, layer + 1);
//        int tmp = i;
//        for(int j=Const::segmentNum - 1;j >= 0;--j, tmp >>= 1){
//            if(tmp % 2) {
//                node->paa_max[j] = paa_max[j];
//                node->paa_min[j] = means[j];
//            }else{
//                node->paa_max[j] = means[j];
//                node->paa_min[j] = paa_min[j];
//            }
//        }
        (*children)[i] = node;
    }

    paa = paas;
    for(int i=0;i<series_number;++i){
        if(i%10000000 == 0) cout << i << endl;
        int nav_id = 0;
        for(int j=0;j<Const::segmentNum;++j, ++paa){
            nav_id <<= 1;
            if(*paa > paa_mu[j])  ++nav_id;
        }
        (*children)[nav_id]->insert(i, true);
    }
}


static int gid = 0;

void IPGNode::partition(IPGNode* root, vector<vector<int>> *g){
    if(root == nullptr || root->children == nullptr || root->size <= Const::th)
        return;
//    if(scoreGraph(root->children) < Const::dense_threshold)
//        IPGPartition::bubble(root->children);
//    else
    if(root->layer == 0)
        IPGPartition::growingGraph(root->children, g, Const::filling_factor_1st, Const::segmentNum);
    else
        IPGPartition::growingGraph(root->children, g, Const::filling_factor, Const::segmentNum);
    if(++gid % 1000 == 0)   cout << gid <<endl;
    for(auto& iter:*root->children)
        partition(iter.second, g);
}

//double IPGNode::scoreGraph(unordered_map<int, IPGNode*>*nodes_map){
//    int node_number = 0, total_size = 0;
//    for(auto& iter:*nodes_map){
//        if(iter.second != nullptr && iter.second->children == nullptr)
//            node_number++, total_size += iter.second->size;
//    }
//    return Const::node_number_weight * node_number / Const::vertexNum + Const::node_size_weight * total_size / (Const::th * node_number);
//}

struct CAND{
    int id;
    double score;

    CAND(int _i, double _score){ id = _i; score = _score;}

    static bool order(CAND &a, CAND &b){return a.score < b.score;}
};

void IPGNode::fuzzyNodeInFirstLayer(IPGNode* node) const{
    float *paa, range;
    IPGNode* temp_node;
    int new_id;
    vector<CAND>candidates;
    auto& node_offsets = node->offsets;
    for(int j=0;j<node->actual_size;++j){
        paa = IPGNode::paas + node->offsets[j] * Const::segmentNum;
        for(int i=0;i<Const::segmentNum;++i){
            new_id = node->id ^ mask_16[i];
            if(children->find(new_id) != children->end()){
                temp_node = (*children)[new_id];
                // temp node must not be an internal node, this is the RULE.
                // 1. pid == -1  => internal node
                // 2. same partition leaf node
                // 3. full leaf node partition
//                if(temp_node->partitionId == -1 || temp_node->partitionId == node->partitionId || temp_node->dataNode->size >= Const::tardis_threshold)    continue;
                if(temp_node->partitionId != -1 && (temp_node->partitionId == node->partitionId || temp_node->dataNode->size >= Const::th))    continue;

                {
                    //                    if(temp_node->partitionId == -1){
                    ////                        assert(temp_node->dataNode == nullptr);
                    ////                        assert(temp_node->size > Const::tardis_threshold);
                    //                        continue;
                    //                    }
                    //                    // now temp_node is a leaf node, partitionId != -1
                    ////                    assert(temp_node->isLeafNode());
                    ////                    assert(temp_node->dataNode != nullptr);
                    //                    if(temp_node->partitionId == node->partitionId){
                    //                        // node is also a leaf node, they are in same partition
                    ////                        assert(temp_node->partitionId != -1);
                    ////                        assert(temp_node->dataNode == node->dataNode);
                    //                        continue;
                    //                    }
                    //                    if(temp_node->dataNode->size >= Const::tardis_threshold){
                    ////                        assert(temp_node->dataNode->size == Const::tardis_threshold);
                    //                        continue;
                    //                    }
                }

                if(paa[i] > 0){
                    range = node->paa_mu[i] * Const::boundary_1st;
                    if(paa[i] <= range) candidates.emplace_back(new_id, paa[i]);
                }else{
                    range = (-node->paa_mu[i]) * Const::boundary_1st;
                    if(-paa[i] <= range)    candidates.emplace_back(new_id, -paa[i]);
                }
            }
        }
        sort(candidates.begin(), candidates.end(), CAND::order);
        int n= 0;
        for(int i = 0; i< candidates.size() && n < Const::max_replica; ++i){
            temp_node = (*children)[candidates[i].id];
//                assert(temp_node->isLeafNode());
//                assert(temp_node->dataNode != nullptr);
            if(temp_node->partitionId != -1 && temp_node->dataNode->size >= Const::th)
                continue;
            temp_node->insert(node_offsets[j], false);
            ++n;
        }
        candidates.clear();
        fuzzy_num += n;
    }
}

void IPGNode::fuzzyFirstLayer(IPGNode* root){
    int count = 0;
    for(auto & iter:*root->children){
        if(++count % 10000 == 0) cout <<count << endl;
        root->fuzzyNodeInFirstLayer(iter.second);
    }
}

void IPGNode::fuzzyNode(IPGNode* node, int chosen_num, vector<int>&chosen_segments) const{
    float *paa;
    double range, lb, ub;
    IPGNode*temp_node;
    vector<CAND>candidates; // max is chosen num
    int new_id, seg;
    for(int j=0;j<node->actual_size;++j){
        paa = IPGNode::paas + node->offsets[j] * Const::segmentNum;
        for(int i=0;i<chosen_num;++i){
            seg = chosen_segments[i];
            new_id = node->id ^ mask_16[16 - chosen_num + i];
            if(children->find(new_id) != children->end()){
                temp_node = (*children)[new_id];
                // leaf node full or same, continue
                if(temp_node->size <= Const::th){
                    assert(temp_node->dataNode != nullptr);
                    if(temp_node->dataNode->size >= Const::th || temp_node->partitionId == node->partitionId)
                        continue;
                }
                if(node->sax[seg] == (1 << node->bits_cardinality[seg]) - 1)
                    SaxUtil::getValueRange(node->sax[seg] - 1, node->bits_cardinality[seg], &lb, &ub);
                else if(node->sax[seg] == 0)
                    SaxUtil::getValueRange(1, node->bits_cardinality[seg], &lb, &ub);
                else
                    SaxUtil::getValueRange(node->sax[seg] , node->bits_cardinality[seg], &lb, &ub);
                range = (ub -lb) * Const::boundary;
                if(node->sax[seg] == (1 << node->bits_cardinality[seg]) - 1){
                    if(paa[seg] - ub <= range) candidates.emplace_back(new_id, paa[seg]);
                }else if(node->sax[seg] == 0){
                    if(lb - paa[seg] <= range)  candidates.emplace_back(new_id, -paa[seg]);
                }else if(node->sax[seg] % 2 == 0){
                    if(abs(ub - paa[seg]) <= range) candidates.emplace_back(new_id, paa[seg]);
                }else{
                    if(abs(paa[seg] - lb) <= range)    candidates.emplace_back(new_id, -paa[seg]);
                }
            }
        }
        sort(candidates.begin(), candidates.end(), CAND::order);
        int n= 0;
        for(int i = 0; i< candidates.size() && n < Const::max_replica; ++i){
            temp_node = (*children)[candidates[i].id];
            if(temp_node->partitionId != -1){   // temp node is a leaf node
                assert(temp_node->dataNode != nullptr);
                if(temp_node->dataNode->size < Const::th){
                    temp_node->insert(node->offsets[j], false);
                    ++n;
                }
            }else{  // temp node is an internal node
                assert(temp_node->actual_size > Const::th);
                temp_node->insert(node->offsets[j], false);
                ++n;
            }
        }
        candidates.clear();
        fuzzy_num +=n;
//            int id=0;
//            for(;id<candidates.size() && id < Const::max_replica;++id)
//                segments[id] = candidates[id].id;
//            num += id;
//            // todo: permutation
//            for(int k=0;k<id;++k)
//                childrens[id ^ mask_16[16 - chosen_num + segments[k]]]->insert(node->offsets[j], false);
    }

}

void IPGNode::fuzzy(IPGNode* parent, vector<int>&chosen_segments){
    int chosen_num = chosen_segments.size(), seg;
    auto childrens = *parent->children;
    for(auto & iter:childrens){
//        if(iter.second->size <= Const::low_filling_factor * Const::tardis_threshold)    continue;
        parent->fuzzyNode(iter.second, chosen_num, chosen_segments);
    }

}

IPGDataNode* IPGNode::buildDataNode(const vector<IPGNode*>& nodes){
    if(nodes.size() < 2)    return buildDataNode(nodes[0]);
    auto* datanode = new IPGDataNode();
    for(IPGNode* node:nodes){
        datanode->size += node->size;
        node->dataNode = datanode;
    }
    return datanode;
}

IPGDataNode* IPGNode::buildDataNode(IPGNode* node){
    auto* datanode = new IPGDataNode();
    node->dataNode = datanode;
    if(node->size > Const::th) {
        datanode->size = Const::th;
        node->size = Const::th;
        node->offsets.resize(Const::th);
    }
    else
        datanode->size = node->size;
    return datanode;
}

vector<OffsetDist *> *IPGNode::preparePqUsingSaxLbNew(double bsf, const float *queryPaa) const {
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

extern int _search_num;
void IPGNode::exactSearchKnn(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap) const {
    assert(isLeafNode());
    double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
    FILE *f = fopen(rowDataFileName.c_str(), "rb");
    float *ts;
    for(int offset:offsets){
        ++_search_num;
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

IPGNode* IPGNode::route(const unsigned short *asax) const{
    if(children == nullptr) return nullptr;
    int * new_bc;
    for(auto &iter:*children)
        if(iter.second != nullptr){
            new_bc = iter.second->bits_cardinality;
            break;
        }
    auto old_bc = bits_cardinality;
    int res = 0;
    for(int i=0;i<Const::segmentNum;++i){
        if(new_bc[i] > old_bc[i]){
            int t = (asax[i] >> (Const::bitsCardinality - new_bc[i])) % 2; // take the new_bc[id] bit of a single sax
            res = (res << 1) + t;
        }
    }
    return (*children)[res];
}

IPGNode* IPGNode::route(const float *paa) const{
    if(children == nullptr) return nullptr;
    int new_id = getDynamicId(paa);
    if(children->find(new_id) == children->end())   return nullptr;
    return (*children)[new_id];
}

IPGNode* IPGNode::routePaaMu(const float *paa) const{
    if(children == nullptr) return nullptr;
    int * new_bc;
    for(auto &iter:*children)
        if(iter.second != nullptr){
            new_bc = iter.second->bits_cardinality;
            break;
        }
    auto old_bc = bits_cardinality;
    int res = 0;
    for(int i=0;i<Const::segmentNum;++i){
        if(new_bc[i] > old_bc[i]){
            res = res << 1;
            if(paa[i] > paa_mu[i]) res++;
        }
    }
    return (*children)[res];
}

IPGNode* IPGNode::routeIn1stLayer(const float *paa) const{
    int nav_id = 0;
    const float *work = paa;
    for(int j=0;j<Const::segmentNum;++j, ++work){
        nav_id <<= 1;
//        if(*work > (paa_max[j]+paa_min[j]) / 2.0)  ++nav_id;
        if(*work > paa_mu[j])   ++nav_id;
    }
    if(children->find(nav_id) == children->end())   return nullptr;
    return (*children)[nav_id];
}

IPGNode* IPGNode::routeBelow1stLayer(const float *paa) const{
    int nav_id = 0;
    for(int j=0;j<chosenSegments.size();++j){
        nav_id <<= 1;
//        if(paa[chosenSegments[j]] > (paa_max[chosenSegments[j]]+paa_min[chosenSegments[j]]) / 2.0)  ++nav_id;
        if(paa[chosenSegments[j]] > paa_mu[chosenSegments[j]])  ++nav_id;
    }
    if(children->find(nav_id) == children->end())   return nullptr;
    return (*children)[nav_id];
}

int IPGNode::loadSax(const string & saxfn){
    long f_size = FileUtil::getFileSize(saxfn.c_str()), series_num = f_size / (sizeof(unsigned short ) * Const::segmentNum);
    saxes = new unsigned short [f_size / sizeof(unsigned short)];
    FILE *f = fopen(saxfn.c_str(), "rb");
    fread(saxes, sizeof(unsigned short ), f_size / sizeof(unsigned short), f);
    fclose(f);
    cout << "Finish loading sax"<<endl;
    return series_num;
}

void IPGNode::loadPaa(const string & paafn){
    long f_size = FileUtil::getFileSize(paafn.c_str());
    paas = new float [f_size / sizeof(float )];
    FILE *f = fopen(paafn.c_str(), "rb");
    fread(paas, sizeof(float ), f_size / sizeof(float ), f);
    fclose(f);
    cout << "Finish loading paa"<<endl;
}

IPGNode::IPGNode() {;

}

IPGNode *IPGNode::loadFromDisk(const string &saxfn, const string &paafn, const string &idxfn, bool loadFull) {
    if(loadFull) {
        loadSax(saxfn);
        loadPaa(paafn);
    }
    ifstream ifs(idxfn.c_str(), ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    auto *g = new IPGNode();
    ia >> (*g);
    ifs.close();
    return g;
}

IPGNode *IPGNode::loadFromDisk(const string &saxfn, const string &idxfn) {
    loadSax(saxfn);
    ifstream ifs(idxfn.c_str(), ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    auto *g = new IPGNode();
    ia >> (*g);
    ifs.close();
    return g;
}

bool IPGNode::isLeafNode() const {
    return children == nullptr;
}

int IPGNode::getHeight() const{
    if(children == nullptr) return 0;
    int tmp = 0;
    for(auto& iter: *children){
        if(iter.second!= nullptr)
            tmp = max(tmp, iter.second->getHeight());
    }
    return tmp + 1;
}

int IPGNode::getFakeLeafNodeNumber() const{
    if(children == nullptr) return 1;
    int sum = 0;
    for(auto &iter: *children){
        if(iter.second != nullptr)
            sum += iter.second->getFakeLeafNodeNumber();
    }
    return sum;
}

long IPGNode::generateSaxAndPaaTbl(){
    string fn = Const::datafn;
    long fs = FileUtil::getFileSize(fn.c_str());
    long series_num = fs / Const::tsLengthBytes;
    cout << "Total Series Number is "<<series_num <<endl;
    float ts[Const::tsLength];
    saxes = new unsigned short[series_num * Const::segmentNum];
    paas = new float[series_num * Const::segmentNum];
    long rest = series_num, cur = 0;
    FILE *f = fopen(fn.c_str(), "rb");

    while(rest > 0){
        int num;
        if(rest > Const::fbl_series_num)    num = Const::fbl_series_num;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

//        auto start = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
//        auto end = chrono::system_clock::now();
//        SAX_PAA_READ_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();

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
//        start = chrono::system_clock::now();
//        SAX_PAA_CPU_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();
    }

    fclose(f);
    return series_num;
}

