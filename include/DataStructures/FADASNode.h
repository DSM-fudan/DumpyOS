//
// Created by wzy on 2021/11/8.
//

#ifndef MULGIFT_FADASNODE_H
#define MULGIFT_FADASNODE_H
#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <unordered_set>
#include "../Const.h"
#include "TimeSeries.h"
#include "PqItemSeries.h"

using namespace std;

struct partUnit{
    int id;
    int size;
    int pid;

    static bool comp_size(partUnit*x, partUnit*y){
        return x->size > y->size;
    }
};

struct PAA_INFO{
    double paa_variance[Const::segmentNum]{};
    int paa_up_size[Const::segmentNum]{},
    paa_below_size[Const::segmentNum]{};
};

struct SAX_INFO{
    double sax_variance[Const::segmentNum]{};
    int up_size[Const::segmentNum]{},
            below_size[Const::segmentNum]{};
};

struct FBL_UNIT{
    int size{0};
    float *buffer{nullptr};
    int pos{0};
};

struct SAX_BUF_UNIT{
    int size{0};
    unsigned short *buffer{nullptr};
};

struct LBL_UNIT{
    vector<int>offsets;
    vector<float*>buffer;
};

struct NODE_RECORDER;

class FADASNode {
    explicit FADASNode(int _layer) {
        layer = _layer;
    }
    FADASNode(int _layer, int pid) {
        layer = _layer;
        partition_id = pid;
        bits_cardinality[0] = -1;
    }
    // only for node in 1st layer
    FADASNode(int _layer, int _size, int _id) {
        layer = _layer;
        size = _size;
        id = _id;
        file_id = to_string(id);
        chosenSegments.clear();
        for(auto &i:sax)    i=0;
        for(auto &i:bits_cardinality)   i = 1;
    }
    FADASNode(const FADASNode* parent, int _size, int _id){
        layer = parent->layer + 1;
        id = _id;
        size = _size;
        file_id  = parent->file_id + "_" + to_string(id);
        chosenSegments.clear();
    }
    FADASNode(const FADASNode* parent, int pid){
        layer = parent->layer + 1;
        partition_id = pid;
        file_id  = parent->file_id + "_" + to_string(partition_id);
        bits_cardinality[0] = -1;
    }
    FADASNode(const FADASNode* parent, int _id, bool is_insert){
        layer = parent->layer + 1;
        file_id  = parent->file_id + "_" + to_string(_id);
        size = 1;
        id = _id;
    }

    static void loadPaa(const string & paafn);
    static int loadSax(const string & saxfn);
    static long generateSaxAndPaaTbl();
    PAA_INFO* statPaa();
    void chooseSegment(PAA_INFO *paa, int chosen_num);
    static int chooseOneSegment(PAA_INFO *node);
    void generateSaxAndCardIn1stLayer(int new_id);
    void generateSaxAndCardinality(FADASNode *node, int new_id);
    void generateSaxAndCardIn1stLayer4LeafNode(int new_id);
    void generateSaxAndCardinality4LeafNode(FADASNode *node, int new_id);
    static int partition1stLayer(partUnit *nodes_map, vector<vector<int>> *g,double filling_factor);
    static int partition(partUnit *nodes_map, int chosen_num);
    void growIndex();
    void growIndexLessPack();
    void growIndexFuzzy(unordered_map<FADASNode *, NODE_RECORDER> &navigating_tbl, vector<vector<int>> *g);

    void collectSAXwords(unsigned short *node_saxes, int *cur, vector<string> &leaf_files, vector<string> &sax_files);
    void deleteSubtree();

    void fuzzySeriesInPartUnit(partUnit *part_units, int actual_size, int chosen_num, vector<int> &node_offsets,
                               vector<int> &series_index_list,
                               unordered_map<FADASNode *, NODE_RECORDER> &navigating_tbl, int _id) const;
    void fuzzy(partUnit *part_units, vector<int> &actual_sizes, vector<vector<int>> &node_offsets,
               vector<vector<int>> &series_index_list, int chosen_num,
               unordered_map<FADASNode *, NODE_RECORDER> &navigating_tbl) const;
    void fuzzySeriesInPartUnitInFirstLayer(partUnit *part_units, vector<int> &node_offsets, int _id,
                                           unordered_map<FADASNode *, NODE_RECORDER> &navigating_tbl,
                                           vector<vector<double >> &paa_mu_part_units) const;
    void fuzzyFirstLayer(partUnit *part_units, const int *nav_ids,
                         unordered_map<FADASNode *, NODE_RECORDER> &navigating_tbl,
                         vector<vector<double>> &paa_mu_part_units) const;

    int getMaxHeight();
    int getNodeNum();
    int getTotalSize();
    int getSumHeight();
    int get1stLayerInterNodesNo();
    int get1stLayerInterNodeSeriesNo();
    int get1stLayerNodesNo();


public:
    const static int power_2[];
    static int* mask;
    static unsigned short *saxes;
    static float *paas;
    static int a,b,c;

    vector<int>pos_cache;
    unsigned short sax[Const::segmentNum]{};
    vector<int> chosenSegments{};
    vector<FADASNode*>children;
    int size = 0;
    int leaf_num = 0;
    int id = -1;
    string file_id{};   // to identify a particular leaf node
    int layer = 0;
    int bits_cardinality[Const::segmentNum]{};
    int partition_id = -1;
    vector<int> offsets{};

    FADASNode *route(const unsigned short *_sax);
    static FADASNode *BuildIndex(string &datafn, string &saxfn);
    static FADASNode *BuildIndexLessPack(string& datafn, string &saxfn, string &paafn, vector<vector<int>> *g);
    static FADASNode *BuildIndexPos(string &datafn, string &saxfn, string &paafn, vector<vector<int>> *g);

    [[nodiscard]] string getFileName() const{
        if(layer == 1)  return "1_" + to_string(partition_id);
        return to_string(layer) + "-" + file_id;
    }
    [[nodiscard]] string getFileNamePack() const{
        if(layer == 1) {
            if(isLeafPack())    return "1_P_" + to_string(partition_id);
            if(isLeafNode())    return "1_" + file_id;
        }
        if(isLeafPack())    return "P-" + to_string(layer) + "-" + file_id;
        if(isLeafNode())    return to_string(layer) + "-" + file_id;
    }
    void getFileNameInsert(const string &index_dir, string &sax_file, string &data_file) const;
    [[nodiscard]] bool isLeafNode() const   {return size <= Const::th && partition_id == -1;}
    [[nodiscard]] bool isLeafPack() const {return size <= Const::th && partition_id != -1;}
    [[nodiscard]] bool isInternalNode() const {return  size > Const::th;}
    void search(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir,
                float *query_reordered, int *ordering) const;
    void search(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir) const;
    void search(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir,
                std::unordered_set<float *, createhash, isEqual> *hash_set) const;

    void search_offset(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir) const;

    static FADASNode*  BuildIndexFuzzy(const string & datafn, const string & saxfn, const string &paafn, vector<vector<int>>* g);

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & partition_id; ar & layer;
        ar & file_id; ar & id;
        ar & size;
        ar & sax;   ar & bits_cardinality;
        ar & children;
        ar & chosenSegments;
    }
    void save2Disk(const string &output) {
        ofstream ofs(output.c_str(), ios::binary);
        boost::archive::binary_oarchive oa(ofs);
        oa << (*this);
        ofs.close();
    }
    static FADASNode *loadFromDisk(const string &saxfn, const string &idxfn, bool need_sax);

    void getIndexStats();
    int getLeafNodeNum();
    int assignLeafNum();
    int getBiasLeafNodeNum();

    FADASNode *route1step(const unsigned short *_sax);


    SAX_INFO *statSAX();

    static int chooseOneSegment(SAX_INFO *node);

    void chooseSegment(SAX_INFO *sax_info, int chosen_num);

    static FADASNode *BuildIndexWOPack(string &datafn, string &saxfn, string &paafn, vector<vector<int>> *g);

    void growIndexWOPack();

    static int partitionLessPack(partUnit *nodes_map, int chosen_segment_number);

    static int partitionNew(partUnit* nodes_map, int chosen_segment_number);

    void searchLessPack(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir) const;

    void determineSegments();

    void determineFanout(int *lambda_min, int *lambda_max) const;

    double compute_score(vector<int> &node_sizes, int *plan, int lambda, vector<double> &data_seg_stdev) const;

    void visitPlanFromBaseTable(unordered_set<int> &visited, int cur_lambda, const int *plan, vector<int> &base_tbl,
                                double *max_score, vector<int> &best_plan, int lambda_min, int mask_code,
                                vector<double> &data_seg_stdev, double base_score);

    void determineSegmentsCluster();

    double
    compute_score_cluster(vector<int> &node_sizes, vector<vector<double>> &plan_node_seg_mean, vector<double> &seg_mean,
                          int *plan, int lambda) const;

    void
    visitPlanFromBaseTableCluster(unordered_set<int> &visited, int cur_lambda, const int *plan,
                                  vector<int> &base_tbl_size,
                                  vector<vector<double>> &base_tbl_seg_sum, double *max_score, vector<int> &best_plan,
                                  int lambda_min, int mask_code, vector<double> &data_seg_mean);

    double compute_score_cluster_weak(vector<int> &node_sizes, vector<vector<double>> &plan_node_seg_mean,
                                      vector<double> &seg_mean, int *plan, int lambda) const;

    void visitPlanFromBaseTableWeakCluster(unordered_set<int> &visited, int cur_lambda, const int *plan,
                                           vector<int> &base_tbl_size, vector<vector<double>> &base_tbl_seg_sum,
                                           double *max_score, vector<int> &best_plan, int lambda_min, int mask_code,
                                           vector<double> &data_seg_mean);

    void determineSegmentsWeakCluster();

    void determineSegmentsAvgVariance();

    void determineSegmentsNaive();

    void searchDTW(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, const string &index_dir) const;

    static long generateSaxTbl();

    FADASNode(){;}

    static void generateSaxTbl(const float *tss, int series_num);

    void routeDuringInsertion(const unsigned short *_sax, int pos);

    void determineSegments(unsigned short *node_saxes);

    void growIndex(unsigned short *node_saxes, bool need_free);

    void reorganize(float *tss, FADASNode *parent);

    void insertBatch(float *tss, int batch_size);
};

struct NODE_RECORDER{
    int actual_size{0};
    vector<int>series_index_list{};

    explicit NODE_RECORDER(int act_size, FADASNode *node) {
        actual_size = act_size;
        // internal node stores all series index list
        series_index_list.resize(node->size);
        for(int i=0;i<node->size;++i)
            series_index_list[i] = i;
    }

    NODE_RECORDER(int size, vector<int>&index_list){
        actual_size = size;
        series_index_list.resize(index_list.size());
        copy(index_list.begin(),  index_list.end(), series_index_list.begin());
    }

    NODE_RECORDER(){actual_size = 0;}
};


#endif //MULGIFT_FADASNODE_H
