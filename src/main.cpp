#include <iostream>
#include <cstdio>
#include <fstream>
#include <vector>
#include <thread>
#include <chrono>
#include "../include/DataStructures/IPGNode.h"
#include "../include/DataStructures/FADASNode.h"
#include "../include/TAR/TARGNode.h"
#include "../include/DataStructures/GraphConstruction.h"
#include "../include/Expr/DNATranslator.h"
#include "../include/Expr/ECGParser.h"
#include "../include/DSTree/DSTreeNodeConstruction.h"
#include "../include/Utils/FileUtil.h"
#include "../include/Utils/SaxUtil.h"
#include "../include/Expr/Recall.h"
#include "../include/Expr/DataDistribution.h"
#include "../include/Expr/RandDataGenerator.h"
#include "../include/Searchers/DSTreeApproxSearcher.h"
#include "../include/Utils/TimeSeriesUtil.h"
#include "../include/Searchers/DSTreeExactSearcher.h"
#include "../include/Searchers/ExactSearcher.h"
#include "../include/Tardis/TardisTreeNode.h"
#include "../include/DataStructures/IPGNode.h"
#include "../include/DataStructures/iSAXNode.h"
#include "../include/Searchers/iSAXSearcher.h"

using namespace std;

string getIdxfn(){
    string idxfn = Const::idxfn;
    int pos = idxfn.rfind(".bin");
    idxfn = idxfn.insert(pos, "_" + Const::method);
    cout << "real idx fn is "<<idxfn<<endl;
    return idxfn;
}

string getOtherIndexfn(const string& name){
    string idxfn = Const::idxfn;
    int pos = idxfn.rfind("IPG");
    string ret = idxfn.substr(0, pos) + name + idxfn.substr(pos+3);
    cout << "real idx fn is "<<ret<<endl;
    return ret;
}

string getIdxfnFuzzy(){
    string idxfn = Const::idxfn;
    int pos = idxfn.rfind(".bin");
    int fuzzy1st = Const::boundary_1st * 100, fuzzy = Const::boundary * 100;
    idxfn = idxfn.insert(pos, "_fuzzy" + to_string(Const::max_replica) + "-" + to_string(fuzzy1st) + to_string(fuzzy) + "_" + Const::method);
    cout << "real idx fn is "<<idxfn<<endl;
    return idxfn;
}

vector<vector<int>>* loadGraphSkeleton(){
    int vd = 0;
    for(int i=1; i<=Const::bitsReserve;++i)
        vd += MathUtil::nChooseK(Const::segmentNum, i);
    auto nnList = new vector<vector<int>>(Const::vertexNum, vector<int>(vd, -1));

    if(!FileUtil::checkFileExists(Const::graphfn.c_str())){
        cout << "File not exists!" << Const::graphfn << endl;
        exit(-1);
    }
    FILE *f = fopen(Const::graphfn.c_str(), "rb");

    for(int i=1;i<Const::vertexNum;++i)
        fread(&((*nnList)[i][0]), sizeof(int), vd, f);

    return nnList;

}

void constructGraph(){
    GraphConstruction::buildAndSave2Disk();
}

void DSTreeExpr(){
    Const::dstreefn += to_string(Const::segmentNum)+ "_" + to_string(Const::th) + "/";
    cout << Const::dstreefn << endl;
    DSTreeNode*root = DSTreeNode::loadFromFile();
    Recall::progressiveSearchInMeomoryDSTree(root);
}

//void generateQueryFile(){
//    string fn = "../data/generator/deep1b-96-100m.bin_le";
//    FileUtil::generateQueryFile(fn, 2000);
//}
//void buildGraph(){
//    auto s1 = chrono::system_clock::now();
//    Graph g;
//    g.loadGraphSkeletonFromFile();
//    auto s2 = chrono::system_clock::now();
//    cout << "Graph Skeleton has been loaded successfully." << endl;
//    g.materializeGraphWithDataFile(Const::rawDataFileName);
//    auto s3 = chrono::system_clock::now();
//    cout << "save graph to disk." << endl;
//    g.save2Disk();
//    auto s4 = chrono::system_clock::now();
//
//    cout << "1. Loading graph skeleton from file costs " << chrono::duration_cast<chrono::microseconds>(s2 - s1).count() << " us." << endl;
//    cout << "2. Building graph costs " << chrono::duration_cast<chrono::microseconds>(s3 - s2).count() << " us."<<endl;
//    cout << "3. Saving to disk costs " << chrono::duration_cast<chrono::microseconds>(s4 - s3).count() << " us."<<endl;
//
//
//}
//
//void buildGraphDebug(){
//    auto s1 = chrono::system_clock::now();
//    Graph g;
//    g.loadGraphSkeletonFromFile();
//    auto s2 = chrono::system_clock::now();
//    cout << "Graph Skeleton has been loaded successfully." << endl;
//    g.materializeGraphWithDataFile(Const::rawDataFileName);
//    auto s3 = chrono::system_clock::now();
//
//    cout << "1. Loading graph skeleton from file costs " << chrono::duration_cast<chrono::microseconds>(s2 - s1).count() << " us." << endl;
//    cout << "2. Building graph costs " << chrono::duration_cast<chrono::microseconds>(s3 - s2).count() << " us."<<endl;
//}
//void buildGraphNew(){
//    auto s1 = chrono::system_clock::now();
//    Graph g;
//    g.loadGraphSkeletonFromFile();
//    auto s2 = chrono::system_clock::now();
//    cout << "Graph Skeleton has been loaded successfully." << endl;
//    g.buildGraphNew(Const::rawDataFileName);
//    auto s3 = chrono::system_clock::now();
//    g.save2Disk();
//    auto s4 = chrono::system_clock::now();
//
//    cout << "1. Loading graph skeleton from file costs " << chrono::duration_cast<chrono::microseconds>(s2 - s1).count() << " us." << endl;
//    cout << "2. Building graph costs " << chrono::duration_cast<chrono::microseconds>(s3 - s2).count() << " us."<<endl;
//    cout << "3. save to disk costs " << chrono::duration_cast<chrono::microseconds>(s4 - s3).count() << " us."<<endl;
//
//}
//Graph* readGraphIndex(){
//    string fn = Const::rawDataFileName;
//    Graph* g = Graph::loadFromFileNew(fn);
//    cout << "load finish, start search"<<endl;
//    return g;
//}
//void loadGraphAndShow(){
//    auto *g = new Graph();
//    g->loadGraphSkeletonFromFile();
//    g->showGraph(100);
//}
//void recallExpr(){
//    Graph * g = readGraphIndex();
//    Recall::doExpr(g);
//}
//
//void recallExprRes(){
//    Graph * g = readGraphIndex();
//    Recall::doExprWithRes(g);
//}
//
//void recallExprDebug(){
//    Graph * g = readGraphIndex();
//    Recall::doExprDebug(g);
//}
//void recallExprPerformance(){
//    Graph * g = readGraphIndex();
//    cout << "load finish, start search"<<endl;
//    Recall::doExprPerformance(g);
//}
//void buildPartitionInvSAXTree(){
//    string saxfn = "../data/sax/Series_gaussian_randomwalk256_100000000_" + to_string(Const::segmentNum) + ".bin_le";
//    auto start = chrono::system_clock::now();
//    TardisTreeNode* root = TardisTreeNode::BuildTardisTree(saxfn);
//    auto end  =chrono::system_clock::now();
//    cout << "build tree structure finish, need " << chrono::duration_cast<chrono::microseconds>(end- start).count() << "us."<<endl;
//    auto g = Graph::loadGraphSkeletonFromFileSimple();
//    auto end2 = chrono::system_clock::now();
//    cout << "start partition, load graph skeleton cost " << chrono::duration_cast<chrono::microseconds>(end2- end).count() << "us."<< endl;
//    root->partition(root, g);
//    auto end3 = chrono::system_clock::now();
//    cout << "partition finish, cost " << chrono::duration_cast<chrono::microseconds>(end3- end2).count() << "us."<< endl;
//    root->save2Disk("../data/tardis/rand_10k_full_model_2.bin");
//    auto end4 = chrono::system_clock::now();
//    cout << "save to disk, cost " << chrono::duration_cast<chrono::microseconds>(end4- end3).count() << "us."<< endl;
//    cout << "totally " << chrono::duration_cast<chrono::microseconds>(end4- start).count() << "us."<<endl;
//}
//void generateGroundTruth2(){
//    string fn = "../data/generator/Series_gaussian_randomwalk256_100000000.bin_le";
//    string result = "../data/results/seismic-in.bin_le";
//    FILE *f = fopen("../data/generator/seismic-256-100m_query_query_in.bin", "rb");
//    auto *queries = new float [2000 * TimeSeries::tsLength];
//    fread(queries, sizeof(float ), 2000 * TimeSeries::tsLength, f);
//    fclose(f);
//    f = fopen(result.c_str(), "wb");
//    Graph * g = readGraphIndex();
//    for(int i=0;i<2000;++i){
//        cout << i << endl;
//        vector<PqItemSeries*> *exactKnn = ExactSearcher::exactKnnSearch(g, queries + i * TimeSeries::tsLength, 500, nullptr);
//        for(int j =0;j<500;++j){
//            fwrite((*exactKnn)[j]->ts, sizeof(float ), 256, f);
//        }
//    }
//}

void buildDSTree(){
    Const::dstreefn += to_string(Const::segmentNum)+ "_" + to_string(Const::th) + "/";
    cout << Const::dstreefn << endl;
    DSTreeNode* root = DSTreeNodeConstruction::buildIndex(Const::series_num);
    root->save2Disk();
}

void generateGroundTruth(){
    ExactSearcher::groundTruthKnnSearch(Const::datafn, Const::queryfn, Const::maxK, Const::query_num, Const::resfn, Const::series_num);
}

void generateGroundTruthDTW(){
    ExactSearcher::groundTruthKnnSearchDTW(Const::datafn, Const::queryfn, Const::maxK, Const::query_num, Const::dtwresfn, Const::series_num);
}

void testExactResDTW(){
    string result = Const::dtwresfn;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");

//    FILE *f = fopen("/mnt/c/Series4Similarity_Search/rand-256-1k.bin", "rb");
    FILE *resf = fopen(result.c_str(), "rb");
    auto *ts1 = new float[Const::tsLength], *ts2 = new float [Const::tsLength];
    for(int i=0;i<Const::query_num;++i){
        cout << i << " : ";
        fread(ts1, sizeof(float ), Const::tsLength, f);
        for(int j=0;j<Const::maxK;++j){
            fread(ts2, sizeof(float ), Const::tsLength, resf);
//            cout << TimeSeriesUtil::euclideanDist(ts1, ts2, Const::tsLength) << " , ";
            cout << TimeSeriesUtil::dtw(ts1, ts2, Const::tsLength, Const::dtw_window_size, numeric_limits<double>::max()) << " , ";

        }
        cout <<endl;
    }
}

void testExactRes(){
    string result = Const::resfn;
    FILE *f = fopen(Const::queryfn.c_str(), "rb");

//    FILE *f = fopen("/mnt/c/Series4Similarity_Search/rand-256-1k.bin", "rb");
    FILE *resf = fopen(result.c_str(), "rb");
    auto *ts1 = new float[Const::tsLength], *ts2 = new float [Const::tsLength];
    for(int i=0;i<Const::query_num;++i){
        cout << i << " : ";
        fread(ts1, sizeof(float ), Const::tsLength, f);
        for(int j=0;j<Const::maxK;++j){
            fread(ts2, sizeof(float ), Const::tsLength, resf);
            cout << TimeSeriesUtil::euclideanDist(ts1, ts2, Const::tsLength) << " , ";
//            cout << TimeSeriesUtil::dtw(ts1, ts2, Const::tsLength, Const::dtw_window_size, numeric_limits<double>::max()) << " , ";

        }
        cout <<endl;
    }
}

//void sortExactRes(){
//    string result = "../data/results/deep.bin_le";
//    FILE *f = fopen("../data/generator/deep1b-96-100m_query.bin", "rb");
//    FILE *resf = fopen(result.c_str(), "rb");
//    FILE *outf = fopen("../data/results/deep_new", "wb");
//    auto *query = new float[Const::tsLength];
//    vector<PqItemSeries>res(500, PqItemSeries());
//    for(int i=0;i<250;++i){
//        cout << i << " : ";
//        fread(query, sizeof(float ), Const::tsLength, f);
//        for(int j=0;j<500;++j){
//            auto *ts2 = new float [Const::tsLength];
//            fread(ts2, sizeof(float ), Const::tsLength, resf);
//            res[j].ts = ts2;
//            res[j].dist = TimeSeriesUtil::euclideanDist(query, ts2, Const::tsLength);
//        }
//        sort(res.begin(),  res.end(), PqItemSeriesMaxHeap2());
//        for(auto &ts:res) fwrite(ts.ts, sizeof(float ), Const::tsLength, outf);
//        cout <<endl;
//    }
//}

void testReadSax(){
    string saxfn = "/mnt/c/Series4Similarity_Search/rand/sax/rand-256-100m-16.bin";
    FILE *f = fopen(saxfn.c_str(), "rb");
    long f_size = FileUtil::getFileSize(saxfn.c_str()), series_num = f_size / (4 * Const::segmentNum), offset;
    unsigned short saxes[100000];
    fread(saxes, sizeof(unsigned short ), 100000, f);
    for(int i=0;i<100;++i){
        for(int j=0;j<Const::segmentNum;++j){
            cout << saxes[i*Const::segmentNum + j] << ",";
        }
        cout << endl;
    }
    fclose(f);
}

void generateSax(){
    SaxUtil::generateSaxFile(Const::datafn, Const::saxfn);
}

void generatePaa(){
    SaxUtil::generatePaaFile(Const::datafn, Const::paafn);
}

void generateSaxDistribution(){
    string fn = Const::saxfn;
    string output = "../ecg.res";
    DataDistribution::getInvSaxHeadDistribution(fn, output);

}

void generatePrecisionLostDistribution(){
    string qf = "../data/generator/Series_gaussian_randomwalk256_100000000_query.bin";
    string rf = "../data/results/Series_gaussian_randomwalk256_100000000.bin_le";
    DataDistribution::getPrecisionLostDistribution(qf, rf);

}

void buildInMemoryIndexFadas(){
    TardisTreeNode* root = TardisTreeNode::BuildTardisTree(Const::saxfn, false);
    cout << "build finish" << endl;
    root->getIndexStats();
    root->save2Disk(getOtherIndexfn("in-memory"));
}

void progressiveSearchExprResInMemory(){
    TardisTreeNode* root = TardisTreeNode::loadFromDisk(Const::saxfn, getOtherIndexfn("in-memory"), false);
    auto g = loadGraphSkeleton();
    cout << "load finish." << endl;
    Recall::progressiveSearchInMeomoryFADAS(root,g);
}

void statMemoryFadas(){
    TardisTreeNode* root = TardisTreeNode::loadFromDisk(Const::saxfn, getOtherIndexfn("in-memory"), true);
    cout << "load finish" << endl;
    root->getIndexStats();

}

void buildTardisTree(){
    TardisTreeNode* root = TardisTreeNode::BuildTardisTree(Const::saxfn, true);
    TardisTreeNode::mergingTardis(root);
    cout << "build finish" << endl;
    root->getIndexStats();
    root->save2Disk(getOtherIndexfn("tardis"));
}

void buildTardisTreeMat(){
    auto start = chrono::system_clock::now();
    TardisTreeNode* root = TardisTreeNode::BuildTardisTree(Const::saxfn, true);
    TardisTreeNode::mergingTardis(root);
    auto end = chrono::system_clock::now();
    long cpu_build_time = chrono::duration_cast<chrono::microseconds>(end - start).count();
    cout << "skeleton build finish" << endl;
    cout <<"CPU building time = " << cpu_build_time / 1000.0 << "ms" <<endl;
    root->materialize(Const::tardisfn);
    root->getIndexStats();
    root->save2Disk(getOtherIndexfn("tardis"));
}

//void mergeTardisTreeNodes(){
//    TardisTreeNode* root = TardisTreeNode::loadFromDisk("../data/tardis/rand_10k.bin");
//    cout << "load finish." << endl;
//    TardisTreeNode::mergingTardis(root);
//    root->save2Disk("../data/tardis/rand_10k_tardis.bin");
//}

void recallExprResTardis(){
    TardisTreeNode* root = TardisTreeNode::loadFromDisk(Const::saxfn, getOtherIndexfn("tardis"), true);
    cout << "load finish." << endl;
    Recall::doExprWithResTardis(root, Const::resfn, Const::queryfn);
}

void recallExprResIncTardis(){
    TardisTreeNode* root = TardisTreeNode::loadFromDisk(Const::saxfn, getOtherIndexfn("tardis"), true);
    cout << "load finish." << endl;
    Recall::doExprWithResIncTardis(root, Const::resfn, Const::queryfn);
}


void statTardis(){
    TardisTreeNode* root = TardisTreeNode::loadFromDisk(Const::saxfn, getOtherIndexfn("tardis"), true);
    cout << "load finish." << endl;
    root->getIndexStats();
}

//void partitionGraphTree(){
//    string saxfn = "../data/sax/Series_gaussian_randomwalk256_100000000_" + to_string(Const::segmentNum) + ".bin_le";
//    TardisTreeNode* root = TardisTreeNode::loadFromDisk(saxfn, "../data/tardis/rand_10k.bin");
//    cout << "load finish." << endl;
//    int i=0;
//    TardisTreeNode *target;
//    for(auto x:*root->children){
//        if(x.second!= nullptr && !x.second->isLeafNode())
//        {
//            ++i;
//            if(i >=5) {
//                target = x.second;
//                break;
//            }
//        }
//    }
//    Partition::bubble(target->children);
////    cout << "Leaf Node Number = " << root->getLeafNodeNumber(root) << endl;
////    cout << "Tree height = " << root->getHeight(root) << endl;
////    cout << "Graph number = " << root->getGraphNum(root) << endl;
////    auto *g = new Graph();
////    g->loadGraphSkeletonFromFile();
////    Partition::growingGraph(root->children, g);
////    Partition::savePartition(root->children, "../data/partition/rand-10k-2.bin");
////    auto *gp = new GraphPartitioner(g);
////    int total_size = 0;
////    for(auto &iter: *root->children){
////        if(iter.second != nullptr && iter.second->isLeafNode()){
////            g->VERTEXES[iter.first]->cnt = iter.second->size;
////            gp->addVertex(g->VERTEXES[iter.first]);
////            total_size += iter.second->size;
////        }
////    }
////    cout << "Start Partition, total size = " << total_size << ", partition number = " <<total_size / Const::tardis_threshold + 1 <<endl;
////    auto start = chrono::system_clock::now();
////    auto res = gp->doPartition(total_size / Const::tardis_threshold + 1);
////    auto end  =chrono::system_clock::now();
////    cout << "Partition need " << chrono::duration_cast<chrono::microseconds>(end- start).count() << "us."<<endl;
////    GraphPartitioner::savePartition("../data/partition/rand-10k-2-ECO.bin", res);
////    GraphPartitioner::outputPartition("../data/partition_new2.out", res);
//}

static float *t_paa;
static bool comp_tardis(const TardisTreeNode* x, const TardisTreeNode* y){
    if(x == nullptr)    return false;
    if(y == nullptr)    return true;
    return SaxUtil::LowerBound_Paa_iSax(t_paa, x->sax, x->layer) < SaxUtil::LowerBound_Paa_iSax(t_paa, y->sax, y->layer);
}

void buildiSAX(){
    iSAXRoot* root = iSAXRoot::buildIndex(Const::saxfn, Const::series_num);
    root->save2Disk(getOtherIndexfn("iSAX"));
}

void recallExprResiSAX(){
    iSAXRoot* root = iSAXRoot::loadFromDisk(Const::saxfn, getOtherIndexfn("iSAX"), Const::series_num);
    Recall::progressiveSearchInMeomoryiSAX(root);
}

void buildiSAXFuzzy(){
    iSAXRoot* root = iSAXRoot::buildIndexFuzzy(Const::saxfn, Const::paafn);
    root->save2Disk(getOtherIndexfn("fuzzy-iSAX"));
}

void recallExprResiSAXFuzzy(){
    iSAXRoot* root = iSAXRoot::loadFromDisk(Const::saxfn, getOtherIndexfn("fuzzy-iSAX"), -1);
    Recall::doExprWithResiSAX(root, Const::resfn, Const::queryfn);
}

void buildIPGFuzzy(){
    auto g = loadGraphSkeleton();
    IPGNode* root = IPGNode::BuildIPGFuzzy(Const::saxfn, Const::paafn, g, Const::method);
    root->save2Disk(getIdxfnFuzzy());
}

void buildIPGandPartition(){
    auto g = loadGraphSkeleton();
    IPGNode* root = IPGNode::BuildIPG(Const::saxfn, Const::paafn, g, Const::method);
    root->save2Disk(getIdxfn());
}

void statIPG(){
    string paafn = "../data/paa/Series_gaussian_randomwalk256_100m_" + to_string(Const::segmentNum) + ".bin_le";
    string saxfn = "../data/sax/rand-256-100m_16.bin_le";
    string idxfn = "../data/IPG/rw256_100m_16.bin_le";
    IPGNode* root = IPGNode::loadFromDisk(saxfn, paafn, idxfn, false);
    cout << "height : " << root->getHeight() << endl;
    cout << "leaf node number (fake) : "<< root->getFakeLeafNodeNumber() << endl;
}

void recallExprResIPG(){
    IPGNode* root = IPGNode::loadFromDisk(Const::saxfn, getIdxfn());
    IPGNode::rowDataFileName = Const::datafn;
    Recall::doExprWithResFADASNonMat(root, Const::queryfn, Const::resfn);
}

void recallExprResIPGFuzzy(){
    IPGNode* root = IPGNode::loadFromDisk(Const::saxfn, getIdxfnFuzzy());
    IPGNode::rowDataFileName = Const::datafn;
    Recall::doExprWithResFADASNonMat(root, Const::queryfn, Const::resfn);
}

void buildIPGDynamic(){
    auto g = loadGraphSkeleton();
    IPGNode* root = IPGNode::BuildDynamicIPG(Const::saxfn, Const::paafn, g);
    root->save2Disk(getOtherIndexfn("dynamic") + "range2");
}

void recallExprResIPGDynamic(){
    IPGNode* root = IPGNode::loadFromDisk(Const::saxfn, getOtherIndexfn("dynamic") + "range2");
    IPGNode::rowDataFileName = Const::datafn;
    Recall::doExprWithResDynamicIPG(root, Const::queryfn, Const::resfn);
}

void buildIPGGrid(){
    auto g = loadGraphSkeleton();
    IPGNode* root = IPGNode::BuildIPGGrid(Const::saxfn, Const::paafn, g);
    root->save2Disk(getOtherIndexfn("grid"));
}

void recallExprResIPGGrid(){
    IPGNode* root = IPGNode::loadFromDisk(Const::saxfn, getOtherIndexfn("grid"));
    IPGNode::rowDataFileName = Const::datafn;
    Recall::doExprWithResFADASNonMat(root, Const::queryfn, Const::resfn);
}

void buildIPGCluster(){
    auto g = loadGraphSkeleton();
    IPGNode* root = IPGNode::BuildIPGCluster(Const::saxfn, Const::paafn, g);
    root->save2Disk(getOtherIndexfn("cluster"));
}

void recallExprResCluster(){
    IPGNode* root = IPGNode::loadFromDisk(Const::saxfn, getOtherIndexfn("cluster"));
    IPGNode::rowDataFileName = Const::datafn;
    Recall::doExprWithResFADASNonMat(root, Const::queryfn, Const::resfn);
}

void buildIPGClusterDynamic(){
    auto g = loadGraphSkeleton();
    IPGNode* root = IPGNode::BuildIPGClusterDynamic(Const::saxfn, Const::paafn, g);
    root->save2Disk(getOtherIndexfn("cluster-dynamic"));
}

void recallExprResClusterDynamic(){
    IPGNode* root = IPGNode::loadFromDisk(Const::saxfn, getOtherIndexfn("cluster-dynamic"));
    IPGNode::rowDataFileName = Const::datafn;
    Recall::doExprWithResIPGDynamicPaaMu(root, Const::queryfn, Const::resfn);
}

void test3(){
    string fn = "../data/generator/Series_gaussian_randomwalk256_100000000.bin_le";
    string saxfn = "../data/sax/rand-256-100m_16.bin_le";
    string idxfn = "../data/IPG/rw256_100m_16.bin_le";
    IPGNode* root = IPGNode::loadFromDisk(saxfn, idxfn);
    IPGNode::rowDataFileName = fn;
    struct NODE{
        IPGNode* node{};
        double dist{};
        NODE(IPGNode *_node, double _dist):node(_node),dist(_dist){;}
        static bool sort(NODE a, NODE b){
            return a.dist < b.dist;
        }

    };
    vector<NODE>candidates(1, NODE((*root->children)[55672], 0));
    for(auto &iter:*root->children){
        if(iter.second!= nullptr && iter.second->isLeafNode() && iter.second->id != 55672){
            candidates.emplace_back(iter.second, TimeSeriesUtil::distanceEstimation(candidates[0].node, iter.second));
        }
    }
    sort(candidates.begin(), candidates.end(), NODE::sort);
    return;

}

//extern int level = 0;
//void testQueryLayer(){
//    string fn = "/mnt/c/Series4Similarity_Search/rand-256-100m.bin";
//    string saxfn = "../data/sax/rand-256-100m_16.bin_le";
//    string idxfn = "../data/IPG/rw-256-100m_16_fuzzy4-1010_sigma.bin_le";
//    string queryfn = "../data/generator/rand-256-2k.bin";
//    string resfn = "../data/results/rand.bin_le";
//    IPGNode* root = IPGNode::loadFromDisk(saxfn, idxfn);
//    IPGNode::rowDataFileName = fn;
////    Recall::doExprWithResFADASNonMat(root, queryfn, resfn);
//}

void generateRandQuery(){
    string fn = "../data/generator/rand-256-2k.bin";
    RandDataGenerator::generate_random_timeseries(256, 2000, fn.c_str());

}

void translateDNA(){
    string output = "/mnt/c/Series4Similarity_Search/dna/dna-1024.bin";
    DNATranslator::translateNonOverlap(output, 1024);
}

void readSIFT(){
    DNATranslator::readSIFT();
}

void transferPaa2Float(){
    string paafn = "../data/paa/seismic-256-100m-16.bin_le";
    long f_size = FileUtil::getFileSize(paafn.c_str());
    auto *paas = new double [f_size / sizeof(double)];
    FILE *f = fopen(paafn.c_str(), "rb");
    fread(paas, sizeof(double ), f_size / sizeof(double), f);
    fclose(f);
    Const::logPrint( "Finish loading paa");
    string newfn = "../data/paa/seismic-256-100m-16.bin";
    auto *new_paas = new float[f_size / sizeof(double)];
    for(int i=0;i<f_size / sizeof(double);++i)
        new_paas[i] = paas[i];
    delete[]paas;
    f = fopen(newfn.c_str(), "wb");
    fwrite(new_paas, sizeof(float ), f_size / sizeof(double), f);
    fclose(f);
}

void transferSax2Short(){
    string saxfn = "../data/sax/rand-256-100m_16.bin_le";
    long f_size = FileUtil::getFileSize(saxfn.c_str()), series_num = f_size / (4 * Const::segmentNum);
    int *saxes = new int[f_size / 4];
    FILE *f = fopen(saxfn.c_str(), "rb");
    fread(saxes, sizeof(int), f_size / 4, f);
    fclose(f);
    Const::logPrint("Finish loading sax");
    string newfn = "../data/sax/rand-256-100m_16.bin";
    auto *new_saxs = new unsigned short[f_size / 4];
    for(int i=0;i<f_size / 4;++i)
        new_saxs[i] = saxes[i];
    delete[]saxes;
    f = fopen(newfn.c_str(), "wb");
    fwrite(new_saxs, sizeof(unsigned short ), f_size / 4, f);
    fclose(f);
}

void buildFADASPos(){
    auto g = loadGraphSkeleton();
    FADASNode* root = FADASNode::BuildIndexPos(Const::datafn, Const::saxfn, Const::paafn, g);
    root->save2Disk(Const::posidxfn + "root.idx");
}

void exactExprFADASPos(){
    FADASNode* root = FADASNode::loadFromDisk(Const::saxfn, Const::posidxfn + "root.idx", true);
    auto *g = loadGraphSkeleton();
    Recall::exactSearchFADASPos(root,g);
}

void recallExprResFADASPos(){
    FADASNode* root = FADASNode::loadFromDisk(Const::saxfn, Const::posidxfn + "root.idx", true);
    auto *g = loadGraphSkeleton();
    Recall::doExprWithResFADASPos(root,g, Const::posidxfn);
}

void buildFADAS(){
    FADASNode* root = FADASNode::BuildIndex(Const::datafn, Const::saxfn);
    root->save2Disk(Const::fidxfn + "root.idx");
}

void buildTARDISORIGIN(){
    TARGNode* root = TARGNode::buildIndex();
    root->save2Disk(Const::tardisfn + "root.idx");
}

void recallExprResFADAS(){
    FADASNode* root = FADASNode::loadFromDisk(Const::saxfn, Const::fidxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    Recall::doExprWithResFADAS(root, g, Const::fidxfn);
}

void approxSearchTARDISORIGIN(){
    TARGNode* root = TARGNode::loadFromDisk(Const::tardisfn + "root.idx");
    Recall::approxTARDISORIGIN(root);
}

void approxSearchDTWTARDISORIGIN(){
    TARGNode* root = TARGNode::loadFromDisk(Const::tardisfn + "root.idx");
    Recall::approxDTWTARDISORIGIN(root);
}

void approxIncSearchTARDISORIGIN(){
    TARGNode* root = TARGNode::loadFromDisk(Const::tardisfn + "root.idx");
    Recall::approxIncSearchTARDISORIGIN(root);
}

void exactSearchTARDISORIGIN(){
    TARGNode* root = TARGNode::loadFromDisk(Const::tardisfn + "root.idx");
    Recall::exactSearchTARDISORIGIN(root);
}

void statTARDISORIGIN(){
    TARGNode* root = TARGNode::loadFromDisk(Const::tardisfn + "root.idx");
    root->getIndexStats();
}

void exactSearchDTWTARDISORIGIN(){
    TARGNode* root = TARGNode::loadFromDisk(Const::tardisfn + "root.idx");
    Recall::exactSearchTARDISORIGINDTW(root);
}


void recallExprResFADASDTW(){
    FADASNode* root = FADASNode::loadFromDisk(Const::saxfn, Const::fidxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    Recall::doExprWithResFADASDTW(root, g, Const::fidxfn);
}

void recallExprResIncFADASDTW(){
    FADASNode* root = FADASNode::loadFromDisk(Const::saxfn, Const::fidxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    Recall::doExprWithResIncFADASDTW(root, g, Const::fidxfn);
}

void recallExprResIncFADAS(){
    FADASNode* root = FADASNode::loadFromDisk(Const::saxfn, Const::fidxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    Recall::doExprWithResIncFADAS(root, g, Const::fidxfn);
}

void buildFADASFuzzy(){
    auto g = loadGraphSkeleton();
    int bound = Const::boundary * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::max_replica) + "/";
    FADASNode* root = FADASNode::BuildIndexFuzzy(Const::datafn, Const::saxfn, Const::paafn, g);
    root->save2Disk(Const::fuzzyidxfn + "root.idx");
}

void recallExprResFADASFuzzy(){
    int bound = Const::boundary * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::max_replica) + "/";
    FADASNode* root = FADASNode::loadFromDisk(Const::saxfn, Const::fuzzyidxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    Recall::doExprWithResFADAS(root, g, Const::fuzzyidxfn);
}

void recallExprResIncFADASFuzzy(){
    int bound = Const::boundary * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::max_replica) + "/";
    FADASNode* root = FADASNode::loadFromDisk(Const::saxfn, Const::fuzzyidxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    Recall::doExprWithResIncFADASFuzzy(root, g, Const::fuzzyidxfn);
}

void recallExprResFADASFuzzyDTW(){
    int bound = Const::boundary * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::max_replica) + "/";
    FADASNode* root = FADASNode::loadFromDisk(Const::saxfn, Const::fuzzyidxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    Recall::doExprWithResFADASDTW(root, g, Const::fuzzyidxfn);
}

void exactExprFADAS(){
    FADASNode* root = FADASNode::loadFromDisk(Const::saxfn, Const::fidxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    Recall::exactSearchFADAS(root,g);
}

void exactExprFADASDTW(){
    FADASNode* root = FADASNode::loadFromDisk(Const::saxfn, Const::fidxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    Recall::exactSearchFADASDTW(root,g);
}

void exactSearchFADAS(){
    FADASNode* root = FADASNode::loadFromDisk(Const::saxfn, Const::fidxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    Recall::exactSearchFADASNoExpr(root,g);
}

void exactSearchFADASPos(){
    FADASNode* root = FADASNode::loadFromDisk(Const::saxfn, Const::posidxfn + "root.idx", true);
    auto *g = loadGraphSkeleton();
    Recall::exactSearchFADASPosNoExpr(root,g);
}

void statIndexFADAS(){
    FADASNode* root = FADASNode::loadFromDisk(Const::saxfn, Const::fidxfn + "root.idx",false);
    root->getIndexStats();
}

void statIndexFADASFuzzy(){
    int bound = Const::boundary * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::max_replica) + "/";
    FADASNode* root = FADASNode::loadFromDisk(Const::saxfn, Const::fuzzyidxfn + "root.idx",false);
    root->getIndexStats();
}

int main() {
    Const::readConfig();

    switch (Const::index) {
        case 0:
            if(Const::ops == 0) buildInMemoryIndexFadas();
            else if(Const::ops == 1) progressiveSearchExprResInMemory();
            else if(Const::ops == 4)    statMemoryFadas();
            exit(0);
        case 1:
            switch (Const::ops) {
                case 0:
                    if (Const::materialized == 0)
                        buildIPGandPartition();
                    else
                        buildFADAS();
                    break;
                case 1:
                    if(Const::materialized == 0)
                        recallExprResIPG();
                    else
                        recallExprResFADAS();
                    break;
                case 2:
                    if(Const::materialized == 1)
                        exactExprFADAS();
                    break;
                case 3:
                    if(Const::materialized == 1)
                        exactSearchFADAS();
                    break;
                case 4:
                    if(Const::materialized == 1)
                        statIndexFADAS();
                    break;
                case 5:
                    if(Const::materialized == 1)
                        recallExprResIncFADAS();
                    break;
                case 6:
                    if(Const::materialized == 1)
                        recallExprResFADASDTW();
                    break;
                case 7:
                    if(Const::materialized == 1)
                        exactExprFADASDTW();
                    break;
                case 8:
                    if(Const::materialized == 1)
                        recallExprResIncFADASDTW();
                    break;
                default:
                    break;
            }
            exit(0);
        case 2:
            if(Const::ops == 0){
                if(Const::materialized == 0)    buildIPGFuzzy();
                else buildFADASFuzzy();
            }
            else if(Const::ops == 1){
                if(Const::materialized == 0)    recallExprResIPGFuzzy();
                else recallExprResFADASFuzzy();
            }else if(Const::ops == 4){
                if(Const::materialized == 1)    statIndexFADASFuzzy();
            }else if(Const::ops == 5){
                if(Const::materialized == 1)    recallExprResIncFADASFuzzy();
            }else if(Const::ops == 6){
                if(Const::materialized == 1)    recallExprResFADASFuzzyDTW();
            }
            exit(0);
        case 3:
            if(Const::materialized == 1){
                switch (Const::ops) {
                    case 0:
                        buildFADASPos();
                        break;
                    case 1:
                        recallExprResFADASPos();
                        break;
                    case 2:
                        exactExprFADASPos();
                        break;
                    case 3:
                        exactSearchFADASPos();
                        break;
                    default:
                        break;
                }
            }
            break;
        case 4:
            if(Const::ops == 0) buildIPGDynamic();
            else if(Const::ops == 1)    recallExprResIPGDynamic();
            break;
        case 5:
            if(Const::ops == 0) buildIPGGrid();
            else if(Const::ops == 1)    recallExprResIPGGrid();
            break;
        case 6:
            if(Const::ops == 0) buildIPGCluster();
            else if(Const::ops == 1)    recallExprResCluster();
            break;
        case 7:
            if(Const::ops == 0) buildIPGClusterDynamic();
            else if(Const::ops == 1)    recallExprResClusterDynamic();
            break;
        case 8:
            if(Const::ops == 0) buildiSAX();
            else if(Const::ops == 1) recallExprResiSAX();
            exit(0);
        case 9:
            if(Const::ops == 0) {
                if(Const::materialized == 0)
                    buildTardisTree();
                else
                    buildTardisTreeMat();
            }
            else if (Const::ops == 1) recallExprResTardis();
            else if(Const::ops == 4)    statTardis();
            else if(Const::ops == 5){
                if(Const::materialized == 0)
                    recallExprResIncTardis();
            }
            exit(0);
        case 10:
            if(Const::ops == 0) {
                if (Const::materialized == 1)
                    buildTARDISORIGIN();
            }else if(Const::ops == 1){
                if(Const::materialized == 1)
                    approxSearchTARDISORIGIN();
            }else if(Const::ops == 2){
                if(Const::materialized == 1)
                    exactSearchTARDISORIGIN();
            }else if(Const::ops == 4){
                if(Const::materialized == 1)
                    statTARDISORIGIN();
            }
            else if(Const::ops == 5){
                if(Const::materialized == 1)
                    approxIncSearchTARDISORIGIN();
            }else if(Const::ops == 6){
                if(Const::materialized == 1)
                    approxSearchDTWTARDISORIGIN();
            } else if(Const::ops == 7){
                if(Const::materialized == 1)
                    exactSearchDTWTARDISORIGIN();
            }

            exit(0);
        case 11:
            if(Const::ops == 0) buildDSTree();
            else if(Const::ops == 1)    DSTreeExpr();
            exit(0);
        default:    break;
    }

    {//generateQueryFile();
//generateGroundTruth();
//generateGroundTruthDTW();
//generateSax();
//generatePaa();
//testExactRes();
//testExactRes2();
//sortExactRes();
//    buildGraphDebug();
//buildGraphNew();
//recallExprRes();
//translateDNA();
//readSIFT();
//generateRandQuery();

//constructGraph();

//test();
//test2();
//testReadSax();
//generateSaxDistribution();
//generatePrecisionLostDistribution();

//buildNonClusteredIndex();
//recallExprResNonClustered();

//buildTardisTree();
//mergeTardisTreeNodes();
//recallExprResTardis();
//partitionGraphTree();
//testCoverage();
//partitionGraph();

//buildPartitionInvSAXTree();

//buildiSAX();
//recallExprResiSAX();

//buildiSAXFuzzy();
//recallExprResiSAXFuzzy();

//buildIPGandPartition();
//buildIPGFuzzy();
//statIPG();
//recallExprResIPG();
//test3();

//        transferPaa2Float();
//        transferSax2Short();
//testgetFiles();
//string data_dir = "/mnt/c/Series4Similarity_Search/ECG";
//string output = "/mnt/c/Series4Similarity_Search/ECG-320.bin";
//ECGParser::generateECG(data_dir, output, 320);
//FILE *f = fopen(Const::queryfn.c_str(), "rb");
//FILE *outf = fopen("/mnt/c/Series4Similarity_Search/ecg-320-20k.bin_clean", "wb");
//int series_num = FileUtil::getFileSize(Const::queryfn.c_str()) / Const::tsLengthBytes;
//float ts[Const::tsLength];
//for(int i=0;i<series_num;++i){
//    fread(ts, sizeof(float ), Const::tsLength, f);
//    int j=0;
//    for(;j<Const::tsLength;++j)
//        if(isnan(ts[j]))
//            break;
//    if(j == Const::tsLength)
//        fwrite(ts, sizeof(float ), Const::tsLength, outf);
//}
//        fclose(f);
//        fclose(outf);

//        float ts[]{1.04835,0.985939,0.923533,0.985939,1.01714,1.04835,1.11075,1.17316,1.11075,1.17316,1.11075,1.07955,1.01714,0.954736,0.923533,0.861128,0.798722,0.861128,0.892331,0.829925,0.892331,0.954736,0.923533,0.892331,0.829925,0.767519,0.829925,0.892331,0.954736,0.892331,0.923533,0.985939,1.01714,1.07955,1.14195,1.17316,1.20436,1.14195,1.20436,1.23556,1.26677,1.23556,1.26677,1.29797,1.32917,1.36037,1.39158,1.42278,1.45398,1.39158,1.42278,1.36037,1.39158,1.45398,1.42278,1.36037,1.29797,1.36037,1.29797,1.23556,1.17316,1.14195,1.17316,1.14195,1.20436,1.17316,1.20436,1.14195,1.11075,1.07955,1.14195,1.07955,1.01714,0.954736,1.01714,1.07955,1.14195,1.17316,1.23556,1.26677,1.20436,1.14195,1.07955,1.04835,0.985939,0.923533,0.892331,0.923533,0.954736,0.923533,0.861128,0.798722,0.767519,0.705113,0.67391,0.736316,0.67391,0.611504,0.549098,0.517895,0.549098,0.486693,0.549098,0.611504,0.580301,0.549098,0.580301,0.517895,0.580301,0.642707,0.67391,0.736316,0.67391,0.611504,0.549098,0.486693,0.45549,0.486693,0.45549,0.424287,0.45549,0.393084,0.330678,0.268272,0.205866,0.14346,0.0810545,0.0498515,0.0186486,-0.0437572,-0.106163,-0.168569,-0.230975,-0.168569,-0.230975,-0.168569,-0.230975,-0.293381,-0.262178,-0.324584,-0.262178,-0.324584,-0.262178,-0.293381,-0.355787,-0.418192,-0.449395,-0.511801,-0.543004,-0.574207,-0.636613,-0.699019,-0.730222,-0.792627,-0.730222,-0.792627,-0.730222,-0.792627,-0.855033,-0.917439,-0.855033,-0.886236,-0.948642,-0.979845,-1.01105,-0.979845,-0.948642,-1.01105,-1.07345,-1.13586,-1.19827,-1.22947,-1.26067,-1.32308,-1.38548,-1.44789,-1.41669,-1.38548,-1.32308,-1.29187,-1.32308,-1.35428,-1.41669,-1.38548,-1.41669,-1.35428,-1.41669,-1.35428,-1.38548,-1.44789,-1.51029,-1.5727,-1.63511,-1.5727,-1.63511,-1.5727,-1.51029,-1.47909,-1.41669,-1.47909,-1.41669,-1.35428,-1.32308,-1.29187,-1.22947,-1.16706,-1.13586,-1.07345,-1.04225,-0.979845,-1.04225,-1.01105,-0.948642,-0.979845,-1.04225,-1.01105,-1.07345,-1.01105,-1.07345,-1.01105,-0.948642,-0.979845,-0.948642,-0.886236,-0.948642,-0.979845,-1.01105,-0.948642,-0.886236,-0.948642,-1.01105,-1.04225,-0.979845,-0.948642,-0.886236,-0.82383,-0.886236,-0.948642,-0.979845,-0.948642,-0.917439,-0.948642,-0.917439,-0.886236,-0.948642,-0.886236,-0.82383,-0.792627,-0.730222,-0.667816,-0.636613,-0.60541,-0.667816,-0.60541,-0.574207,-0.636613};
//        unsigned short sax[16];
//        SaxUtil::saxFromTs(ts, sax, Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
//        for(auto &a:sax)
//            cout << a <<endl;
//        cout << SaxUtil::invSaxHeadFromSax(sax, Const::bitsCardinality, Const::segmentNum);
    }



}
