//
// Created by wzy on 2021/8/11.
//

#ifndef MULGIFT_RECALL_H
#define MULGIFT_RECALL_H
#include "../Tardis/TardisTreeNode.h"
#include "../DataStructures/IPGNode.h"
#include "../DataStructures/iSAXNode.h"
#include "../DataStructures/FADASNode.h"
#include "../DSTree/DSTreeNode.h"
#include "../TAR/TARGNode.h"

class Recall {

public:
//    static void doExpr(Graph *g);
//
//    static void doExprDebug(Graph *g);
//
//    static void doExprPerformance(Graph *g);
//
//    static void doExprWithRes(Graph *g);


    static vector<float *> *getResult(const string &fn, int queryNo, int k);

    static void doExprWithResTardis(TardisTreeNode *root, const string &resFile, const string &queryFile);

    static void doExprWithResFADASNonMat(IPGNode *root, const string &queryFile, const string &resFile);

    static void doExprWithResIPGDynamicPaaMu(IPGNode *root, const string &queryFile, const string &resFile);

    static void doExprWithResiSAX(iSAXRoot *root, const string &resFile, const string &queryFile);

    static void doExprWithResInMemory(TardisTreeNode *root, const string &resFile, const string &queryFile);

    static void doExprWithResIPGWO1st(IPGNode *root, const string &queryFile, const string &resFile);

    static void doExprWithResDynamicIPG(IPGNode *root, const string &queryFile, const string &resFile);

    static void doExprWithResFADAS(FADASNode *root, vector<vector<int>> *g, const string &index_dir);

    static void ngSearchDumpy(FADASNode *root, vector<vector<int>> *g);

    static void exactSearchFADAS(FADASNode *root, vector<vector<int>> *g);

    static void exactSearchFADASDTW(FADASNode *root, vector<vector<int>>*g);

    static void exactSearchFADASNoExpr(FADASNode *root, vector<vector<int>> *g);

    static void exactSearchFADASPos(FADASNode *root, vector<vector<int>> *g);

    static void doExprWithResFADASPos(FADASNode *root, vector<vector<int>> *g, const string &index_dir);

    static void exactSearchFADASPosNoExpr(FADASNode *root, vector<vector<int>> *g);

    static void approxIncSearchTARDISORIGIN(TARGNode *root);

    static void progressiveSearchInMeomoryFADAS(TardisTreeNode *root, vector<vector<int>> *g);

    static void progressiveSearchInMeomoryDSTree(DSTreeNode *root);

    static void doExprWithResIncFADAS(FADASNode *root, vector<vector<int>> *g, const string &index_dir);

    static void doExprWithResFADASDTW(FADASNode *root, vector<vector<int>> *g, const string &index_dir);

    static void doExprWithResIncFADASDTW(FADASNode *root, vector<vector<int>> *g, const string &index_dir);

    static void progressiveSearchInMeomoryiSAX(iSAXRoot *root);

    static void doExprWithResIncFADASFuzzy(FADASNode *root, vector<vector<int>> *g, const string &index_dir);

    static void doExprWithResIncTardis(TardisTreeNode *root, const string &resFile, const string &queryFile);

    static void exactSearchTARDISORIGIN(TARGNode *root);

    void doExprWithResIncFADASFuzzyDTW(FADASNode *root, vector<vector<int>> *g, const string &index_dir);

    static void approxTARDISORIGIN(TARGNode *root);

    static void approxDTWTARDISORIGIN(TARGNode *root);

    static void exactSearchTARDISORIGINDTW(TARGNode *root);

    static void completeWorkload();

    static void ngSearchDumpyFuzzy(FADASNode *root);

    static void ngSearchTARDISORIGIN(TARGNode *root);
};


#endif //MULGIFT_RECALL_H
