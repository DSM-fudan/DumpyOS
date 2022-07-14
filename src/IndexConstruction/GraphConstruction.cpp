//
// Created by wzy on 2021/8/9.
//

#include <cstdio>
#include <vector>
#include <iostream>
#include <thread>
#include "../../include/DataStructures/GraphConstruction.h"
#include "../../include/Utils/MathUtil.h"
#include "../../include/Utils/FileUtil.h"
int segmentNum  = 16, bitsReserve =7;
string graphFileName =  "../data/RowGraph_" + to_string(segmentNum) + "_" + to_string(bitsReserve) + ".bin";

void func(int start, int end, int neighborNum, int neighborBits, int segmentNum, const string & graphFileName){
    int graphSkeleton[neighborNum];
    FILE *f = fopen((graphFileName + to_string(start)).c_str(), "wb");
    int a = MathUtil::nChooseK(segmentNum, 1), b = MathUtil::nChooseK(segmentNum, 2), c = MathUtil::nChooseK(segmentNum, 3);
    for(int i=start; i< end; ++i)
    {
        if((i - start) % 10000 == 0) cout <<start << ": " << (i - start) <<endl;
        int s = 0;
        if(neighborBits >= 1) {
            MathUtil::get1BitChangedNums(i, segmentNum, graphSkeleton, s);
            s += a;
        }
        if(neighborBits >= 2){
            MathUtil::get2BitsChangedNums(i, segmentNum, graphSkeleton, s);
            s += b;
        }
        if(neighborBits >= 3){
            MathUtil::get3BitsChangedNums(i, segmentNum, graphSkeleton, s);
            s += c;
        }
        if(neighborBits >= 4){
            MathUtil::get4BitsChangedNums(i, segmentNum, graphSkeleton, s);
        }
        fwrite(graphSkeleton, sizeof(int), neighborNum, f);
    }
    fclose(f);

}

void GraphConstruction::buildAndSave2Disk() {

    int arrayLength = 1 << segmentNum, neighborBits = bitsReserve;
    int neighborNum = 0;
    for(int i=1;i<=neighborBits;++i)
        neighborNum += MathUtil::nChooseK(segmentNum, i);
    cout << "neighbor number = " << neighborNum << endl;
    int threadNo = 20, chunk_size = arrayLength / threadNo;
    thread threads[threadNo];
    for(int i=0;i<threadNo - 1;++i)
        threads[i] = thread(func, i * chunk_size, (i + 1) * chunk_size, neighborNum, neighborBits, segmentNum, graphFileName);
    threads[threadNo - 1] = thread(func, (threadNo - 1) * chunk_size, arrayLength, neighborNum, neighborBits, segmentNum, graphFileName);

    for(int i=0;i<threadNo;++i)
        threads[i].join();

    string sources[threadNo];
    for(int i=0;i<threadNo;++i)
        sources[i] = (graphFileName + to_string(i * chunk_size));
    FileUtil::mergeFiles(sources, graphFileName, threadNo);


}

