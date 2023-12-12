//
// Created by pengwang5 on 2022/2/13.
//

#ifndef MULGIFT_ECGPARSER_H
#define MULGIFT_ECGPARSER_H
#include <string>
using std::string ;
class ECGParser {
public:
    static void generateECG(const string &data_dir, const string &output, int ts_length);

};


#endif //MULGIFT_ECGPARSER_H
