//
// Created by wzy on 2021/9/21.
//

#ifndef MULGIFT_DATADISTRIBUTION_H
#define MULGIFT_DATADISTRIBUTION_H


#include <string>

class DataDistribution {
public:
    static void
    getInvSaxHeadDistribution(const std::__cxx11::basic_string<char> &fn,
                              const std::string &outfn);

    static void getPrecisionLostDistribution(const std::string &query_file, const std::string &result_file);
};


#endif //MULGIFT_DATADISTRIBUTION_H
