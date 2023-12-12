//
// Created by caucher on 2021/12/3.
//

#ifndef MULGIFT_DNATRANSLATOR_H
#define MULGIFT_DNATRANSLATOR_H


#include <string>
using namespace  std;

class DNATranslator {
public:
    static void translate(const string &output, int len);

    static void z_normalize_and_save(const float *ts, int size, FILE *f);

    static void readSIFT();

    static void z_normalize(float *ts, int size);

    static void translateNonOverlap(const string &output, int len);
};


#endif //MULGIFT_DNATRANSLATOR_H
