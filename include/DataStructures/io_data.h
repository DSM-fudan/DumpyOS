//
// Created by zeyuwang on 2022/12/21.
//

#ifndef FADAS_IO_DATA_H
#define FADAS_IO_DATA_H

struct io_data{
    int node_size{};
    float lb_dist{};
    int fd{};
    FILE*f{};
    float* tss{};
};

#endif //FADAS_IO_DATA_H
