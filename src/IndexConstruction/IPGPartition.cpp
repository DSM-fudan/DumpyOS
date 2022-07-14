//
// Created by wzy on 2021/11/9.
//

#include "../../include/DataStructures/IPGPartition.h"
#include "../../include/DataStructures/IPGNode.h"
#include "../../include/Utils/MathUtil.h"
#include <iostream>
#include <numeric>
#include <deque>
#include <unordered_set>
static bool comp_size(IPGNode* x, IPGNode* y){
    if(x== nullptr) return false;
    if(y== nullptr) return true;
    return x->size > y->size;
}
struct failNode{
    IPGNode *node{};
    int neighbor_size{};

    failNode(IPGNode* a, int b){ node = a; neighbor_size = b;}
};
static bool comp_fail_node(failNode *x, failNode *y){
    return x->neighbor_size > y->neighbor_size;
}
int IPGPartition::growingGraph(unordered_map<int, IPGNode *> *nodes_map, vector<vector<int>> *g,
                               double filling_factor, int seg_num) {
    if(nodes_map->size() <=2)   return 0;
    vector<IPGNode*>nodes;
    int total_size = 0, node_number;
    for(auto &iter: *nodes_map) {
        if(iter.second->isLeafNode() && iter.second->size < Const::th){
            nodes.push_back(iter.second);
            total_size += iter.second->size;
        }
    }
    node_number = nodes.size();
    if(node_number < 2) return 0;
    if(node_number == 2 && nodes[0]->size + nodes[1]->size > Const::th) return 0;
//    cout << "node number = " << node_number <<", total size = " << total_size << endl;
    sort(nodes.begin(),  nodes.end(), comp_size);
    int pid = 0;
    vector<int>partition_size;
    int k = total_size / Const::th + 1;
    partition_size.reserve(k);

    int id, a = MathUtil::nChooseK(seg_num, 1), b=MathUtil::nChooseK(seg_num, 2), c=MathUtil::nChooseK(seg_num, 3),
    finish_num = 0, finish_size = 0, temp_size = 0, fail_partition_size = 0;
    IPGNode *temp_node;
    vector<IPGNode*>temp, candidates;
    // first loop, build 2-clique
    {
        for(int i=0;i<node_number;++i){
            if(nodes[i]->partitionId != -1)   continue;
            temp_size = nodes[i]->size;  temp.clear(); candidates.clear();  temp.push_back(nodes[i]);
            id = nodes[i]->id;
            for(int j=0; j<a; ++j)
            {
                if((*nodes_map).find((*g)[id][j]) == (*nodes_map).end())    continue;
                temp_node = (*nodes_map)[(*g)[id][j]];
                if(temp_node!= nullptr && temp_node->isLeafNode() && temp_node->partitionId == -1 && temp_node->size < Const::th)
                    candidates.push_back(temp_node);
            }

            sort(candidates.begin(),candidates.end(), comp_size);
            for(IPGNode* node:candidates){
                if(node->size + temp_size > Const::th) continue;
                temp.push_back(node);
                temp_size += node->size;
            }

            // fulfill the partition requirement
            if(temp_size >= filling_factor * Const::th){
                nodes[i]->partitionId = pid;
                for(IPGNode* cur: temp) { cur->partitionId = pid; IPGNode::buildDataNode(temp);}
                partition_size.push_back(temp_size);
                finish_num += temp.size();
                finish_size += temp_size;
                ++pid;
            }
        }
    }

    if(finish_num >= node_number)   return pid;
//    cout<< "After first loop, finish nodes number = " << finish_num <<", finish series number = " << finish_size <<endl
//    <<"Built " << pid <<" partitions, with avg size = " << accumulate(partition_size.begin(), partition_size.end(), 0.0) / (double)pid<<endl;

    vector<IPGNode*>candidates2;
    // second loop, build 4-clique
    {
        for(int i=0;i<node_number;++i){
            if(nodes[i]->partitionId != -1)   continue;
            id = nodes[i]->id;
            temp_size = nodes[i]->size;  temp.clear(); temp.push_back(nodes[i]);  candidates2.clear();
            for(int j=0; j<a; ++j) {
                if((*nodes_map).find((*g)[id][j]) == (*nodes_map).end())    continue;
                temp_node = (*nodes_map)[(*g)[id][j]];
                if(temp_node!= nullptr && temp_node->partitionId == -1 && temp_node->size + temp_size <= Const::th){
                    temp.push_back(temp_node);
                    temp_size += temp_node->size;
                }
            }
            for(int j=a;j<a+b;++j){
                if((*nodes_map).find((*g)[id][j]) == (*nodes_map).end())    continue;
                temp_node = (*nodes_map)[(*g)[id][j]];
                if(temp_node!= nullptr && temp_node->partitionId == -1 && temp_node->size < Const::th){
                    candidates2.push_back(temp_node);
                }
            }

            sort(candidates2.begin(),candidates2.end(), comp_size);
            for(IPGNode* node:candidates2){
                if(node->size + temp_size > Const::th) continue;
                temp.push_back(node);
                temp_size += node->size;
            }

            // fulfill the partition requirement
            if(temp_size >= filling_factor * Const::th){
                nodes[i]->partitionId = pid;
                for(IPGNode* cur: temp) { cur->partitionId = pid; IPGNode::buildDataNode(temp);}
                partition_size.push_back(temp_size);
                finish_num += temp.size();
                finish_size += temp_size;
                ++pid;
            }
        }
    }

    if(finish_num >= node_number)   return pid;
//    cout<< "After second loop, finish nodes number = " << finish_num <<", finish series number = " << finish_size <<endl
//    <<"Built " << pid <<" partitions, with avg size = " << accumulate(partition_size.begin() + one_loop_pid, partition_size.end(), 0.0) / (double)(pid-one_loop_pid)<<endl;
    vector<failNode*>fail_partition_node; fail_partition_node.reserve(Const::vertexNum - finish_num);

    {
        int no_part;
        for(int i=0;i<node_number;++i){
            if(nodes[i]->partitionId != -1)   continue;
            no_part = 0;
            id = nodes[i]->id;
            for(int j=0;j<a+b+c;++j){
                if((*nodes_map).find((*g)[id][j]) == (*nodes_map).end())    continue;
                temp_node = (*nodes_map)[(*g)[id][j]];
                if(temp_node == nullptr || !temp_node->isLeafNode() || temp_node->size > Const::th)    continue;
                if(temp_node->partitionId == -1){
                    no_part+= temp_node->size;
                }
            }
            {
                //            int max_par_size = 0, max_id = -1;

//            for(int j=a;j<a+b;++j){
//                if((*nodes_map).find((*g)[id][j]) == (*nodes_map).end())    continue;
//                temp_node = (*nodes_map)[(*g)[id][j]];
//                if(temp_node == nullptr || !temp_node->isLeafNode() || temp_node->size < Const::tardis_threshold)    continue;
//                if(temp_node->partitionId == -1){
//                    no_part+=temp_node->size;
//                    continue;
//                }
                //                int par_size = partition_size[temp_node->partitionId];
                //                if(par_size + nodes[id]->size > Const::tardis_threshold) {
                //                    continue;
                //                }
                //                if(par_size > max_par_size){
                //                    max_par_size = par_size;
                //                    max_id = temp_node->partitionId;
                //                }
//            }
                //            if(max_id != -1){
                //                nodes[id]->partitionId = max_id;
                //                partition_size[max_id] += nodes[id]->size;
                //                finish_num++; finish_size+=nodes[id]->size;
                //                continue;
                //            }

//            for(int j=a+b;j<a+b+c;++j){
//                if((*nodes_map).find((*g)[id][j]) == (*nodes_map).end())    continue;
//                temp_node = (*nodes_map)[(*g)[id][j]];
//                if(temp_node == nullptr || !temp_node->isLeafNode() || temp_node->size < Const::tardis_threshold)    continue;
//                if(temp_node->partitionId == -1){
//                    no_part+= temp_node->size;
//                    continue;
//                }
                //                int par_size = partition_size[temp_node->partitionId];
                //                if(par_size + nodes[id]->size > Const::tardis_threshold) {
                //                    continue;
                //                }
                //                if(par_size > max_par_size){
                //                    max_par_size = par_size;
                //                    max_id = temp_node->partitionId;
                //                }
//            }
                //            if(max_id != -1){
                //                nodes[id]->partitionId = max_id;
                //                partition_size[max_id] += nodes[id]->size;
                //                finish_num++; finish_size+=nodes[id]->size;
                //                continue;
                //            }
                        }

            fail_partition_node.push_back(new failNode(nodes[i], no_part + nodes[i]->size));
            fail_partition_size += nodes[i]->size;
        }
    }
//    cout<< "After filling stage, finish nodes number = " << finish_num <<", finish series number = " << finish_size <<endl
//    <<"Built " << pid <<" partitions, with avg size = " << accumulate(partition_size.begin(), partition_size.end(), 0.0) / (double)pid<<endl;
//    cout <<"Fail partition node number = " << fail_partition_node.size() << " , total size = " << fail_partition_size << endl;

    sort(fail_partition_node.begin(),  fail_partition_node.end(), comp_fail_node);
    vector<bool>flag(Const::vertexNum, false);
    for(failNode* node:fail_partition_node){
        if(flag[node->node->id])    continue;
        temp_size = node->node->size;   temp.clear();   temp.push_back(node->node);
        node->node->partitionId = pid;  flag[node->node->id] = true;
        finish_num++; finish_size+=node->node->size;
        id = node->node->id;

        int j;
        for(j=0; j<a +b +c && temp_size < Const::th; ++j)
        {
            if((*nodes_map).find((*g)[id][j]) == (*nodes_map).end())    continue;
            temp_node = (*nodes_map)[(*g)[id][j]];
            if(temp_node!= nullptr && temp_node->isLeafNode() && temp_node->partitionId == -1&& temp_node->size + temp_size <= Const::th)
                temp_node->partitionId = pid, flag[temp_node->id] = true, finish_num++, finish_size += temp_node->size, temp_size+=temp_node->size, temp.push_back(temp_node);
        }

        if(temp_size >= filling_factor * Const::th) {
            partition_size.push_back(temp_size);
            IPGNode::buildDataNode(temp);
            ++pid; continue;
        }
        for(int i=4; i <= Const::max_radius&& temp_size < Const::th; ++i){
            int mask = (1U << i) - 1, r, last1;
            do{
                {
                    int cur = mask ^ node->node->id;
                    if((*nodes_map).find(cur) != (*nodes_map).end())    temp_node = (*nodes_map)[cur];
                    else temp_node = nullptr;
                    if(temp_node!= nullptr && temp_node->isLeafNode() && temp_node->partitionId == -1 && temp_node->size + temp_size < Const::th)
                        temp_node->partitionId = pid, flag[temp_node->id] = true, finish_num++, finish_size += temp_node->size, temp_size+=temp_node->size, temp.push_back(temp_node);
                    if(temp_size >= Const::th) break;
                }
                last1 = mask & -mask;
                r = mask + last1;
                mask = (((r ^ mask) >> 2) / last1) | r;
            } while (r < (1 << Const::segmentNum));
        }
        partition_size.push_back(temp_size);
        IPGNode::buildDataNode(temp);
        ++pid;
    }

//    cout<< "After stretch stage, finish nodes number = " << finish_num <<", finish series number = " << finish_size <<endl
//    <<"Built " << pid <<" partitions, with avg size = " << accumulate(partition_size.begin(), partition_size.end(), 0.0) / (double)pid<<endl;

    for(auto* _:fail_partition_node)    delete _;
    return pid;
}

int IPGPartition::findFirstLE(vector<IPGNode*>&nodes, int start, int end, int target){
    int mid, nn = end + 1;
    while (start < end){
        mid = (start + end) / 2;
        if(nodes[mid]->size > target){
            start = mid + 1;
        }else{
            end = mid;
        }
    }
    if(start == end && nodes[start]->size <= target)    return  start;
    return  nn;
}


int IPGPartition::balanceMatch(unordered_map<int, IPGNode *> *nodes_map){
    vector<IPGNode*>nodes;
    int total_size = 0, node_number;
    for(auto &iter: *nodes_map) {
        if(iter.second->isLeafNode() && iter.second->size < Const::th){
            nodes.push_back(iter.second);
            total_size += iter.second->size;
        }
    }

    node_number = nodes.size();
    if(node_number < 2) return 0;
    if(node_number == 2 && nodes[0]->size + nodes[1]->size > Const::th) return 0;
//    cout << "node number = " << node_number <<", total size = " << total_size << endl;
    sort(nodes.begin(),  nodes.end(), comp_size);
    int pid = 0;
    vector<IPGNode*>tmp;
    for(int i=0;i<node_number;++i){
        if(nodes[i]->partitionId != -1) continue;
        tmp.clear();
        int target = Const::th - nodes[i]->size;
        int start = findFirstLE(nodes, i + 1, node_number - 1, target);
        nodes[i]->partitionId = pid;    tmp.push_back(nodes[i]);
        for(int j=start; j < node_number; ++j){
            if(nodes[j]->partitionId != -1 || nodes[j]->size > target) continue;
            nodes[j]->partitionId = pid;
            target -= nodes[j]->size;
            tmp.push_back(nodes[j]);
        }
        IPGNode::buildDataNode(tmp);
        ++pid;
    }
    return pid;
}



static bool comp_node(IPGNode *x, IPGNode *y){
    return x->size > y->size;
}
struct Candidate{
    int id{};
    int dist{};

    Candidate(){id = -1; dist = -1;}
    Candidate(int _id, int _dist){id = _id; dist = _dist;}
};
static bool cand_min_heap(Candidate * a, Candidate *b){
    return a->dist > b->dist;
}
static bool cand_comp(Candidate& a, Candidate& b){
    return a.dist < b.dist;
}

//void IPGPartition::bubble(unordered_map<int, IPGNode *> *nodes_map) {
//    vector<IPGNode*>nodes;nodes.reserve(Const::vertexNum);
//    int node_number=0, total_size = 0;
//    for(auto &iter:*nodes_map) {
//        if (iter.second != nullptr && iter.second->isLeafNode()) {
//            nodes.push_back(iter.second);
//            total_size += iter.second->size;
//            ++node_number;
//        }
//    }
//    if(node_number <= 2) return;
//    cout << "node number = " << node_number <<", total size = " << total_size << endl;
//    sort(nodes.begin(),  nodes.end(), comp_node);
//    int k = total_size / Const::th + 1;
//    int pid = 1;
//    unordered_map<int, vector<Candidate>>dist_matrix;
//    vector<Candidate*>seeds(1, new Candidate(0, nodes[0]->size));
//    nodes[0]->partitionId = 0;
//    vector<long>dist_product(node_number, 1);
//
//    // seeds pre-selection for large nodes
//    // compute distance for the biggest size node, O(n)
//    dist_matrix[0].resize(node_number);
//    int sum_dist = 0;
//    for(int i=0;i<node_number;++i){
//        int dist = MathUtil::bitDiffNum(nodes[0]->id, nodes[i]->id);
//        sum_dist += dist;
//        dist_matrix[0][i] = Candidate(i, dist);
//        dist_product[i] *= dist;
//    }
//    double avg_dist = (double)sum_dist / (node_number - 1);
//    for(int i=1;i<node_number && nodes[i]->size > Const::th * Const::large_node_size_ratio && pid < k; ++i){
//        int j=0;
//        for(;j<pid;++j)
//            if(MathUtil::bitDiffNum(nodes[seeds[j]->id]->id, nodes[i]->id) < avg_dist)  break;
//            if(j < pid) continue;
//            seeds.push_back(new Candidate(i, nodes[i]->size));
//            nodes[i]->partitionId = pid;
//            ++pid;
//    }
//
//    // compute distance matrix for seeds
//    for(int i=1;i<pid;++i){
//        int seed_id = seeds[i]->id;
//        dist_matrix[seed_id].resize(node_number);
//        vector<Candidate>&tmp_dist_list = dist_matrix[seed_id];
//
//        for(int j=0;j<node_number;++j){
//            int dist = MathUtil::bitDiffNum(nodes[seed_id]->id, nodes[j]->id);
//            tmp_dist_list[j] = Candidate(j, dist);
//            dist_product[j] *= dist;
//        }
//    }
//
//    // seed selection   O(2kn)
//    for(;pid < k;++pid){
//        // find the largest dist sum(if same choose larger size) node, O(n)
//        long max_sum_dist = 0;
//        int max_id = -1;
//        for(int i=0;i<node_number;++i){
//            if(nodes[i]->partitionId != -1) continue;
//            if(dist_product[i] > max_sum_dist ||
//            (dist_product[i] == max_sum_dist && nodes[i]->size > nodes[max_id]->size))
//                max_sum_dist = dist_product[i], max_id = i;
//        }
//
//        // determine seed
//        seeds.push_back(new Candidate(max_id, nodes[max_id]->size));
//        nodes[max_id]->partitionId = pid;
//
//        dist_matrix[max_id].resize(node_number);
//        vector<Candidate>&tmp_dist_list = dist_matrix[max_id];
//
//        // compute distance according to chosen node, O(n)
//        for(int i=0;i<node_number;++i){
//            int dist = MathUtil::bitDiffNum(nodes[max_id]->id, nodes[i]->id);
//            tmp_dist_list[i] = Candidate(i, dist);
//            dist_product[i] *= dist;
//        }
//    }
//
//    // sort dist matrix, O(knlogn)
//    for(auto &dist_list:dist_matrix)
//        sort(dist_list.second.begin(),  dist_list.second.end(), cand_comp);
//
//    // little size first bubble, O(kn)
//    // first assign pid for node neighborhood
//    int pos[k];
//    for(int &_:pos) _ = 1;
//    int par_size[k];
//    int finish_num = pid;
//    int seed_ids[k];
//    for(int i=0;i<k;++i)    seed_ids[i] = seeds[i]->id;
//    for(Candidate* seed:seeds) {
//        int seed_id = seed->id, par_id = nodes[seed_id]->partitionId, tmp_size= seed->dist;
//        int work = pos[par_id];
//        vector<Candidate>& dist_list = dist_matrix[seed_id];
//        for(;work < node_number && tmp_size < Const::th && dist_list[work].dist <= Const::bitsReserve; ++work){
//            nodes[dist_list[work].id]->partitionId = par_id;
//            tmp_size += nodes[dist_list[work].id]->size;
//            ++finish_num;
//        }
//        pos[par_id] = work;
//        par_size[par_id] = tmp_size;
//    }
//    // little size first bubble
//    for(; finish_num < node_number && !seeds.empty(); )
//    {
//        make_heap(seeds.begin(), seeds.end(), cand_min_heap);
//        int cur_id = seeds[0]->id, cur_par = nodes[cur_id]->partitionId;
//        int work = pos[cur_par];
//        for(;work< node_number;++work){
//            int neighbor_id = dist_matrix[cur_id][work].id;
//            if(nodes[neighbor_id]->partitionId != -1 || nodes[neighbor_id]->size + par_size[cur_par] > Const::th)   continue;
//            nodes[neighbor_id]->partitionId = cur_par;
//            par_size[cur_par] += nodes[neighbor_id]->size;
//            seeds[0]->dist = par_size[cur_par];
//            ++finish_num;
//            break;
//        }
//        if(work < node_number)  pos[cur_par] = work;
//        else {
//            delete seeds[0];
//            seeds.erase(seeds.begin());
//        }
//    }
//
//    cout << "After first bubble, finish series number = " << finish_num << endl;
//    // resolve rest nodes, assign partition as near as possible
//    if(finish_num < node_number){
//        for(int i=0;i<node_number;++i){
//            if(nodes[i]->partitionId != -1) continue;
//            int min_dist = Const::segmentNum, min_seed = -1;
//            for(int seed_id:seed_ids){
//                int dist = MathUtil::bitDiffNum(nodes[seed_id]->id, nodes[i]->id);
//                if(dist < min_dist || (dist == min_dist && par_size[nodes[seed_id]->partitionId] < par_size[nodes[min_seed]->partitionId])){
//                    min_dist = dist;    min_seed = seed_id;
//                }
//            }
//            nodes[i]->partitionId = nodes[min_seed]->partitionId;
//            par_size[nodes[i]->partitionId] += nodes[i]->size;
//        }
//    }
//
//
//}

void IPGPartition::zOrderCurve(unordered_map<int, IPGNode *> *nodes_map){
    cout << "z-order-partition"<<endl;
    if(nodes_map->size() <=2)   return;
    vector<IPGNode*>nodes;
    int total_size = 0, node_number;
    for(auto &iter: *nodes_map) {
        if(iter.second->isLeafNode() && iter.second->size < Const::th){
            nodes.push_back(iter.second);
            total_size += iter.second->size;
        }
    }
    node_number = nodes.size();
    if(node_number < 2) return;
    if(node_number == 2 && nodes[0]->size + nodes[1]->size > Const::th) return;

    int pid=0, tmp_size = 0;
    for(int i=0;i<(1<<Const::segmentNum);++i){
        if(nodes_map->find(i) == nodes_map->end())  continue;
        auto& cur = (*nodes_map)[i];
        if(cur->isLeafNode() && cur->size < Const::th){
            if(tmp_size + cur->size <= Const::th){
                cur->partitionId = pid;
                tmp_size += cur->size;
            }else{
                ++pid;
                cur->partitionId = pid;
                tmp_size = cur->size;
            }
        }
    }

}

void IPGPartition::savePartition(unordered_map<int, IPGNode *> *nodes_map, const string& output){
    FILE *f = fopen(output.c_str(), "wb");
    for(auto& iter: *nodes_map){
        if(iter.second != nullptr && iter.second->isLeafNode())
        {
            fwrite(&iter.first, sizeof(int), 1, f);
            fwrite(&iter.second->partitionId, sizeof(int),1,f);
        }
    }

    fclose(f);
}

void IPGPartition::outputPartition(unordered_map<int, IPGNode *> *nodes_map, const string& output, int k){
    ofstream f(output.c_str(),ios::out);
    vector<vector<int>>graph(k,vector<int>());
    vector<int>node_number(k,0);
    vector<int>sum_size(k,0);
    for(auto &iter:*nodes_map){
        if(iter.second!= nullptr && iter.second->isLeafNode()){
            int pid = iter.second->partitionId;
            graph[pid].push_back(iter.second->id);
            node_number[pid]++;
            sum_size[pid] += iter.second->size;
        }
    }
    for(int i=0;i<k;++i){
        f << "IPGPartition " << i << ": Total node number = " << node_number[i] << ", Total size = "<< sum_size[i] <<endl;
        for(int id:graph[i]){
            f << id << ": " << (*nodes_map)[id]->size << ", ";
        }
        f<<endl;
    }
    f.close();
}

