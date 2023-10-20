//
// Created by andrea on 20/12/22.
//#

#ifndef DECREMENTAL_TOPK_H
#define DECREMENTAL_TOPK_H
#include <vector>
#include <tuple>
#include <sys/time.h>
#include <stdint.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <queue>
#include <set>
#include <limits>
#include "progressBar.h"
#include "networkit/graph/Graph.hpp"
#include <string>
#include "mytimer.h"

using vertex = uint32_t;
using dist = uint8_t;
using edge_id = uint32_t;
using label_ds = std::vector<std::pair<vertex,dist>>;
using indexed_paths = std::vector<std::vector<std::vector<vertex>>>;

class DecrementalTopK{
        struct index_t{
            label_ds label_offset;
            std::vector<indexed_paths> p_array;
            index_t(){
                label_offset.clear();
                p_array.clear();
            }
            ~index_t(){
                label_offset.clear();
                p_array.clear();
            }
        };

public:
    DecrementalTopK(NetworKit::Graph*, dist, bool, dist, bool);

    ~DecrementalTopK();

    void build();

    // int query(vertex, vertex, dist, std::vector<dist> &);
    void query(vertex, vertex, std::vector<dist> &);
    // int query(vertex, vertex, dist);
    // int query(vertex, vertex);
    

    vertex x;
    vertex y;
    void update_loops(bool);
    void update_lengths();
    void incremental_lengths();



    double n_reached_nodes();
    double n_reached_nodes_mbfs();
    // void mod_bfs(vertex, vertex, std::vector<dist>&);


    double loops_time;
    double lengths_time;
    uint64_t total_bits;
    vertex aff_hubs;
    vertex aff_cycles;
    vertex total_pruning;
    void deallocate_aux();

private:
    std::vector<std::vector<vertex>> dists;
    std::vector<vertex> updated;
    std::queue<std::pair<std::vector<vertex>, bool>> new_labels;
    std::set<vertex> vertices_to_update;
    bool is_from_scratch_only;
    std::vector<dist> tmp_v;
    void verify_sizes();
    void pruned_bfs (vertex, bool);
    void reset_temp_vars(vertex, bool);
    void reset_temp_update_vars(vertex, bool);
    void set_temp_vars(vertex, bool);
    void set_temp_update_vars(vertex, bool);
    size_t prune(vertex,  dist, bool);
    void compute_loop_entries(vertex);
    void resume_pbfs(vertex, bool);
    void incremental_resume_pbfs(vertex, vertex, std::vector<vertex> &, bool);
    void allocate_label(vertex, vertex, std::vector<vertex>&, bool);
    void extend_label(vertex, vertex, std::vector<vertex>&, bool, size_t);
    void extend_label_repair(vertex, vertex, std::vector<vertex>&, bool);
    void incremental_extend_label_repair(vertex, vertex, std::vector<vertex>&, bool);

    static const dist null_distance;
    static const vertex null_vertex;
    std::queue<std::pair<vertex, std::vector<vertex>>> * node_que;
    std::vector<std::pair<vertex,vertex>> deleted;
    vertex* ordering;
    vertex* reverse_ordering; 
    std::pair<double,vertex>* ordering_rank;
    vertex total;
    NetworKit::Graph * graph;
    dist K;
    bool directed;
    dist ordering_type;

    label_ds old_label_a;
    label_ds old_label_b;
    std::vector<indexed_paths> old_paths_a;
    std::vector<indexed_paths> old_paths_b;

    std::vector<vertex> visited_in_update_lengths_from_a;
    std::vector<vertex> visited_in_update_lengths_from_b;

    std::vector<vertex> reached_nodes;
    std::vector<vertex> reached_mbfs;

    indexed_paths* loop_labels;
    index_t ** length_labels;

    bool* tmp_pruned;
    dist* tmp_offset;
    vertex* tmp_count;
    std::vector<dist>*  tmp_dist_count;
    dist*  tmp_s_offset;
    std::vector<dist>* tmp_s_count;
    dist* visited_in_update_loops;
    std::set<vertex> union_of_reached_nodes;
    std::vector<std::pair<dist,bool>>* update_temp_dist_flag;
    std::vector<dist> auxiliary_prune;
};
#endif //DECREMENTAL_TOPK_H
