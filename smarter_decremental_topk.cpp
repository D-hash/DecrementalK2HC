//
// created by anonym 20/06/23
//#
#include "smarter_decremental_topk.h"
#include <algorithm>
#include <cassert>
#include <climits>
#include <networkit/centrality/EstimateBetweenness.hpp>
#include <networkit/centrality/KPathCentrality.hpp>
#include "networkit/centrality/DegreeCentrality.hpp"

using namespace std;

const vertex DecrementalTopK::null_vertex = round(std::numeric_limits<vertex>::max()/2);
const dist DecrementalTopK::null_distance = round(std::numeric_limits<dist>::max()/2);

DecrementalTopK::DecrementalTopK(NetworKit::Graph* g, dist k_value, bool dir, dist ord, bool is_from_scratch){

    this->graph = g;
    this->K = k_value;
    this->directed = dir;
    this->ordering_type = ord;
    this->is_from_scratch_only = is_from_scratch;


    
    this->loops_time = 0.0;
    this->lengths_time = 0.0;
    
    this->ordering_rank = new std::pair<double,vertex>[graph->numberOfNodes()];
    this->ordering = new vertex[graph->numberOfNodes()];
    this->reverse_ordering = new vertex[graph->numberOfNodes()];

    this->graph->parallelForNodes([&] (vertex i){
        assert(graph->hasNode(i));
        this->ordering[i]=null_vertex;
        this->reverse_ordering[i]=null_vertex;
        this->ordering_rank[i] = {null_vertex,null_vertex};


    });
  
    
    double centr_time = 0.0;
    this->total_bits = 0;

    if(this->ordering_type==0){

        INFO("BY DEGREE");       
        mytimer local_constr_timer;
        local_constr_timer.restart();
        const NetworKit::Graph& hand = *graph;
        NetworKit::DegreeCentrality* rank = new NetworKit::DegreeCentrality(hand);
        rank->run();
        this->graph->forNodes([&] (vertex i){
            assert(graph->hasNode(i));
            this->ordering_rank[i]=std::make_pair(rank->score(i),i);

        });
        delete rank;
        centr_time = local_constr_timer.elapsed();


    }
    else if(this->ordering_type==1){
        INFO("BY APX BETW");
        mytimer local_constr_timer;
        double max_time = 30.0;
        double cumulative_time = 0.0;
        double fract = 0.66;
        double n_samples =  round(std::pow((double)graph->numberOfNodes(),fract));

        const NetworKit::Graph& hand = *graph;
        while(cumulative_time<max_time && n_samples<(double)graph->numberOfNodes()){
            local_constr_timer.restart();

            std::cout<<"fract: "<<fract<<" "<<n_samples<<" SAMPLES\n";
            NetworKit::EstimateBetweenness* rank = new NetworKit::EstimateBetweenness(hand,n_samples,false,true);

            
            rank->run();
            
            this->graph->forNodes([&] (vertex i){
                assert(graph->hasNode(i));
                assert(i<graph->numberOfNodes());
                this->ordering_rank[i]=std::make_pair(rank->score(i),i);

            });
            delete rank;
            cumulative_time+=local_constr_timer.elapsed();
            n_samples*=2;
        }
        centr_time = cumulative_time;


    }
    else{
        assert(this->ordering_type==2);
        INFO("BY kPATH");
        mytimer local_constr_timer;
        local_constr_timer.restart();
        const NetworKit::Graph& hand = *graph;

        NetworKit::KPathCentrality* rank = new NetworKit::KPathCentrality(hand,0.0,round(std::pow((double)graph->numberOfNodes(),0.3)));

            
        rank->run();
            
        this->graph->forNodes([&] (vertex i){
            assert(graph->hasNode(i));
            assert(i<graph->numberOfNodes());
            this->ordering_rank[i]=std::make_pair(rank->score(i),i);

        });
        delete rank;

        centr_time = local_constr_timer.elapsed();
        //TODO : USE IT


    }
        

    std::sort(this->ordering_rank, this->ordering_rank+graph->numberOfNodes(), [](const std::pair<double,vertex>  &a, const std::pair<double,vertex>  &b) {
        if(a.first == b.first)
            return a.second > b.second;
        else{
            return a.first > b.first;
        }
    });
    

    for(size_t count = 0; count < graph->numberOfNodes();count++){
        this->reverse_ordering[count]=this->ordering_rank[count].second;
        this->ordering[this->ordering_rank[count].second]=count;
    }

    for(size_t count = 0; count < 10 ;count++)
        std::cout<<"In position "<<count<<" we have vertex "<<this->reverse_ordering[count]<<" rank "<<this->ordering_rank[count].first<<std::endl;

    
    delete[] this->ordering_rank;

};
DecrementalTopK::~DecrementalTopK(){
    delete[] this->ordering;
    delete[] this->reverse_ordering;
    delete[] this->loop_labels;
    


    for (size_t dir = 0; dir < 1 + directed; dir++){
        delete[] this->length_labels[dir];

    }
    delete[] this->length_labels;


};

void DecrementalTopK::deallocate_aux() {
    
    delete[] this->tmp_pruned;
    delete[] this->tmp_offset;
    delete[] this->tmp_count;
    delete[] this->tmp_s_offset;
    delete[] this->tmp_s_count;
    delete[] this->update_temp_dist_flag;
    if(!this->is_from_scratch_only){
        delete[] this->visited_in_update_loops;
    }
    delete[] this->tmp_dist_count;



}


void DecrementalTopK::build(){
    total_pruning = 0;
    this->tmp_pruned = new bool[graph->numberOfNodes()];
    this->tmp_offset = new dist[graph->numberOfNodes()];
    this->tmp_count = new vertex[graph->numberOfNodes()];
    auxiliary_prune.resize(50);
    this->tmp_dist_count = new std::vector<dist>[2];
    for (size_t j = 0; j < 2; j++){
        this->tmp_dist_count[j].resize(graph->numberOfNodes(), 0);
    }

    this->tmp_s_offset = new dist[graph->numberOfNodes()+1];
    this->tmp_s_count = new std::vector<dist>[graph->numberOfNodes()];
    this->update_temp_dist_flag = new std::vector<std::pair<dist,bool>>[graph->numberOfNodes()];
    this->graph->parallelForNodes([&] (vertex i){
        assert(graph->hasNode(i));
        this->tmp_pruned[i]=false;        
        this->tmp_offset[i]=null_distance;
        this->tmp_count[i]=0;

        this->tmp_s_offset[i] = null_distance;
        this->tmp_s_count[i].clear();


    });
    this->tmp_s_offset[graph->numberOfNodes()] = 0;


    
    
    if(!this->is_from_scratch_only){
        this->visited_in_update_loops = new dist[graph->numberOfNodes()];
        this->graph->parallelForNodes([&] (vertex i){
            assert(graph->hasNode(i));
            this->visited_in_update_loops[i]=null_distance;
        });
    }

    this->reached_mbfs.clear();

    this->updated.clear();

    this->length_labels = new index_t*[2];

    this->loop_labels = new indexed_paths [graph->numberOfNodes()];

    for (size_t dir = 0; dir < 1 + directed; dir++){
        this->length_labels[dir] = new index_t[graph->numberOfNodes()];

        for (size_t v = 0; v < graph->numberOfNodes(); v++){
            this->loop_labels[v].clear();
            this->length_labels[dir][v].label_offset.clear();
            this->length_labels[dir][v].p_array.clear();
        }
    }

    #ifndef NDEBUG
        for (size_t v = 0; v < graph->numberOfNodes(); v++)
        assert(loop_labels[v].size()==0);
    #endif
    
    
    mytimer build_timer;
    build_timer.restart();
    ProgressStream loop_bar(graph->numberOfNodes());

    loop_bar.label() << "Loops construction";
    
    for(size_t v = 0; v < graph->numberOfNodes(); v++){
        this->compute_loop_entries(v);
        ++loop_bar;
    }


    loops_time = build_timer.elapsed();

    build_timer.restart();
    ProgressStream index_bar(graph->numberOfNodes());
    index_bar.label() << "Index construction";
    
    for(size_t v = 0; v < graph->numberOfNodes(); v++){

        this->pruned_bfs(v, false);
        ++index_bar;
        if (directed){
            this->pruned_bfs(v, true);
        }
    }
    lengths_time = build_timer.elapsed();

//    #ifndef NDEBUG
//    #verify_sizes();
//    #endif

}

inline void DecrementalTopK::compute_loop_entries(vertex s){

    this->node_que = new std::queue<std::pair<vertex, std::vector<vertex>>>;
    std::vector<vertex> p = {s};
    this->node_que->push(std::make_pair(s,p));
    this->updated.push_back(s);
    vertex v, to_vert;
    dist distance;
    std::pair<vertex,std::vector<vertex>> pvv;
    vertex num_reached = 0;
    size_t count = 0;
    while (!this->node_que->empty()) {
        pvv = this->node_que->front();
        v = pvv.first;
        distance = pvv.second.size()-1;
        this->node_que->pop();
        if (v == s){
            if(pvv.second.size()-1 >= loop_labels[s].size()){
                loop_labels[s].resize(pvv.second.size(), {});
            }
            this->loop_labels[s][pvv.second.size()-1].push_back({pvv.second});
            this->total_bits += 32*pvv.second.size();
            count ++;
            if(count == K) break;
        }
        num_reached++;
        for(vertex u : graph->neighborRange(reverse_ordering[v])){

            to_vert = ordering[u];
            if (this->tmp_count[to_vert] == 0){
                // assert(std::find(this->updated.begin(),this->updated.end(),to_v)==this->updated.end());
                this->updated.push_back(to_vert);
            }

            if (to_vert >= s && this->tmp_count[to_vert] < this->K){
                this->tmp_count[to_vert] += 1;
                std::vector<vertex> temp(pvv.second.begin(),pvv.second.end());
                temp.push_back(to_vert);
                this->node_que->push(std::make_pair(to_vert,temp));
            }
        }
    }
    delete node_que;
//    for(size_t i = 1; i < loop_labels[s].size(); i++){
//        loop_labels[s][i] += loop_labels[s][i-1];
//    }
//    assert(loop_labels[s][0] == 1);
    
    this->reached_mbfs.push_back(num_reached);
    
    reset_temp_vars(s, false);
}


//inline void DecrementalTopK::pruned_bfs(vertex s, bool reversed){
//
//    set_temp_vars(s, reversed);
//
//    vertex curr = 0;
//    vertex next = 1;
//    dist   distance = 0;
//
//    this->node_que = new std::queue<std::pair<vertex, std::vector<vertex>>>[2];
//    this->updated.clear();
//    std::vector<vertex> p = {s};
//    this->node_que[curr].push(std::make_pair(s, p));
//    this->tmp_dist_count[curr][s] = 1;
//    this->updated.push_back(s);
//    std::pair<vertex, std::vector<vertex>> pvv;
//    vertex v, to_vert;
//    dist c;
//
//    for (;;){
//
//        while (!this->node_que[curr].empty()){
//
//            pvv = this->node_que[curr].front();
//            v = pvv.first;
//            this->node_que[curr].pop();
//            c = this->tmp_dist_count[curr][v];
//            this->tmp_dist_count[curr][v] = 0;
//
//            if(c == 0 || tmp_pruned[v]){
//                continue;
//            }
//            this->tmp_pruned[v] = prune(v, distance, reversed);
//
//            if(this->tmp_pruned[v]){
//                continue;
//            }
//
//            if(this->tmp_offset[v] == null_distance){
//                this->tmp_offset[v] = distance;
//                allocate_label(v, s, pvv.second, reversed);
//            }
//            else{
//                extend_label(v, s, pvv.second, reversed, 0);
//            }
//
//            for(vertex u : graph->neighborRange(reverse_ordering[v])){
//                to_vert = ordering[u];
//                if(this->tmp_count[to_vert] == 0){
//                    this->updated.push_back(to_vert);
//                }
//
//                if(to_vert > s && this->tmp_count[to_vert] < K){
//                    this->tmp_count[to_vert] += c;
//                    std::vector<vertex> temp(pvv.second.begin(), pvv.second.end());
//                    temp.push_back(to_vert);
//                    this->node_que[next].push(std::make_pair(to_vert, temp));
//                    this->tmp_dist_count[next][to_vert] += c;
//                }
//            }
//        }
//
//        if (this->node_que[next].empty()){
//            break;
//            }
//        std::swap(curr, next);
//        distance++;
//    }
//    delete[] node_que;
//
//    reset_temp_vars(s, reversed);
//};


inline void DecrementalTopK::pruned_bfs(vertex s, bool reversed){
    this->updated.clear();
    this->updated.push_back(s);
    std::vector<vertex> p = {s};
    this->node_que = new std::queue<std::pair<vertex, std::vector<vertex>>>;
    allocate_label(s, s, p, reversed);
    vertex v, to_vert;
    std::pair<vertex, std::vector<vertex>> pvv;
    pvv.first = s;
    pvv.second = p;
    for(vertex u : graph->neighborRange(reverse_ordering[s])){
        to_vert = ordering[u];
        if(to_vert > s && this->tmp_count[to_vert] < K){
            if(this->tmp_offset[to_vert] == null_distance){
                this->updated.push_back(to_vert);
            }
            std::vector<vertex> temp(pvv.second.begin(), pvv.second.end());
            temp.push_back(to_vert);
            this->node_que->push(std::make_pair(to_vert, temp));
        }
    }
    set_temp_vars(s, reversed);

    dist distance;
     while (!this->node_que->empty()) {
         pvv = this->node_que->front();
         v = pvv.first;
         distance = pvv.second.size()-1;
         this->node_que->pop();
         if(tmp_pruned[v]){
             continue;
         }
         tmp_pruned[v] = prune(v, distance, reversed) >= K;
//         dists.clear();
//
//         query(reverse_ordering[s], reverse_ordering[v], dists);
//         tmp_pruned[v] = dists.size() >= K && dists.rbegin()->size() - 1 <= distance;
         if(tmp_pruned[v]) {
             total_pruning += 1;
             continue;
         }
         if(this->tmp_offset[v] == null_distance){
             this->tmp_offset[v] = distance;
             allocate_label(v, s, pvv.second, reversed);
             this->tmp_count[v] = 1;
         }
         else{
             extend_label(v, s, pvv.second, reversed, 0);
             this->tmp_count[v] += 1;
         }
         for(vertex u : graph->neighborRange(reverse_ordering[v])){
             to_vert = ordering[u];
             if(this->tmp_offset[to_vert] == null_distance){
                 this->updated.push_back(to_vert);
             }

             if(to_vert > s && this->tmp_count[to_vert] < K){
                 std::vector<vertex> temp(pvv.second.begin(), pvv.second.end());
                 temp.push_back(to_vert);
                 this->node_que->push(std::make_pair(to_vert, temp));
             }
         }
     }
    delete node_que;

    reset_temp_vars(s, reversed);
};


void DecrementalTopK::query(vertex s, vertex t, std::vector<dist> & container){
    container.clear();

    s = ordering[s];
    t = ordering[t];
    size_t pos1 = 0;
    size_t pos2 = 0;

    std::vector<dist> count(30, 0);
    const index_t &ids = length_labels[directed][s];
    const index_t &idt = length_labels[0][t];
    vertex W;
    dist d_tmp, c_tmp;
    std::vector<vertex> combination;
    for (;;){
        if (pos1 >= ids.label_offset.size()){
            break;
        }
        if (pos2 >= idt.label_offset.size()){
            break;
        }
        if (ids.label_offset[pos1].first == idt.label_offset[pos2].first){
            if(ids.label_offset[pos1].second >= null_distance || idt.label_offset[pos2].second >= null_distance){
                pos1++;
                pos2++;
                continue;
            }
            W = ids.label_offset[pos1].first;
            for(size_t i = 0; i < ids.p_array[pos1].size(); i++){
                for(size_t j = 0; j < idt.p_array[pos2].size(); j++){
                    for(size_t m = 0; m < loop_labels[W].size(); m++){
                        d_tmp = ids.label_offset[pos1].second + idt.label_offset[pos2].second + i + j + m;
                        c_tmp = loop_labels[W][m].size();
                        if (count.size() <= d_tmp) {
                            count.resize(d_tmp + 1, 0);
                        }
                        count[d_tmp] += (vertex)ids.p_array[pos1][i].size() * idt.p_array[pos2][j].size() * c_tmp;
                    }
                }
            }

            pos1++;
            pos2++;
        }
        else {
            if (ids.label_offset[pos1].first < idt.label_offset[pos2].first){
                pos1++;
            }
            else {
                pos2++;
            }
        }
    }

    for (size_t i = 0; i < count.size(); i++){
        while (container.size() < this->K && count[i]-- > 0){
            container.push_back(i);
        }
    }


    //return container.size() < this->K ? INT_MAX : 0;
}

//inline void DecrementalTopK::verify_sizes(){
//
//    vertex sz = 0;
//    for(vertex i=0;i<graph->numberOfNodes();i++){
//        assert(graph->hasNode(i));
//        sz += loop_labels[i].size(); // loopcount
//
//    }
//    //assert(sz==loop_entries);
//    sz = 0;
//    for (size_t dir = 0; dir < 1 + directed; dir++){
//        for(vertex i=0;i<graph->numberOfNodes();i++){
//            assert(graph->hasNode(i));
//            for(vertex j = 0; j < length_labels[dir][i].label.size(); j++){
//                sz += 1;
//                for(size_t w = 0; w < length_labels[dir][i].label[j].second.size(); w++){
//                    sz += length_labels[dir][i].label[j].second[w].size();
//                }
//            }
//        }
//    }
//    // std::cout<<sz<<" "<<length_entries<<"\n";
//    //assert(sz==length_entries);
//
//}

void DecrementalTopK::update_loops(bool decremental) {
    if(decremental) {
        graph->addEdge(this->x, this->y);
    }
    // todo improve with global variable distance and only one updated vertex set
    this->vertices_to_update.clear();
    std::set<vertex> reset_visited;
    std::queue<vertex> *q = new std::queue<vertex>();

    vertex dequeued_v;
    vertex to_v;
    q->push(ordering[this->x]);
    this->visited_in_update_loops[ordering[this->x]] = 0;
    
    while(!q->empty()){
        dequeued_v = q->front();
        q->pop();
        if(this->visited_in_update_loops[dequeued_v] > K){
            continue;
        }
        if(decremental) {
            for(size_t i = 0; i < this->loop_labels[dequeued_v].size(); i++) {
                for (const auto &arr: this->loop_labels[dequeued_v][i]) {
                    for (vertex w = 0; w < arr.size() - 1; w++) {
                        if ((arr[w] == ordering[this->x] && arr[w + 1] == ordering[this->y])
                            || (arr[w] == ordering[this->y] && arr[w + 1] == ordering[this->x])) {
                            this->vertices_to_update.insert(dequeued_v);
                            goto to_update_from_x;
                        }
                    }
                }
            }
            to_update_from_x:
            {};
        }
        else {
            this->vertices_to_update.insert(dequeued_v);
        }
        for(vertex u : graph->neighborRange(reverse_ordering[dequeued_v])){
            to_v = ordering[u];
            if(this->visited_in_update_loops[to_v] == null_distance && to_v != ordering[this->y]){
                q->push(to_v);
                this->visited_in_update_loops[to_v] = this->visited_in_update_loops[dequeued_v] + 1;
                if(reset_visited.size() != graph->numberOfNodes()){
                    reset_visited.insert(to_v);
                }
            }
        }
    }
    assert(q->empty());
    
    
    
    q->push(ordering[this->y]);

    this->visited_in_update_loops[ordering[this->y]] = 0;

    while(!q->empty()){
        dequeued_v = q->front(); 
        q->pop();
        if(this->visited_in_update_loops[dequeued_v] > K){
            continue;
        }
        if(decremental) {
            for(size_t i = 0; i < this->loop_labels[dequeued_v].size(); i++) {
                for (const auto &arr: this->loop_labels[dequeued_v][i]) {
                    for (vertex w = 0; w < arr.size() - 1; w++) {
                        if ((arr[w] == ordering[this->x] && arr[w + 1] == ordering[this->y])
                            ||
                            (arr[w] == ordering[this->y] && arr[w + 1] == ordering[this->x])) {
                            this->vertices_to_update.insert(dequeued_v);
                            goto to_update_from_y;
                        }
                    }
                }
            }
            to_update_from_y:
            {};
        }
        else{
            this->vertices_to_update.insert(dequeued_v);
        }
        for(vertex u : graph->neighborRange(reverse_ordering[dequeued_v])){

            to_v = ordering[u];
            if(this->visited_in_update_loops[to_v] == null_distance){
                q->push(to_v);
                this->visited_in_update_loops[to_v] = this->visited_in_update_loops[dequeued_v] + 1;
                if(reset_visited.size() != graph->numberOfNodes()){
                    reset_visited.insert(to_v);
                }
            }
        }
    }
    delete q;
    if(decremental) {
        assert(graph->hasEdge(this->x,this->y));
        graph->removeEdge(this->x, this->y);
    }
    else{
        graph->addEdge(this->x, this->y);
    }

    
    aff_cycles = this->vertices_to_update.size();
    
    this->reached_mbfs.clear(); // tracing nodes visited while updating loops

    
    vertex ordered_degree = 0;
    vertex u;
    std::set<vertex>::iterator it;

    for(it=this->vertices_to_update.begin();it!=this->vertices_to_update.end();it++){
    // for(vertex u: to_update){
        u = *it;

        for(const auto& dist_arr: this->loop_labels[u])
            for(const auto& arr: dist_arr)
                total_bits -= 32*arr.size();
        this->loop_labels[u].clear();
        // this->loop_labels[u].shrink_to_fit();
        this->compute_loop_entries(u);
    }
    for(it=reset_visited.begin();it!=reset_visited.end();it++){
        this->visited_in_update_loops[*it] = null_distance;
    }
    reset_visited.clear();
    
    this->visited_in_update_loops[ordering[this->x]] = null_distance;
    this->visited_in_update_loops[ordering[this->y]] = null_distance;
    #ifndef NDEBUG
    this->graph->parallelForNodes([&] (vertex i){
            assert(graph->hasNode(i));
            assert(this->visited_in_update_loops[i]==null_distance);
        });
    #endif
}

void DecrementalTopK::incremental_lengths() {
    this->aff_hubs = 0;
    this->reached_nodes.clear();

    const index_t &idva = length_labels[0][ordering[this->x]];
    const index_t &idvb = length_labels[0][ordering[this->y]];


    this->old_label_a.clear();
    this->old_label_b.clear();
    this->old_paths_a.clear();
    this->old_paths_b.clear();

    vertex w_a, w_b;

    for(auto& elem: idva.label_offset){
        this->old_label_a.push_back(elem);
    }

    for(auto& elem: idvb.label_offset){
        this->old_label_b.push_back(elem);
    }
    this->old_paths_a.resize(idva.p_array.size());
    this->old_paths_b.resize(idvb.p_array.size());
    for(size_t i = 0; i < idva.p_array.size(); i++){
        for(size_t j = 0; j < idva.p_array[i].size(); j++){
            this->old_paths_a[i].push_back(idva.p_array[i][j]);
        }
    }
    for(size_t i = 0; i < idvb.p_array.size(); i++){
        for(size_t j = 0; j < idvb.p_array[i].size(); j++){
            this->old_paths_b[i].push_back(idvb.p_array[i][j]);
        }
    }
    size_t pos_a = 0;
    size_t pos_b = 0;
    while (pos_a != this->old_label_a.size() && pos_b != this->old_label_b.size()){

        assert(this->new_labels.empty());

        w_a = pos_a < this->old_label_a.size() ? this->old_label_a[pos_a].first : graph->numberOfNodes();
        w_b = pos_b < this->old_label_b.size() ? this->old_label_b[pos_b].first : graph->numberOfNodes();

        if(w_a < w_b){
            if(w_a < ordering[this->y]){
                aff_hubs++;
                for(size_t i = 0; i < this->old_paths_a[pos_a].size(); i++){
                    for(size_t j = 0; j < this->old_paths_a[pos_a][i].size(); j++) {
                        std::vector<vertex> temp(this->old_paths_a[pos_a][i][j].begin(),
                                                 this->old_paths_a[pos_a][i][j].end());
                        temp.push_back(ordering[this->y]);
                        incremental_resume_pbfs(w_a, ordering[this->y], temp, false);
                    }
                }
                assert(this->updated.empty());
                reset_temp_vars(w_a, directed);
            }
            pos_a++;

        }
        else if(w_b < w_a){
            if(w_b < ordering[this->x]){
                aff_hubs++;
                for(size_t i = 0; i < this->old_paths_b[pos_b].size(); i++){
                    for(size_t j = 0; j < this->old_paths_b[pos_b][i].size(); j++) {
                        std::vector<vertex> temp(this->old_paths_b[pos_b][i][j].begin(),
                                                 this->old_paths_b[pos_b][i][j].end());
                        temp.push_back(ordering[this->x]);
                        incremental_resume_pbfs(w_b, ordering[this->x], temp, false);
                    }
                }
                assert(this->updated.empty());
                reset_temp_vars(w_b, directed);
            }
            pos_b++;

        }
        else {
            aff_hubs++;

            if(w_a < ordering[this->y]){
                for(size_t i = 0; i < this->old_paths_a[pos_a].size(); i++){
                    for(size_t j = 0; j < this->old_paths_a[pos_a][i].size(); j++) {
                        std::vector<vertex> temp(this->old_paths_a[pos_a][i][j].begin(),
                                                 this->old_paths_a[pos_a][i][j].end());
                        temp.push_back(ordering[this->y]);
                        incremental_resume_pbfs(w_a, ordering[this->y], temp, false);
                    }
                }
                assert(this->updated.empty());
                reset_temp_vars(w_a, directed);
            }
            if(w_b < ordering[this->x]){
                for(size_t i = 0; i < this->old_paths_b[pos_b].size(); i++){
                    for(size_t j = 0; j < this->old_paths_b[pos_b][i].size(); j++) {
                        std::vector<vertex> temp(this->old_paths_b[pos_b][i][j].begin(),
                                                 this->old_paths_b[pos_b][i][j].end());
                        temp.push_back(ordering[this->x]);
                        incremental_resume_pbfs(w_b, ordering[this->x], temp, false);
                    }
                }
                assert(this->updated.empty());
                reset_temp_vars(w_b, directed);
            }
            pos_a++;
            pos_b++;
        }

        while(!this->new_labels.empty()){
            incremental_extend_label_repair(*this->new_labels.front().first.rbegin(),this->new_labels.front().first[0],
                                this->new_labels.front().first, this->new_labels.front().second);
            this->new_labels.pop();
        }
    }
    // std::cout<<"done!\n";
}



void DecrementalTopK::update_lengths() {
    graph->removeEdge(this->x,this->y);
    deleted.emplace_back(x,y);
    this->aff_hubs = 0;
    this->reached_nodes.clear();

    const index_t &idva = length_labels[0][ordering[this->x]];
    const index_t &idvb = length_labels[0][ordering[this->y]];

    
//    this->old_label_a.clear();
//    this->old_label_b.clear();
//
//    vertex w_a, w_b;
//
//    for(auto& elem: idva.label){
//        this->old_label_a.push_back(elem);
//    }
//
//    for(auto& elem: idvb.label){
//        this->old_label_b.push_back(elem);
//    }
    // DISCOVER AFFECTED NODES
    set_temp_update_vars(ordering[this->y],directed);
    this->visited_in_update_lengths_from_a.clear();
    this->visited_in_update_lengths_from_b.clear();
    this->reached_mbfs.clear();

    std::queue<vertex> *q = new std::queue<vertex>();

    vertex dequeued_v;
    vertex to_v;
    q->push(ordering[this->x]);
    this->visited_in_update_lengths_from_a.push_back(ordering[this->x]);

    vertex min_order = min(ordering[this->x], ordering[this->y]);
    this->union_of_reached_nodes.clear();
    std::vector<std::vector<vertex>> query_results;
    std::cout << "Affected nodes from x started..\n";
    while(!q->empty()){
        dequeued_v = q->front();
        bool jump = true;
        q->pop();
        if(this->tmp_pruned[dequeued_v]) continue;

        this->tmp_pruned[dequeued_v] = true;
        this->visited_in_update_lengths_from_a.push_back(dequeued_v);
        //std::cout << "Number of visited vertices " << this->visited_in_update_lengths_from_a.size() << "\n";

        std::vector<std::pair<dist,bool>> reconstructed_k_paths;
        for(size_t i=0; i < this->length_labels[directed][dequeued_v].label_offset.size(); i++){
            if(update_temp_dist_flag[this->length_labels[directed][dequeued_v].label_offset[i].first].empty())
                continue;
            for(auto & p : this->length_labels[directed][dequeued_v].p_array[i]){
                jump = false;
                if(p.empty()) continue;
                for(auto & t : p){
                    for(size_t e = 0; e < t.size() - 1; e++)
                    if ((t[e] == ordering[this->x] && t[e+1] == ordering[this->y])
                        || (t[e] == ordering[this->y] && t[e+1] == ordering[this->x])){
                        jump = true;
                        break;
                    }
                }
                for(const auto & prefix: update_temp_dist_flag[this->length_labels[directed][dequeued_v].label_offset[i].first])
                    reconstructed_k_paths.emplace_back(p[0].size()-1 + prefix.first, jump || prefix.second);
                if(reconstructed_k_paths.size() > K){
                    std::sort(reconstructed_k_paths.begin(), reconstructed_k_paths.end());
                }
                while(reconstructed_k_paths.size() > K)
                    reconstructed_k_paths.pop_back();
            }
        }
        for(const auto & df: reconstructed_k_paths){
            if(df.second){
                goto toresumefromx;
            }
        }
        reconstructed_k_paths.clear();
        continue;
        toresumefromx: {};
        reconstructed_k_paths.clear();
//        bool check = false;
//         query(this->y, reverse_ordering[dequeued_v], query_results);
//
//        for(size_t j = 0; j < query_results.size(); j++){
//            for(size_t w = 0; w < query_results[j].size() - 1; w++){
//                assert(w+1 < query_results[j].size());
//                if((query_results[j][w] == this->x && query_results[j][w+1] == this->y) ||
//                        (query_results[j][w] == this->y && query_results[j][w+1] == this->x)){
//
//                    check = true;
//
//                }
//            }
//        }
//        assert(check);
//        query_results.clear();
        union_of_reached_nodes.insert(dequeued_v);
        for(vertex u : graph->neighborRange(reverse_ordering[dequeued_v])){
            to_v = ordering[u];
            if(this->tmp_pruned[to_v] || to_v == ordering[this->y]){
                continue;
            }
            q->push(to_v);
        }
    }
    assert(q->empty());
    for(const auto& vert: this->visited_in_update_lengths_from_a){
        this->tmp_pruned[vert] = false;
    }
    reset_temp_update_vars(ordering[this->y], directed);
    std::cout << "Affected nodes from x ended..\n";
    std::cout << "Affected nodes from y started..\n";

    set_temp_update_vars(ordering[this->x], directed);
    q->push(ordering[this->y]);
    this->visited_in_update_lengths_from_b.push_back(ordering[this->y]);

    while(!q->empty()){
        dequeued_v = q->front();
        bool jump = true;
        q->pop();
        if(this->tmp_pruned[dequeued_v]) continue;
        this->tmp_pruned[dequeued_v] = true;
        this->visited_in_update_lengths_from_b.push_back(dequeued_v);
        // std::cout << "Number of visited vertices " << this->visited_in_update_lengths_from_b.size() << "\n";
//        query(this->x, reverse_ordering[dequeued_v], query_results);
//        for(size_t j = 0; j < query_results.size(); j++){
//            for(size_t w = 0; w < query_results[j].size()-1; w++){
//                assert(w+1 < query_results[j].size());
//                if((query_results[j][w] == this->x && query_results[j][w+1] == this->y) ||
//                   (query_results[j][w] == this->y && query_results[j][w+1] == this->x)){
//                    union_of_reached_nodes.insert(dequeued_v);
//                    jump = false;
//                }
//            }
//        }
//        query_results.clear();
        std::vector<std::pair<dist,bool>> reconstructed_k_paths;
        for(size_t i = 0; i < this->length_labels[directed][dequeued_v].label_offset.size(); i++){
            if(update_temp_dist_flag[this->length_labels[directed][dequeued_v].label_offset[i].first].empty())
                continue;
            for(auto & p : this->length_labels[directed][dequeued_v].p_array[i]){
                if(p.empty()) continue;
                jump = false;
                for(auto & t : p){
                    for(size_t e = 0; e < t.size() -1; e++) {
                        if ((t[e] == ordering[this->x] && t[e+1] == ordering[this->y])
                            || (t[e] == ordering[this->y] && t[e + 1] == ordering[this->x])) {
                            jump = true;
                            break;
                        }
                    }
                }
                for(const auto & prefix: update_temp_dist_flag[this->length_labels[directed][dequeued_v].label_offset[i].first])
                    reconstructed_k_paths.emplace_back(p[0].size() - 1 + prefix.first, jump || prefix.second); // todo controllare se ce ne sono di stessa lunghezza con flag false
                if(reconstructed_k_paths.size() > K){
                    std::sort(reconstructed_k_paths.begin(), reconstructed_k_paths.end());
                }
                while(reconstructed_k_paths.size() > K)
                    reconstructed_k_paths.pop_back();
            }
        }
        for(const auto & df: reconstructed_k_paths){
            if(df.second){
                goto toresumefromy;
            }
        }
        reconstructed_k_paths.clear();
        continue;
        toresumefromy: {};
        reconstructed_k_paths.clear();
//        bool check = false;
//        query(this->x, reverse_ordering[dequeued_v], query_results);
//
//        for(size_t j = 0; j < query_results.size(); j++){
//            for(size_t w = 0; w < query_results[j].size() - 1; w++){
//                assert(w+1 < query_results[j].size());
//                if((query_results[j][w] == this->x && query_results[j][w+1] == this->y) ||
//                   (query_results[j][w] == this->y && query_results[j][w+1] == this->x)){
//
//                    check = true;
//
//                }
//            }
//        }
//        assert(check);
//        query_results.clear();
        union_of_reached_nodes.insert(dequeued_v);

        for(vertex u : graph->neighborRange(reverse_ordering[dequeued_v])){
            to_v = ordering[u];
            if(this->tmp_pruned[to_v] || to_v == ordering[this->x]){
                continue;
            }
            q->push(to_v);
        }
    }
    assert(q->empty());
    for(const auto& vert: this->visited_in_update_lengths_from_b){
        this->tmp_pruned[vert] = false;
    }
    reset_temp_update_vars(ordering[this->x], directed);
    std::cout << "Affected nodes from y ended..\n";

    // LOOPS UPDATE
    this->update_loops(true);

    // REMOVE OBSOLETE ENTRIES
    std::cout << "Obsolete entries removal started..\n";
    std::set<vertex> to_add;
    for(const vertex& resume_hub: union_of_reached_nodes){
        for(size_t i = 0; i < this->length_labels[0][resume_hub].p_array.size(); i++){
            std::set<vertex> indices_to_remove;
            for(size_t j = 0; j < this->length_labels[0][resume_hub].p_array[i].size(); j++){
                for(size_t w = 0; w < this->length_labels[0][resume_hub].p_array[i][j].size(); w++){
                    for(size_t e = 0; e < this->length_labels[0][resume_hub].p_array[i][j][w].size(); e++)
                    if(
                            (this->length_labels[0][resume_hub].p_array[i][j][w][e] == ordering[this->x]
                            && this->length_labels[0][resume_hub].p_array[i][j][w][e+1] == ordering[this->y]) ||
                            (this->length_labels[0][resume_hub].p_array[i][j][w][e] == ordering[this->y]
                            && this->length_labels[0][resume_hub].p_array[i][j][w][e+1] == ordering[this->x])
                    ){
                        indices_to_remove.insert(w);
                        to_add.insert(this->length_labels[0][resume_hub].label_offset[i].first);
                    }
                }
                for(auto f = indices_to_remove.rbegin(); f != indices_to_remove.rend(); f++){
                    this->length_labels[0][resume_hub].p_array[i][j].erase(this->length_labels[0][resume_hub].p_array[i][j].begin()+*f);
                }
                indices_to_remove.clear();
            }
            size_t first_not_empty = null_vertex;
            for(size_t j = 0; j < this->length_labels[0][resume_hub].p_array[i].size(); j++){
                if(this->length_labels[0][resume_hub].p_array[i][j].size() > 0){
                    first_not_empty = j;
                    break;
                }
            }
            if(first_not_empty == 0) continue;
            if(first_not_empty == null_vertex){ // no more paths, remove label?
                this->length_labels[0][resume_hub].label_offset[i].second = null_distance;
            }
            else{
                for(size_t c = first_not_empty; c < this->length_labels[0][resume_hub].p_array[i].size(); c++){
                    this->length_labels[0][resume_hub].p_array[i][c-first_not_empty] = this->length_labels[0][resume_hub].p_array[i][c];
                    this->length_labels[0][resume_hub].p_array[i][c].clear();
                }
                this->length_labels[0][resume_hub].label_offset[i].second = this->length_labels[0][resume_hub].p_array[i][0][0].size()-1;
            }
        }
    }
    for(auto & v: to_add) union_of_reached_nodes.insert(v);
    std::cout << "Obsolete entries removal ended..\n";
    for(vertex v = 0; v < graph->numberOfNodes(); v++){
        for(const auto & l: this->length_labels[0][v].p_array){
            for(const auto & j: l){
                for(const auto & i: j){
                    for(vertex h = 0; h < i.size() -1; h++){
                        if((i[h] == ordering[this->x] && i[h+1] == ordering[this->y]) ||
                           (i[h] == ordering[this->y] && i[h+1] == ordering[this->x])){
                            throw new std::runtime_error("edge still present in length");
                        }
                    }
                }
            }
        }
        for(const auto & vec: this->loop_labels[v]){
            for(const auto & i: vec){
                for(vertex h= 0; h < i.size()-1; h++){
                    if((i[h] == ordering[this->x] && i[h+1] == ordering[this->y]) ||
                       (i[h] == ordering[this->y] && i[h+1] == ordering[this->x])){
                        throw new std::runtime_error("edge still present in cycle");
                    }
                }
            }
        }

    }
    // RESUME kBFS
    std::cout << "Resumed kBFS started for " << union_of_reached_nodes.size() << " nodes out of " << this->graph->numberOfNodes() << " total vertices..\n";
    for(const vertex& resume_hub: union_of_reached_nodes){
        resume_pbfs(resume_hub, false);
    }
    visited_in_update_lengths_from_a.clear();
    visited_in_update_lengths_from_b.clear();
    std::cout << "Resumed kBFS ended\n";
}

inline size_t DecrementalTopK::prune(vertex v,  dist d, bool rev){
    const index_t &idv = length_labels[rev][v];

    size_t pcount = 0;
    vertex w = 0;
    // cerr << "prune start" << endl;
    for (size_t pos = 0; pos < idv.label_offset.size(); pos++){
        w = idv.label_offset[pos].first;
        if(idv.label_offset[pos].second == null_distance) // todo empty label due to path removal - should be deleted as well?
            continue;
        if (tmp_s_offset[w] == null_distance) continue;

        const vector<dist> &dcs = tmp_s_count[w];

        int l = dcs.size() - 1;
        int c = d - tmp_s_offset[w] - idv.label_offset[pos].second;

        for (int i = 0; i <= c; i++){
                pcount += (int)dcs[std::min(c - i, l)] * idv.p_array[pos][i].size();
        }
        if (pcount >= K) return pcount;
    }
    return pcount;

}
inline void DecrementalTopK::incremental_resume_pbfs(vertex s, vertex res, std::vector<vertex> & p, bool rev) {
    this->updated.clear();

    assert(*p.rbegin() == res);
    this->node_que = new std::queue<std::pair<vertex, std::vector<vertex>>>;
    node_que->push(make_pair(res,p));
    vertex v, to_vert;
    std::pair<vertex, std::vector<vertex>> pvv;
    set_temp_vars(s, rev);

    dist distance;
    while (!this->node_que->empty()) {
        pvv = this->node_que->front();
        v = pvv.first;

        distance = pvv.second.size()-1;
        assert(*pvv.second.rbegin() == v);
        assert(*pvv.second.begin() == s);
        this->node_que->pop();
        if(tmp_pruned[v]){
            continue;
        }
        this->updated.push_back(v);

        tmp_pruned[v] = prune(v, distance, rev) >= K;
//         dists.clear();
//
//         query(reverse_ordering[s], reverse_ordering[v], dists);
//         tmp_pruned[v] = dists.size() >= K && dists.rbegin()->size() - 1 <= distance;
        if(tmp_pruned[v]) {
            total_pruning += 1;
            continue;
        }
        this->new_labels.push(make_pair(pvv.second, rev));
        for(vertex u : graph->neighborRange(reverse_ordering[v])){
            to_vert = ordering[u];

            if(to_vert > s && this->tmp_count[to_vert] < K){
                std::vector<vertex> temp(pvv.second.begin(), pvv.second.end());
                temp.push_back(to_vert);
                this->node_que->push(std::make_pair(to_vert, temp));
            }
        }
    }
    delete node_que;

    reset_temp_vars(s, rev);
}
inline void DecrementalTopK::resume_pbfs(vertex s, bool rev) {
    set_temp_vars(s, rev);
    this->updated.clear();

    this->node_que = new std::queue<std::pair<vertex, std::vector<vertex>>>;
    std::vector<vertex> p = {s};
    this->node_que->push(std::make_pair(s,p));
    std::pair<vertex, std::vector<vertex>> pvv;
    vertex v, to_vert;
    dist distance;
    vertex currently_reached_nodes = 0;
    this->dists.clear();
    this->updated.push_back(s);
    while (!this->node_que->empty()) {
        pvv = this->node_que->front();
        v = pvv.first;
        distance = pvv.second.size()-1;
        this->node_que->pop();
        if(tmp_pruned[v] || tmp_count[v] == K){
            continue;
        }
        currently_reached_nodes += 1;
        // TODO il codice commentato dovrebbe essere sfruttato in qualche modo (naive: fai query in update_lengths anziché guardare solo length_labels)
//        if(union_of_reached_nodes.find(v) == union_of_reached_nodes.end()){
//            goto no_query;
//        }
        // if path in labeling, avoid query and insertion
        for(size_t t = 0; t < this->length_labels[rev][v].label_offset.size(); t++){
            if(this->length_labels[rev][v].label_offset[t].first < s){
                continue;
            }
            if(this->length_labels[rev][v].label_offset[t].first > s){
                break;
            }
            if(this->length_labels[rev][v].label_offset[t].second >= null_distance
                ||
                this->length_labels[rev][v].p_array[t].size() <= pvv.second.size() - 1 - this->length_labels[rev][v].label_offset[t].second){
                break;
            }

            for(const auto& stored_path: this->length_labels[rev][v].p_array[t][pvv.second.size() - 1 - this->length_labels[rev][v].label_offset[t].second]){
                if(stored_path == pvv.second){
                    goto no_query;
                }
            }
        }
//        if(
//                std::find(visited_in_update_lengths_from_b.begin(), visited_in_update_lengths_from_b.end(), v) != visited_in_update_lengths_from_b.end()
//                ||
//                std::find(visited_in_update_lengths_from_a.begin(), visited_in_update_lengths_from_a.end(), v) != visited_in_update_lengths_from_a.end()
//            ) {
            // this->query(reverse_ordering[s], reverse_ordering[v], dists);

            this->tmp_pruned[v] = prune(v,distance,rev) >= K;

            if (this->tmp_pruned[v]) {
                continue;
            }
            extend_label_repair(v, s, pvv.second, rev);
        //}
        no_query: {};
        tmp_count[v] += 1;

        for(vertex u : graph->neighborRange(reverse_ordering[v])){
            to_vert = ordering[u];
            if(this->tmp_count[to_vert] == 0){
                this->updated.push_back(to_vert);
            }

            if(to_vert > s && this->tmp_count[to_vert] < K){
                std::vector<vertex> temp(pvv.second.begin(), pvv.second.end());
                temp.push_back(to_vert);
                this->node_que->push(std::make_pair(to_vert, temp));
            }
        }
    }

    this->reached_nodes.push_back(currently_reached_nodes);
    delete node_que;

    reset_temp_vars(s, rev);
}

inline void DecrementalTopK::set_temp_vars(vertex s, bool rev){

    const index_t &ids = length_labels[directed && !rev][s];

    vertex w;

    for(size_t pos = 0; pos < ids.label_offset.size(); pos++){
        w = ids.label_offset[pos].first;
        if(ids.label_offset[pos].second == null_vertex){
            continue;
        }
        this->tmp_s_offset[w] = ids.label_offset[pos].second;

        this->tmp_v.clear();
        for (size_t i = 0; i < ids.p_array[pos].size(); i++) {
            this->tmp_v.push_back(ids.p_array[pos][i].size());
        }

        this->tmp_s_count[w].resize(tmp_v.size() + loop_labels[w].size() - 1, 0);
        // this->tmp_s_count[w].shrink_to_fit();
        for(size_t i = 0; i < tmp_v.size(); i++){
            for(size_t j = 0; j < loop_labels[w].size(); j++){
                this->tmp_s_count[w][i+j] += this->tmp_v[i] * this->loop_labels[w][j].size();
            }
        }
    }
}

inline void DecrementalTopK::reset_temp_vars(vertex s, bool rev){

    const index_t &ids = length_labels[directed && !rev][s];
    vertex w;
    for(size_t pos = 0; pos < ids.label_offset.size(); pos++){
        w = ids.label_offset[pos].first;
        this->tmp_s_offset[w] = null_distance;
        this->tmp_s_count[w].clear();
        // this->tmp_s_count[w].shrink_to_fit();
    }

    for(std::vector<vertex>::iterator i=this->updated.begin(); i != this->updated.end(); i++){
        w = *i;
        this->tmp_count[w] = 0;
        this->tmp_offset[w] = null_distance;
        this->tmp_pruned[w] = false;
        for(size_t j = 0; j < 2; j++){
            this->tmp_dist_count[j][w] = 0;
        }
    }

    this->updated.clear();
}

inline void DecrementalTopK::set_temp_update_vars(vertex s, bool rev){
    const index_t &ids = length_labels[directed && !rev][s];

    vertex w;

    for(size_t pos = 0; pos < ids.label_offset.size(); pos++){
        std::vector<bool> flags;
        w = ids.label_offset[pos].first;
        if(ids.label_offset[pos].second == null_vertex){
            continue;
        }
        tmp_v.clear();
        for(size_t t = 0; t < ids.p_array[pos].size(); t++) {
            for (const auto &path: ids.p_array[pos][t]) {
                tmp_v.push_back(path.size() - 1);
                bool flag;
                flag = false;
                for (size_t i = 0; i < path.size() - 1; i++) {
                    if ((path[i] == ordering[this->x] && path[i + 1] == ordering[this->y])
                        || (path[i] == ordering[this->y] && path[i + 1] == ordering[this->x])) {
                        flag = true;
                        break;
                    }
                }
                flags.push_back(flag);
            }
        }

        this->update_temp_dist_flag[w].clear();
        // this->tmp_s_count[w].shrink_to_fit();
        for(size_t i = 0; i < tmp_v.size(); i++){
            for(size_t j = 0; j < loop_labels[w].size(); j++){
                bool flag = false;
                for(const auto & l: loop_labels[w][j]){
                    for(size_t e = 0; e < l.size() -1; e++) {
                        if ((l[e] == ordering[this->x] && l[e + 1] == ordering[this->y])
                            || (l[e] == ordering[this->y] &&
                                l[e + 1] == ordering[this->x])) {
                            flag = true;
                            break;
                        }
                    }
                    assert(j == l.size()-1);
                    this->update_temp_dist_flag[w].emplace_back(this->tmp_v[i] + j, flags[i] || flag);
                    if(this->update_temp_dist_flag[w].size() > this->K){
                        std::sort(this->update_temp_dist_flag[w].begin(), this->update_temp_dist_flag[w].end());
                        this->update_temp_dist_flag[w].pop_back();
                    }
                }

            }
        }
    }
}

inline void DecrementalTopK::reset_temp_update_vars(vertex s, bool rev){

    const index_t &ids = length_labels[directed && !rev][s];
    vertex w;
    for(size_t pos = 0; pos < ids.label_offset.size(); pos++){
        w = ids.label_offset[pos].first;
        update_temp_dist_flag[w].clear();
        // this->tmp_s_count[w].shrink_to_fit();
    }

    this->updated.clear();
}

inline void DecrementalTopK::allocate_label(vertex v, vertex start, std::vector<vertex>& path, bool dir){
    index_t &idv = length_labels[dir][v];
    
    idv.label_offset.emplace_back(start,path.size()-1);
    idv.p_array.push_back({{path}});
    total_bits += 32; // new hub
    total_bits += 8; // offset distance from new hub
    total_bits += 32*path.size(); // path contribution
}

double DecrementalTopK::n_reached_nodes(){
    double sum = 0.0; 
    for(auto& element:reached_nodes){
        sum+=(double)element;
    }
    return reached_nodes.size()>0 ? sum / (double)reached_nodes.size() : 0.0;
}

 double DecrementalTopK::n_reached_nodes_mbfs(){
     double sum = 0.0;
     for(auto& element:reached_mbfs){
         sum+=(double)element;
     }
     return reached_mbfs.size()>0 ? sum / (double) reached_mbfs.size() : 0.0;
 }

inline void DecrementalTopK::extend_label(vertex v, vertex start, std::vector<vertex>& path, bool dir, size_t pos){
    index_t &idv = length_labels[dir][v];
    if(pos == 0){
        for(; pos < idv.label_offset.size(); pos++){
            if (idv.label_offset[pos].first == start){
                break;
            }
        }
    }
    if(idv.label_offset[pos].second >= null_distance){
        assert(idv.p_array[pos][0].empty());
        idv.label_offset[pos].second = path.size()-1;
        idv.p_array[pos][0].push_back(path);
    }
    else {
        dist offset = path.size() - 1 - idv.label_offset[pos].second;
        if (offset < idv.p_array[pos].size()) {
            idv.p_array[pos][offset].push_back(path);
        } else {
            idv.p_array[pos].resize(offset + 1, {});
            idv.p_array[pos][offset].push_back(path);
        }
    }
    total_bits += 32 * path.size();
}

inline void DecrementalTopK::incremental_extend_label_repair(vertex v, vertex start, std::vector<vertex>& path, bool dir){
    index_t &idv = length_labels[dir][v];

    size_t last = 0;
    for(; last < idv.label_offset.size(); last++) {
        if (idv.label_offset[last].first == start) {
            break;
        }
        if (idv.label_offset[last].first > start){
            break;
        }
    }
    if (last == idv.label_offset.size()){
        allocate_label(v, start, path, dir);
    }
    else if(idv.label_offset[last].first == start){
        extend_label(v, start, path, dir, last);
    }
    else {
        idv.label_offset.insert(idv.label_offset.begin()+last, std::make_pair(start,path.size()-1));
        idv.p_array.insert(idv.p_array.begin()+last, {{path}});
        total_bits += 32; // new hub
        total_bits += 8; // offset distance from new hub
        total_bits += 32*path.size(); // path contribution
    }

}

inline void DecrementalTopK::extend_label_repair(vertex v, vertex start, std::vector<vertex>& path, bool dir){
    index_t &idv = length_labels[dir][v];

    size_t last = 0;
    for(; last < idv.label_offset.size(); last++) {
        if (idv.label_offset[last].first == start) {
            break;
        }
        if (idv.label_offset[last].first > start){
            break;
        }
    }
    if (last == idv.label_offset.size()){
        allocate_label(v, start, path, dir);
        return;
    }
    else if(idv.label_offset[last].first == start){
        extend_label(v, start, path, dir, last);
    }
    else {
        idv.label_offset.insert(idv.label_offset.begin()+last, std::make_pair(start,path.size()-1));
        idv.p_array.insert(idv.p_array.begin()+last, {{path}});
        total_bits += 32; // new hub
        total_bits += 8; // offset distance from new hub
        total_bits += 32*path.size(); // path contribution
    }

}

// void DecrementalTopK::mod_bfs(vertex s, vertex t, std::vector<dist> &ret) {
//     s = ordering[s];
//     t = ordering[t];
//     std:vector<std::vector<dist> > distanze(graph->numberOfNodes());
//     std::priority_queue<pair<dist, vertex > > que;

//     que.push(make_pair(0, s));
//     vertex v;
//     vertex dist_negative;
//     while (!que.empty()) {
//         v = que.top().second;
//         dist_negative = -que.top().first;
//         que.pop();

//         if(distanze[t].size() >= this->K){
//             break;
//         }
//         if (distanze[v].size() >= this->K){
//             continue;
//         }

//         distanze[v].push_back(dist_negative);
//         for(vertex u : graph->neighborRange(reverse_ordering[v])){

//         // graph->forNeighborsOf(reverse_ordering[v], [&](NetworKit::node u) {
            
//             que.push(make_pair(-(1 + dist_negative), ordering[u]));
//         }
//     }
//     ret = distanze[t];
// }
