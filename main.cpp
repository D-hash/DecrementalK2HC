// Created by andrea
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <cassert>
#include <map>
#include <algorithm>
#include "decremental_topk.h"
#include <string>
#include <bits/stdc++.h>
#include <boost/program_options.hpp>
#include "networkit/components/ConnectedComponents.hpp"
#include <networkit/graph/GraphTools.hpp>
#include "mytimer.h"
#include "GraphManager.hpp"
#include <networkit/distance/Diameter.hpp>

using namespace std;
double median(std::vector<double>& arr) { //SORTS
    size_t n = arr.size() / 2;
    if (n % 2 == 0) {
        std::nth_element(arr.begin(),arr.begin() + n/2,arr.end());
        std::nth_element(arr.begin(),arr.begin() + (n - 1) / 2,arr.end());
        return (double) (arr[(n-1)/2]+ arr[n/2])/2.0;
    }

    else{
        std::nth_element(arr.begin(),arr.begin() + n / 2,arr.end());
        return (double) arr[n/2];
    }
    assert(false);
}



double average(std::vector<double> & arr) {

    auto const count = static_cast<double>(arr.size());
    double sum = 0;
    for(double value: arr) sum += value;
    return sum / count;
}


int main(int argc, char **argv) {
    srand (time(NULL));
    
    //declare supported options
	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	
	desc.add_options()
	("graph_location,g", po::value<std::string>(), "Input Graph File Location")
	("k_paths,k", po::value<int>(), "Number of Top Paths to Compute")
	("num_insertions,n", po::value<int>(), "Number of Insertions to Be Performed")
    ("num_queries,q", po::value<int>(), "Number of Queries to Be Performed")
	("directed,d",po::value<int>(), "[FALSE(0) TRUE(1)]")
	("ordering,o",po::value<int>(), "Type of Node Ordering [DEGREE(0) APPROX-BETWEENESS(1) k-PATH(2)]")
    ("experiment,e",po::value<int>(), "Type of Experiment [RANDOM(0) SEMI-REALISTIC(1) TEMPORAL(2)]")
    ;

    po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);
    if(vm.empty()){
		std::cout << desc << "\n";
		throw std::runtime_error("empty argument array");
	}
	std::string graph_location;
	
    if(vm.count("graph_location")){
		graph_location = vm["graph_location"].as<std::string>();
    }

    
    int K = -1;

	if (vm.count("k_paths")){
		K = vm["k_paths"].as<int>();
    }

	if (K < 2){
		std::cout << desc << "\n";
		throw std::runtime_error("K must be at least 2");
	}

    int ordering = -1; 

	if (vm.count("ordering")){
		ordering = vm["ordering"].as<int>();
	}
	if(ordering != 0 && ordering != 1 && ordering != 2){
		std::cout << desc << "\n";
		throw std::runtime_error("Wrong ordering selection (0 or 1 or 2)");
	}

    int num_insertions = -1;
	
    if (vm.count("num_insertions")){
		num_insertions = vm["num_insertions"].as<int>();
	}
	if(num_insertions < 2){
		std::cout << desc << "\n";
		throw std::runtime_error("Wrong num_insertions");
	}

    int num_queries = -1;
	if (vm.count("num_queries")){
		num_queries = vm["num_queries"].as<int>();
	}
	if(num_queries < 2){
		std::cout << desc << "\n";
		throw std::runtime_error("Wrong num queries");
	}

	int directed = -1;
	if (vm.count("directed")){
		directed = vm["directed"].as<int>();
	}
	if(directed != 0 && directed !=1){
		std::cout << desc << "\n";
		throw std::runtime_error("Wrong directed selection");
	}
    
    if(directed==0){
        std::cout << "Graph is undirected\n";
    }
    else{ 
        throw std::runtime_error("not yet implemented");

    }

    int experiment = -1;
    if(vm.count("experiment")){
        experiment = vm["experiment"].as<int>();
    }

    if(experiment != 0 && experiment != 1 && experiment != 2){
        std::cout << desc << "\n";
        throw std::runtime_error("Wrong experiment type");
    }

    std::cout << "Reading " << graph_location << " building with K = " << K << " num queries = "<<num_queries<<" ordering "<<ordering<<" experiment " << experiment << "\n";


    // LETTURA
    NetworKit::Graph *graph = new NetworKit::Graph();

    std::vector<std::pair<uint64_t, std::pair<uint32_t , uint32_t>> > edges; // used to sort temporal edges

	if(graph_location.find(".hist") != std::string::npos || graph_location.find(".temporal") != std::string::npos){
		GraphManager::read_hist(graph_location,&graph,&edges);
        std::cout << "Sorting temporal edges...\n";
        sort(edges.begin(), edges.end());
        std::cout << "Sorting done!\n";
	}
	else if(graph_location.find(".nde") != std::string::npos || graph_location.find(".el") != std::string::npos){
		GraphManager::read_nde(graph_location,&graph);

	}
	else{
		throw std::runtime_error("Unsupported graph file format");
	}

	*graph = NetworKit::GraphTools::toUnweighted(*graph);
	*graph = NetworKit::GraphTools::toUndirected(*graph);

    std::vector<std::pair<int,std::pair<vertex,vertex>>> edge_updates;

    if(num_insertions>=graph->numberOfEdges()){
        num_insertions=graph->numberOfEdges();
    }
    edge_id original_n_edges = graph->numberOfEdges();

    if(experiment == 2){
        long long int ni = 0;
        for(size_t i = edges.size()-1; ni < num_insertions; i--) {
            if (find(edge_updates.begin(), edge_updates.end(), make_pair(1,make_pair(edges[i].second.first, edges[i].second.second))) !=
                        edge_updates.end() ||
                find(edge_updates.begin(), edge_updates.end(), make_pair(1,make_pair(edges[i].second.second, edges[i].second.first))) !=
                        edge_updates.end())
                continue;
            if (edges[i].second.first == edges[i].second.second) continue;
            edge_updates.emplace_back(1, make_pair(edges[i].second.first, edges[i].second.second));
            graph->removeEdge(edges[i].second.first, edges[i].second.second);
            ni++;
        }
        std::cout << "Edges after removal " << graph->numberOfEdges() << '\n';
        assert(graph->numberOfEdges()==original_n_edges-edge_updates.size());
        size_t num_removal = num_insertions / 10;
        std::cout << "Finding " << num_removal << " edges for future removal..\n";

        for (size_t j = 0; j < num_removal; j++) {

            vertex a = NetworKit::GraphTools::randomNode(*graph);
            vertex b = NetworKit::GraphTools::randomNeighbor(*graph, a);
            while(find(edge_updates.begin(), edge_updates.end(), make_pair(-1,make_pair(a, b))) !=
                  edge_updates.end() ||
                  find(edge_updates.begin(), edge_updates.end(), make_pair(-1,make_pair(b, a))) !=
                  edge_updates.end() ||
                    find(edge_updates.begin(), edge_updates.end(), make_pair(1,make_pair(a, b))) !=
                    edge_updates.end() ||
                    find(edge_updates.begin(), edge_updates.end(), make_pair(1,make_pair(b, a))) !=
                    edge_updates.end()){
                a = NetworKit::GraphTools::randomNode(*graph);
                b = NetworKit::GraphTools::randomNeighbor(*graph, a);
            }
            edge_updates.insert(edge_updates.begin()+10*j, make_pair(-1, make_pair(a,b)));
        }
    }
    else {
        const NetworKit::Graph &graph_handle = *graph;
        NetworKit::ConnectedComponents *cc = new NetworKit::ConnectedComponents(graph_handle);
        *graph = cc->extractLargestConnectedComponent(graph_handle, true);
        graph->shrinkToFit();
        graph->indexEdges();
        std::cout << "Graph after CC has " << graph->numberOfNodes() << " vertices and " << graph->numberOfEdges()
                  << " edges\n";
        double density = NetworKit::GraphTools::density(*graph);
        NetworKit::Diameter *dm = new NetworKit::Diameter(graph_handle);
        dm->run();


        double diameter = dm->getDiameter().first;
        delete dm;

        std::cout << "Density: " << density << "\n";
        std::cout << "Diameter: " << diameter << "\n";

        std::cout << "Number of insertions " << num_insertions << "\n";

        if (experiment == 1) {
            std::cout << "Finding " << num_insertions << " edges for future removal..\n";

            uint16_t attempts = 0;
            for (size_t i = 0; i < num_insertions; i++) {

                vertex a = NetworKit::GraphTools::randomNode(*graph);
                vertex b = NetworKit::GraphTools::randomNeighbor(*graph, a);

                graph->removeEdge(a, b);
                attempts = 0;
                cc->run();
                while (cc->numberOfComponents() > 1) {
                    attempts += 1;
                    graph->addEdge(a, b);
                    a = NetworKit::GraphTools::randomNode(*graph);
                    b = NetworKit::GraphTools::randomNeighbor(*graph, a);
                    graph->removeEdge(a, b);
                    cc->run();
                    if (attempts > graph->numberOfEdges())
                        throw new std::runtime_error("experiment fails, too many removals, cc cannot be preserved");
                }
                edge_updates.emplace_back(-1, make_pair(a, b));
            }
            for(const auto& edge: edge_updates){
                graph->addEdge(edge.second.first, edge.second.second);
            }
        }
        else {
            uint16_t attempts = 0;
            for (int i = 0; i < num_insertions; i++) {
                vertex a = NetworKit::GraphTools::randomNode(graph_handle);
                vertex b = NetworKit::GraphTools::randomNeighbor(*graph, a);
                graph->removeEdge(a, b);
                attempts = 0;
                cc->run();
                while (cc->numberOfComponents() > 1) {
                    attempts += 1;
                    graph->addEdge(a, b);
                    a = NetworKit::GraphTools::randomNode(*graph);
                    b = NetworKit::GraphTools::randomNeighbor(*graph, a);
                    graph->removeEdge(a, b);
                    cc->run();
                    if (attempts > graph->numberOfEdges())
                        throw new std::runtime_error("experiment fails, too many removals, cc cannot be preserved");
                }
                edge_updates.emplace_back(1, make_pair(a, b));
            }
            assert(edge_updates.size() == num_insertions);
            size_t num_removal = num_insertions / 10;
            std::cout << "Finding " << num_removal << " edges for future removal..\n";

            for (size_t j = 0; j < num_removal; j++) {

                vertex a = NetworKit::GraphTools::randomNode(*graph);
                vertex b = NetworKit::GraphTools::randomNeighbor(*graph, a);
                while(find(edge_updates.begin(), edge_updates.end(), make_pair(-1,make_pair(a, b))) !=
                      edge_updates.end() ||
                      find(edge_updates.begin(), edge_updates.end(), make_pair(-1,make_pair(b, a))) !=
                      edge_updates.end() ||
                      find(edge_updates.begin(), edge_updates.end(), make_pair(1,make_pair(a, b))) !=
                      edge_updates.end() ||
                      find(edge_updates.begin(), edge_updates.end(), make_pair(1,make_pair(b, a))) !=
                      edge_updates.end()){
                    a = NetworKit::GraphTools::randomNode(*graph);
                    b = NetworKit::GraphTools::randomNeighbor(*graph, a);
                }
                edge_updates.insert(edge_updates.begin()+10*j, make_pair(-1, make_pair(a,b)));
            }
        }
    }
    edges.clear();

    //edge_updates = {{1,{0,3}}};
    DecrementalTopK* kpll = new DecrementalTopK(graph, K, directed,ordering, false);

    kpll->build();

    std::cout << "First Labeling Loop time: " << kpll->loops_time << " s | First Labeling Indexing time: " << kpll->lengths_time<< " s\n";
    std::cout << "First Labeling Total MBs: " << kpll->compute_index_size() / 1000000 << "\n";
    std::cout << "Number of Vertices: " << graph->numberOfNodes() << "\n";
    std::cout << "Total prunings: " << kpll->total_pruning << "\n";
    std::string order_string;

    switch(ordering){
        case (0):
            order_string = "DEG";
            break;
        case (1):
            order_string = "BET";
            break;
        case (2):
            order_string = "KPT";
            break;
        default:
            throw new std::runtime_error("problem");

    
    }

    std::string experiment_string;

    switch(experiment){
        case (0):
            experiment_string = "RND";
            break;
        case (1):
            experiment_string = "SRL";
            break;
        case (2):
            experiment_string = "TMP";
            break;
        default:
            throw new std::runtime_error("problem on experiment string");


    }

    //timestring
    std::time_t rawtime;
    std::tm* timeinfo;
    char buffer [80];

    std::time(&rawtime);
    timeinfo = std::localtime(&rawtime);
    std::strftime(buffer,80,"%d-%m-%Y-%H-%M-%S",timeinfo);
    std::string tmp_time = buffer;   
	std::string shortenedName = graph_location.substr(0, 16);
    stringstream string_object_name;
    string_object_name<<K;
    std::string k_string,n_insert_string;
    string_object_name>>k_string;
    string_object_name<<num_insertions;
    string_object_name>>n_insert_string;
	std::string timestampstring = "memory_"+shortenedName+"_"+k_string+"_"+n_insert_string+"_DECREM_"+order_string+"_"+experiment_string+"_"+tmp_time;

	std::string logFile = timestampstring +".csv";
    

    std::ofstream ofs(logFile);


    ofs << "G,V,E,K,Insertions,x,y,UTLoops,UTLengths,UTotalBits,UAvgQT,UMedQT,affhubs,reached,affcycles,reachedMbfs\n";
    ofs << graph_location << "," << graph->numberOfNodes() << "," << graph->numberOfEdges() << "," << K << "," << 0 << ","
        << 0 << "," << 0 << ","  << kpll->loops_time << "," << kpll->lengths_time << "," << kpll->total_bits << ","
        << 0 << "," << 0 << "," << 0 << "," << 0 << "," << 0 << "," << 0 << "\n";
    

    std::vector<double> update_loops_times;
    std::vector<double> update_lengths_times;


    std::vector<size_t> index_total_bits;


    std::vector<vertex> affected_hubs;
    std::vector<vertex> affected_nodes;
    std::vector<vertex> affected_cycles;
    std::vector<double> reached_nodes;
    std::vector<double> reached_nodes_mbfs;
    std::vector<pair<vertex, vertex>> added_edges;


    mytimer time_counter;

    for(size_t t=0;t<edge_updates.size();t++){
        int update_type = edge_updates[t].first;
        kpll->x = edge_updates[t].second.first;
        kpll->y = edge_updates[t].second.second;
        if(update_type == 1) std::cout << "Inserted ";
        else std::cout << "Removed ";
        std::cout << "edge " << kpll->x << " " << kpll->y << "\n" << std::flush;
        if(update_type == 1) {
            std::cout << "Updating loops...";
            time_counter.restart();
            kpll->update_loops(false);
            update_loops_times.push_back(time_counter.elapsed());
            std::cout << "done! \nUpdate loops time: " << time_counter.elapsed() << "\n" << std::flush;
            affected_cycles.push_back(kpll->aff_cycles);
            reached_nodes_mbfs.push_back(kpll->n_reached_nodes_mbfs());
            std::cout << "Updating lengths...";
            time_counter.restart();
            kpll->incremental_lengths();
            update_lengths_times.push_back(time_counter.elapsed());
            std::cout << "done! \nUpdate lengths time: " << time_counter.elapsed()<<"\n"<<std::flush;
            std::cout << t+1 << "-th update (insertion) done!" << "\n";
        }
        else{
            update_loops_times.push_back(0);
            std::cout << "Updating lengths...";
            time_counter.restart();
            kpll->update_lengths();
            update_lengths_times.push_back(time_counter.elapsed());
            affected_cycles.push_back(kpll->up_loops());
            reached_nodes_mbfs.push_back(kpll->n_reached_nodes_mbfs());
            std::cout << "done! \nUpdate lengths time: " << time_counter.elapsed()<<"\n"<<std::flush;
            std::cout << t+1 << "-th update (deletion) done!" << "\n";
        }
        index_total_bits.push_back(kpll->total_bits);

        affected_hubs.push_back(kpll->affected_hubs());
        affected_nodes.push_back(kpll->affected());
        reached_nodes.push_back(kpll->n_reached_nodes());
        added_edges.push_back(edge_updates[t].second);

    }
    kpll->deallocate_aux();

    assert(added_edges.size()==edge_updates.size());
    std::vector<double> sl_time;
    std::vector<double> khl_time;

    std::vector<std::pair<vertex,vertex>> queries;
    std::vector<std::vector<std::vector<vertex>>> dists;
    ProgressStream query_bar(num_queries);
    vertex u,v;
    query_bar.label() << "Queries Generation and DynKPLL";
    for(uint64_t j=0; j<num_queries; j++){
        u = NetworKit::GraphTools::randomNode(*graph);
        v = NetworKit::GraphTools::randomNode(*graph);
        std::vector<std::vector<vertex>> distances;
        time_counter.restart();
        kpll->query(u, v, distances);
        khl_time.push_back(time_counter.elapsed());
        queries.emplace_back(u,v);
        dists.emplace_back(distances);
        ++query_bar;
    }
    uint64_t final_total_bits = kpll->compute_index_size(), final_aff_hubs = kpll->aff_hubs, final_aff_cycles = kpll->aff_cycles;
    double final_reached = kpll->n_reached_nodes();
    double final_reached_mbfs = kpll->n_reached_nodes_mbfs();
    delete kpll;

    DecrementalTopK* scratch_kpll = new DecrementalTopK(graph, K, directed, ordering, true);

    scratch_kpll->build();
    //OK, no updates
    scratch_kpll->deallocate_aux();
    uint64_t scratch_total_bits = scratch_kpll->compute_index_size();
    std::cout << "From Scratch Loop time: " << scratch_kpll->loops_time << " s | From Scratch Indexing time: "<< scratch_kpll->lengths_time<< " s\n";
    std::cout << "From Scratch Total MBs: " << scratch_total_bits / 1000000 << "\n";

    assert(queries.size()==num_queries);
    assert(dists.size()==num_queries);

    ProgressStream second_query_bar(num_queries);


    second_query_bar.label() << "Queries KPLL and Comparison";
    for(uint64_t j=0; j<num_queries; j++){
        u = queries[j].first;
        v = queries[j].second;
        time_counter.restart();

        std::vector<std::vector<vertex>> distances;
        scratch_kpll->query(u, v, distances);
        sl_time.push_back(time_counter.elapsed());

        if(dists[j].size() != distances.size()){
            std::cout << dists[j].size() << " DYNkPLL paths\n";
            for(const auto& path: dists[j]){
                for(const auto& vert: path){
                    std::cout << vert << "-";
                }
                std::cout << "\n";
            }
            std::cout << distances.size() << " Scratch-kPLL paths\n";
            for(const auto& path: distances){
                for(const auto& vert: path){
                    std::cout << vert << "-";
                }
                std::cout << "\n";
            }
            throw new std::runtime_error("cardinality problem");
        }
        for(size_t l=0; l < dists[j].size(); l++){
            if(dists[j][l].size() != distances[l].size()){

                std::cout << "Error bw " << u << "-" << v << "\n";
                std::cout << "Updated labeling distance: " << dists[j][l].size() << "\n";
                std::cout << "Scratch labeling distance: " << distances[l].size() << "\n";
                for(size_t id=0; id < dists[j].size(); id++){
                    std:: cout << "Up " << dists[j][id].size() << " | Scratch " << distances[id].size() << "\n";
                }
                std::cout << dists[j].size() << " DYNkPLL paths\n";
                for(const auto& path: dists[j]){
                    for(const auto& vert: path){
                        std::cout << vert << "-";
                    }
                    std::cout << "\n";
                }
                std::cout << distances.size() << " Scratch-kPLL paths\n";
                for(const auto& path: distances){
                    for(const auto& vert: path){
                        std::cout << vert << "-";
                    }
                    std::cout << "\n";
                }
                throw new std::runtime_error("correctness problem");
            }
        }
        ++second_query_bar;
    }

    std::cout << "Writing CSV file...";
    for(size_t j = 0; j < num_insertions; j++) {
        ofs << graph_location << "," << graph->numberOfNodes() << "," << graph->numberOfEdges() << "," << K << "," << j << ","
            << added_edges[j].first << "," << added_edges[j].second << ","  << update_loops_times[j] << "," << update_lengths_times[j] << "," << index_total_bits[j]  << ","
            << 0 << "," << 0 << "," << affected_hubs[j] << "," << reached_nodes[j] << "," << affected_cycles[j] << "," << reached_nodes_mbfs[j] <<"\n";
    }
    ofs << graph_location << ","
        << graph->numberOfNodes() << ","
        << graph->numberOfEdges() << ","
        << K << ","
        << edge_updates.size()+1 << ","
        << "none" << "," << "none" << ","  << "final" << "," << "final" << "," << final_total_bits << ","
        << average(khl_time) << ","
        << median(khl_time) << ","
        << final_aff_hubs << ","
        << final_reached << ","
        << final_aff_cycles << ","
        << final_reached_mbfs <<"\n";

    ofs << graph_location << ","
        << graph->numberOfNodes() << ","
        << graph->numberOfEdges() << ","
        << K << ","
        << num_insertions << ","
        << "scratch" << "," << "scratch" << ","  << scratch_kpll->loops_time << "," << scratch_kpll->lengths_time << "," << scratch_total_bits << ","
        << average(sl_time) << ","
        << median(sl_time) << ",scratch,scratch,scratch,scratch\n";
    std::cout << "done!\n";
    ofs.close();
    delete scratch_kpll;
    delete graph;


    exit(EXIT_SUCCESS);
}
