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

    if(experiment != 0 && experiment != 1 && experiment != 2 && experiment != 3){
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
        case (3):
            experiment_string = "CML";
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
	std::string timestampstring = "cumulative_"+shortenedName+"_"+k_string+"_"+n_insert_string+"_DECREM_"+order_string+"_"+experiment_string+"_"+tmp_time;

	std::string logFile = timestampstring +".csv";
    

    std::ofstream ofs(logFile);

    double cumulative_time_hc = 0;

    cumulative_time_hc += kpll->loops_time;
    cumulative_time_hc += kpll->lengths_time;
    ofs << "G,V,E,K,Insertions,x,y,UTLoops,UTLengths,UTotalBits,HCCumulative,kBFS,AvgkBFS,MedkBFS\n";
    ofs << graph_location << "," << graph->numberOfNodes() << "," << graph->numberOfEdges() << "," << K << "," << 0 << ","
        << 0 << "," << 0 << ","  << kpll->loops_time << "," << kpll->lengths_time << "," << kpll->total_bits << ","
        <<  cumulative_time_hc << "," << 0 << "," << 0 << "," << 0 << "\n";


    std::vector<double> update_loops_times;
    std::vector<double> update_lengths_times;


    std::vector<size_t> index_total_bits;


    std::vector<pair<vertex, vertex>> added_edges;


    mytimer time_counter;
    double cumulative_time_kbfs = 0;

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
            cumulative_time_hc += *update_loops_times.rbegin();
            std::cout << "done! \nUpdate loops time: " << time_counter.elapsed() << "\n" << std::flush;
            std::cout << "Updating lengths...";
            time_counter.restart();
            kpll->incremental_lengths();
            update_lengths_times.push_back(time_counter.elapsed());
            cumulative_time_hc += *update_lengths_times.rbegin();
            std::cout << "done! \nUpdate lengths time: " << time_counter.elapsed()<<"\n"<<std::flush;
            std::cout << t+1 << "-th update (insertion) done!" << "\n";
        }
        else{
            update_loops_times.push_back(0);
            std::cout << "Updating lengths...";
            time_counter.restart();
            kpll->update_lengths();
            update_lengths_times.push_back(time_counter.elapsed());
            cumulative_time_hc += *update_lengths_times.rbegin();
            std::cout << "done! \nUpdate lengths time: " << time_counter.elapsed()<<"\n"<<std::flush;
            std::cout << t+1 << "-th update (deletion) done!" << "\n";
        }
        index_total_bits.push_back(kpll->total_bits);

        added_edges.push_back(edge_updates[t].second);

    }
    kpll->deallocate_aux();

    assert(added_edges.size()==edge_updates.size());
    std::vector<double> kbfs_time;
    std::vector<double> khl_time;

    ProgressStream query_bar(num_queries);
    vertex u,v;
    query_bar.label() << "Queries Generation and DynKPLL";
    for(uint64_t j=0; j<num_queries; j++){
        u = NetworKit::GraphTools::randomNode(*graph);
        v = NetworKit::GraphTools::randomNode(*graph);
        std::vector<std::vector<vertex>> distances;
        std::vector<std::vector<vertex>> kbfsdistances;
        time_counter.restart();
        kpll->query(u, v, distances);
        khl_time.push_back(time_counter.elapsed());
        cumulative_time_hc += *khl_time.rbegin();
        time_counter.restart();
        kpll->kbfs(u,v,kbfsdistances);
        kbfs_time.push_back(time_counter.elapsed());
        cumulative_time_kbfs += *kbfs_time.rbegin();
        if(distances.size() != kbfsdistances.size()){
            std::cout << distances.size() << " DYNkPLL paths\n";
            for(const auto& path: distances){
                for(const auto& vert: path){
                    std::cout << vert << "-";
                }
                std::cout << "\n";
            }
            std::cout << kbfsdistances.size() << " kbfs paths\n";
            for(const auto& path: kbfsdistances){
                for(const auto& vert: path){
                    std::cout << vert << "-";
                }
                std::cout << "\n";
            }
            throw new std::runtime_error("cardinality problem");
        }

        for(size_t l=0; l < distances.size(); l++){
            if(distances[l].size() != kbfsdistances[l].size()){

                std::cout << "Error bw " << u << "-" << v << "\n";
                std::cout << "Updated labeling distance: " << distances[l].size() << "\n";
                std::cout << "kBFS distance: " << kbfsdistances[l].size() << "\n";
                for(size_t id=0; id < distances.size(); id++){
                    std:: cout << "Up " << distances[id].size() << " | kBFS " << kbfsdistances[id].size() << "\n";
                }
                std::cout << distances[l].size() << " DYNkPLL paths\n";
                for(const auto& path: distances){
                    for(const auto& vert: path){
                        std::cout << vert << "-";
                    }
                    std::cout << "\n";
                }
                std::cout << kbfsdistances.size() << " kBFS paths\n";
                for(const auto& path: kbfsdistances){
                    for(const auto& vert: path){
                        std::cout << vert << "-";
                    }
                    std::cout << "\n";
                }
                throw new std::runtime_error("correctness problem");
            }
        }
        ++query_bar;
    }
    uint64_t final_total_bits = kpll->compute_index_size(), final_aff_hubs = kpll->aff_hubs, final_aff_cycles = kpll->aff_cycles;
    double final_reached = kpll->n_reached_nodes();
    double final_reached_mbfs = kpll->n_reached_nodes_mbfs();
    std::cout << "Cumulative Dyn-KPLL " << cumulative_time_hc << " | Cumulative kBFS " << cumulative_time_kbfs << "\n";
    std::cout << "Mean kBFS " << average(kbfs_time) << " | Median kBFS " << median(kbfs_time) << "\n";
    delete kpll;


    std::cout << "Writing CSV file...";
    ofs << graph_location << ","
        << graph->numberOfNodes() << ","
        << graph->numberOfEdges() << ","
        << K << ","
        << edge_updates.size()+1 << ","
        << "none" << "," << "none" << ","  << "final" << "," << "final" << "," << final_total_bits << ","
        << cumulative_time_hc << ","
        << cumulative_time_kbfs << ","
        << average(kbfs_time) << ","
        << median(kbfs_time) <<"\n";

    std::cout << "done!\n";
    ofs.close();
    delete graph;


    exit(EXIT_SUCCESS);
}