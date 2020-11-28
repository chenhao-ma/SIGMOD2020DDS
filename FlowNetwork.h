//
// Created by MA Chenhao on 7/5/2019.
//

#ifndef DIRECTEDDENSESTSUBGRAPH_FLOWNETWORK_H
#define DIRECTEDDENSESTSUBGRAPH_FLOWNETWORK_H

#include <vector>
#include <queue>
#include "EdgeFN.h"
#include <cmath>

class FlowNetwork {
public:
    int n;
    std::vector <std::vector <EdgeFN>> adj;
    std::vector <double> excess;
    std::vector <int> dist, count;
    std::vector <bool> active;
    std::vector <std::vector <int>> B;
    std::vector<int> nums;
    std::vector<int> ori_id;
    int b;
    std::queue <int> Q;

    FlowNetwork(int max_n, std::vector<std::pair<int, int>> &edges, double g, double sqrt_a,
                    std::vector<double> &weights);

    void add_edge(int from, int to, double cap);

    void enqueue (int v);

    void push (EdgeFN &e);

    void gap (int k);

    void relabel (int v);

    void discharge (int v);

    double get_maxflow (int s, int t);

    double get_mincut(int s, int t, std::vector<int> &S, std::vector<int> &T, double &edge_num,
                          std::vector<double> &weights);

};


#endif //DIRECTEDDENSESTSUBGRAPH_FLOWNETWORK_H
