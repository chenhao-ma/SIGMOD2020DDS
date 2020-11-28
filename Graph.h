//
// Created by MA Chenhao on 27/4/2019.
//

#ifndef DIRECTEDDENSESTSUBGRAPH_GRAPH_H
#define DIRECTEDDENSESTSUBGRAPH_GRAPH_H

#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cmath>
#include <set>
#include <map>
#include <utility>
#include <boost/heap/fibonacci_heap.hpp>
#include <ctime>
#include "FlowNetwork.h"

class Graph {
public:
    explicit Graph(FILE* file, bool weighted);

    explicit Graph(int node_num, const std::vector<std::pair<int, int> > &edges);

    Graph(const Graph& g);

    ~Graph();

    double approx_core_ds();

    double approx_ds();

    double approx_para_ds(double delta, double epsilon);

    double approx_icalp_ds();

    double approx_ficalp_ds();

    double exact_core_ds();

    double exact_ds();

    double exact_ad_core_ds();

    void output_ds();


private:
    int n;
    int m;
    std::vector<std::vector<int> > *adj;
    std::vector<int> *deg;
    std::vector<int> *bin;
    std::vector<int> *vert;
    std::vector<int> *pos;
    std::vector<int> max_deg;
    std::vector<int> flag;
    std::vector<int> mapping;
    std::vector<int> S;
    std::vector<int> T;
    std::vector<double> weights;
    double min_weight;
    double density;
    double r_max;
    bool weighted;
    bool size_reported;

    int skyline_core_num(int cur, int alpha, int beta, bool reduced = false);

    double gen_alpha_beta_core(int alpha, int beta);

    void add_edge(int u, int v);

    void gen_supporting_structures();

    void sort_adj_list();

    void init_adj_deg();

    void dec_deg(int cur, int t);

    std::pair<int, int> get_max_core_num_pair();

    int get_delta();

    void exact_ds_a(double l, double r, double bias, double sqrt_left, double sqrt_a, double sqrt_right,
                        std::pair<int, int> &ratio_b);

    void divide_conquer(double left, double right);

    void peeling(double c, double epsilon);

    void peel(int k, int cur);
};



struct compare_pair
{
    bool operator()(const std::pair<int, int>& n1, const std::pair<int, int>& n2) const
    {
        return n1.first > n2.first;
    }
};


#endif //DIRECTEDDENSESTSUBGRAPH_GRAPH_H
