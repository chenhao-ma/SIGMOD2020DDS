//
// Created by MA Chenhao on 7/5/2019.
//


#include "FlowNetwork.h"

void FlowNetwork::enqueue (int v) {
    if (!active[v] && excess[v] > 0 && dist[v] < n) {
        active[v] = true;
        B[dist[v]].push_back(v);
        b = std::max(b, dist[v]);
    }
}

void FlowNetwork::push (EdgeFN &e) {
    double amt = std::min(excess[e.from], e.cap - e.flow);
    if (dist[e.from] == dist[e.to] + 1 && amt > 0) {
        e.flow += amt;
        adj[e.to][e.index].flow -= amt;
        excess[e.to] += amt;
        excess[e.from] -= amt;
        enqueue(e.to);
    }
}

void FlowNetwork::gap (int k) {
    for (int v = 0; v < n; v++) {
        if (dist[v] >= k) {
            count[dist[v]]--;
            dist[v] = std::max(dist[v], n);
            count[dist[v]]++;
            enqueue(v);
        }
    }
}

void FlowNetwork::relabel(int v) {
    count[dist[v]]--;
    dist[v] = n;
    for (auto e : adj[v]) if (e.cap - e.flow > 0) {
        dist[v] = std::min(dist[v], dist[e.to] + 1);
    }
    count[dist[v]]++;
    enqueue(v);
}

void FlowNetwork::discharge(int v) {
    for (auto &e : adj[v]) {
        if (excess[v] > 0) {
            push(e);
        } else {
            break;
        }
    }

    if (excess[v] > 0) {
        if (count[dist[v]] == 1) {
            gap(dist[v]);
        } else {
            relabel(v);
        }
    }
}

double FlowNetwork::get_maxflow(int s, int t) {
    dist =  std::vector <int>(n, 0);
    excess = std::vector <double>(n, 0);
    count = std::vector <int>(n + 1, 0);
    active =  std::vector <bool>(n, false);
    B = std::vector <std::vector <int>>(n);
    b = 0;
//    dist[s] = n;

    for (auto &e: adj[s]) {
        excess[s] += e.cap;
    }

    count[0] = n;
    enqueue(s);
    active[t] = true;

    while (b >= 0) {
        if (!B[b].empty()) {
            int v = B[b].back();
            B[b].pop_back();
            active[v] = false;
            discharge(v);
        } else {
            b--;
        }
    }

    return excess[t];
}

double FlowNetwork::get_mincut(int s, int t, std::vector<int> &S, std::vector<int> &T, double &edge_num,
                               std::vector<double> &weights) {
    auto ret = get_maxflow(s, t);
    S.clear();
    T.clear();
    edge_num = 0;
    for (int v = 0; v < n; v++) {
        if (dist[v] >= n) {
            if (v == s)
                continue;
            if (v > 0 && v <= nums[0]) {
                S.push_back(ori_id[v]);
                for (auto edge : adj[v]) {
                    if (edge.to > nums[0] && edge.to <= nums[1] && dist[edge.to] >= n) {
                        edge_num += weights[ori_id[edge.to]];
                    }
                }
            }
            else
                T.push_back(ori_id[v]);
        }
    }
    return ret;
}

FlowNetwork::FlowNetwork(int max_n, std::vector<std::pair<int, int>> &edges, double g, double sqrt_a,
                         std::vector<double> &weights) {
    std::vector<bool> flag[2];
    std::vector<int> mapping[2];
    for (int cur = 0; cur < 2; cur++) {
        flag[cur].resize(static_cast<unsigned long>(max_n), false);
        mapping[cur].resize(static_cast<unsigned long>(max_n));
    }
    for (auto &edge : edges) {
        flag[0][edge.first] = true;
        flag[1][edge.second] = true;
    }
    n = 0;
    ori_id.emplace_back();
    for (int cur = 0; cur < 2; cur++) {
        for (int i = 0; i < max_n; i++) {
            if (flag[cur][i]) {
                mapping[cur][i] = ++n;
                ori_id.push_back(i);
            }
        }
        nums.push_back(n);
    }
    ++n;
    ++n;
    std::vector<int> deg(n, 0);
    adj.resize(n);
    double m = 0;
    for (auto &edge : edges) {
        m += weights[edge.second];
    }
    for (auto &edge : edges) {
        int u = mapping[1][edge.second];
        int v = mapping[0][edge.first];
        add_edge(u, v, 2 * weights[edge.second]);
        deg[u] += weights[edge.second];
    }
    for (int i = 1; i < n - 1; i++) {
        add_edge(0, i, m);
        if (i <= nums[0]) {
            add_edge(i, n - 1, m + g / sqrt_a);
        } else {
            add_edge(i, n - 1, m + g * sqrt_a - 2 * deg[i]);
        }
    }
}

void FlowNetwork::add_edge(int from, int to, double cap) {
    adj[from].push_back(EdgeFN(from, to, cap, 0, adj[to].size()));
    if (from == to) {
        adj[from].back().index++;
    }
    adj[to].push_back(EdgeFN(to, from, 0, 0, adj[from].size() - 1));
}
