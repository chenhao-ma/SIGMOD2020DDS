//
// Created by MA Chenhao on 27/4/2019.
//

//#include <tkDecls.h>
#include "Graph.h"


Graph::Graph(FILE* file, bool weighted) {
    clock_t begin = clock();
    int edge_num;
    fscanf(file, "%d%d", &n, &edge_num);
    size_reported = false;
    init_adj_deg();

    flag.resize(n);
    mapping.resize(n);
    weights.resize(n);

    m = 0;
    min_weight = 100;
    r_max = 0;
    this->weighted = weighted;
    for (int i = 0; i < edge_num; i++) {
        int u, v;
        fscanf(file, "%d%d", &u, &v);
        if (weighted){
            fscanf(file, "%lf", &weights[v]);
        } else {
            weights[v] = 1;
        }
        r_max += weights[v];
        min_weight = std::min(min_weight, weights[v]);
        add_edge(u, v);
    }

    sort_adj_list();

    gen_supporting_structures();

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("graph construction: %.4f\n", elapsed_secs);
}

Graph::Graph(const int node_num, const std::vector<std::pair<int, int> > &edges) {
    n = node_num;
    init_adj_deg();

    for (auto& edge : edges) {
        add_edge(edge.first, edge.second);
    }

    gen_supporting_structures();
//    printf("ok reduced\n");
}

/**
 * Copy constructor of Graph
 * @param g the graph to be copied
 */
Graph::Graph(const Graph &g) {
    n = g.n;
    m = g.m;
    max_deg = g.max_deg;
    adj = new std::vector<std::vector<int> >[2];
    deg = new std::vector<int>[2];
    bin = new std::vector<int>[2];
    vert = new std::vector<int>[2];
    pos = new std::vector<int>[2];
    for (int i = 0; i < 2; i++) {
        adj[i] = g.adj[i];
        deg[i] = g.deg[i];
        bin[i] = g.bin[i];
        vert[i] = g.vert[i];
        pos[i] = g.pos[i];
    }
    flag = g.flag;
    mapping = g.mapping;
}

void Graph::init_adj_deg() {
    max_deg.resize(2, 0);
    adj = new std::vector<std::vector<int> >[2]; // 0 out edges, 1 in edges
    deg = new std::vector<int>[2];
    m = 0;
    for (int i = 0; i < 2; i++) {
        adj[i].resize(static_cast<unsigned long>(n));
        deg[i].resize(static_cast<unsigned long>(n), 0);
    }
}

void Graph::gen_supporting_structures() {
    bin = new std::vector<int>[2];
    vert = new std::vector<int>[2];
    pos = new std::vector<int>[2];
    for (int i = 0; i < 2; i++){
        bin[i].resize(static_cast<unsigned long>(max_deg[i] + 1), 0);
        vert[i].resize(static_cast<unsigned long>(n), 0);
        pos[i].resize(static_cast<unsigned long>(n), 0);
        for (int v = 0; v < n; v++) {
            ++bin[i][deg[i][v]];
        }
        int start = 0;
        for (int d = 0; d <= max_deg[i]; d++) {
            int num = bin[i][d];
            bin[i][d] = start;
            start += num;
        }
        for (int v = 0; v < n; v++) {
            pos[i][v] = bin[i][deg[i][v]]++;
            vert[i][pos[i][v]] = v;
        }
        for (int d = max_deg[i]; d > 0; d--) {
            bin[i][d] = bin[i][d - 1];
        }
        bin[i][0] = 0;
    }
}

Graph::~Graph() {
    for (int i = 0; i < 2; i++) {
//        for (int v = 0; v < n; i++) {
//            adj[i][v].clear();
//        }
        adj[i].clear();
        deg[i].clear();
        bin[i].clear();
        vert[i].clear();
        pos[i].clear();
    }
    delete[] adj;
    delete[] deg;
    delete[] bin;
    delete[] vert;
    delete[] pos;
    max_deg.clear();
}

void Graph::sort_adj_list() {
    for (int i = 0; i < 2; i++) {
        for (int v = 0; v < n; v++) {
            std::sort(adj[i][v].begin(), adj[i][v].end(),
                      [&](const int& a, const int& b) -> bool
                      {
                          return deg[1 - i][a] > deg[1 - i][b];
                      });
        }
    }
}

inline void Graph::add_edge(int u, int v) {
//    printf("%d %d\n", u, v);
    adj[0][u].push_back(v);
    deg[0][u] += 1;
    adj[1][v].push_back(u);
    deg[1][v] += 1;
    max_deg[0] = std::max(max_deg[0], deg[0][u]);
    max_deg[1] = std::max(max_deg[1], deg[1][v]);
    ++m;
}

double Graph::approx_core_ds() {
    std::pair<int, int> core_nums = get_max_core_num_pair();
    density = gen_alpha_beta_core(core_nums.first, core_nums.second);
    return density;
}

/**
 * Find the maximum core number pair (alpha, beta) when alpha is fixed
 * @param cur S side (0) or T side (1)
 * @param alpha the fixed alpha value
 * @param beta the lower bound for beta to beat the current maximum core number pair
 * @param reduced whether the graph is already reduced
 * @return the maximum value for beta
 */

int Graph::skyline_core_num(int cur, int alpha, int beta, bool reduced) {
//    printf("cur %d alpha %d beta %d reduced %d\n", cur, alpha, beta, reduced);
    if (beta > std::upper_bound(bin[1 - cur].begin(), bin[1 - cur].end(), n - alpha) - bin[1 - cur].begin() - 1)
        return 0;

    if (reduced) {
        std::vector<bool> flag(static_cast<unsigned long>(n), false);
        int ret = 0;
        for (int i = 0; i < n; i++) {
            int u = vert[1 - cur][i];
            ret = std::max(ret, deg[1 - cur][u]);
            for (auto v : adj[1 - cur][u]) {
                if (flag[v]) continue;
                if (--deg[cur][v] < alpha) {
                    flag[v] = true;
                    for (auto t : adj[cur][v]) {
                        if (t != u && deg[1 - cur][t] > deg[1 - cur][u]) {
                            dec_deg(1 - cur, t);
                        }
                    }
                }
            }
        }
        return ret;

    } else {
//        clock_t begin = clock();
        std::vector<std::pair<int, int> > edges;
        for (int i = bin[cur][alpha]; i < n; i++) {
            int u = vert[cur][i];
            int v = adj[cur][u][alpha - 1];
            if (deg[1 - cur][v] < beta) continue;
            flag[u] = alpha;
            for (auto t :adj[cur][u]) {
                if (deg[1 - cur][t] < beta)
                    break;
                flag[t] = alpha;
                edges.emplace_back(u, t);
            }
        }
        int nodes_num = n;
        if (edges.size() < m / 5) {
            nodes_num = 0;
            for (int i = 0; i < n; i++) {
                if (flag[i] == alpha)
                    mapping[i] = nodes_num++;
            }
            for (auto &edge : edges) {
                edge.first = mapping[edge.first];
                edge.second = mapping[edge.second];
            }
        }
//        printf("new n %d new m %lu\n", nodes_num, edges.size());
        Graph reduced_g(nodes_num, edges);
//        clock_t end = clock();
//        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//        printf("graph construction: %.4f\n", elapsed_secs);

        auto ret = reduced_g.skyline_core_num(0, alpha, beta, true);
//        clock_t end1 = clock();
//        elapsed_secs = double(end1 - end) / CLOCKS_PER_SEC;
//        printf("calculation time: %.4f\n", elapsed_secs);
        return ret;
    }
}

double Graph::gen_alpha_beta_core(int alpha, int beta) {
    int i = 0, j = 0;
    while (deg[0][vert[0][i]] < alpha || deg[1][vert[1][j]] < beta) {
        for (; deg[0][vert[0][i]] < alpha; i++) {
            int u = vert[0][i];
            for (auto v : adj[0][u]) {
                if (deg[1][v] >= beta) {
                    dec_deg(1, v);
                }
            }
        }

        for (; deg[1][vert[1][j]] < beta; j++) {
            int u = vert[1][j];
            for (auto v : adj[1][u]) {
                if (deg[0][v] >= alpha) {
                    dec_deg(0, v);
                }
            }
        }
    }

    int edge_num = 0;
    int s = n - i, t = n - j;
    for (; i < n; i++) {
        int u = vert[0][i];
        for (auto v : adj[0][u]) {
            if (deg[1][v] >= beta) {
                ++edge_num;
            }
        }
    }
    printf("s %d, t %d, e %d\n", s, t, edge_num);
    return edge_num / sqrt(s * t);
}

inline void Graph::dec_deg(int cur, int t) {
    int dt = deg[cur][t], pt = pos[cur][t];
    int pw = bin[cur][dt], w = vert[cur][pw];
    if (t != w) {
        pos[cur][t] = pw; vert[cur][pt] = w;
        pos[cur][w] = pt; vert[cur][pw] = t;
    }
    ++bin[cur][dt];
    --deg[cur][t];
}

using Heap = boost::heap::fibonacci_heap<std::pair<int, int>, boost::heap::compare<compare_pair> >;
double Graph::approx_ds() {
    int s_num = n - bin[0][1];
    int t_num = n - bin[1][1];
    std::vector<Heap::handle_type> handles[2];
    std::vector<bool> flags[2];
    for (int i = 0; i < 2; i++) {
        handles[i].resize(n);
        flags[i].resize(n);
    }
    double density = 0, ratio = 0;
    int s_size = 0, t_size = 0;
    for (int s_test = 1; s_test < s_num; ++s_test) {
        for (int t_test = 1; t_test < t_num; ++t_test) {
            double c = s_test * 1.0 / t_test;
            Heap heaps[2];
            for (int i = 0; i < 2; i++) {
                for (int j = bin[i][1]; j < n; j++) {
                    int u = vert[i][j];
                    handles[i][u] = heaps[i].push(std::make_pair(deg[i][u], u));
                    flags[i][u] = false;
                }
            }
            int ss = s_num, tt = t_num;
            double e = m;
            auto remove_node = [&](int cur) {
                int u = heaps[cur].top().second;
                heaps[cur].pop();
                flags[cur][u] = true;
                for (auto v : adj[cur][u]) {
                    if (!flags[1 - cur][v]) {
                        (*handles[1 - cur][v]).first--;
                        heaps[1 - cur].increase(handles[1 - cur][v]);
                        e -= 1;
                    }
                }
            };
            while (ss > 0 && tt > 0) {
                if (e / sqrt(ss * tt) > density) {
                    density = e / sqrt(ss * tt);
                    ratio = c;
                    s_size = ss;
                    t_size = tt;
                }
                if (heaps[0].top().first <= heaps[1].top().first / c) {
                    remove_node(0);
                    --ss;
                } else {
                    remove_node(1);
                    --tt;
                }
            }
        }
    }
    printf("s %d, t %d, c %.2f\n", s_size, t_size, ratio);
    return density;
}

double Graph::exact_core_ds() {
    double l = 0;
    if (!weighted) {
        std::pair<int, int> max_core_nums = get_max_core_num_pair();
        l = sqrt((double) max_core_nums.first * max_core_nums.second);
        r_max = l * 2;
    }
    double a_lower = l / (2 * max_deg[0]); a_lower *= a_lower;
    double a_upper = 2 * max_deg[1] / l; a_upper *= a_upper;
    auto cmp_larger = [](const std::pair <int, int> &lhs, const std::pair <int, int> &rhs) {
        return (long long) lhs.first * rhs.second > (long long) rhs.first * lhs.second;
    };
    auto cmp_smaller = [](const std::pair <int, int> &lhs, const std::pair <int, int> &rhs) {
        return (long long) lhs.first * rhs.second < (long long) rhs.first * lhs.second;
    };
    std::priority_queue <std::pair <int, int>, std::vector <std::pair <int, int>>, decltype(cmp_larger) > q_lower(cmp_larger);

    for (int t = 1; t <= n - bin[1][1]; t++) {
        auto s = static_cast<int>(ceil(t * a_lower));
        if (s > n - bin[0][1] || s > t * a_upper) continue;
        q_lower.push(std::make_pair(s, t));
    }
    l = 0;

    density = 0;

//    bool flag = true;

    while (!q_lower.empty()) {
        std::pair<int, int> ratio;
        ratio = q_lower.top();
        q_lower.pop();
        if ((double) ratio.first > ratio.second * a_upper)
            break;
        printf("a %d/%d %.2f\n", ratio.first, ratio.second, ratio.first * 1.0 / ratio.second);
        if (ratio.first > n - bin[0][1])
            continue;
        if ((double) ratio.first < ratio.second * a_lower) {
            q_lower.push(std::make_pair(static_cast<int>(ceil(a_lower * ratio.second)), ratio.second));
            continue;
        }
        q_lower.push(std::make_pair(ratio.first + 1, ratio.second));
        if ((long long) ratio.first * q_lower.top().second == (long long) ratio.second * q_lower.top().first)
            continue;

        double sqrt_a = sqrt((double) ratio.first / ratio.second);
        auto alpha = std::max(static_cast<int>(ceil(l / 2 / sqrt_a)), 1);
        auto beta = std::max(static_cast<int>(ceil(sqrt_a * l / 2)), 1);

        double r = r_max;
//        l = 0;
        double bias = min_weight / sqrt((double) (n - bin[0][alpha]) * (n - bin[1][beta]) - std::min(n - bin[0][alpha], n - bin[1][beta]))
                      - min_weight / sqrt((double) (n - bin[0][alpha]) * (n - bin[1][beta]));
        std::pair<int, int> ratio_b;

        exact_ds_a(l, r, bias, sqrt_a, sqrt_a, sqrt_a, ratio_b);

        printf("g %.4f density %.4f b %d/%d \n", l, density, ratio_b.first, ratio_b.second);

        //update a_lower and a_upper
        l = density;
        a_lower = l / (2 * max_deg[0]); a_lower *= a_lower;
        a_upper = 2 * max_deg[1] / l; a_upper *= a_upper;
    }

    return density;
}

double Graph::exact_ad_core_ds() {
    if (!weighted) {
        std::pair<int, int> max_core_nums = get_max_core_num_pair();
        r_max = 2 * sqrt((double) max_core_nums.first * max_core_nums.second);
    }
    density = 0;
    double bias = 1 / sqrt((double) (n - bin[0][1]) * (n - bin[1][1]) - std::min(n - bin[0][1], n - bin[1][1]))
                  - 1 / sqrt((double) (n - bin[0][1]) * (n - bin[1][1]));
    std::pair<int, int> a_lower_pair = std::make_pair(1, n - bin[1][1]);
    std::pair<int, int> a_upper_pair = std::make_pair(n - bin[0][1], 1);
    double left_cover = (double) a_lower_pair.first / a_lower_pair.second;
    double right_cover = (double) a_upper_pair.first / a_upper_pair.second;
    divide_conquer(left_cover, right_cover);

    return density;
}

double Graph::exact_ds() {

    auto cmp = [](std::pair <int, int> &lhs, std::pair <int, int> &rhs) {
        return (double) lhs.first / lhs.second > (double) rhs.first / rhs.second;
    };
    std::priority_queue <std::pair <int, int>, std::vector <std::pair <int, int>>, decltype(cmp) > q(cmp);

    for (int t = 1; t <= n; t++) {
        q.push(std::make_pair(1, t));
    }

    double density = 0;
    std::vector<int> S;
    std::vector<int> T;

    while (!q.empty()) {
        auto ratio = q.top();
        q.pop();
        if (ratio.first > n)
            continue;
        q.push(std::make_pair(ratio.first + 1, ratio.second));
        if ((long long) ratio.first * q.top().second == (long long) ratio.second * q.top().first)
            continue;
        double sqrt_a = sqrt((double) ratio.first / ratio.second);

        double l = 0;
        double r = m;

        double bias = min_weight / sqrt((double) (n - bin[0][1]) * (n - bin[1][1]) - std::min(n - bin[0][1], n - bin[1][1]))
                      - min_weight / sqrt((double) (n - bin[0][1]) * (n - bin[1][1]));

        while (r - l > bias) {
            double mid = (r + l) / 2;

            std::vector<std::pair<int, int> > edges;
            int cur = 0;
            for (int i = bin[cur][1]; i < n; i++) {
                int u = vert[cur][i];
                for (auto t :adj[cur][u]) {
                    edges.emplace_back(u, t);
                }
            }

            //construct flow network and compute min cut
            std::pair<int, int> ratio_tmp = std::make_pair(0, 0);
            FlowNetwork fn = FlowNetwork(n, edges, mid, sqrt_a, weights);
            std::vector<int> tmp_S, tmp_T;
            double edge_num;
            fn.get_mincut(0, fn.n - 1, tmp_S, tmp_T, edge_num, weights);

            if (!tmp_S.empty() && !tmp_T.empty()) {
                l = mid;
                if (density < edge_num / sqrt(tmp_S.size() * tmp_T.size())) {
                    density = edge_num / sqrt(tmp_S.size() * tmp_T.size());
                    S = tmp_S;
                    T = tmp_T;
                }
            } else {
                r = mid;
            }
        }
        printf("a %d/%d g %.4f density %.4f b %lu/%lu\n", ratio.first, ratio.second, l, density, S.size(), T.size());
    }
    return density;
}

std::pair<int, int> Graph::get_max_core_num_pair() {
    long long max_prod = 0;
    int core_nums[2];
    int cur = (max_deg[0] < max_deg[1]) ? 0 : 1;
    core_nums[cur] = 1;
    core_nums[1 - cur] = max_deg[1 - cur];
    max_prod = max_deg[1 - cur];

    Graph g(*this);
    int delta = g.get_delta();
    printf("delta %d\n", delta);
    printf("max deg %d %d\n", max_deg[0], max_deg[1]);
    for (int i = 0; i < 2; i++) {
        for (int d = 1; d <= delta; d++) {
//            printf("d %d bin[d] %d \n", d, bin[i][d]);
            if (((long long) d) * max_deg[1 - i] <= max_prod)
                continue;
//            if (d < max_deg[i] && bin[i][d] == bin[i][d + 1])
//                continue;
            int d_opst;
            d_opst = skyline_core_num(i, d, static_cast<int>(max_prod / d + 1));
//            printf("%d %d\n", d, d_opst);
//            if (d_opst < d)
//                break;
            if (((long long) d) * d_opst > max_prod) {
                max_prod = d * d_opst;
                core_nums[i] = d;
                core_nums[1 - i] = d_opst;
            }
        }
    }
    printf("alpha %d beta %d\n", core_nums[0], core_nums[1]);
    return std::make_pair(core_nums[0], core_nums[1]);
}

int Graph::get_delta() {
    int i = bin[0][1], j = bin[1][1];
    int delta = 0;

//    for (int i = 0; i < n; i++) {
//        printf("i %d u %d deg %d\n", i, vert[0][i], deg[0][vert[0][i]]);
//    }

    while (i < n && j < n) {
        int cur, u;

        if (deg[0][vert[0][i]] <= deg[1][vert[1][j]]) {
            cur = 0;
            u = vert[0][i];
            i++;
        } else {
            cur = 1;
            u = vert[1][j];
            j++;
        }
        delta = std::max(delta, deg[cur][u]);
        for (auto v : adj[cur][u]) {
            if (deg[1 - cur][v] > deg[cur][u]) {
                dec_deg(1 - cur, v);
            }
        }
//        printf(" %d deg %d cur %d i %d u %d\n", delta, deg[cur][u], cur, i, u);
    }
    return delta;
}

void Graph::exact_ds_a(double l, double r, double bias, double sqrt_left, double sqrt_a, double sqrt_right,
                       std::pair<int, int> &ratio_b) {
    int alpha, beta;
    alpha = std::max(static_cast<int>(ceil(l / 2 / sqrt_right)), 1);
    beta = std::max(static_cast<int>(ceil(sqrt_left * l / 2)), 1);
    int it = 0;
    while (r - l > bias) {
        double mid = (l + r) / 2;

        //select eligible edges
        std::vector<std::pair<int, int>> edges;
        int cur = 0;
        if (alpha > max_deg[cur]) {
            r = mid;
            continue;
        }
        for (int i = bin[cur][alpha]; i < n; i++) {
            int u = vert[cur][i];
            int v = adj[cur][u][std::max(alpha - 1, 0)];
            if (deg[1 - cur][v] < beta) continue;
            for (auto t :adj[cur][u]) {
                if (deg[1 - cur][t] < beta)
                    break;
                edges.emplace_back(u, t);
            }
        }
        if (!size_reported) {
            printf("it %d, %lu\n", it++, edges.size());
        }
        //construct flow network and compute min cut
        FlowNetwork fn = FlowNetwork(n, edges, mid, sqrt_a, weights);
        std::vector<int> tmp_S, tmp_T;
        double edge_num;
        fn.get_mincut(0, fn.n - 1, tmp_S, tmp_T, edge_num, weights);

        //update l r

        if (!tmp_S.empty() && !tmp_T.empty()) {
            l = mid;
            ratio_b = std::make_pair(tmp_S.size(), tmp_T.size());
            alpha = std::max(static_cast<int>(ceil(l / 2 / sqrt_right)), 1);
            beta = std::max(static_cast<int>(ceil(sqrt_left * l / 2)), 1);
            if (density < edge_num / sqrt(tmp_S.size() * tmp_T.size())) {
                density = edge_num / sqrt(tmp_S.size() * tmp_T.size());
                S = tmp_S;
                T = tmp_T;
            }
        } else {
            r = mid;
        }
    }
//    size_reported = true;
}

void Graph::divide_conquer(double left, double right) {
    printf("lcover %.4f rcover %.4f\n", left, right);
    if (right < left + 1.0 / n)
        return;
    double mid = (left + right) / 2;

    std::pair<int, int> ratio_b;
    double bias = min_weight / sqrt((double) (n - bin[0][1]) * (n - bin[1][1]) - std::min(n - bin[0][1], n - bin[1][1]))
                  - min_weight / sqrt((double) (n - bin[0][1]) * (n - bin[1][1]));

    double sqrt_left = sqrt(left);
    double sqrt_right = sqrt(right);
    double a = mid;
    double sqrt_a = sqrt(mid);
//    printf("pre Ratio b : %d/%d\n", ratio_b.first, ratio_b.second);
    exact_ds_a(0, r_max, bias, sqrt_left, sqrt_a, sqrt_right, ratio_b);
//    exact_ds_a(0, r_max, bias, sqrt_a, sqrt_a, sqrt_a, ratio_b);
    double b = mid;
    if (ratio_b.first != 0 && ratio_b.second != 0)
        b = (double) ratio_b.first / ratio_b.second;
    printf("Ratio b : %d/%d\n", ratio_b.first, ratio_b.second);
    double c = sqrt(b) / sqrt_a + sqrt_a / sqrt(b);
    if (b > mid) {
        //1/2 (-2 a + a c^2 - a Sqrt[-4 c^2 + c^4])
        double mid_cover = (-2 * a + a * c * c - a * sqrt(pow(c, 4) - 4 * c * c)) / 2;
        divide_conquer(left, mid_cover);
        if (b < right) {
            divide_conquer(b, right);
        }
    } else {
        //1/2 (-2 a + a c^2 + a Sqrt[-4 c^2 + c^4])
        double mid_cover = (-2 * a + a * c * c + a * sqrt(pow(c, 4) - 4 * c * c)) / 2;
        divide_conquer(mid_cover, right);
        if (left < b) {
            divide_conquer(left, b);
        }
    }
}

double Graph::approx_para_ds(double delta, double epsilon) {
    double c = delta;
    density = 0;
    while (c >= 1.0 / n) {
        peeling(c, epsilon);
        c /= delta;
    }

    c = delta * delta;
    while (c <= n) {
        peeling(c, epsilon);
        c *= delta;
    }

    return density;
}

void Graph::peeling(double c, double epsilon) {
    Graph tmp_g(*this);
    S = tmp_g.vert[0];
    T = tmp_g.vert[1];
    int e = m;
    int s = static_cast<int>(S.size());
    int t = static_cast<int>(T.size());
    density = e / sqrt((double) s * t);
    while (s > 0 && t > 0) {
        if (s >= t * c) {
            std::vector<int> tmp_s;
            int e_delta = 0;
            for (auto u : tmp_g.vert[0]) {
                if (tmp_g.deg[0][u] >= (1 + epsilon) * e / s) {
                    tmp_s.push_back(u);
                } else {
                    for (auto v : adj[0][u]) {
                        --tmp_g.deg[1][v];
                    }
                    e_delta += tmp_g.deg[0][u];
                }
            }
            tmp_g.vert[0] = tmp_s;
            s = static_cast<int>(tmp_g.vert[0].size());
            e -= e_delta;
        } else {
            std::vector<int> tmp_t;
            int e_delta = 0;
            for (auto u : tmp_g.vert[1]) {
                if (tmp_g.deg[1][u] >= (1 + epsilon) * e / t) {
                    tmp_t.push_back(u);
                } else {
                    for (auto v : adj[1][u]) {
                        --tmp_g.deg[0][v];
                    }
                    e_delta += tmp_g.deg[1][u];
                }
            }
            tmp_g.vert[1] = tmp_t;
            t = static_cast<int>(tmp_g.vert[1].size());
            e -= e_delta;
        }
        if (e > density * sqrt((double) s * t)) {
            density = e / sqrt((double) s * t);
            S = tmp_g.vert[0];
            T = tmp_g.vert[1];
        }
    }
}

double Graph::approx_icalp_ds() {
    // 0 for out edges, for in edges
    density = (double) m / n;
    int best_pos_in = 0;
    int best_pos_out = 0;
    int e = m;
    int pos_in = 0;
    int pos_out = 0;

    while (pos_in < n && pos_out < n) {
        int u = vert[1][pos_in];
        int v = vert[0][pos_out];
        if (deg[1][u] <= deg[0][v]) {
            e -= deg[1][u];
            for (auto w : adj[1][u]) {
                if (pos[0][w] >= pos_out) {
                    dec_deg(0, w);
                }
            }
            pos_in++;
        } else {
            e -= deg[0][v];
            for (auto w : adj[0][v]) {
                if (pos[1][w] >= pos_in) {
                    dec_deg(1, w);
                }
            }
            pos_out++;
        }
        if (e / sqrt((double) (n - pos_in) * (n - pos_out)) > density) {
            density = e / sqrt((double) (n - pos_in) * (n - pos_out));
            best_pos_in = pos_in;
            best_pos_out = pos_out;
        }
    }

    S.clear();
    T.clear();
    for (int i = best_pos_out; i < vert[0].size(); i++) {
        S.push_back(vert[0][i]);
    }
    for (int i = best_pos_in; i < vert[1].size(); i++) {
        T.push_back(vert[1][i]);
    }

    return density;
}

double Graph::approx_ficalp_ds() {
    density = 0;
    for (int i = 0; i < n; i++) {
        Graph g(*this);
        g.peel(i, 0);
        if (g.density > density) {
            density = g.density;
            S = g.S;
            T = g.T;
        }
    }
    for (int i = 0; i < n; i++) {
        Graph g(*this);
        g.peel(i, 1);
        if (g.density > density) {
            density = g.density;
            S = g.T;
            T = g.S;
        }
    }
    return density;
}

void Graph::peel(int k, int cur) {
    int s = n, t = n, m = this->m;
    for (int i = 0; i < k; i++){
        int u = vert[cur][i];
        for (auto v : adj[cur][u]) {
            dec_deg(1 - cur, v);
            dec_deg(cur, u);
            m--;
        }
    }
    s -= k;
    for (int i = 0; i < n; i++) {
        if (density < m / sqrt((double) s * t)) {
            density = m / sqrt((double) s * t);
            S.clear();
            T.clear();
            for (int t = 0; t < n; t++) {
                if (deg[0][t] != 0) {
                    S.push_back(t);
                }
                if (deg[1][t] != 0) {
                    T.push_back(t);
                }
            }
        }
        int u = vert[1-cur][i];
        for (auto v : adj[1-cur][u]) {
            if (deg[cur][v] > 0) {
                dec_deg(cur, v);
                m--;
                if (deg[cur][v] == 0) {
                    s--;
                }
            }
        }
        t--;
    }
}



void Graph::output_ds() {
//    printf("S:\n");
    FILE* dsFile = fopen("ds.txt", "w");
    fprintf(dsFile, "s %lu, t %lu\n", S.size(), T.size());
    fprintf(dsFile, "S:\n");
    for (auto u : S){
        fprintf(dsFile, "%d\n", u);
    }

    fprintf(dsFile, "T:\n");
    for (auto u : T){
        fprintf(dsFile, "%d\n", u);
    }

    fclose(dsFile);

}













