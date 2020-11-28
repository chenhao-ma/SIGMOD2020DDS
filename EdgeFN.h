//
// Created by MA Chenhao on 7/5/2019.
//

#ifndef DIRECTEDDENSESTSUBGRAPH_EDGEFN_H
#define DIRECTEDDENSESTSUBGRAPH_EDGEFN_H


class EdgeFN {
public:
    int from, to, index;
    double cap, flow;

    EdgeFN(int from, int to, double cap, double flow, double index);
};


#endif //DIRECTEDDENSESTSUBGRAPH_EDGEFN_H
