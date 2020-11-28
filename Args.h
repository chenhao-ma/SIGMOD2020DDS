//
// Created by MA Chenhao on 29/4/2019.
//

#ifndef DIRECTEDDENSESTSUBGRAPH_ARGS_H
#define DIRECTEDDENSESTSUBGRAPH_ARGS_H

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <getopt.h>

#define USAGE_TXT							   \
    "usage: \n"                                \
    "\t[-g input directed graph file]\n"      \
    "\t[-m method b for baseline, p for parameterized, c for core, a for advanced, i for icalp]\n" \
    "\t[-a accuracy e for exact, a for approximate]\n" \
    "\t[-d delta, used in parameterized VLDB'12]\n" \
    "\t[-e epsilon, used in parameterized VLDB'12]\n" \
    "\t[-w weighted]\n"

class Args {
public:
    char *address;
    bool baseline;
    bool exact;
    bool advanced;
    bool parameterized;
    bool icalp;
    bool weighted;
    bool ficalp;
    double delta;
    double epsilon;
    Args();
    Args(const Args& args);
    ~Args();
    void usage(char *msg, int exit_status);
    void parse_args(int argc, char *argv[]);


};


#endif //DIRECTEDDENSESTSUBGRAPH_ARGS_H
