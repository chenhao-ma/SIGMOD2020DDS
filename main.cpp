#include <iostream>
#include "Graph.h"
#include "Args.h"

int main(int argc, char *argv[]) {
    clock_t begin = clock();
    Args *args = new Args();
    args->parse_args(argc, argv);
    FILE* dFile = fopen(args->address, "r");
    Graph g = Graph(dFile, args->weighted);
    clock_t io_end = clock();
    if (args->baseline) {
        if (args->exact) {
            printf("density %.2f\n", g.exact_ds());
        } else {
            printf("density %.2f\n", g.approx_ds());
        }
    } else {
        if (args->exact) {
            if (args->advanced) {
                printf("density %.2f\n", g.exact_ad_core_ds());
            } else {
                printf("density %.2f\n", g.exact_core_ds());
            }
        } else {
            if (args->parameterized) {
                printf("density %.2f\n", g.approx_para_ds(args->delta, args->epsilon));
            } else if (args->icalp) {
                printf("density %.2f\n", g.approx_icalp_ds());
            } else if (args->ficalp) {
                printf("density %.2f\n", g.approx_ficalp_ds());
            } else {
		printf("density %.2f\n", g.approx_core_ds());
	    }

        }
    }

    clock_t end = clock();
    double io_secs = double(io_end - begin) / CLOCKS_PER_SEC;
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    g.output_ds();
    printf("io time: %.4f, total time: %.4f\n", io_secs, elapsed_secs);

    return 0;
}
