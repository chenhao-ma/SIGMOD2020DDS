//
// Created by MA Chenhao on 29/4/2019.
//

#include "Args.h"

Args::Args() {
    address = nullptr;
}

Args::Args(const Args &args) {
    address = strdup(args.address);
}

Args::~Args() {
    delete address;
}

void Args::usage(char *msg, int exit_status) {
    fprintf(exit_status == 0 ? stdout : stderr, "%s", USAGE_TXT);

    if (msg) {
        fprintf(exit_status == 0 ? stdout : stderr, "\n%s\n", msg);
    }
    exit(exit_status);
}

void Args::parse_args(int argc, char *argv[]) {
    int c;
    opterr = 0;
    if (argc < 2) {
        usage(nullptr, 0);
    }
    weighted = false;
    while ((c = getopt(argc, argv, "g:m:a:d:e:w")) != -1) {
        switch (c) {
            case 'g': {
                printf("%s\n", optarg);
                address = strdup(optarg);
                break;
            }

            case 'm': {
                char* tmp;
                tmp = strdup(optarg);
                baseline = tmp[0] == 'b';
                advanced = tmp[0] == 'a';
                parameterized = tmp[0] == 'p';
                icalp = tmp[0] == 'i';
                ficalp = tmp[0] == 'f';
                delete tmp;
                break;
            }

            case 'a': {
                char* tmp;
                tmp = strdup(optarg);
                exact = tmp[0] == 'e';
                delete tmp;
                break;
            }

            case 'd': {
                sscanf(optarg, "%lf", &delta);
                printf("delta %.2f\n", delta);
                break;
            }

            case 'e': {
                sscanf(optarg, "%lf", &epsilon);
                printf("epsilon %.2f\n", epsilon);
                break;
            }

            case 'w': {
                weighted = true;
                break;
            }

            default:
                break;
        }
    }

}