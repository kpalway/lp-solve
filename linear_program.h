#ifndef __LINEAR_PROGRAM_H__
#define __LINEAR_PROGRAM_H__

#include "matrix.h"

struct lp;

typedef struct lp * LP;

Matrix solve(LP p);

LP read_LP();

void print_LP(LP P);

void destroy_LP();

void vis_init(int *argc, char **argv);

void visualize(LP P);

#endif
