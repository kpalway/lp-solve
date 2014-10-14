#ifndef __MATRIX_H__
#define __MATRIX_H__

struct matrix {
   int n, m;
   double **data;
};

typedef struct matrix * Matrix;

/*
   create_matrix
   PRE: n, m > 0
   POST: produces a Matrix
   create_matrix(n, m) creates an n*m zero matrix
*/
Matrix create_matrix(int n, int m);

void set_data(Matrix mtr, double **data);

void transpose(Matrix m);

int rref(Matrix m);

void invert(Matrix m);

Matrix create_identity(int n);

void matrix_set_entry(Matrix m, int i, int j, double x);

double matrix_get_entry(Matrix m, int i, int j);

Matrix copy_matrix(Matrix m);

void destroy_matrix(Matrix m);

void join_right(Matrix to, Matrix from);

void join_bottom(Matrix to, Matrix from);

void print_matrix(Matrix m);

Matrix read_matrix();

Matrix mult_new(Matrix mtr1, Matrix mtr2);

void multiply_scalar(Matrix mtr, double s);

void add_matrix(Matrix to, Matrix from);

Matrix take_columns(Matrix mtr, int *cols, int m);

Matrix take_rows(Matrix mtr, int *rows, int n);

double single_to_num(Matrix mtr);

#endif
