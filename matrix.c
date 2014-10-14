#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

void destroy_data(double **d, int n);
double **create_data(int n, int m);

Matrix create_matrix(int n, int m) {
   Matrix mtr = malloc(sizeof(struct matrix));
   mtr->data = create_data(n, m);
   mtr->n = n;
   mtr->m = m;
   return mtr;
}

double **create_data(int n, int m) {
   double **d = malloc(sizeof(double*) * n);
   for(int i=0; i < n; i++) {
      d[i] = malloc(sizeof(double) * m);
      for(int j=0; j < m; j++) {
         d[i][j] = 0;
      }
   }
   return d;
}

void destroy_data(double **d, int n) {
   for(int i=0; i < n; i++) {
      free(d[i]);
   }
   free(d);
}

void set_data(Matrix mtr, double **data) {
   mtr->data = data;
}

void destroy_matrix(Matrix mtr) {
   assert(mtr != NULL);
   destroy_data(mtr->data, mtr->n);
   free(mtr);
   mtr = NULL;
}

Matrix copy_matrix(Matrix mtr) {
   Matrix c = create_matrix(mtr->n, mtr->m);
   for(int i=0; i < mtr->n; i++) {
      for(int j=0; j < mtr->m; j++) {
         c->data[i][j] = mtr->data[i][j];
      }
   }
   return c;
}

void matrix_set_entry(Matrix m, int i, int j, double x) {
   m->data[i][j] = x;
}

double matrix_get_entry(Matrix m, int i, int j) {
   return m->data[i][j];
}

void transpose(Matrix mtr) {
   double **nd = create_data(mtr->m, mtr->n);
   for(int i=0; i < mtr->n; i++) {
      for(int j=0; j < mtr->m; j++) {
         nd[j][i] = mtr->data[i][j];
      }
   }
   destroy_data(mtr->data, mtr->n);
   mtr->data = nd;
   int temp = mtr->m;
   mtr->m = mtr->n;
   mtr->n = temp;
}

void mult_row(double *row, int m, double factor) {
   for(int i=0; i < m; i++) {
      row[i] *= factor;
   }
}

void add_row(double *row1, double *row2, int m, double factor) {
   for(int i=0; i < m; i++) {
      row1[i] += factor * row2[i];
   }
}

void swap_rows(double *row1, double *row2, int m) {
   int temp;
   for(int i=0; i < m; i++) {
      temp = row1[i];
      row1[i] = row2[i];
      row2[i] = temp;
   }
}

// returns rank(mtr)
int rref(Matrix mtr) {
   int cols = 0;
   int leading_ones = 0;
   int total_cols = mtr->m < mtr->n ? mtr->m : mtr->n;
   
   while(cols < total_cols) {
      if(mtr->data[leading_ones][cols] == 0) {
         int row2;
         for(row2 = leading_ones; row2 < mtr->n; row2++) {
            if(mtr->data[row2][cols] != 0) break;
         }
         if(row2 < mtr->n)
            swap_rows(mtr->data[leading_ones], mtr->data[row2], mtr->m);
         else {
            cols++;
            continue;
         }
      }
      
      if(mtr->data[leading_ones][cols] != 1) {
         mult_row(mtr->data[leading_ones], mtr->m,
            1 / mtr->data[leading_ones][cols]);
      }
      
      for(int i=0; i < mtr->n; i++) {
         if(i != leading_ones && mtr->data[i][cols] != 0) {
            add_row(mtr->data[i], mtr->data[leading_ones], mtr->m,
               -mtr->data[i][cols]);
         }
      }
      cols++;
      leading_ones++;
   }
   return leading_ones;
}

void join_right(Matrix to, Matrix from) {
   assert(to->n == from->n);
   double **data = create_data(to->n, to->m + from->m);
   for(int i=0; i < to->n; i++) {
      for(int j=0; j < to->m; j++) {
         data[i][j] = to->data[i][j];
      }
      for(int j=0; j < from->m; j++) {
         data[i][to->m + j] = from->data[i][j];
      }
   }
   destroy_data(to->data, to->n);
   to->data = data;
   to->m += from->m;
}

// not implemented
double determinant(Matrix mtr) {
   assert(mtr->n == mtr->m);
   if(mtr->n==2) {
      return mtr->data[0][0] * mtr->data[1][1]
         - mtr->data[0][1] * mtr->data[1][0];
   }
   else {
      return 0;
   }
}

void invert(Matrix mtr) {
   assert(mtr->m == mtr->n);
   Matrix im = create_identity(mtr->n);
   join_right(mtr, im);
   destroy_matrix(im);
   rref(mtr);
   int *cols = malloc(sizeof(int)*mtr->n);
   for(int i=0; i < mtr->n; i++) {
      cols[i] = i + mtr->n;
   }
   Matrix inv = take_columns(mtr, cols, mtr->n);
   mtr->data = inv->data;
   mtr->m = inv->m;
   free(cols);
   free(inv);
}

Matrix mult_new(const Matrix mtr1, const Matrix mtr2) {
   assert(mtr1->m == mtr2->n);
   Matrix mtr = malloc(sizeof(struct matrix));
   mtr->data = create_data(mtr1->n, mtr2->m);
   mtr->n = mtr1->n;
   mtr->m = mtr2->m;
   for(int i=0; i < mtr1->n; i++) {
      for(int j=0; j < mtr2->m; j++) {
         for(int p=0; p < mtr1->m; p++) {
            mtr->data[i][j] += mtr1->data[i][p]*mtr2->data[p][j];
         }
      }
   }
   return mtr;
}

void add_matrix(Matrix to, Matrix from) {
   assert(to->n == from->n);
   assert(to->m == from->m);
   for(int i=0; i < to->n; i++) {
      for(int j=0; j < to->m; j++) {
         to->data[i][j] += from->data[i][j];
      }
   }
}

void multiply_scalar(Matrix mtr, double s) {
   for(int i=0; i < mtr->n; i++) {
      for(int j=0; j < mtr->m; j++) {
         mtr->data[i][j] *= s;
      }
   }
}

Matrix create_identity(int n) {
   Matrix im = malloc(sizeof(struct matrix));
   im->data = create_data(n, n);
   im->n = im->m = n;
   for(int i=0; i < n; i++) {
      im->data[i][i] = 1;
   }
   return im;
}

Matrix take_columns(Matrix mtr, int *cols, int m) {
   for(int i=0; i < m; i++) {
      assert(cols[i] < mtr->m);
   }
   Matrix res = malloc(sizeof(struct matrix));
   res->data = create_data(mtr->n, m);
   res->n = mtr->n;
   res->m = m;
   for(int i=0; i < mtr->n; i++) {
      for(int j=0; j < m; j++) {
         res->data[i][j] = mtr->data[i][cols[j]];
      }
   }
   return res;
}

Matrix take_rows(Matrix mtr, int *rows, int n) {
   Matrix res = malloc(sizeof(struct matrix));
   res->data = create_data(n, mtr->m);
   res->n = n;
   res->m = mtr->m;
   for(int i=0; i < n; i++) {
      for(int j=0; j < mtr->m; j++) {
         res->data[i][j] = mtr->data[rows[i]][j];
      }
   }
   return res;
}

void print_matrix(Matrix mtr) {
   for(int i=0; i < mtr->n; i++) {
      printf("[ ");
      for(int j=0; j < mtr->m; j++) {
         printf("%5.1f ", mtr->data[i][j]);
      }
      printf("]\n");
   }
   printf("\n");
}

Matrix read_matrix() {
   int n, m;
   scanf("%d", &n);
   scanf("%d", &m);
   Matrix mtr = create_matrix(n, m);
   for(int i=0; i < n; i++) {
      for(int j=0; j < m; j++) {
         scanf("%lf", &mtr->data[i][j]);
      }
   }
   return mtr;
}

double single_to_num(Matrix mtr) {
   assert(mtr->m == 1 && mtr->n == 1);
   return mtr->data[0][0];
}
