#include "linear_program.h"
#include <float.h>
#include <stdbool.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <GL/glew.h>
#include <GL/freeglut.h>

typedef enum { OPTIMAL, UNBOUNDED, INFEASIBLE } simp;

typedef enum { EQ, LE, GE } eq;

struct constraint {
   double *coefficients;
   double constant;
   eq type;
};

typedef struct constraint * Constraint;

struct lp {
   bool is_SEF;
   
   // non-SEF data
   Constraint *cons;
   int num_cons;
   int variables;
   bool *ge0;
   int num_not_ge0;
   bool obj_max;

   // SEF data
   Matrix A;
   Matrix b;
   Matrix c;
   double z;
};

void vis_reshape(int w, int h);
void vis_display();
void vis_timer(int);

int vis_current_height = 600;
int vis_current_width = 800;

unsigned int frame_count = 0;

const char *WINDOW_TITLE = "LP Visualization";

void vis_init(int *argc, char **argv) {
  glutInit(argc, argv);
  //glutInitContextVersion(4, 0);
  //glutInitContextFlags(GLUT_FORWARD_COMPATIBLE);
  glutInitContextProfile(GLUT_CORE_PROFILE);
  glutSetOption(
    GLUT_ACTION_ON_WINDOW_CLOSE,
    GLUT_ACTION_GLUTMAINLOOP_RETURNS
  );
  glutInitWindowSize(vis_current_width, vis_current_height);
  glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
}

void visualize(LP P) {
  assert(P->variables == 2);
  int vis_window = 0;
  vis_window = glutCreateWindow("LP Visualization");
  assert(vis_window);

  glewExperimental = GL_TRUE;
  GLenum err = glewInit();
  if(err != GLEW_OK) {
    fprintf(stderr, "Error: %s", glewGetErrorString(err));
  }

  glutReshapeFunc(vis_reshape);
  glutDisplayFunc(vis_display);
  glutTimerFunc(0, vis_timer, 0);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, vis_current_width, vis_current_height, 0.0, 1.0, -1.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0, 0, 0, 0);
 
  glBegin(GL_QUADS);
    glColor3f(1, 0, 0); glVertex2f(-50, -50);
    glColor3f(1, 1, 0); glVertex2f(50, -50);
    glColor3f(0, 1, 0); glVertex2f(50, 50);
    glColor3f(0, 0, 1); glVertex2f(-50, 50);
  glEnd();

  glutMainLoop();
}

void vis_reshape(int w, int h) {
  vis_current_width = w;
  vis_current_height = h;
  glViewport(0, 0, w, h);
}

void vis_display() {
  glClear(GL_COLOR_BUFFER_BIT);
  
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glTranslatef(vis_current_width/2, vis_current_height/2, 0);
  
  frame_count++;
}

void vis_timer(int v) {
  if(0 != v) {
    char *tmp = malloc(sizeof(char)*40 +
      sizeof(char)*strlen(WINDOW_TITLE));
    sprintf(tmp, "%s: %d fps, %d x %d", WINDOW_TITLE,
      frame_count*4, vis_current_width, vis_current_height);
    glutSetWindowTitle(tmp);
    free(tmp);
  }
  frame_count = 0;
  glutTimerFunc(250, vis_timer, 1);
}

Constraint create_constraint(int vars) {
   Constraint con = malloc(sizeof(struct constraint));
   con->coefficients = malloc(sizeof(double)*vars);
   return con;
}

void destroy_constraint(Constraint con) {
   free(con->coefficients);
   free(con);
}

void destroy_LP(LP P) {
   for(int i=0; i < P->num_cons; i++) {
      destroy_constraint(P->cons[i]);
   }
   free(P->cons);
   if(P->is_SEF) {
      destroy_matrix(P->A);
      destroy_matrix(P->b);
   }
   destroy_matrix(P->c);
}

LP copy_LP(LP P) {
   LP Q = malloc(sizeof(struct lp));
   Q->A = copy_matrix(P->A);
   Q->b = copy_matrix(P->b);
   Q->c = copy_matrix(P->c);
   Q->z = P->z;
   Q->obj_max = P->obj_max;
   return Q;
}

void canonical_form(LP P, int *B) {
   Matrix A_B = take_columns(P->A, B, P->A->n);
   invert(A_B);
   Matrix A_BT = copy_matrix(A_B);
   transpose(A_BT);
   Matrix c_B = take_rows(P->c, B, P->A->n);
   Matrix y = mult_new(A_BT, c_B);
   destroy_matrix(A_BT);
   destroy_matrix(c_B);
   
   transpose(y);
   
   Matrix z_diff = mult_new(y, P->b);
   P->z += single_to_num(z_diff);
   destroy_matrix(z_diff);
   
   
   Matrix c_diff = mult_new(y, P->A);
   
   multiply_scalar(c_diff, -1);
   transpose(c_diff);
   add_matrix(P->c, c_diff);
   destroy_matrix(c_diff);
   destroy_matrix(y);
   
   Matrix new_A = mult_new(A_B, P->A);
   Matrix new_b = mult_new(A_B, P->b);
   destroy_matrix(A_B);
   destroy_matrix(P->A);
   destroy_matrix(P->b);
   
   P->A = new_A;
   P->b = new_b;
}

bool check_basis(Matrix mtr, int *B) {
   Matrix A_B = take_columns(mtr, B, mtr->n);
   bool res = rref(A_B) == mtr->n;
   destroy_matrix(A_B);
   return res;
}

// subsets of size n out of numbers 0..m-1
// s : starting index in B
bool basis(Matrix mtr, int *B, int sn, int sm, int n, int m) {
   if(n==0) {
      return check_basis(mtr, B);
   }
   for(int i=sm; i < m-n+1; i++) {
      B[sn] = i;
      if(basis(mtr, B, sn+1, i, n-1, m)) {
         return true;
      }
   }
   return false;
}

simp simplex(LP P, int *B) {
   canonical_form(P, B);
   
   int k = -1;
   for(int i=0; k == -1 && i < P->c->n; i++) {
      if(P->c->data[i][0] > 0)
         k = i;
   }
   if(k == -1) return OPTIMAL;
   
   Matrix A_k = take_columns(P->A, &k, 1);
   bool unbounded = true;
   int r;
   double t = DBL_MAX;
   for(int i=0; i < A_k->n; i++) {
      if(A_k->data[i][0] > 0) {
         unbounded = false;
         int temp = P->b->data[i][0] / A_k->data[i][0];
         if(temp < t) {
            t = temp;
            r = i;
         }
      }
   }
   if(unbounded) return UNBOUNDED;
   
   B[r] = k;
   return simplex(P, B);
}

void SEF(LP P) {
   assert(!P->is_SEF);
   if(!P->obj_max) {
      P->z *= -1;
      multiply_scalar(P->c, -1);
      P->obj_max = true;
   }
   
   int ineqs = 0;
   for(int i=0; i < P->num_cons; i++) {
      if(P->cons[i]->type != EQ) {
         ineqs++;
      }
   }
   P->A = create_matrix(P->num_cons, P->variables+ P->num_not_ge0 + ineqs);
   P->b = create_matrix(P->num_cons, 1);
   Matrix new_c = create_matrix(P->variables + P->num_not_ge0 + ineqs, 1);
   for(int i=0, v=0; i < P->variables; v++, i++) {
      new_c->data[i][0] = P->c->data[v][0];
      
      if(!P->ge0[v]) {
         i++;
         new_c->data[i][0] = -1 * P->c->data[v][0];
      }
      
      for(int j=0; j < P->num_cons; j++) {
         P->A->data[j][i] = P->cons[j]->coefficients[v];
         if(!P->ge0[v])
            P->A->data[j][i] = -1 * P->cons[j]->coefficients[v];
      }
   }
   
   destroy_matrix(P->c);
   P->c = new_c;
   
   int count_ineq=0;
   for(int i=0; i < P->num_cons; i++) {
      P->b->data[i][0] = P->cons[i]->constant;
      if(P->cons[i]->type != EQ) {
         bool le = P->cons[i]->type == LE;
         P->A->data[i][P->variables*2+count_ineq] = le ? 1 : -1;
         count_ineq++;
      }
   }
   
   P->is_SEF = true;
}

// P must be in canonical form for basis B
Matrix basic_solution(LP P, int *B) {
   Matrix x = create_matrix(P->A->m, 1);
   for(int i=0; i < P->A->n; i++) {
      x->data[B[i]][0] = P->b->data[i][0];
   }
   return x;
}

Matrix solve(LP P) {
   assert(!P->is_SEF);
   SEF(P);
   LP Q = copy_LP(P);
   Matrix im = create_identity(P->num_cons);
   join_right(Q->A, im);
   destroy_matrix(im);
   
   destroy_matrix(Q->c);
   Q->c = create_matrix(Q->A->m, 1);
   
   int *QB = malloc(sizeof(int)*P->A->n);
   for(int i=0; i < P->A->n; i++) {
      Q->c->data[i+P->A->m][0] = -1;
      QB[i] = P->A->m + i;
   }
   
   
   simplex(Q, QB);
   Matrix x = basic_solution(Q, QB);
   bool feasible = true;
   for(int i=0; i < P->A->n; i++) {
      if(x->data[i+P->A->m][0] != 0)
         feasible = false;
   }
   destroy_matrix(x);
   if(!feasible) return NULL;
   
   simp final = simplex(P, QB);
   if(final == UNBOUNDED) return NULL;
   
   x = basic_solution(P, QB);
   Matrix rx = create_matrix(P->variables, 1);
   for(int v=0, i=0; v < P->variables; v++, i++) {
      if(P->ge0[v]) {
         rx->data[v][0] = x->data[i][0];
      }
      else {
         rx->data[v][0] = x->data[i][0] - x->data[i+1][0];
         i++;
      }
   }
   destroy_matrix(x);
   return rx;
}

LP read_LP() {
   LP P = malloc(sizeof(struct lp));
   
   P->is_SEF = false;
   P->z = 0;
   P->num_not_ge0 = 0;
   
   // max or min
   int maxmin;
   scanf("%d", &maxmin);
   if(maxmin==0) P->obj_max = false;
   else P->obj_max = true;
   
   // num vars
   scanf("%d", &P->variables);
   
   // read in c
   P->c = create_matrix(P->variables, 1);
   for(int i=0; i < P->variables; i++) {
      scanf("%lf", &P->c->data[i][0]);
   }
   
   // num constraints
   scanf("%d", &P->num_cons);
   
   P->cons = malloc(sizeof(struct constraint *)*P->num_cons);
   // read in constraints
   for(int i=0; i < P->num_cons; i++) {
      P->cons[i] = create_constraint(P->variables);
      for(int v=0; v < P->variables; v++) {
         scanf("%lf", &P->cons[i]->coefficients[v]);
      }
      char t[3];
      scanf("%s", t);
      if(t[0]=='=') P->cons[i]->type = EQ;
      else if(t[0]=='>') P->cons[i]->type = GE;
      else P->cons[i]->type = LE;
      scanf("%lf", &P->cons[i]->constant);
   }
   
   // read in which vars >= 0
   P->ge0 = malloc(sizeof(bool)*P->variables);
   for(int i=0; i < P->variables; i++) {
      int g0;
      scanf("%d", &g0);
      if(g0==0) {
         P->ge0[i] = false;
         P->num_not_ge0++;
      }
      else
         P->ge0[i] = true;
   }
   
   return P;
}

void print_LP(LP P) {
   if(P->obj_max)
      printf("max");
   else
      printf("min");
   printf(" %5.1f + [ ", P->z);
   for(int i=0; i < P->c->n; i++) {
      printf("%5.1f ", P->c->data[i][0]);
   }
   printf("]x\n");
   for(int i=0; i < P->A->n; i++) {
      printf("[ ");
      for(int j=0; j < P->A->m; j++) {
         printf("%5.1f ", P->A->data[i][j]);
      }
      if(i==P->A->n/2)
         printf("] x = ");
      else
         printf("]     ");
      printf("[ %5.1f ]\n", P->b->data[i][0]);
   }
   printf("\n");
}
