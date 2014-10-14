#include "matrix.h"
#include "linear_program.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
   
   LP P = read_LP();
   Matrix mtr = solve(P);
   
   if(mtr != NULL) {
      print_matrix(mtr);
      destroy_matrix(mtr);
   }
   else
      printf("NULL\n");
      
   destroy_LP(P);
   return 0;
}
