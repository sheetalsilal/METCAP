
#include<R.h>

//*********************************************************************
  // EQ - derivative function 
void EQ(int *L, int *N, double *oldeq, double *transit, int *transitionsiu1, int *transitionsiu2, int *transitionsiv1, int *transitionsiv2, double *eq){
  int i,iu1, iu2, iv1, iv2;
  
  for(i=0;i<*N;i++){
    eq[i]=oldeq[i];
  }
  
  for(i=0;i < *L;i++){
    
    iu1 = transitionsiu1[i]-1;
    iv1 = transitionsiv1[i];
    iu2 = transitionsiu2[i]-1;
    iv2 = transitionsiv2[i];
    
    eq[iu1] = eq[iu1] + transit[i]*iv1;
    eq[iu2] = eq[iu2] + transit[i]*iv2;
    
  }  
}
