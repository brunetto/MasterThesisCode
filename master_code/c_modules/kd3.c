#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define POWER 2
#define LEAFSIZE = 400

/*
 Points
 */

typedef float xtp;

typedef struct{
  xtp x;
  xtp y;
  xtp z;
} point;

int particles_allocation(point **set, int n_point){
  printf("Allocating particles...");
  *set = (point*) malloc(sizeof(point)*n_point);
  assert(*set != NULL);
  return 0;
}

void particles_deallocation(point *set){
  free(set);
}

int minkowsky_distance(point *a , point *b, int *result, POWER){
  for(int i=0; i<0; ++i){
    xtp dx = a->x - b->x;
    xtp dy = a->y - b->y;
    xtp dz = a->z - b->z;
    result[i] = sqrt(dx**2 + dy**2 + dz**2);
  }
  return 0;
}

/*
 Rectangles
 */

typedef struct{
  xtp *maxes, *mins;
} rectangle;

int rect_max_dist(){
  //FIXME
  return 0;
}

int rect_min_dist(){
  //FIXME
  return 0;
}

/*
 Tree
 */

typedef struct{
  /*FIXME e se dentro ci mettessi direttamente il vettore di struct???*/
//   int *indexes;
  int n_elements;
  int split_dim;
  int tag;
  node *sup;
  node *inf;
} node;

int max(point *data, int length, dim){
  /* Find the maximum of an array of float */
  int max = 0;
  switch(dim){
    case 0:{
      for(int i=0; i<length; ++i){
        max = (max > data[i].x) ? max : data[i].x;
      }
      break;
    }
    case 1:{
      for(int i=0; i<length; ++i){
        max = (max > data[i].y) ? max : data[i].y;
      }
      break;
    }
    case 2:{
      for(int i=0; i<length; ++i){
        max = (max > data[i].z) ? max : data[i].z;
      }
      break;
    }
  }
  return max;
}

int min(point *data, int length, dim){
  /* Find the minimum of an array of float */
  int min =0;
  switch(dim){
    case 0:{
      for(int i=0; i<length; ++i){
        min = (min < data[i].x) ? min : data[i].x;
      }
      break;
    }
    case 1:{
      for(int i=0; i<length; ++i){
        min = (min < data[i].y) ? min : data[i].y;
      }
      break;
    }
    case 2:{
      for(int i=0; i<length; ++i){
        min = (max < data[i].z) ? min : data[i].z;
      }
      break;
    }
  }
  return min;
}

int max_id(point *data, int length, dim){
  /* Find the index of the maximum of an array of float */
  int max_id =0;
  switch(dim){
    case 0:{
      for(int i=0; i<length; ++i){
        max_id = (data[max_id].x > data[i].x) ? max_id : i;
      }
      break;
    }
    case 1:{
      for(int i=0; i<length; ++i){
        max_id = (data[max_id].y > data[i].y) ? max_id : i;
      }
      break;
    }
    case 2:{
      for(int i=0; i<length; ++i){
        max_id = (data[max_id].z > data[i].z) ? max_id : i;
      }
      break;
    }
  }
  return max_id;
}

int min_id(point *data, int length, dim){
  /* Find the id of the minimum of an array of float */
  int min_id =0;
  switch(dim){
    case 0:{
      for(int i=0; i<length; ++i){
        min_id = (data[min_id].x < data[i].x) ? min_id : i;
      }
      break;
    }
    case 1:{
      for(int i=0; i<length; ++i){
        min_id = (data[min_id].y > data[i].y) ? min_id : i;
      }
      break;
    }
    case 2:{
      for(int i=0; i<length; ++i){
        min_id = (data[min_id].z > data[i].z) ? min_id : i;
      }
      break;
    }
  }
  return min_id;
}

int tree_build(point *data, int length, node **tree){
  //FIXME
  /* Build the tree as a series of nested node structs.*/
  static tag; /* static because we want tag to be 
		unique across the whole tree */
  if(length <= LEAFSIZE){
    *tree = (node*) malloc(sizeof(node));
    assert(*tree != NULL);
    node mynode;/*FIXME controllare se si può 
		  fare in modo migliore*/
    *tree = mynode;
    *tree->n_elements = length;
    /* Find coord with max range */
    max_x
    
    
    *tree->split_dim = //FIXME;
    *tree->tag = tag++; /*Copy node tag and increment it */
    /* FIXME scoprire se si può fare così */
    /*node *sup = */tree_build(point *sup, n_sup, &(tree->sup));
    /*node *inf = */tree_build(point *inf, n_inf, &(tree->inf));
  }
  else{
    
  }
  return 0;
}

int count_neighbors(){
  //FIXME
  return 0;
}



