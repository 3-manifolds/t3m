#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "vertex.h"

#define PRIME ((unsigned int)0x7fffffff)

static void print_vertices(vertex_stack_t *stack, int dimension);
static PyObject *build_vertex_list(vertex_stack_t *stack, int dimension);
static void no_memory(void);
static void add_block(reservoir_t *reservoir);

static void no_memory(void){
  fprintf(stderr, "No memory.\n");
  exit(ENOMEM);
}


reservoir_t *new_reservoir(int dimension){
  reservoir_t *result = calloc(1, sizeof(reservoir_t));

  if (result == NULL)
    no_memory();

  result->dimension = dimension;
  return result;
}


void destroy_reservoir(reservoir_t *reservoir){
  block_t *first_block, *next_block;

  if (reservoir == NULL)
    return;

  for(first_block = reservoir->blocklist; first_block != NULL; ) {
    next_block = first_block->next;
    free(first_block);
    first_block = next_block;
  }

  free(reservoir);
}


static void add_block(reservoir_t *reservoir){
  block_t *new;
  int i, vertex_size = sizeof(vertex_t) + reservoir->dimension*sizeof(int);
  char *p;
  vertex_t *v;

  new = calloc(1, sizeof(block_t) + BLOCKSIZE*vertex_size );
  if (new == NULL) no_memory();
  
  p = (char *)(new->vertices);
  v = (vertex_t *)p;
  for (i=0; i<BLOCKSIZE ; i++) {
    v->next = (vertex_t *)(p += vertex_size);
    v = v->next;
  }
  reservoir->first_vertex = new->vertices;
  new->next = reservoir->blocklist;
  reservoir->blocklist = new;
}


vertex_t *new_vertex(reservoir_t *reservoir){
  vertex_t *result = reservoir->first_vertex;

  if (result != NULL){
    reservoir->first_vertex = result->next;
    memset(result, 0, sizeof(vertex_t) + reservoir->dimension*sizeof(int));
    return result;
  }
  
  else {
    add_block(reservoir);
    return new_vertex(reservoir);
  }
}

void push_vertex(vertex_t *vertex, vertex_stack_t *stack){
  vertex->next = *stack;
  *stack = vertex;
}

vertex_t *pop_vertex(vertex_stack_t *stack){
  vertex_t *result = *stack;

  if (result != NULL) {
    *stack = result->next;
    result->next = NULL;
  }
  return result;
}


void recycle_vertices(vertex_stack_t *stack, reservoir_t *reservoir){
  vertex_t *v = *stack;

  if (v == NULL)
    return;

  while (v->next !=NULL)
    v = v->next;

  v->next = reservoir->first_vertex;
  reservoir->first_vertex = *stack;
  *stack = NULL;
}

void set_support(unsigned int index, support_t *support){
  int i;

  if (index & 0x01 != 0)
    i = 2 + (index >> 6);
  else 
    i = (index >> 6);
  support->supp[i] |= (1 << (0x1f & (index >> 1)) );
}

vertex_t *unit_vertex(unsigned int index, reservoir_t *reservoir){
  vertex_t *result;
  int i;

  result = new_vertex(reservoir);
  result->vector[index] = 1;
  set_support(index, &result->support);

  return result;
}

inline int gcd(int x, int y){
  int r;
  if (x == 0)
    return y;
  if (x < 0)  x = -x;
  if (y < 0)  y = -y;
  while (y != 0) {
    r = x%y;
    x = y;
    y = r;
  }
  return x;
}

void reduce(int dimension, vertex_t *v){
  int *a = v->vector;
  int i = dimension, j = i;
  int x, d = a[i - 1];

  while (i > 0) {
    --i;
    x = a[i];
    d = gcd(x,d);
  }

  if (d == 0) return;

  while (j > 0) {
    --j;
    x = a[j];
    a[j] = x/d;
  }
}

void evaluate(matrix_t *A, int row, vertex_t *v){
  int size = A->columns;

  if ( dot(size, A->matrix + size*row, v->vector, &v->value) !=0 )
       fprintf(stderr, "Overflow in dot product!\n");
}

void support_union(support_t *x, support_t *y, support_t *result){
  int i;
  for (i=0; i <4; i++) {
    result->supp[i] = x->supp[i] | y->supp[i];
  }
}

matrix_t *new_matrix(int rows, int columns){
  matrix_t *result;

  result = (matrix_t *)calloc(1, sizeof(matrix_t) + (rows*columns + 1)*sizeof(int));
  result->rows = rows;
  result->columns = columns;
  return result;
};

void *destroy_matrix(matrix_t *matrix){
  if (matrix != NULL)
    free(matrix);
};

static void print_matrix(matrix_t *matrix){
  int i,j, *p = matrix->matrix;

  for (i=0; i<matrix->rows; i++){ 
    for(j=0; j<matrix->columns; j++)
      printf("%d ", *p++);
    printf("\n");
  }
  printf("\n");
}

int test_corank(matrix_t *M, int threshold){
  int i, j, k, a, b, d;
  int corank = 0;
  int numcols = M->columns, numrows = M->rows;
  int *p, *A[numrows];

  // build a permuted row index
  p = M->matrix;
  for (i=0; i<numrows; i++){
    A[i] = p;
    p += numcols;
  }

  i = 0;
  for (j = 0; j < numcols; j++){

    // search for a non-zero entry in column j
    for(k=i; k<numrows; k++)
      if (A[k][j] != 0) break;

    // if there are none, increase the corank and continue
    if (k == numrows) {
      ++corank;
      if (corank > threshold)
	return -1;
    }

    else{
      // transpose rows
      p = A[k];
      A[k] = A[i];
      A[i] = p;

      // zero out the rest of the column
      a = A[i][j];
      for (k=i+1; k<numrows; k++){
	b = A[k][j];
	if (b != 0) {
	  d = gcd(a,b);
	  b = -b/d;
	  if ( ax_plus_by(numcols-j, b, a/d, A[i]+j, A[k]+j) != 0)
	    printf("Overflow in ax_plus_by!\n");
	}
      }

      // go onto the next row
      ++i;

      // but if we reached the last row then all rows are independent
      if (i == numrows) {
	corank = numcols - numrows;
	if (corank > threshold)
	  return -1;
	else
	  return corank;
      }

    } // end else
  } // end for j
  return corank;
}


// This computation is done mod p = 2^31 -1
int test_corank_mod_p(matrix_t *M, int threshold){
  int i, j, k, a, b;
  int corank = 0;
  int numcols = M->columns, numrows = M->rows;
  int *p, *A[numrows];

  // build a permuted row index
  p = M->matrix;
  for (i=0; i<numrows; i++){
    A[i] = p;
    p += numcols;
  }

  i = 0;
  for (j = 0; j < numcols; j++){

    // search for a non-zero entry in column j
    for(k=i; k<numrows; k++)
      if (A[k][j] != 0) break;

    // if there are none, increase the corank and continue
    if (k == numrows) {
      ++corank;
      if (corank > threshold)
	return -1;
    }

    else{
      // transpose rows
      p = A[k];
      A[k] = A[i];
      A[i] = p;

      // zero out the rest of the column
      a = A[i][j];
      for (k=i+1; k<numrows; k++){
	b = A[k][j];
	if (b != 0) {
	  b = PRIME - b;
	  ax_plus_by_mod_p(numcols-j, b, a, A[i]+j, A[k]+j);
	}
      }

      // go onto the next row
      ++i;

      // but if we reached the last row then all rows are independent
      if (i == numrows) {
	corank = numcols - numrows;
	if (corank > threshold)
	  return -1;
	else
	  return corank;
      }

    } // end else
  } // end for j
  return corank;
}

filter_list_t *embedded_filter(int tets){
  filter_list_t *result;
  int i,size = 3*tets;

  if (tets > 42){
    printf("Too many tetrahedra!\n");
    exit(-1);
  }

  // NOTE: this assumes that calloc will give us 8-byte aligned memory.
  // gcc says it will do that.
  // http://www.gnu.org/manual/glibc-2.2.3/html_node/libc_30.html
  //  result = (filter_list_t *)calloc(1,sizeof(filter_list_t)+size*sizeof(support_t));
  result = (filter_list_t *)calloc(1,8+sizeof(filter_list_t)+size*sizeof(support_t));
  if (result == NULL) no_memory();

  result->size = size;

  for (i=0; i<size; i+=3) {
    set_support(i,   result->filter + i);
    set_support(i+1, result->filter + i);
    set_support(i,   result->filter + i+1);
    set_support(i+2, result->filter + i+1);
    set_support(i+1, result->filter + i+2);
    set_support(i+2, result->filter + i+2);
  }
  return result;
}

void *destroy_filter_list(filter_list_t *filterlist){
  if (filterlist != NULL)
    free(filterlist);
};


PyObject *find_vertices(matrix_t *matrix, filter_list_t *filter_list){
  PyObject *result;
  int i, x, dimension = matrix->columns, size = dimension*matrix->rows;
  int slice, filtered = 0, interior = 0;
  int numpos, numneg, numzero;
  int *p;
  vertex_stack_t positives = NULL, negatives = NULL, zeros = NULL, current = NULL;
  reservoir_t *reservoir = new_reservoir(dimension);
  vertex_t *vertex = NULL, *P, *N;
  matrix_t *temp_matrix = new_matrix(matrix->rows, matrix->columns);

  for (i = 0; i < dimension; i++) {
    push_vertex( unit_vertex(i, reservoir), &current);
  }

  for (slice = 0; slice < matrix->rows; slice++) {
    printf("slice %3d : ", slice);
    numpos = 0;
    numneg = 0;
    numzero = 0;
    while ( vertex = pop_vertex(&current) ){
      evaluate(matrix, slice, vertex);
      if (vertex->value > 0) {
	push_vertex(vertex, &positives);
	++numpos;
      }
      if (vertex->value < 0) {
	push_vertex(vertex, &negatives);
	++numneg;
      }
      if (vertex->value == 0) {
	push_vertex(vertex, &zeros);
	++numzero;
      }
    }

    printf(" %5d positive %5d negative %5d zero\n",numpos, numneg, numzero);

    for (P = positives; P != NULL; P = P->next) {
      for (N = negatives; N != NULL; N = N->next) {
	if (vertex == NULL)
	  vertex = new_vertex(reservoir);

	support_union(&P->support, &N->support, &vertex->support);

	if (filter(vertex, filter_list) != 0) {
	  if (extract_matrix(matrix, slice+1, &vertex->support, temp_matrix) == 1
              && test_corank(temp_matrix, 1) == 1){
	    for ( i=0; i<dimension; i++ )
	      vertex->vector[i] = N->vector[i];
	    ax_plus_by(dimension, -(N->value), P->value, P->vector, vertex->vector);
	    reduce(dimension, vertex);
	    push_vertex(vertex, &zeros);
	    vertex = NULL;
	  }
	  else ++interior;
	}
	else ++filtered;
      }
    }

    current = zeros;
    zeros = NULL;
    recycle_vertices(&positives, reservoir);
    recycle_vertices(&negatives, reservoir);
  }  

  printf("DONE.  %d vertices were filtered;   %d were interior. \n",
	 filtered, interior);
  result = build_vertex_list(&current, dimension);
  recycle_vertices(&current, reservoir);
  destroy_reservoir(reservoir);
  reservoir = NULL;
  destroy_matrix(temp_matrix);
  return result;
}

PyObject *find_vertices_mod_p(matrix_t *matrix, filter_list_t *filter_list){
  PyObject *result;
  int i, x, dimension = matrix->columns, size = dimension*matrix->rows;
  int slice, filtered = 0, interior = 0;
  int numpos, numneg, numzero;
  int *p;
  vertex_stack_t positives = NULL, negatives = NULL, zeros = NULL, current = NULL;
  reservoir_t *reservoir = new_reservoir(dimension);
  vertex_t *vertex = NULL, *P, *N;
  matrix_t *temp_matrix = new_matrix(matrix->rows, matrix->columns);
  matrix_t *mod_p_matrix = new_matrix(matrix->rows, matrix->columns);

  for (i = 0; i < size; i++) {
    x = matrix->matrix[i];
    if (x < 0)
      x = PRIME + x;
    mod_p_matrix->matrix[i] = x;
  }

  for (i = 0; i < dimension; i++) {
    push_vertex( unit_vertex(i, reservoir), &current);
  }

  for (slice = 0; slice < matrix->rows; slice++) {
    printf("slice %3d : ", slice);
    numpos = 0;
    numneg = 0;
    numzero = 0;
    while ( vertex = pop_vertex(&current) ){
      evaluate(matrix, slice, vertex);
      if (vertex->value > 0) {
	push_vertex(vertex, &positives);
	++numpos;
      }
      if (vertex->value < 0) {
	push_vertex(vertex, &negatives);
	++numneg;
      }
      if (vertex->value == 0) {
	push_vertex(vertex, &zeros);
	++numzero;
      }
    }

    printf(" %5d positive %5d negative %5d zero\n",numpos, numneg, numzero);

    for (P = positives; P != NULL; P = P->next) {
      for (N = negatives; N != NULL; N = N->next) {
	if (vertex == NULL)
	  vertex = new_vertex(reservoir);

	support_union(&P->support, &N->support, &vertex->support);

	if (filter(vertex, filter_list) != 0) {
        if (extract_matrix(mod_p_matrix, slice+1, &vertex->support, temp_matrix) == 1
	    && test_corank_mod_p(temp_matrix, 1) == 1){
	    for ( i=0; i<dimension; i++ )
	      vertex->vector[i] = N->vector[i];
	    ax_plus_by(dimension, -(N->value), P->value, P->vector, vertex->vector);
	    reduce(dimension, vertex);
	    push_vertex(vertex, &zeros);
	    vertex = NULL;
	  }
	  else ++interior;
	}
	else ++filtered; 
      }
    }

    current = zeros;
    zeros = NULL;
    recycle_vertices(&positives, reservoir);
    recycle_vertices(&negatives, reservoir);
  }  

  printf("DONE.  %d vertices were filtered;   %d were interior. \n",
	 filtered, interior);
  result = build_vertex_list(&current, dimension);
  recycle_vertices(&current, reservoir);
  destroy_reservoir(reservoir);
  reservoir = NULL;
  destroy_matrix(temp_matrix);
  destroy_matrix(mod_p_matrix);
  return result;
}



static void print_vertices(vertex_stack_t *stack, int dimension){
  vertex_t *V = *stack;
  int i;

  for (; V != NULL; V = V->next ) {
    printf("[ ");
    for (i=0; i<dimension ; i++)
      printf("%d ", V->vector[i]);
    printf("]\n");
  }
}


static PyObject *build_vertex_list(vertex_stack_t *stack, int dimension){
  PyObject *result, *coeff;
  vertex_t *V = *stack;
  int i;

  result = PyList_New(0);
  if (result != NULL){
    for (; V != NULL; V = V->next ) {
      coeff = PyTuple_New(dimension);
      for (i=0; i<dimension ; i++)
	PyTuple_SetItem(coeff, i, PyInt_FromLong((long)V->vector[i]));
      PyList_Append(result, coeff);
    }
  }
  return result;
}