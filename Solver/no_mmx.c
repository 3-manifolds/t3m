#include "vertex.h"

// This computation is done mod p = 2^31 -1

void ax_plus_by_mod_p(int size, int a, int b, int *x, int *y){

  register int A = a, B = b;
  register unsigned int S;
  register int *X = x, *Y = y;
  register long long prod;

  while (size--) {
    prod = ((long long)(*Y))*((long long)B);
    S = (prod & 0x7fffffff) + (prod >> 31);
    if (S >= (unsigned int)PRIME) S -= PRIME;
    *Y = S;
    prod = ((long long)(*X++))*((long long)A);
    S = (prod & 0x7fffffff) + (prod >> 31);
    if (S >= (unsigned int)PRIME) S -= PRIME;
    S += *Y;
    if (S >= (unsigned int)PRIME) S -= PRIME;
    *Y++ = S;
  }
}

// This returns a non-zero value if the value of any coordinate of
// ax+by overflows a 32 bit word.  Note that the products ax and by
// and their sum are computed as 64 bit numbers.  So the products can
// overflow without causing an overflow of the final linear combination.

int ax_plus_by(int size, int a, int b, int *x, int *y){

  register int A = a, B = b;
  register int *X = x, *Y = y;
  register long long prod;

  int result = 0;
  while (size--) {
    prod = ((long long)(*Y))*((long long)B) + ((long long)(*X++))*((long long)A);
    *Y++ = (prod & 0xffffffff);
    result |= (((prod >> 32) + 1) >> 1);
  }
  return result;
}

// This returns a non-zero value if any partial sum of the dot product
// overflows a 32 bit integer.

int dot(int size, int *x, int *y, int *dotprod){
  register int result = 0;
  register int *X = x, *Y = y;
  register long long accumulator = 0;
  while (size--){
    accumulator += ((long long)*X++)*((long long)*Y++);
    result |= (((accumulator >> 32) + 1) >> 1);
  }
  *dotprod = accumulator & 0xffffffff;
  return result;
}

// This extracts those columns of the input matrix specified by the
// support vector.  If dimension considerations show that the the
// resulting matrix could not possibly have co-rank 1, the extraction
// is aborted and 0 is returned (leaving garbage in the output
// matrix). Otherwise the return value is 1.
//
// WARNING: This may write one int past the end of the array.  Allow
// extra space in your output matrix!!!
//
// To avoid unpredictable branching, we just overwrite each of the
// columns that corresponds to a 0 in the support vector.

int extract_matrix(matrix_t *in, int rows, support_t *support, matrix_t *out) {
  register int *in_coeff = in->matrix;
  register int *out_coeff = out->matrix;
  register int supp1, supp2, temp;
  int count1, count2, count, columns_out = 0;

  if (in->columns > 64){
    count1 = 64;
    count2 = in->columns - 64;
  }
  else{
    count1 = in->columns;
    count2 = 0;
  }

  out->rows = rows;
  // We have to look at the first row to count how many columns the
  // output matrix will have.
  supp1 = support->supp[0];
  supp2 = support->supp[2];
  count = count1;
  while (count--){
    *out_coeff = *in_coeff++;
    out_coeff += (supp1 & 0x1);
    temp = supp1 >> 1;
    supp1 = supp2;
    supp2 = temp;
  }
  supp1 = support->supp[1];
  supp2 = support->supp[3];
  count = count2;
  while (count--){
    *out_coeff = *in_coeff++;
    out_coeff += (supp1 & 0x1);
    temp = supp1 >> 1;
    supp1 = supp2;
    supp2 = temp;
  }
  columns_out = out_coeff - out->matrix;
  // Bail out if there aren't enough rows to have co-rank 1.
  if (rows < columns_out - 1)
    return 0;
  out->columns = columns_out;
  // Otherwise extract the rest of the rows.
  while (--rows) {
    supp1 = support->supp[0];
    supp2 = support->supp[2];
    count = count1;
    while (count--){
      *out_coeff = *in_coeff++;
      out_coeff += (supp1 & 0x1);
      temp = supp1 >> 1;
      supp1 = supp2;
      supp2 = temp;
    }
    supp1 = support->supp[1];
    supp2 = support->supp[3];
    count = count2;
    while (count--){
      *out_coeff = *in_coeff++;
      out_coeff += (supp1 & 0x1);
      temp = supp1 >> 1;
      supp1 = supp2;
      supp2 = temp;
    }
  }
  return 1;
}
