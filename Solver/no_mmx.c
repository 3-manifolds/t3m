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
  int result = 0;
  
  if (size != 0)
    asm volatile ("
        movd  %%eax, %%mm3   ;
        xorl  %%eax, %%eax   ;
        pxor  %%mm0, %%mm0   ;
        pxor  %%mm1, %%mm1   ;
        xorl  %%esi, %%esi   ;
        dotloop:
        movl  (%%ebx), %%eax ;
        imull (%%ecx)        ;
        addl  $0x4, %%ebx    ;
        addl  $0x4, %%ecx    ;
        movd  %%edx, %%mm2   ;
        movd  %%mm0, %%edx   ;
        addl  %%edx, %%eax   ;
        movd  %%eax, %%mm0   ;
        movd  %%mm1, %%eax   ;
        movd  %%mm2, %%edx   ;
        adcl  %%edx, %%eax   ;
        movd  %%eax, %%mm1   ;
        addl  $1, %%eax      ;
        shrl  $1, %%eax      ;
        or    %%eax, %%esi   ;
        decl  %%edi          ;
        jnz    dotloop       ;
        movd  %%mm3, %%eax   ;
        movd  %%mm0, (%%eax) ;
        movl  %%esi, %%eax   ;
        emms ;"
	:"=a"(result) 
	:"D"(size),"b"(x),"c"(y),"a"(dotprod)
	:"%esi","%edx"
	);
  return result;
}

// WARNING: This may write one int past the end of the array.  Allow extra
// space in your matrix.
// To avoid unpredictable branching, we just overwrite each of the columns that
// corresponds to a 0 in the support vector.  

int extract_matrix(matrix_t *in, int rows, support_t *support, matrix_t *out) {
  int columns;
  int *in_coeff = in->matrix;
  int *out_coeff = out->matrix;
  int columns_out = 0;

  out->rows = rows;
  columns = in->columns;
  asm volatile ("
        movq    (%%esi), %%mm0;
        movq    %%mm0, %%mm2
        movq    0x8(%%esi), %%mm1;
        movq    %%mm1, %%mm3;
        firstrow:
        movl    (%%ebx), %%esi;
        movl    %%esi, (%%ecx);
        movd    %%mm0, %%esi;
        addl    $4, %%ebx;
        psrlq   $1, %%mm0;
        movq    %%mm0, %%mm4;
        movq    %%mm1, %%mm0;
        movq    %%mm4, %%mm1;
        andl    $1, %%esi;
        addl    %%esi, %%eax;
        addl    %%esi, %%esi;
        addl    %%esi, %%esi;
        addl    %%esi, %%ecx;
        decl    %%edi;
        jnz     firstrow;"
	:"=a"(columns_out),"=S"(support),"=D"(columns),"=b"(in_coeff),"=c"(out_coeff)
	:"a"(columns_out),"S"(support),"D"(columns),"b"(in_coeff),"c"(out_coeff)
	);
  if (rows < columns_out - 1)
    return 0;
  out->columns = columns_out;
  columns = in->columns;
  if (--rows)
    asm volatile ("
        nextrow:
        movq    %%mm2, %%mm0;
        movq    %%mm3, %%mm1;
        movl    %%edi, %%edx;
        rowloop:
        movl    (%%ebx), %%esi;
        movl    %%esi, (%%ecx);
        movd    %%mm0, %%esi;
        addl    $4, %%ebx;
        psrlq   $1, %%mm0;
        movq    %%mm0, %%mm4;
        movq    %%mm1, %%mm0;
        movq    %%mm4, %%mm1;
        andl    $1, %%esi;
        addl    %%esi, %%esi;
        addl    %%esi, %%esi
        addl    %%esi, %%ecx;
        decl    %%edi;
        jnz     rowloop;
        movl    %%edx, %%edi;
        decl    %%eax;
        jnz     nextrow;
        emms ;"
	:
	:"a"(rows),"S"(support),"D"(columns),"b"(in_coeff),"c"(out_coeff)
        :"%edx"
	);
  return 1;
}
