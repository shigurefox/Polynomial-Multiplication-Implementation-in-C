/*
 * poly.h
 *
 *  Created on: Aug 27, 2019
 *      Author: shigurefox
 */

#ifndef POLY_H_
#define POLY_H_


typedef struct Poly {
	int degree;
	int* coefficients;
	int allocated_num_coeffs;
} Poly;

// ========================================
//   Helper functions
// ========================================

Poly* poly_init(int value);
Poly* poly_gen();
void poly_free(Poly*);
void poly_assign(Poly*, const Poly*);
void poly_read(Poly*, const char*);
void poly_assert_degree(Poly*, const int);
void poly_slice(Poly* , Poly*, int, int);
void poly_print(const Poly*);
void poly_rshift(Poly*, Poly*, int);
void poly_lshift(Poly*, Poly*, int);
int poly_compare(const Poly*, const Poly*);

// ========================================
//   Arithmetic operations
//   results are stored in the first argument
// ========================================

void poly_add(Poly*, Poly*, Poly*);
void poly_add_mod(Poly*, Poly*, Poly*);
void poly_add_int(Poly*, Poly*, const int);
void poly_add_mod_int(Poly*, Poly*, const int);

void poly_sub(Poly*, Poly*, Poly*);
void poly_sub_mod(Poly*, Poly*, Poly*);
void poly_sub_int(Poly*, Poly*, const int);
void poly_sub_mod_int(Poly*, Poly*, const int);

void poly_mul(Poly*, Poly*, Poly*);
void poly_mul_mod(Poly*, Poly*, Poly*);
void poly_mul_int(Poly*, Poly*, const int);
void poly_mul_mod_int(Poly*, Poly*, const int);

void poly_div_mod(Poly*, Poly*, const int);
void poly_div_mod_int(Poly*, Poly*, const int);

// ========================================
//   Miscellaneous
// ========================================

void poly_eval_print(Poly* p);
int poly_eval_mod(Poly*, int);

#endif /* POLY_H_ */
