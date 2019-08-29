/*
 ============================================================================
 Name        : example.c
 Author      : Shigurefox
 Copyright   : Copyright Shigurefox 2019
 Description : Karatsuba implementation in C
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "macros.h"
#include "poly.h"
#include "poly_karatsuba.h"

int main(void) {
	srand(time(NULL));

	clock_t t0, t1, t2;
	Poly *x = poly_gen();
	Poly *y = poly_gen();
	Poly *result1 = poly_init(42), *result2 = poly_init(42), *result3 = poly_init(42);

	t0 = clock();
	poly_mul(result1, x, y);
	t1 = clock();
	poly_karatsuba(result2, x, y);
	t2 = clock();

	printf("\nTest case: ");
	poly_eval_print(x); printf(" * "); poly_eval_print(y);
	printf("\n\n");

	printf("Result using schoolbook is: \n");
	poly_eval_print(result1);
	printf("\nTime spent: %f ms\n\n", (double)(t1-t0)*1000/CLOCKS_PER_SEC);
	printf("Result using karatsuba is: \n");
	poly_eval_print(result2);
	printf("\nTime spent: %f ms\n\n", (double)(t2-t1)*1000/CLOCKS_PER_SEC);

	poly_free(x);poly_free(y);
	poly_free(result1);poly_free(result2);poly_free(result3);

	return EXIT_SUCCESS;
}
