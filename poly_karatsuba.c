/*
 * karatsuba.c
 *
 *  Created on: Aug 22, 2019
 *      Author: shigurefox
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "macros.h"
#include "poly.h"
#include "poly_karatsuba.h"

void poly_karatsuba(Poly* r, Poly* p, Poly* q) {
	if (DEBUG) { printf("Entering function poly_karatsuba(): \n"); }
	Poly* result = poly_init(0);
	int res_degree = p->degree + q->degree;

	if ( p->degree < SCHOOLBOOK_CUTOFF || q->degree < SCHOOLBOOK_CUTOFF) {
		if (DEBUG) { printf("degree smaller than cutoff, calling poly_mul()...\n"); }
		poly_mul(result, p, q);
	} else {
		/*
		 * p = p0 + p1*x^m
		 * q = q0 + q1*x^m
	 	 * p*q = p0*q0 + [(p0+p1)*(q0+q1)-p0*q0-p1*q1]*x^m + p1*q1*x^2m
	 	 *     = r0 + [s0*s1 - r0 - r2]*x^m + r2*x^2m
	 	 *     = r0 + r1*x^m + r2*x^2m
	 	 */
		int greater_degree = MAX(p->degree, q->degree);
		int m = (greater_degree + 1)/2; // ceil
		Poly *p0 = poly_init(0), *p1 = poly_init(0);
		Poly *q0 = poly_init(0), *q1 = poly_init(0);
		Poly *s0 = poly_init(0), *s1 = poly_init(0), *s2 = poly_init(0), *s3 = poly_init(0);
		Poly *s4 = poly_init(0), *s5 = poly_init(0), *s6 = poly_init(0);
		Poly *r0 = poly_init(0), *r1 = poly_init(0), *r2 = poly_init(0);
		poly_slice(p0, p, 0, m - 1);
		poly_slice(p1, p, m, greater_degree - m);
		poly_slice(q0, q, 0, m - 1);
		poly_slice(q1, q, m, greater_degree - m);

		if (DEBUG) {
			printf("\tBefore slicing, p is: \n"); poly_print(p);
			printf("\tBefore slicing, q is: \n"); poly_print(q);
			printf("\tAfter slicing, p0 is: \n"); poly_print(p0);
			printf("\tAfter slicing, p1 is: \n"); poly_print(p1);
			printf("\tAfter slicing, q0 is: \n"); poly_print(q0);
			printf("\tAfter slicing, q1 is: \n"); poly_print(q1);
		}

		poly_mul(r0, p0, q0);        // r0 = p0 * q0
		poly_mul(r2, p1, q1);        // r2 = p1 * q1
		poly_add(s0, p0, p1);
		poly_add(s1, q0, q1);
		poly_mul(s2, s0, s1);
		poly_sub(s3, s2, r0);
		poly_sub(r1, s3, r2);        // r1 = s0 * s1 - r0 - r2
		poly_rshift(s4, r1, m);
		poly_rshift(s5, r2, 2*m);
		poly_add(s6, r0, s4);
		poly_add(result, s5, s6);


		result->degree = res_degree;

		if (DEBUG) {
			printf("Alignment check for r0 to r2:\n");
			printf("r0: \n"); poly_print(r0);
			printf("r1: \n"); poly_print(s4);
			printf("r2: \n"); poly_print(s5);
			printf("r0+r1*x: \n"); poly_print(s6);
			printf("Result r = p * q is: \n"); poly_print(result);

			printf("Exiting poly_karatsuba():\n");
		}

		// free memory
		poly_free(p0);poly_free(p1);poly_free(q0);poly_free(q1);
		poly_free(s0);poly_free(s1);poly_free(s2);poly_free(s3);
		poly_free(s4);poly_free(s5);poly_free(s6);
		poly_free(r0);poly_free(r1);poly_free(r2);
	}
	poly_assign(r, result);
	poly_free(result);

}

void poly_karatsuba_mod(Poly* r, Poly* p, Poly* q) {
	if (DEBUG) { printf("Entering function poly_karatsuba_mod(): \n"); }
	Poly* result = poly_init(0);
	int res_degree = p->degree + q->degree;

	if ( p->degree < SCHOOLBOOK_CUTOFF || q->degree < SCHOOLBOOK_CUTOFF) {
		if (DEBUG) { printf("degree smaller than cutoff, calling poly_mul_mod()...\n"); }
		poly_mul_mod(result, p, q);
	} else {
		/*
		 * p = p0 + p1*x^m
		 * q = q0 + q1*x^m
	 	 * p*q = p0*q0 + [(p0+p1)*(q0+q1)-p0*q0-p1*q1]*x^m + p1*q1*x^2m
	 	 *     = r0 + [s0*s1 - r0 - r2]*x^m + r2*x^2m
	 	 *     = r0 + r1*x^m + r2*x^2m
	 	 */
		int greater_degree = MAX(p->degree, q->degree);
		int m = (greater_degree + 1)/2; // ceil
		Poly *p0 = poly_init(0), *p1 = poly_init(0);
		Poly *q0 = poly_init(0), *q1 = poly_init(0);
		Poly *s0 = poly_init(0), *s1 = poly_init(0), *s2 = poly_init(0), *s3 = poly_init(0);
		Poly *s4 = poly_init(0), *s5 = poly_init(0), *s6 = poly_init(0);
		Poly *r0 = poly_init(0), *r1 = poly_init(0), *r2 = poly_init(0);
		poly_slice(p0, p, 0, m - 1);
		poly_slice(p1, p, m, greater_degree - m);
		poly_slice(q0, q, 0, m - 1);
		poly_slice(q1, q, m, greater_degree - m);

		if (DEBUG) {
			printf("\tBefore slicing, p is: \n"); poly_print(p);
			printf("\tBefore slicing, q is: \n"); poly_print(q);
			printf("\tAfter slicing, p0 is: \n"); poly_print(p0);
			printf("\tAfter slicing, p1 is: \n"); poly_print(p1);
			printf("\tAfter slicing, q0 is: \n"); poly_print(q0);
			printf("\tAfter slicing, q1 is: \n"); poly_print(q1);
		}

		poly_mul_mod(r0, p0, q0);        // r0 = p0 * q0
		poly_mul_mod(r2, p1, q1);        // r2 = p1 * q1
		poly_add_mod(s0, p0, p1);
		poly_add_mod(s1, q0, q1);
		poly_mul_mod(s2, s0, s1);
		poly_sub_mod(s3, s2, r0);
		poly_sub_mod(r1, s3, r2);        // r1 = s0 * s1 - r0 - r2
		poly_rshift(s4, r1, m);
		poly_rshift(s5, r2, 2*m);
		poly_add_mod(s6, r0, s4);
		poly_add_mod(result, s5, s6);


		result->degree = res_degree;

		if (DEBUG) {
			printf("Alignment check for r0 to r2:\n");
			printf("r0: \n"); poly_print(r0);
			printf("r1: \n"); poly_print(s4);
			printf("r2: \n"); poly_print(s5);
			printf("r0+r1*x: \n"); poly_print(s6);
			printf("Result r = p * q is: \n"); poly_print(result);

			printf("Exiting poly_karatsuba_mod():\n");
		}

		// free memory
		poly_free(p0);poly_free(p1);poly_free(q0);poly_free(q1);
		poly_free(s0);poly_free(s1);poly_free(s2);poly_free(s3);
		poly_free(s4);poly_free(s5);poly_free(s6);
		poly_free(r0);poly_free(r1);poly_free(r2);
	}
	poly_assign(r, result);
	poly_free(result);

}
