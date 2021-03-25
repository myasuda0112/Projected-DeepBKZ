#ifndef LATTICE_H
#define LATTICE_H

#include <stdio.h>
#include <iostream>
#include <NTL/LLL.h>
#include <vector>
// #include <math.h>
#include <cmath>
#include <float.h>
#include <algorithm>
#include <Eigen/Dense>
// #include <Eigen/QR>
// #include <Eigen/LU>

using namespace std;
using namespace Eigen;
using namespace NTL; 

// typedef int INT; 
typedef long int INT;
// typedef double FLOAT; 
typedef long double FLOAT;
// typedef quad_float FLOAT; 

typedef Matrix<INT, 1, Dynamic> VectorINT;
typedef Matrix<FLOAT, 1, Dynamic> VectorFLOAT; 
typedef Matrix<INT, Dynamic, Dynamic> MatrixINT;
typedef Matrix<FLOAT, Dynamic, Dynamic> MatrixFLOAT; 

/**************************
Class for a lattice basis
**************************/
class Lattice {
public:
	MatrixINT basis; 	/* lattice basis */
	MatrixFLOAT r; 		/* GNR information */
	MatrixFLOAT mu; 	/* Gram-Schmidt coefficients */
	VectorFLOAT B; 		/* squared Gram-Schmidt lengths */
	
	/* lattice.h */
	bool GSO(int n);
	bool MGSO(int n);
	
	/* DeepLLL.h */
	void Size_Reduce(int k);
	bool DeepLLL(int start, int end, FLOAT alpha, int gamma, int stage);
	
	/* Enum.h */
	bool ENUM(VectorXi &result, int g, int h); 
	
	/* DeepBKZ.h */
	void Babai(VectorINT &sv, VectorXi v, int j, int r);
	void Add_SizeReduce(VectorXi &v, int j, int k); 
	bool DeepBKZ(int start, int end, int block, FLOAT alpha, int gamma, FLOAT sigma, int l);
	bool ExDeepBKZ(int start, int block, FLOAT alpha, int gamma, FLOAT sigma); 
};


/*********************************
Compute Gram-Schmidt information
**********************************/
inline bool Lattice::GSO(int n)
{
	int i, j, k; 
	MatrixFLOAT s(n, n); 
	r.setZero(n, n); mu.setZero(n, n); B.setZero(n); 
	
	for (i=0; i<n; ++i) {
		for (j=0; j<i; ++j) {
			r(i, j) = basis.row(i).dot(basis.row(j));		
			for (k=0; k<j; ++k) {
				r(i, j) -= mu(j, k)*r(i, k); 
			}
			mu(i, j) = r(i, j)/r(j, j); 
		}
		mu(i, i) = 1.0;
		
		s(i, 0) = basis.row(i).dot(basis.row(i));
		for (j=1; j<=i; ++j) {
			s(i, j) = s(i, j-1) - mu(i, j-1)*r(i, j-1); 
		}
		if (s(i, i) < DBL_EPSILON ) {
			cout << "GSO error at i = " << i << endl; 
			return false; 
		}
		B(i) = s(i, i); 
		r(i, i) = s(i, i);
	}
	return true; 
}

/****************************************
Modified Gram-Schmidt Orthogonalization 
****************************************/
inline bool Lattice::MGSO(int n)
{
	int i, j, k;
	MatrixFLOAT a, q; a.setZero(n, n), q.setZero(n, n);
	r.setZero(n, n); mu.setZero(n, n); B.setZero(n);

	for (k=0; k<n; ++k) {
		q.col(k) = basis.row(k).cast<FLOAT>(); 
		
		for (j=0; j<k; ++j) {
			r(j, k) = q.col(j).dot(q.col(k));
			q.col(k) -= r(j, k)*q.col(j);
		}
		r(k, k) = sqrtl(q.col(k).dot(q.col(k)));
		q.col(k) /= r(k, k);
	}
	
	for (i=0; i<n; ++i) {
		mu(i, i) = 1.0L;
		B(i) = r(i, i)*r(i, i);
		if (B(i) < LDBL_EPSILON) {
			cout << "Length error in MGSO" << endl; 
			return false; 
		}
		for (j=0; j<i; ++j) {
			mu(i, j) = r(j, i)/r(j, j); 
		}
	}
	a = r.transpose();
	r = a;
	
	return true;
}

#endif //LATTICE_H
