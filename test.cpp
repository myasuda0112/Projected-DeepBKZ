#include <iostream>
#include <fstream>
#include <stdio.h> 
#include <string>
#include <iomanip>
// #include <Eigen/Dense>

#include "lattice.h"
#include "DeepLLL.h"
#include "Enum.h"
#include "DeepBKZ.h"

int main(int argc, char *argv[1]){
	int block, i, j, k, l, n = 50, d = atoi(argv[1]), q = 2503; 
	FLOAT alpha = 0.010L;
	FLOAT sigma = alpha*q, bound;
	string line, token;
	stringstream ss, ss1;
	bool res;
	Lattice L; L.basis.setZero(d+1, d+1);
	VectorXi b; b.setZero(d);
	MatrixXi A; A.setZero(d, n);
	mat_ZZ C, D; C.SetDims(d+n, d), D.SetDims(d+1, d+1);
	ZZ det;
	
	// setNbThreads(1);
	cout << "#Threads = " << nbThreads() << endl;
	
	/* Read files */
	ifstream File1("n50a0010/target.txt");
	if (!File1) {
		cout << "File1 open error" << endl; 
		return -1; 
	}
	getline(File1, line);
	ss << line; i = 0;
	for (i=0; i<d; ++i) {
		getline(ss, token, ',');
		b(i) = atoi(token.c_str());
	}
	
	ifstream File2("n50a0010/matrix.txt");
	if (!File2) {
		cout << "File2 open error" << endl; 
		return -1; 
	}
	for (i=0; i<d; ++i) {
		getline(File2, line);
		ss1 << line;
		for (j=0; j<n; ++j) {
			getline(ss1, token, ',');
			A(i, j) = atoi(token.c_str());
		}
		ss1.clear(stringstream::goodbit); 
	}
	
	/* Construction of matrices */
	for (i=1; i<=d; i++) {C(i, i) = q;}
	for (i=1; i<=n; i++) {
		for (j=1; j<=d; j++) {C(d+i, j) = A(j-1, i-1);}
	}
	
	/* Remove the linear dependency of C by LLL */
	LLL(det, C, 0);
	for (i=1; i<=d; i++) {
		for (j=1; j<=d; j++) { D(i, j) = C(i+n, j); }
	}
	for (j=1; j<=d; ++j) { D(d+1, j) = b(j-1); }
	D(d+1, d+1) = 1;
	
	/* BKZ-20 (NTL-library) */
	BKZ_FP(D, 0.99, 20, 0, 0, 0);
	
	/* Switch to Eigen-library */
	for (i=0; i<d+1; ++i) {
		for (j=0; j<d+1; ++j) {
			L.basis(i, j) = to_int(D(i+1, j+1)); 
		}
	}
	
#if 1
	res = L.MGSO(d+1);
	long double vol, GH, tmp, tmp1;
	int h; 
	
	for (j=0; j<=40; ++j) {
		tmp1 = 1.0;
		h = 0; 
		for (l=0; l<=d; ++l) {
			bound = sqrtl(1.0*(d-l+1))*sigma;
			vol = 1.0; 
			for (i=l; i<=d; ++i) { vol *= L.B(i); }
			GH = pow(vol, 1.0/(2.0*(d-l+1)));
			GH *= sqrtl((d-l+1)/(2*3.1416*2.7182));
			
			tmp = GH/(0.3*bound);
			tmp = pow(tmp, 1.0/(d-l+1));
			// cout << "tmp = " << tmp << endl; 
			if (tmp >= tmp1) {
				tmp1 = tmp; 
				h = l;
			}
			
			// if (bound >= GH) {
			if (sqrtl(L.B(l)) < 2.0*bound) {
				// cout << "position = " << l << endl;
				// cout << "beta = " << d-l+1 << endl; 
				break;
			}
		}
		cout << "h = " << h << ", l = " << l << ", d = " << d << endl;
		cout << "tmp = " << tmp << endl;
		
		block = 20 + j;
		res = L.DeepBKZ(h, d, block, 0.99L, block, sigma, l);
		if (res != true) {
			cout << "h = " << h << ", l = " << l << ", d = " << d << endl; 
			cout << "tmp = " << tmp << endl;
			cout << "block = " << block << endl; 
			return 0; 
		}
	}
#endif

	return 0; 
}