#ifndef ENUM_H
#define ENUM_H

#include "lattice.h"

// using namespace NTL;
// extern void set_approximate_extreme_pruning(double* pf, quad_float* c, int bs, double prob, double delta);

/*********************
Enumeration 
**********************/
inline bool Lattice::ENUM(VectorXi &result, int g, int h)
{
	int i, j, k, k_1, last_nonzero, n = basis.cols(), flag;
	static MatrixXd sigma; sigma.setZero(n+1, n); 
	static VectorXi rr; rr.setZero(n+2); 
	static VectorXi v; v.setZero(n); 
	static VectorXd c; c.setZero(n); 
	static VectorXd w; w.setZero(n); 
	static VectorXd rho; rho.setZero(n+1);
	double dtmp, R;

#if 0
	/* Step 0 prepare extreme pruning bounding function */
	int n_sub = h-g+1;
	double* pruning_func = new double[n_sub];
	double* normal_pruning_func = new double[n_sub+1];
	quad_float* c_tmp = new quad_float[n_sub+1];
	
	if (n_sub <= 45) {
		/* Full Enumeration Setting */
		for(i=0; i<n_sub; i++) { 
			normal_pruning_func[i] = 0.99;
			pruning_func[i] = normal_pruning_func[i] * static_cast<double>(B(g));
		}
	} else {
		/* Extreme Pruning Enumeration Setting */
		double prob; // success probability of extreme pruning to find a short vector in the pruned tree 
							 //(within a bound of alpha*GH, we set alpha = 1.05 default here).
		prob = 2.0 / pow(1.05, n_sub);

		for (i = 0; i<n_sub; ++i) {
			c_tmp[i] = static_cast<double>(B(i+g));
		}
		double clim = B(g); // or set clim = alpha*GH(B_sub)
#if 0
		long double tmpp = 1.0; 
		for (i=g; i<=h; ++i) { tmpp *= B(i); }
		tmpp = pow(tmpp, 1.0/n_sub);
		tmpp *= n_sub/(2*M_PI*expl(1));	
		clim = 1.1*tmpp; if (B(g) < clim) { clim = B(g); }
#endif		
		//cout << "debug::: start generating pf..." << endl;
		set_approximate_extreme_pruning(normal_pruning_func, c_tmp, n_sub, prob, 0);
		for (i=0;i<n_sub;i++) {
			pruning_func[i] = normal_pruning_func[i] * clim;
		}
	}
#endif 	
	R = 0.999*static_cast<double>(B(g)); 
	for (i=0; i<=n+1; ++i) { rr(i) = i-1; }
	v(g) = 1; last_nonzero = g;
	
	k=g; flag = 0;
	while (1) {
		// cout << "k =" << k << endl; 
		dtmp = static_cast<double>(v(k)) - c(k); dtmp *= dtmp;  
		rho(k) = rho(k+1) + dtmp*static_cast<double>(B(k));
		// if (rho(k) < pruning_func[k-g]) {
		if (rho(k) <= R) {
			if (k == g) {
				flag += 1; 
				result = v;
				/* pruning bound update */
				// for (i=0; i<n_sub; ++i) { pruning_func[i] = 0.99*normal_pruning_func[i]*rho(g); }
				R = 0.999*rho(g);
			} else {
				--k;
				rr(k) = max(rr(k), rr(k+1));
				for (i=rr(k+1); i>=k+1; --i) {
					sigma(i, k) = sigma(i+1, k) + mu(i, k)*static_cast<double>(v(i)); 
				}
				c(k) = -sigma(k+1, k);
				v(k) = round(c(k));
				w(k) = 1;
			}
		} else {
			++k;
			if (k == h+1) {
				if (flag == 0) {
					// cout << "No solution" << endl;
					// delete[] pruning_func, normal_pruning_func, c_tmp; 				
					return false;
				} else {
					// cout << "Solution found" << endl;
					// delete[] pruning_func, normal_pruning_func, c_tmp;
					return true;
				}
			}
			rr(k) = k;
			if (k >= last_nonzero) {
				last_nonzero = k; 
				v(k) += 1;
			} else {
				if (v(k) > c(k)) {
					v(k) -= w(k); 
				} else {
					v(k) += w(k); 
				}
				++w(k);
			}
		}
	}
}

#endif // ENUM_H
