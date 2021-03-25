#ifndef DBKZ_H
#define DBKZ_H

#include "lattice.h"
#include "DeepLLL.h"
#include "Enum.h"

/***********************************
Babai's lifting
************************************/
inline void Lattice::Babai(VectorINT &sv, VectorXi v, int j, int m)
{
	int i, n = basis.cols(), s, c;
	VectorFLOAT nu; nu.setZero(n); 
	
	for (i=0; i<=m; ++i) {
		nu(i) = static_cast<FLOAT>(v(i));
		for (s=i+1; s<=m; ++s) {
			nu(i) += static_cast<FLOAT>(v(s))*mu(s, i); 
		}
	}
	for (i=j-1; i>=0; --i) {
		if (fabsl(nu(i)) > ETA) {
			c = roundl(nu(i));
			sv -= c*basis.row(i);
			nu -= static_cast<FLOAT>(c)*mu.row(i);
		}
	} 
}

#if 0
/********************************
Additional size reduction
*********************************/
inline void Lattice::Add_SizeReduce(VectorXi &v, int j, int k)
{
	int i, h, f, l, n = basis.cols();
	FLOAT tmp, tmp1, tmp2;
	
	for (f=k+1; f<=n-1; ++f) {
		tmp = 0.0; 
		for (l=j; l<=f; ++l) {
			tmp += mu(f, l)*mu(f, l)*B(l);
		}
		tmp2 = 0.0; 
		for (i=j; i<f; ++i) {
			if (v(i) != 0) {
				tmp1 = 0,0; 
				for (h=j; h<=i; ++h) {
					tmp1 += mu(i, h)*mu(f, h)*B(h); 
				}
				tmp1 *= v(i);
				tmp2 += tmp1; 
			}
		}
		if (fabs(tmp2) > ETA*fabs(tmp)) {
			// v(f) = roundl(-tmp2/tmp);
			v(f) = to_int(floor(-tmp2/tmp + 0.50)); 
			cout << "j = " << j << ", k = " << k; 
			cout << ", f = " << f << ", v(f) = " << v(f) << endl;
		}
	}
}
#endif

/****************************
DeepBKZ reduction algorithm
****************************/
inline bool Lattice::DeepBKZ(int start, int end, int block, FLOAT alpha, int gamma, FLOAT sigma, int ll)
{
	int i, j, k, h, hh, l, m, n = basis.cols(), z, tour, N; 
	bool res;
	FLOAT current, current0, bound, tmp, tmp1, tmp2, tmp3, rho, rho1, GH, vol; 
	VectorXi v; v.setZero(n);
	VectorINT vv;
	mat_ZZ Basis;
	
	/* Compute Gram-Schmidt information  */
	res = MGSO(n);
	if (res != true) {
		cout << "GSO error in DeepBKZ" << endl;
		return false; 
	}
	cout << "\nDeepBKZ-" << block << endl; 
	
	/* DeepLLL reduction */
	res = DeepLLL(start, end, alpha, gamma, start+1);
	if (res != true) {
		cout << "DeepLLL error in DeepBKZ" << endl; 
		return false;
	}
	
	if (basis.row(start).norm() < 1.3*sigma*sqrt(n)) {
		cout << basis.row(start) << endl; 
		cout << "Solution found at l = " << start << endl;
		return false; 
	}
	
	current = B(start);
	cout << "Projected Norm = " << sqrt(current) << endl;
	
	m = end-start+1;
	tmp = 0.0; tmp1 = 0.0; 
	for (i=start; i<=end; ++i) {
		tmp += (i-start+1)*log(B(i)); 
		tmp1 += log(B(i)); 
	}
	tmp2 = 0.50*(m+1); 
	tmp3 = m*(m*m-1)/12.0;
	rho = (tmp - tmp1*tmp2)/tmp3; 
	cout << "rho = " << rho << endl; 
	
	// z = start-1; k = start-1;
	z = ll-1; k = ll-1; 
	tour = 0; N = 0;
	while (z < end) {
		++k;
		l = min(k+block-1, end);
		h = min(l+1, end);
		if (k==end) {
			// k = start;
			k = ll; 
			l = k+block-1; h = l+1;
			
			m = end-start+1; 
			tmp = 0.0; tmp1 = 0.0; 
			for (i=start; i<=end; ++i) {
				tmp += (i-start+1)*log(B(i)); 
				tmp1 += log(B(i));
			}
			tmp2 = 0.50*(m+1);
			tmp3 = m*(m*m-1)/12.0; 
			rho1 = (tmp - tmp1*tmp2)/tmp3;
			if (rho1 >= rho) { N = 0; }
			else {
				++N;
				if (N >= 5 && block >= 45) {
					cout << "Early termination" << endl; 
					return true; 
				}
			}
			rho = rho1;
			++tour;
			if (tour % 100 == 0) {			
				cout << "#Tour = " << tour << ", rho = " << rho << endl;
			}
		}
		/* Enumeration */
		// cout << "k = " << k << ", l = " << l << endl; 
		res = ENUM(v, k, l);
		if (res == false) {
			res = DeepLLL(start, h, alpha, gamma, h-1);
			if (res != true) {
				cout << "DeepLLL error in DeepBKZ" << endl; 
				return false;
			}		
			++z; 
		} else {
			// z = start-1;
			z = ll-1; 
			/* Insertion vector */
			tmp = B(k);
			// Add_SizeReduce(v, k, l); 
			m = 0; 
			for (i=n-1; i>=k; --i) {
				if (v(i) != 0) { m = i; break; }
			}
			vv.setZero(n);
			for (i=k; i<=m; ++i) { vv += v(i)*basis.row(i); }
			if (abs(v(m))==1) {
				for (i=m; i>k; --i) {
					basis.row(i).swap(basis.row(i-1));
				}
				Babai(vv, v, k, m);
				basis.row(k) = vv;
				
				/* DeepLLL reduction */
				res = MGSO(n);
				if (res != true) {
					cout << "GSO error in DeepBKZ" << endl;
					return false;
				}
				if (B(k) >= tmp) {
					cout << "ENUM error in DeepBKZ" << endl; 
					return false;
				}
				res = DeepLLL(start, h, alpha, gamma, k);
				if (res != true) {
					cout << "DeepLLL error in DeepBKZ" << endl; 
					return false; 
				}
			} else {
#if 0
				Basis.SetDims(h+1, n);
				for (i=0; i<k; ++i) {
					for (j=1; j<=n; ++j) { Basis(i+1, j) = basis(i, j-1); }
				}
				for (j=1; j<=n; ++j) {Basis(k+1, j) = vv(j-1); }
				for (i=k+1; i<=h; ++i) {
					for (j=1; j<=n; ++j) { Basis(i+1, j) = basis(i-1, j-1); }
				}
				LLL_FP(Basis, alpha, 0, 0, 0);
				for (i=0; i<h; ++i) {
					for (j=0; j<n; ++j) { basis(i, j) = to_int(Basis(i+2, j+1)); }
				}
				
				res = MGSO(n);
				if (res != true) {
					cout << "GSO error in DeepBKZ" << endl;
					return false;
				}
				res = DeepLLL(start, h, alpha, gamma, k);
				if (res != true) {
					cout << "DeepLLL error in DeepBKZ" << endl; 
					return false; 
				}
#endif 				
			}
		}
		/* Size-reduction */
		for (i=start; i<=end; ++i) {
			Size_Reduce(i);
		}
		if (B(start) <= 0.99*current) {
			current = B(start); 
			cout << "Projected Norm = " << sqrt(current) << endl;
			cout << basis.row(start) << endl;
			if (basis.row(start).norm() < 1.3*sigma*sqrt(n)) {
				cout << "Solution found at l = " << start << endl;
				return false; 
			}
		}
	}
	return true; 
}

#if 0
/************************************
Extended DeepBKZ reduction algorithm
*************************************/
inline bool Lattice::ExDeepBKZ(int start, int block, FLOAT alpha, int gamma, FLOAT sigma)
{
	int b, i, j, k, h, l, m, n = basis.cols(), z, tour, N;
	bool res;
	FLOAT current, current0, bound, tmp, tmp1, tmp2, tmp3, rho, rho1; 
	VectorXi v; v.setZero(n);
	VectorINT vv;
	mat_ZZ Basis;
	
	/* Compute Gram-Schmidt information  */
	res = MGSO(n);
	if (res != true) {
		cout << "GSO error in DeepBKZ" << endl;
		return false; 
	}
	cout << "\nDeepBKZ-" << block << endl;
	current = B(start);
	cout << "Projected Norm = " << sqrt(current) << endl;
	
	/* DepBKZ reduction */
	b = round(block/2); 
	res = DeepBKZ(start, n-1, b, alpha, b, sigma);

	m = n-start;
	tmp = 0.0; tmp1 = 0.0; 
	for (i=start; i< n; ++i) {
		tmp += (i-start+1)*log(B(i)); 
		tmp1 += log(B(i)); 
	}
	tmp2 = 0.50*(m+1); 
	tmp3 = m*(m*m-1)/12.0; 
	rho = (tmp - tmp1*tmp2)/tmp3; 
	cout << "rho = " << rho << endl; 
	
	z = start-1; k = start-1; tour = 0; N = 0; 
	while (z < n-1) {
		++k; 
		l = min(k+block-1, n-1);
		h = min(l+1, n-1);
		if (k==n-1) {
			k = start; l = k+block-1; h = l+1;
			
			m = n-start; 
			tmp = 0.0; tmp1 = 0.0; 
			for (i=start; i< n; ++i) {
				tmp += (i-start+1)*log(B(i)); 
				tmp1 += log(B(i));
			}
			tmp2 = 0.50*(m+1);
			tmp3 = m*(m*m-1)/12.0; 
			rho1 = (tmp - tmp1*tmp2)/tmp3;
			if (rho1 >= rho) { N = 0; }
			else {
				++N;
				if (N >= 5 && block >= 45) {
					cout << "Early termination" << endl; 
					return true; 
				}
			}
			rho = rho1;
			++tour;
			if (tour % 100 == 0) {
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
						break;
					}
				}				
				
				
				
				
				
				cout << "#Tour = " << tour << ", rho = " << rho << endl; 
			}
		}
		/* Enumeration */
		// cout << "k = " << k << ", l = " << l << end
		res = ENUM(v, k, l);
		if (res == false) {
			// res = DeepBKZ(start, h, b, alpha, b, sigma);
			++z; 
		} else {
			z = start-1;
			/* Insertion vector */
			tmp = B(k);
			m = 0; 
			for (i=n-1; i>=k; --i) {
				if (v(i) != 0) { m = i; break; }
			}
			vv.setZero(n);
			for (i=k; i<=m; ++i) { vv += v(i)*basis.row(i); }
			if (abs(v(m))==1) {
				for (i=m; i>k; --i) {
					basis.row(i).swap(basis.row(i-1));
				}
				Babai(vv, v, k, m);
				basis.row(k) = vv;
				
				/* DeepLLL reduction */
				res = MGSO(n);
				if (res != true) {
					cout << "GSO error in DeepBKZ" << endl;
					return false;
				}
				if (B(k) >= tmp) {
					cout << "ENUM error in DeepBKZ" << endl; 
					return false;
				}
				res = DeepBKZ(start, h, b, alpha, b, sigma);
			} else {
				// cout << "point" << endl; 		
			}
		}
		/* Size-reduction */
		for (i=start; i<n; ++i) {
			Size_Reduce(i);
		}
		if (B(start) <= 0.99*current) {
			current = B(start); 
			cout << "Projected Norm = " << sqrt(current) << endl;
			cout << basis.row(start) << endl;
			if (basis.row(start).norm() < 1.3*sigma*sqrt(n)) {
				cout << "Solution found at l = " << start << endl; 
				return false; 
			}
		}	
	}
	return true; 
}
#endif 

#endif // DBKZ_H
