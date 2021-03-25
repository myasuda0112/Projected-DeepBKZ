#ifndef DEEP_H
#define DEEP_H

const FLOAT ETA= 0.51L; 

/**************************************
size-reduce at index k
**************************************/
inline void Lattice::Size_Reduce(int k)
{
	int i;
	INT q;
	for (i=k-1; i>=0; --i) {
		if (fabsl(mu(k, i)) >= ETA ) {
			q = roundl(mu(k, i));
			basis.row(k) -= q*basis.row(i);
			mu.row(k) -= static_cast<FLOAT>(q)*mu.row(i);
		}
	}
}

/*****************************************
DeepLLL basis reduction
*****************************************/
inline bool Lattice::DeepLLL(int start, int end, FLOAT alpha, int gamma, int stage)
{
	int i, j, k, l, s, t, flag, n = basis.cols();
	int ii, s_1, t_1;  
	bool res; 
	FLOAT tmp0, tmp, tmp1;
	static VectorFLOAT D, BB, DD, tmpld, mut; 
	D.setZero(n); BB.setZero(n); DD.setZero(n); 
	tmpld.setZero(n); mut.setZero(n);
		
	k = stage;
	if (k <= start) {k = start+1;}
	while (k<=end) {
		// cout << "k = " << k << endl;
		Size_Reduce(k);
		tmp = basis.row(k).dot(basis.row(k));
		
		i = 0; flag = 0; 
		while (i<k) {
			/* deep exchange condition */
			if (i>=start && tmp < alpha*B(i)) {
				if ( (i-start+1) <= gamma || (k-i+1) <= gamma) {
					tmp0 = B(i);
					flag = 1;
				}
			}
			if (flag == 0) {
				tmp -= (mu(k, i)*mu(k, i))*B(i); 
				++i; 
			} else {
				/* deep insertion */
				for (t = k; t != i; --t) {
					basis.row(t).swap(basis.row(t-1)); 
				}					
				/* Update of Gram-Schmidt information */
				D(k) = B(k);
				BB(k) = B(k); 
				DD(k) = 1.0/D(k); 
				ii = i-1;
				for (l=k-1; l != ii; --l) {
					BB(l) = mu(k, l)*B(l);
					D(l) = D(l+1) + mu(k, l)*BB(l); 
					DD(l) = 1.0/D(l); 
					B(l+1) = D(l+1)*B(l)*DD(l); 
				}
				B(i) = D(i);
				
				/* for debug */
				if (B(i) >= tmp0) {
					cout << "Insertion error at " << i << endl;
					MGSO(n); 
					return true; 
					// return false; 
				}				
				for (t=i; t<=k; ++t) {
					if (B(t) <= DBL_EPSILON) {
						cout << "GSO update error in DLLL" << endl;
						MGSO(n); 
						return true; 
						// return false; 
					}
				}
				
				t_1 = k-1; 
				tmp1 = mu(k, t_1)*DD(k); 
				for (s =k+1; s != n; ++s) {
					tmpld(s) = BB(k)*mu(s, k); 
					mu(s, k) = mu(s, t_1) - tmp1*tmpld(s); 
				}
				for (t = k-1; t != i; --t) {
					t_1 = t-1; 
					tmp1 = mu(k, t_1)*DD(t);
					for (s = k+1; s != n; ++s) {
						tmpld(s) += BB(t)*mu(s, t);
						mu(s, t) = mu(s, t_1) - tmp1*tmpld(s); 
					}
					for (s=k; s != t+1; --s) {
						s_1 = s-1;
						tmpld(s) += BB(t)*mu(s_1, t);
						mu(s, t) = mu(s_1, t_1) - tmp1*tmpld(s); 
					}
					s_1 = t+1;
					tmpld(s_1) = BB(t); 
					mu(s_1, t) = mu(t, t-1) - tmp1*tmpld(s_1); 
				}
				for (s=k+1; s != n; ++s) {
					mu(s, i) = (tmpld(s) + BB(i)*mu(s, i))*DD(i); 
				}
				ii = i+1;
				for (s=k; s != ii; --s) {
					mu(s, i) = (tmpld(s) + BB(i)*mu(s-1, i))*DD(i); 
				}
				mu(ii, i) = BB(i)*DD(i);
				mut = mu.row(k); 
				for (s = k; s != i; --s) {
					for (t = 0; t!=i; ++t) {
						mu(s, t) = mu(s-1, t); 
					}
				}
				for (t = 0; t != i; ++t) { mu(i, t) = mut(t); }

				/* size-reduction */
				if (i != 0) { Size_Reduce(i); }			
				/* Update indices */
				k = i; 
				// k = max(i, start+1)-1;
				break; 
			}
		}
		++k; 
	}
	return true; 
}

#endif // DEEP_H
