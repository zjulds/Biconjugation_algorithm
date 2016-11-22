#include <iostream>
#include <vector>
#include <math.h>
using namespace std;
typedef  struct _sort_type_{
		int pos;
		double value;
	}sort_type;
extern void drop(vector<double> &x,int nmax,double droptol,int j);
extern int nnz(vector<double> a);
extern void Heapfy(sort_type A[],int idx,int max) ;
extern vector<int> buildHeap(vector<double> &A,const int ll) ;