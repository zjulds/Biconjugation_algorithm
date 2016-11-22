#include "drop.h"

void drop(vector<double> &x,int nmax,double droptol,int j){
	/**
	author	:liadashuang 2016-11-3
	x       : vector where the drop will be applied
	nmax    : max nnz entries in the resulting s
	droptol : drop tolerance
	 j      : vector that specifies components that cannot be
              dropped (eg. diagonal)
	*/
	
	//initial the vector
			vector<double> y,result(x.size(),0);
			vector<int> order;
			order.reserve(x.size());
			y.reserve(x.size());
	//
			int num0=nnz(x);
			if(num0<=nmax ){
				return;
			}

			for(int i=0;i<x.size();i++){
				y.push_back(abs(x[i]));
			}
			order=buildHeap(y,y.size());
			//buildHeap(x,9);

			//keep just the nmax biggest entries (the others are dropped)
			int n;
			if(nmax<x.size()){
				n=nmax;			
			}else{
				n=x.size();
			}
			for(int i=0;i<n;i++){
				if(y[i]<=droptol){
					y[i]=0;
				}else{
					y[i]=x[order[i]];
				}		
			}
			for(int i=n;i<y.size();i++){
				y[i]=0;
			}
			//set the new vector result
			for(int i=0;i<y.size();i++){
				result[order[i]]=y[i];
			}

			//restore the original value of the j-th entry
			result[j]=x[j];
			for(int i=0;i<x.size();i++){
				x[i]=result[i];
			}

	}

//this function is used to count the number of nonzeros in vector
int nnz(vector<double> a){
	int count=0;
	int size=a.size();
	for(int i=0;i<size;i++){
		if(a[i]!=0)
			count++;
	}
	return count;
}
//sort function
//descend ! this code is modified according to increase algorithm


void Heapfy(sort_type A[],int idx,int max)      //建立最大堆  
{  
    int left=idx*2+1;  
    int right=left+1;  
  
    int largest=idx;  
  
	//if(left<max&&A[left].value>A[idx].value){largest=left;}  
    if(left<max&&A[left].value<A[idx].value){largest=left;}
   // if(right<max&&A[largest].value<A[right].value){largest=right;}  
    if(right<max&&A[largest].value>A[right].value){largest=right;}  
    if(largest!=idx)  
    {  
        sort_type temp=A[largest];   //较大的节点值将交换到其所在节点的父节点  
        A[largest]=A[idx];  
        A[idx]=temp;  
  
        Heapfy(A,largest,max); //递归遍历  
  
    }  
}  
 // 
vector<int> buildHeap(vector<double> &A,const int ll)  
{  
    int len=ll;  
	vector<int> order;
	order.reserve(ll);
	sort_type *p=new sort_type[ll];
	for(int i=0;i<A.size();i++){
		p[i].value=A[i];
		p[i].pos=i;
	};
    for(int i=len/2-1;i>=0;--i)  
    {  
        Heapfy(p,i,len);     //建立最大堆，将堆中最大的值交换到根节点  
    }  
  //////////////////////////////////////maybe some problems////////////////////
    for(int i=len-1;i>=1;--i)  
    {  
      /*  int temp=A[0];   //将当前堆的根节点交换到堆尾的指定位置  
        A[0]=A[i];  
        A[i]=temp;  
	  */
		sort_type temp=p[0];
		p[0]=p[i];
		p[i]=temp;
        Heapfy(p,0,i);  //建立下一次的最大堆  
    }  
	for(int i=0;i<A.size();i++){
		A[i]=p[i].value;
		order.push_back(p[i].pos);
	};
	return order;
}  

