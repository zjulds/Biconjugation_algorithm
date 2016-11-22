#ifndef __GENREARE_DATA_CPP
#define __GENREARE_DATA_CPP
#include <math.h>
#include "ILU.h"
#include "spars_ma.h"

void generate_data(sparse_block_ILU_CSR &A,matrix_ &x,matrix_ &b){
	const int row=20;
	const int cloumn=20;
	float data[row][cloumn];
	

	for(int i=0;i<row;i++){
		for(int j=0;j<cloumn;j++){
			if(i==j){
				data[i][j]=1;
				cout<<data[i][j]<<" ";
			}else{
				float temp = abs(i-j);
				data[i][j]=sin(temp)/temp;
				cout<<data[i][j]<<" ";
				if(j==cloumn-1)
					cout<<endl;
			}
		}
	}	

	matrix_ *AA ;
	int *JA;
	int *IA;
	int *uptr;

	AA=A.AA;
	JA=A.JA;
	IA=A.IA;
	uptr=A.uptr;
	A.Rows=row;
	A.Cols=cloumn;
	AA = new matrix_[row*cloumn];
	JA = new int[row*cloumn];
	IA = new int[row];
	uptr = new int[row];


	int k=0,l=0;
	bool flag;
	for(int i=0;i<row;i++){
		flag = true;
		for(int j=0;j<cloumn;j++){
			if(data[i][j]!=0){
				AA[k]=matrix_(1,1);
				AA[k].array[0] = data[i][j];
				JA[k]= j;
				if(flag){
					IA[l++]=k;
					flag=false;
				}
				if(i==j){
					uptr[l-1]=k;
				}
				k++;
			}
		}
	}
	ILU_decomp(A,1,1);

	matrix_ X,y;
	X=matrix_(row,row);
	x=matrix_(row,1);
	y=matrix_(row,1);
	b=matrix_(row,1);

	matrix_ *xx,*yy,*bb;
	xx=new matrix_[row];
	yy=new matrix_[row];
	bb=new matrix_[row];
	for(int i=0;i<row;i++)
	{
		xx[i]=matrix_(1,1);
		yy[i]=matrix_(1,1);
		bb[i]=matrix_(1,1);
	}
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<row;j++)
		{
			if(i==j)
			b[j]=1;
			else
			b[j]=0;
			bb[j][0]=b[j];
		}
		forward_substitution(A,yy,bb,1);
		backward_substitution(A,xx,yy,1);
		for(int j=0;j<row;j++)
		X[i*row+j]=xx[j][0];
	}
	foutput_ma(X,"X.txt");
		



	
}

#endif