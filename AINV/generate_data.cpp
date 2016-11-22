#ifndef __GENREARE_DATA_CPP
#define __GENREARE_DATA_CPP
#include <math.h>
#include "ILU.h"
#include "spars_ma.h"

void generate_data(){//sparse_block_ILU_CSR &A,matrix_ &x,matrix_ &b){
	sparse_block_ILU_CSR A;
	matrix_ x;matrix_ b;matrix_ A1;
	const int row=100;
	const int cloumn=100;
	float data[row][cloumn];
	A1=matrix_(row,cloumn);
	for(int i=0;i<row;i++){
		for(int j=0;j<cloumn;j++){
			if(i==j){
				data[i][j]=1;
				A1[i*cloumn+j]=data[i][j];
				cout<<data[i][j]<<" ";
			}else{
				float temp = abs(i-j);
				data[i][j]=sin(temp)/temp;
				A1[i*cloumn+j]=data[i][j];
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

	
	A.Rows=row;
	A.Cols=cloumn;
	AA = new matrix_[row*cloumn];
	JA = new int[row*cloumn];
	IA = new int[row];
	uptr = new int[row];


	int k=0,l=0,count0=0;
	bool flag;
	for(int i=0;i<row;i++){
		flag = true;
		for(int j=0;j<cloumn;j++){
			if(data[i][j]!=0){
				AA[k]=matrix_(1,1);
				AA[k].array[0] = data[i][j];
				JA[k]= j;
				count0++;
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
	IA[row]=count0;
	A.Terms=count0;
	A.AA=AA;
	A.JA=JA;
	A.IA=IA;
	A.uptr=uptr;
	ILU_decomp(A,1,1);

	matrix_ L,U;
	L=matrix_(row,row);
	U=matrix_(row,row);
	copy(L,0);
	copy(U,0);
	for(int i=0;i<A.Rows;i++)
	{
		for(int j0=A.IA[i];j0<A.IA[i+1];j0++)
		{
			int j=A.JA[j0];
			if(i>j)
			L[i*A.Rows+j]=A.AA[j0][0];
			else if(i==j)
			{
				L[i*A.Rows+j]=1.;
				U[i*A.Rows+j]=1./A.AA[j0][0];
				
			}
			else
			U[i*A.Rows+j]=A.AA[j0][0];
		}
	}
	matrix_ A2;
	A2=matrix_(row,row);
	multi(A2,L,U);
	sub(A2,A2,A1);
	double error;
	error=abs_matrix(A2)/abs_matrix(A1);
	cout<<"\nerror="<<error<<'\n';


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
	matrix_ I,I0;
	I=matrix_(row,cloumn);
	I0=matrix_(row,cloumn);
	multi(I,A1,X);
	foutput_ma(I,"I.txt");

	copy(I0,0);
	for(int i=0;i<row;i++)
	{
		I0[i*cloumn+i]=1.;
	}
	sub(I,I,I0);
	double error_I;
	error_I=abs_matrix(I)/abs_matrix(I0);
	cout<<"\nerror of I="<<error_I<<'\n';
		



	
}

#endif