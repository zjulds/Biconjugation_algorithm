#include "iostream"
#include <vector>
#include "drop.h"
#include "Binary_search.h"
using namespace std;


void main(){
	//set values
	const int row =3;
	const int column=3;
	double array[row][column]={{1,2,3},{2,2,1},{3,4,3}};
	double array_t[row][column];
	int t;
    for(int i=0; i<row; i++){
       for(int j=0; j<column; j++)
       {
		array_t[i][j]=array[j][i];
       }
	}
	//generate sparse matrix for calculating 
	vector<double> AA ;vector<int> JA;
	vector<int> IA;vector<int> uptr;
	AA.reserve(row*column);
	JA.reserve(row*column);
	IA.reserve(row);
	uptr.reserve(row);

	int k=0,l=0;
	bool flag;
	for(int i=0;i<row;i++){
		flag = true;
		for(int j=0;j<column;j++){
			if(array[i][j]!=0){
				AA.push_back(array[i][j]);
				JA.push_back(j);
				if(flag){
					IA.push_back(k);
					l++;
					flag=false;
				}
				if(i==j){
					uptr.push_back(k);
				}
				k++;
			}
		}
	}
	IA.push_back(AA.size());
	// print the value of AA,IA,JA
	cout<<"------AA---------"<<endl;
	 vector<double>::iterator it ;
    for(it=AA.begin(); it!=AA.end(); it++)
	{  cout<<*it<<" " ;}
	cout<<endl;
	cout<<"------JA---------"<<endl;
	vector<int>::iterator it1 ;
    for(it1=JA.begin(); it1!=JA.end(); it1++)
	{  cout<<*it1<<" " ;}
	cout<<endl;
	cout<<"------IA---------"<<endl;
	vector<int>::iterator it2 ;
    for(it2=IA.begin(); it2!=IA.end(); it2++)
	{  cout<<*it2<<" " ;}
	cout<<endl;cout<<endl;cout<<endl;cout<<endl;

	vector<double> AA_t ;vector<int> JA_t;
	vector<int> IA_t;vector<int> uptr_t;
	AA_t.reserve(row*column);
	JA_t.reserve(row*column);
	IA_t.reserve(row);
	uptr_t.reserve(row);
	int k_t=0,l_t=0;
	bool flag_t;
	for(int i=0;i<row;i++){
		flag_t = true;
		for(int j=0;j<column;j++){
			if(array[i][j]!=0){
				AA_t.push_back(array_t[i][j]);
				JA_t.push_back(j);
				if(flag_t){
					IA_t.push_back(k_t);
					l_t++;
					flag_t=false;
				}
				if(i==j){
					uptr_t.push_back(k_t);
				}
				k_t++;
			}
		}
	}
	IA_t.push_back(AA_t.size());
	cout<<"------AA_t---------"<<endl;
	 vector<double>::iterator it_t ;
    for(it_t=AA_t.begin(); it_t!=AA_t.end(); it_t++)
	{  cout<<*it_t<<" " ;}
	cout<<endl;
	cout<<"------JA_t---------"<<endl;
	vector<int>::iterator it1_t ;
    for(it1_t=JA_t.begin(); it1_t!=JA_t.end(); it1_t++)
	{  cout<<*it1_t<<" " ;}
	cout<<endl;
	cout<<"------IA_t---------"<<endl;
	vector<int>::iterator it2_t ;
    for(it2_t=IA_t.begin(); it2_t!=IA_t.end(); it2_t++)
	{  cout<<*it2_t<<" " ;}
	cout<<endl;cout<<endl;cout<<endl;cout<<endl;



	
	//initializing the vector Z and vector W
	vector<double> AA_Z ;vector<int> JA_Z;
	vector<int> IA_Z;vector<int> uptr_Z;
	AA_Z.reserve(row*column/2+1);
	JA_Z.reserve(row*column/2+1);
	IA_Z.reserve(row+1);
	uptr_Z.reserve(row);
	vector<double> AA_W ;vector<int> JA_W;
	vector<int> IA_W;vector<int> uptr_W;
	AA_W.reserve(row*column/2+1);
	JA_W.reserve(row*column/2+1);
	IA_W.reserve(row+1);
	uptr_W.reserve(row);
	
	for(int i=0;i<row;i++){
		AA_Z.push_back(1);
		AA_W.push_back(1);
	}
	for(int i=0;i<row;i++){
		JA_Z.push_back(i);
		JA_W.push_back(i);
	}
	for(int i=0;i<=row;i++){
		IA_Z.push_back(i);
		IA_W.push_back(i);
	}




	//////////////////////////////////////////////////////////////////////////////
	/////////////////////////////Biconjugation Algorithm//////////////////////////
	//////////////////////////////////////////////////////////////////////////////

	double p[column];
	double q[column];	
		for (int j=0;j<row;j++){
			p[j]=0;
			q[j]=0;
		}	
	int nmax=2;//this algorithm can be packaged into a function
	int number_non=0;//those arguments can be transfered into function
	double droptol=0.5;
	for (int i=0;i<row;i++){
		for (int j=0;j<row;j++){
			if(array[i][j]!=0)
				number_non++;
	}
	}
	int nmax_col=floor((double)nmax*(double)number_non/(double)row);

	for(int i=0;i<column;i++){
	/////////////////////////////////////////////////////////////////////////////
	///////////////////calculate the vector p and q /////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
		for(int j=i;j<row;j++){
			if(i==0){
				for(int k=0;k<row;k++){
				p[k]=AA[IA[k]];
				q[k]=AA_t[IA_t[k]];
				}
		}else{
		double temp=0,temp_t=0;
		vector<int> temp_j_z_column;
		temp_j_z_column.reserve(IA_Z[j+1]-IA_Z[j-1]);
		for(int l=IA_Z[j];l<IA_Z[j+1];l++){
			temp_j_z_column.push_back(JA_Z[l]);
		}
		for(int k=IA[i];k<IA[i+1];k++){
			int test=Binary_search(temp_j_z_column,JA[k]);
			if(test>=0){
				temp+=AA[k]*AA_Z[IA_Z[j]+test];
			}
		}
		p[j]=temp;


		vector<int> temp_j_w_column;
		temp_j_w_column.reserve(IA_W[j+1]-IA_W[j-1]);
		for(int l=IA_W[j];l<IA_W[j+1];l++){
			temp_j_w_column.push_back(JA_W[l]);
		}
		for(int k=IA_t[i];k<IA_t[i+1];k++){
			int test=Binary_search(temp_j_w_column,JA_t[k]);
			if(test>=0){
				temp_t+=AA_t[k]*AA_W[IA_W[j]+test];
			}
		}
		q[j]=temp_t;
		}
		}
	

	// I guess this is useless for if(i==row-1){**}
		if(i==column-1){
			break;
		} else{
			//update Z and W
			// Zj=Zj-p(j)/p(i)*Zi, j=(i+1):n

			//if(i>0){
		for(int j=i+1;j<column;j++){
					/////////////////////////////////////////////////////////////////////////////
				    //////this part is used to get zero and nonzero elements in column j ////////
					/////////////////////////////////////////////////////////////////////////////
				int count=0;
				vector<double> result_column1(j+1),result_column2(j+1);
				vector<int> column_nozero_Z,column_zero_Z;
				
				//to calculate the number of nonzero elements in j column;
				//then save their index of row into column_nozero
				int count_j_Z=IA_Z[j+1]-IA_Z[j]; 
				column_nozero_Z.reserve(count_j_Z);
				for(int i=IA_Z[j];i<IA_Z[j]+count_j_Z;i++){
					column_nozero_Z.push_back(JA_Z[i]);
					}
				//to calculate the number of zero elements in j column;
				//then save their index of row into column_nozero
				column_zero_Z.reserve(j+1);
				for(int i=0;i<=j;i++){
					column_zero_Z.push_back(i);
				}
				for(int i=0;i<count_j_Z;i++){
				vector<int>::iterator it = column_zero_Z.begin()+column_nozero_Z[i]-i;
					column_zero_Z.erase(it);
				}



				vector<int> column_nozero_W,column_zero_W;
				//result_column1.reserve(j);result_column2.reserve(j);	
				
				//to calculate the number of nonzero elements in j column;
				//then save their index of row into column_nozero
				int count_j_W=IA_W[j+1]-IA_W[j]; 
				column_nozero_W.reserve(count_j_W);
				for(int i=IA_W[j];i<IA_W[j]+count_j_W;i++){
					column_nozero_W.push_back(JA_W[i]);
					}
				//to calculate the number of zero elements in j column;
				//then save their index of row into column_nozero
				column_zero_W.reserve(j+1);
				for(int i=0;i<=j;i++){
					column_zero_W.push_back(i);
				}
				for(int i=0;i<count_j_W;i++){
				vector<int>::iterator it = column_zero_W.begin()+column_nozero_W[i]-i;
					column_zero_W.erase(it);
				}
				//////////////////////////////////////////////////////////////////////////////
				//////this part is used to get zero and nonzero elements in column i /////////
				//////////////////////////////////////////////////////////////////////////////
				int count2=0;
				vector<double> result_column3(i+1),result_column4(i+1);
				vector<int> column_nozero_i_Z,column_zero_i_Z;
				int count_i_Z=IA_Z[i+1]-IA_Z[i]; 
				column_nozero_i_Z.reserve(count_i_Z);
				
				for(int k=IA_Z[i];k<IA_Z[i]+count_i_Z;k++){
					column_nozero_i_Z.push_back(JA_Z[k]);
					}
				
				column_zero_i_Z.reserve(i+1);
				for(int k=0;k<=i;k++){
					column_zero_i_Z.push_back(k);
				}
				for(int i=0;i<count_i_Z;i++){
				vector<int>::iterator it = column_zero_i_Z.begin()+column_nozero_i_Z[i]-i;
					column_zero_i_Z.erase(it);
				}

				vector<int> column_nozero_i_W,column_zero_i_W;
				int count_i_W=IA_W[i+1]-IA_W[i]; 
				column_nozero_i_W.reserve(count_i_W);
				
				for(int k=IA_W[i];k<IA_W[i]+count_i_W;k++){
					column_nozero_i_W.push_back(JA_W[k]);
				}
				
				column_zero_i_W.reserve(i+1);
				for(int k=0;k<=i;k++){
					column_zero_i_W.push_back(k);
				}
				for(int i=0;i<count_i_W;i++){
				vector<int>::iterator it = column_zero_i_W.begin()+column_nozero_i_W[i]-i;
					column_zero_i_W.erase(it);
				}
			/////////////////////////////////////////////////////////////////////////////////////////
			/////////////////////calculate the nozero elements in column j///////////////////////////	
			/////////////////////////////////////////////////////////////////////////////////////////
			 for(int l=0; l<count_j_Z;l++){
					double zki=0;
					int index_i_Z=Binary_search(column_nozero_i_Z,column_nozero_Z[l]);
				if(index_i_Z>=0){
				
					zki=AA_Z[IA_Z[i]+index_i_Z];		
					double zkj=0;
					zkj=AA_Z[IA_Z[j]+l];

					zkj=zkj-p[j]/p[i]*zki;	
					result_column1[column_nozero_Z[l]]=zkj;
					
				}else if(column_nozero_Z[l]<=i){
					double zkj=0;		
					zkj=AA_Z[IA_Z[j]+l];
					zkj=zkj-p[j]/p[i]*zki;	
					result_column1[column_nozero_Z[l]]=zkj;
								

				}else{				
					double zkj=0;				
					zkj=AA_Z[IA_Z[j]+l];
					result_column1[column_nozero_Z[l]]=zkj;
					
					}
				}





			 for(int l=0; l<count_j_W;l++){
					double wki=0;
					int index_i_W=Binary_search(column_nozero_i_W,column_nozero_W[l]);
				if(index_i_W>=0){
					
					wki=AA_W[IA_W[i]+index_i_W];					
					double wkj=0;		
					wkj=AA_W[IA_W[j]+l];			
					wkj=wkj-q[j]/q[i]*wki;
					//transfer the value into result_column1
					result_column2[column_nozero_W[l]]=wkj;

				}else if(column_nozero_W[l]<=i){
					double wkj=0;			
					wkj=AA_W[IA_W[j]+l];//this wkj  need change
					
					wkj=wkj-q[j]/q[i]*wki;
					//transfer the value into result_column1
					result_column2[column_nozero_W[l]]=wkj;				

				}else{				
					double wkj=0;
					wkj=AA_W[IA_W[j]+l];		
					result_column2[column_nozero_W[l]]=wkj;
					}
				}
			

			/////////////////////////////////////////////////////////////////////////////////////////
			/////////////////////calculate the zero elements in column j///////////////////////////	
			/////////////////////////////////////////////////////////////////////////////////////////
				for(int l=0; l<column_zero_Z.size();l++){
					double zkj=0,wkj=0;
				int index_i_Z=Binary_search(column_zero_i_Z,column_zero_Z[l]);
				
				if(index_i_Z>=0){
					double zkj=0;
					double zki=0;
				
					zkj=zkj-p[j]/p[i]*zki;	
					//transfer the value into result_column1
					result_column1[column_zero_Z[l]]=zkj;
					
				
				}else if(column_zero_Z[l]<=i){
					double zkj=0;
					double zki=0;
					int index_i_Z_nozero=0;
					index_i_Z_nozero=Binary_search(column_nozero_i_Z,column_zero_Z[l]);
					zki=AA_Z[IA_Z[i]+index_i_Z_nozero];
					zkj=zkj-p[j]/p[i]*zki;	
				
					//transfer the value into result_column1
					result_column1[column_zero_Z[l]]=zkj;
				}else{					
					double zkj=0;//,wkj=0;
					result_column1[column_zero_Z[l]]=zkj;
					
				}
				}


				for(int l=0; l<column_zero_W.size();l++){
					double wkj=0;
				int index_i_W=Binary_search(column_zero_i_W,column_zero_W[l]);
				
					if(index_i_W>=0){
					double  wkj=0;
					double  wki=0;
				
					wkj=wkj-q[j]/q[i]*wki;
					//transfer the value into result_column1
					result_column2[column_zero_W[l]]=wkj;
				}else if(column_zero_W[l]<=i){
					double wkj=0;
					double wki=0;
					int index_i_W_nozero=0;
					index_i_W_nozero=Binary_search(column_nozero_i_W,column_zero_W[l]);
					wki=AA_W[IA_W[i]+index_i_W_nozero];
					wkj=wkj-q[j]/q[i]*wki;
					//transfer the value into result_column1
					result_column2[column_zero_W[l]]=wkj;			
				}else{					
					double wkj=0;
	
					result_column2[column_zero_W[l]]=wkj;
				}
				}

			////////////////////////////////////////////////////////////////////////
			////////////////////////drop process////////////////////////////////////
			////////////////////////////////////////////////////////////////////////
			vector<double> temp(result_column1);
			vector<double> temp2(result_column2);
			int num=nnz(temp);
		 	int num2=nnz(temp2);
			if(num>nmax_col){
				drop(temp,nmax_col,droptol,j);
			}
			if(num2>nmax_col){
				drop(temp2,nmax_col,droptol,j);
			}
			////////////////////////////////////////////////////////////////////////////
			////////// ////count the number of elements that added in column j//////////
			/////////////////////then update the value of IA,AA,JA//////////////////////
			////////////////////////////////////////////////////////////////////////////
			int delta_number=nnz(temp)-count_j_Z;
			int delta_number2=nnz(temp2)-count_j_W;

			if(j==column-1){
				IA_Z[j+1]=IA_Z[j+1]+delta_number;
			}else{
			for (int i=IA_Z[j+1];i<IA_Z.size();i++){
				IA_Z[i]=IA_Z[i]+delta_number;
			}
			}

			if(j==column-1){
				IA_W[j+1]=IA_W[j+1]+delta_number2;
			}else{
			for (int i=IA_W[j+1];i<IA_W.size();i++){
				IA_W[i]=IA_W[i]+delta_number2;
			}
			}

			//update the value of AA
			for(int k=0;k<count_j_Z;k++){
				vector<double>::iterator it = AA_Z.begin()+IA_Z[j];
					AA_Z.erase(it);
				}
			for(int k=0;k<count_j_Z;k++){
				vector<int>::iterator it = JA_Z.begin()+IA_Z[j];
					JA_Z.erase(it);
				}


			for(int k=0;k<count_j_W;k++){
				vector<double>::iterator it = AA_W.begin()+IA_W[j];
					AA_W.erase(it);
				}
			for(int k=0;k<count_j_W;k++){
				vector<int>::iterator it = JA_W.begin()+IA_W[j];
					JA_W.erase(it);
				}


			int count1_n=0;int count2_n=0;
			vector<double> ::iterator itr=AA_Z.begin()+IA_Z[j];
			for (int i=temp.size()-1;i>=0;i--){	
				if (temp[i]!=0){			 					
						AA_Z.insert(itr,temp[i]);	
						count1_n++;
				}
			}	
			for (int i=temp2.size()-1;i>=0;i--){		
				if (temp2[i]!=0){
					 vector<double> ::iterator itr=AA_W.begin()+IA_W[j];					
						AA_W.insert(itr,temp2[i]);		
						count2_n++;
				}
			}	
			int count3_n=0;int count4_n=0;
			for (int i=temp.size()-1;i>=0;i--){		
				if (temp[i]!=0){	
					 vector<int> ::iterator itr=JA_Z.begin()+IA_Z[j];						
						JA_Z.insert(itr,i);
						count3_n++;
				}
			}
		 for (int i=temp2.size()-1;i>=0;i--){		
				if (temp2[i]!=0){
					 vector<int> ::iterator itr=JA_W.begin()+IA_W[j];					
						JA_W.insert(itr,i);	
							count4_n++;
			  }
			 }	
		   }
		 }
	   }
	//}

// print the value of AA_Z,IA_Z,JA_Z
	cout<<"-----------------------AA_Z----------------------"<<endl;
	 vector<double>::iterator it11 ;
    for(it11=AA_Z.begin(); it11!=AA_Z.end(); it11++)
	{  cout<<*it11<<" " ;}
	cout<<endl;
	cout<<"-----------------------JA_Z----------------------"<<endl;
	vector<int>::iterator it12;
    for(it12=JA_Z.begin(); it12!=JA_Z.end(); it12++)
	{  cout<<*it12<<" " ;}
	cout<<endl;
	cout<<"-------------------------IA_Z--------------------"<<endl;
	vector<int>::iterator it13 ;
    for(it13=IA_Z.begin(); it13!=IA_Z.end(); it13++)
	{  cout<<*it13<<" " ;}
	cout<<endl;cout<<endl;cout<<endl;cout<<endl;


	cout<<"-----------------------AA_W----------------------"<<endl;
	 vector<double>::iterator it21 ;
    for(it21=AA_W.begin(); it21!=AA_W.end(); it21++)
	{  cout<<*it21<<" " ;}
	cout<<endl;
	cout<<"-----------------------JA_W----------------------"<<endl;
	vector<int>::iterator it22 ;
    for(it22=JA_W.begin(); it22!=JA_W.end(); it22++)
	{  cout<<*it22<<" " ;}
	cout<<endl;
	cout<<"-------------------------IA_W--------------------"<<endl;
	vector<int>::iterator it23 ;
    for(it23=IA_W.begin(); it23!=IA_W.end(); it23++)
	{  cout<<*it23<<" " ;}
	cout<<endl;


	cout<<"-------------------------P--------------------"<<endl;
	for (int i=0;i<row;i++){	
		cout<<p[i]<<"  ";
	}
	cout<<endl;cout<<endl;cout<<endl;cout<<endl;
}


