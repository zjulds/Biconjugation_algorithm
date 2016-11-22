#include "Binary_search.h"
/**Binary_search this algorithm is programmed according to <大话数据结构> p299
	author :liaodashuang 
*/
int Binary_search(vector<int> wait_find,double value){
	int low,high,mid;
	low=0;
	high=wait_find.size()-1;
	while (low<=high){
		mid = (high+low)/2;
		if (value<wait_find[mid]){
			high=mid-1;
		}else if(value>wait_find[mid]){
			low = mid+1;
		}else{
			return mid;
		}
	}
	return -1;
}