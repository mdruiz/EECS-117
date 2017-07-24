#include <stdio.h>
#define max 10




int main() { 
 
	int a[10] = { 10, 14, 19, 26, 27, 31, 33, 35, 42, 44 };
   printf("List before sorting\n");
 
	int k = BS(a, 9, 0,28 );
	printf("%d",k);
}

 int BS ( int a[], int high, int low, int search){
	int l = low;
	int h = high;
	while(true){
		int mid = ( l + (h - l)/2  );
		if(search >= a[mid]){
		if(search < a[mid+1]){
				return mid+1;
			}
		else if( mid + 2 > high){
				return high;
			}

			else{
				l = mid;
			}
		}
		
		else{
			if(seach > a[mid-1]){
				return mid;
			}
			else if( mid - 2 < low){
				return low;
			}

			else{
				h = mid;
			}
		}
		
	}
}

	 






