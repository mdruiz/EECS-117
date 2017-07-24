/**
 *  \file parallel-mergesort.cc
 *
 *  \brief Implement your parallel mergesort in this file.
 */
#include<omp.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "sort.hh"

void doMerge(keytype* A,int lower, int higher, int num_threads, int* temp_arr );
void MergeParts(int low1, int high1, int low2, int high2, keytype* array, int* temp_arr , int k);
void merge(int lower, int middle, int high, keytype* array, int num_threads, int* temp_arr );
int binary_search(keytype* a, int high, int low, int search);
void serial_merge( int low, int mid, int high, keytype* a);
//void serial_merge(int low, int mid, int high, keytype* A);


void parallelSort (int N, keytype* A)
	{
	  /* Lucky you, you get to start from scratch */
	  int num_threads;
	  

     #pragma omp parallel
	 #pragma omp single
		{
			//#pragma omp master
			//{
				//number of physical threads
				num_threads = omp_get_num_threads();
				//printf("num: %d", num_threads);
				int* temp_arr;
				temp_arr = new int[N];
			
				doMerge(A, 0, N-1, ( (N-1)/num_threads ), temp_arr );	//calling doMerge with number of logical threads needed
			//}
		}
	  
	}


	void doMerge(keytype* A, int lower, int higher, int num_threads, int* temp_arr ){

    int middle;
		if( lower < higher ){ // exit condition
			
			if( ( higher - lower +1) > num_threads && num_threads!= 0){   // parallelism end condition
			  // if( num_threads > 1){
			  
			  //omp_set_num_threads(num_threads);
			  //	printf("Number of threads: %d  , %d    \n", omp_get_num_threads(),num_threads) ;
        	middle = (lower + (higher-lower)/2);
				
				//int inter = (middle-lower)/num_threads;
        	#pragma omp task
			doMerge(A,lower,middle, (middle-lower)/num_threads, temp_arr ); // seperate task for half of array
					//doMerge(A,lower,middle, num_threads  );
				
        	#pragma omp task
			doMerge( A, middle+1, higher, (higher-middle-1)/num_threads, temp_arr);
					//doMerge( A, middle+1, higher, num_threads   );
          				
        	#pragma omp taskwait    
			
/* 				for (int i = 0; i <= higher; i++){
				temp_arr[i] = A[i];
				//printf("%d,",array[i]);
				}
					 */
				merge(lower,middle,higher, A, num_threads, temp_arr); // merging two sorted subarrays
			}
      
			else{		 // code for serial codce when thread limit is reached
					middle = (lower + (higher - lower)/2);
					//if(num_threads == 0) num_threads =1;
					//printf("%d",(higher-lower+1));
					doMerge(A,lower,middle, num_threads/4, temp_arr);
					//doMerge(A,lower,middle, (middle-lower)/num_threads, temp_arr  );
					
					doMerge( A, middle+1, higher, num_threads/4, temp_arr );
					
					//doMerge( A, middle+1, higher, (higher-middle-1)/num_threads, temp_arr  );
					//#pragma omp single
					//merge(lower,middle,higher, A, 1);
					serial_merge(lower,middle,higher,A);
			}
      
      //merge(lower,middle,higher, A, num_threads);
		}
	}
 
 void serial_merge (int lower, int middle, int high, keytype* array) { // regular mergesort

		int tempar[ high+1];
		
		for(int i = lower; i<= high; i++){
			tempar[i] = array[i];
		} 
		int i = lower;
		int j = middle+1;
		int k = lower;
		
		while( i <= middle && j <= high){
			if( tempar[i] >= tempar[j] ){
				array[k] = tempar[j];
				j++;
			}
			else{
				array[k] = tempar[i];
				i++;
			}
			k++;
		}
		while(i <= middle){
			array[k] = tempar[i];
			i++;
			k++;
		}
	}


void merge(int lower, int middle, int high, keytype* array, int num_threads, int* temp_arr ){ //combining two arrays with parallelism
	 
	int mid_index1 = lower + ((middle - lower)/2);  //middle index of first array
    int mid_value1 = array[mid_index1];          //value of middle index
    int pos = binary_search(array,high,(middle+1),mid_value1);
	
	//#pragma omp parallel
	//#pragma omp single nowait
    int i = 0;
 
    for (i = lower; i <= high; i++){
        temp_arr[i] = array[i];
		//printf("%d,",array[i]);
    } 
    int k = (mid_index1 - lower + 1) + ( pos - (middle+1)) + lower;
	  
	if( num_threads > 10*( high - lower ) ){ //if logical threads available 
				
				int mid_index11 = lower + ((middle - lower)/2);
				int mid_index2 = middle + ((high - middle)/2);		
			
			//	 merge( lower, mid_index11, middle, array , num_threads -1 );
			merge( lower, mid_index11, middle, array , (middle-lower)/num_threads, temp_arr );
				 //merge( middle+1, mid_index2, high, array, num_threads -1 );
			merge( middle+1, mid_index2, high, array, (high-middle-1)/num_threads, temp_arr );
	}

		// parallelizing MergeParts()
		#pragma omp task
			MergeParts(lower, mid_index1, (middle + 1), pos-1, array, temp_arr, lower);
		#pragma omp task
			MergeParts((mid_index1 + 1), middle, pos, high, array, temp_arr, k);
		#pragma omp taskwait
  }
 
 	void MergeParts(int low1, int high1, int low2, int high2, keytype* array, int temp_arr[], int kk){
		
		int i = low1;
		int j = low2;
		//int k = low1;
		int k = kk;
		while( i <= high1 && j <= high2){
			if( temp_arr[i] >= temp_arr[j] ){
				array[k] = temp_arr[j];
				j++;
			}
			else{
				array[k] = temp_arr[i];
				i++;
			}
			k++;
		}
		while(i <= high1){
			array[k] = temp_arr[i];
			i++;
			k++;
		}
		while(j <= high2){ 
			 array[k] = temp_arr[j];
			 j++;
			 k++;
		   }
  }
 
	int binary_search(keytype* a, int r, int p, int x){
		int mid = 0;
		int low = p;
		int high = 0;
		if (p>(r+1)) {
			high = p;
		}
		else{
			high = r+1;
		}
		//int high = r;
		int i = 0;
		while (low< high){
			
			mid = (low+high)/2;
			
			if (x <= a[mid]){	
				high = mid;	
			}
			else{
				low = mid + 1;
			}
		}
		return high; 
	}	


/* eof */
