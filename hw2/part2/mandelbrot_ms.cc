
/**
 *  \file mandelbrot_ms.cc
 *
 *  \brief Implement your parallel mandelbrot set in this file.
 */

#include <iostream>
#include <cstdlib>
#include "render.hh"
#include <mpi.h>
using namespace std;
#define WIDTH 1000
#define HEIGHT 1000

int mandelbrot(double x, double y) {
	int maxit = 511;
	double cx = x;
	double cy = y;
	double newx, newy;
	
	int it = 0;
	for (it = 0; it < maxit && (x*x + y*y) < 4; ++it) {
		newx = x*x - y*y + cx;
		newy = 2*x*y + cy;
		x = newx;
		y = newy;
	}
	return it;
}


int main (int argc, char* argv[]){
	double minX = -2.1;
	double maxX = 0.7;
	double minY = -1.25;
	double maxY = 1.25;

	int height, width;
	if (argc == 3) {
		height = atoi (argv[1]);
		width = atoi (argv[2]);
		assert (height > 0 && width > 0);
	} else {
		fprintf (stderr, "usage: %s <height> <width>\n", argv[0]);
		fprintf (stderr, "where <height> and <width> are the dimensions of the image.\n");
		return -1;
	}

	double it = (maxY - minY)/height;
	double jt = (maxX - minX)/width;
	double x, y;

	//Start parallelizing 
	int rank, size;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;
	int rownum;

	// master process
	if(rank == 0){
		//int Sending[height * width];
		
		//start timing 
		double t1,t2;
		t1 = MPI_Wtime();
		
		rownum = 0;
		int s = size-1;    
		int receiver[width];
		int imagearr[height][width];
          
		//Send every process a row to compute
		while( s != 0 ){
			MPI_Send( &rownum, 1, MPI_INT, s, 1, MPI_COMM_WORLD);
			s--;
			rownum++;
		}    
		
       //Recieve and send till no more rows left
		while( rownum < height ){
			MPI_Recv( receiver, width, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Send( &rownum, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
			for(int i =0; i < width; i++){
				imagearr[status.MPI_TAG][i] = receiver[i];    
			}
			rownum++;
		}

		s = size - 1;
		int temp = -1;
		
		// recieve any remaining rows that haven't been recieved yet
		while( s != 0 ){
			MPI_Recv(receiver, width, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Send( &temp, 1 ,MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
			for(int i = 0; i < width ; i++){
			 imagearr[status.MPI_TAG][i] = receiver[i];
			}
			s--;
		}
		// rendering image
		gil::rgb8_image_t img(height, width);
		auto img_view = gil::view(img);
		for( int i = 0; i < height; i++){
			for( int j = 0; j < width; j++){
				img_view(j, i) = render( imagearr[i][j]/512.0 );
			}
		}
		
		t2 = MPI_Wtime();
		printf("(MS) Time taken: %f  with %d processes \n", (t2-t1),size);
		
		gil::png_write_view("mandelbrot_ms.png", const_view(img)); 
		
		
	}	
  
	//slave processes
	else {
		int sender[width];
		double x, y;
		bool flag = true;
		while( flag == true ){	
			//receive row number to compute
			MPI_Recv(&rownum, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			// if no more work to be done, master will send row number = -1
			if(rownum < 0) 
				break;          
			y = minY + (rownum * it);
			x = minX;
			for (int i = 0; i < width; ++i) {
				sender[i] = mandelbrot(x, y);
				x += jt;
			}
			//send back computed row
			MPI_Send(sender, width, MPI_INT, 0, rownum, MPI_COMM_WORLD);
		}  
	}
	//end parallelization
	MPI_Finalize();

}