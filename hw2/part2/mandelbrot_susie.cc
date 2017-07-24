/**
 *  \file mandelbrot_susie.cc
 *
 *  \brief Implement your parallel mandelbrot set in this file.
 */

#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include "render.hh"
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

int main(int argc, char* argv[]) {
	/* Lucky you, you get to write MPI code */
	
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

	//Start Parallelization
	MPI_Init(&argc, &argv);
	
	int size, rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	double t1, t2;
	
	// start timer on main process
	if (rank == 0) {
		t1 = MPI_Wtime();
	} 
	
	
	//int maxRows = (height / size)+1; 
	int maxRows = (height / size)+1; 
	int dataSize = maxRows * width;
	int sendBuffer[dataSize];  
	int * receiveBuffer = NULL;
	int row, start;
	
	//printf ("I am %d of %d\n", rank, size);

	row = 0;
	y = minY + (rank*it);
	
	//calculating mandelbrot values
	for(int i = rank; i < height; i+=size){
		x = minX;
		for(int j = 0; j< width; ++j){
			sendBuffer[(row * width)+j] = mandelbrot(x,y); 
			x += jt;
		}
		//y = minY + (i+1)*it;
		y += (it * size);
		//row += 1;
		row++;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	//if main proccess
	if(rank == 0){
		receiveBuffer = new int [dataSize*size];
	}
	
	//gathering data from all processes
	MPI_Gather(sendBuffer, dataSize, MPI_INT, receiveBuffer, dataSize, MPI_INT, 0, MPI_COMM_WORLD);
	
	//main process renders image
	if (rank == 0){
		
		gil::rgb8_image_t img(height, width);
		auto img_view = gil::view(img);

		row = 0;
		start = 0;
		for (int i = 0; i < height; ++i) {
			row = i /size;  
			for (int j = 0; j < width; ++j) {
				start = (i%size) * dataSize;
				img_view(j, i) = render(receiveBuffer[start + (row*width) +j]/512.0);
			}
			//row = i /size;
		} 
		
		//calculating run time 
		t2 = MPI_Wtime();
		printf("(Susie) Time taken: %f  with %d processes \n", (t2-t1),size);
		
		gil::png_write_view("mandelbrot_susie.png", const_view(img));
		
	}
	
	//end parallelization
	MPI_Finalize();
	
}

	/* eof */
