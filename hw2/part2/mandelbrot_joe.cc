/**
 *  \file mandelbrot_joe.cc
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

int
mandelbrot(double x, double y) {
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

int
main(int argc, char* argv[]) {
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
	
	//Start parallelisation
	int rank, size;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	double t1,t2;
	
	if (rank == 0) {
		t1 = MPI_Wtime();
	} 
	
	int rowDis = height/size; 				  // number of rows by proc size so that all get even number
	int dataprow = rowDis * width; 
	int SepAr[(width*rowDis)];               //Each proc has seperate array of size (width * rowdis )
	//int twoDAr[rowDis][width];			 // 2d doesnt work. 
	int *recieveBuffer = NULL;
	int LeftOut;
	int *LeftArray;	
	y = minY + ( rank * ( it * rowDis) );
	for (int i = 0; i < rowDis; ++i) {
		x = minX;
		for (int j = 0; j < width; ++j) {
			//img_view(j, i) = render(mandelbrot(x, y)/512.0);
			SepAr[ (i * width) + j] = mandelbrot(x,y);
			//twoDAr[i][j] = mandelbrot(x,y);
			x += jt;
		}
		y += it;
	}		
	
	MPI_Barrier(MPI_COMM_WORLD);
		
	// Trailing shit
	if(rank == 0 ) { 
		recieveBuffer = new int[ size * dataprow];
		LeftOut = ( height - ( size * rowDis) );
		LeftArray = new int[ LeftOut * width ];
	
		if(LeftOut > 0 ){
			y = minY + ( LeftOut * ( it * rowDis) );
			for (int i = 0; i < LeftOut; ++i) {
				x = minX;
				for (int j = 0; j < width; ++j) {
					//img_view(j, i) = render(mandelbrot(x, y)/512.0);
					LeftArray[ ( i * width ) + j] = mandelbrot(x,y);      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					//twoDAr[i][j] = mandelbrot(x,y);
					x += jt;
				}
				
				y += it;
			}
		} 
	}
	
	MPI_Gather(SepAr,dataprow, MPI_INT,recieveBuffer,dataprow, MPI_INT,0, MPI_COMM_WORLD);
	
	if (rank == 0) {

		gil::rgb8_image_t img(height, width);
		auto img_view = gil::view(img);
		for (int i = 0; i < (rowDis * size); ++i) {
			for (int j = 0; j < width; ++j) {
				img_view(j, i) = render( recieveBuffer[ (i * width) + j] / 512.0);
			}
		} 			
		for (int k = 0; k < LeftOut; ++k) {
			for (int p = 0; p < width; ++p) {
				//img_view(p, k + (rowDis * size) ) = render( LeftArray[ (k * width) + p] / 512.0);
			}
		} 
		
			// CAlculating runtime
		t2 = MPI_Wtime ();
		printf("(Joe) Time taken: %f  with %d processes \n", (t2-t1),size);
		
	
		gil::png_write_view("mandelbrot_joe.png", const_view(img));

		

		
	}

	MPI_Finalize();
	return 1;
 }

/* eof */