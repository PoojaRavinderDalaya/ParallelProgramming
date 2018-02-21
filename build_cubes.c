#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FIRST_SCUBE(X,Y,Z, SIDE) (X,Y,Z)
#define SECOND_S

int main (int argc, char*argv[])
{

	float x, y, z;
	float side;
	if (argc!=5) {
		printf("Usage: \n");
		exit(1);
	}
	x = atoi(argv[1]);
	y = atoi(argv[2]);
	z = atoi(argv[3]);
	side = atoi(argv[4]);
	/*
	printf("Input parameters: Origin: %f %f %f  Side: %f\n\n", x, y , z, side);

	printf("Co-ordinate of subcube 1: %f %f %f\n", x, y ,z);
	printf("Co-ordinate of subcube 2: %f %f %f\n", x+side/(float)8, y,z);
	printf("Co-ordinate of subcube 3: %f %f %f\n", x, y+side/(float)8,z);
	printf("Co-ordinate of subcube 4: %f %f %f\n", x+side/8, y+side/8, z);
	printf("Co-ordinate of subcube 5: %f %f %f\n", x,y,z+side/8);
	printf("Co-ordinate of subcube 6: %f %f %f\n", x+side/8, y, z+side/8);
	printf("Co-ordinate of subcube 7: %f %f %f\n", x, y+side/8, z+side/8);
	printf("Co-ordinate of subcube 8: %f %f %f\n", x+side/8, y+side/8, z+side/8);
	*/
	
	printf("%f %f %f\n", 2.00- 3.5,fmaxf(-3.00, -4.05), fmaxf(-4.00, 3.00));
}
