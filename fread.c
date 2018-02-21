#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main (int argc, char*argv[])
{
	FILE *pdbfile;
	char line[82];
	float x, y, z;
	int cutoff=0, num_atoms=0;

	if (argc!=2) {
		printf("Usage: ./nblist input_file_name num_atoms cutoff\n");
		exit(1);
	}

	pdbfile = fopen(argv[1], "r");
	if (!pdbfile)
		return -1;

	printf("Input parameters %s  %d  %d\n\n", argv[1], num_atoms, cutoff);

	while (fgets(line, 82, pdbfile)) {
		if (!(strncmp(line, "ATOM", 4)) || !(strncmp(line, "HETATM", 6))) {
			//printf("%s", line);
			sscanf(&line[31], "%08f", &x);
			sscanf(&line[39], "%08f", &y);
			sscanf(&line[47], "%08f", &z);
			printf("%.3f %.3f %.3f\n",x, y, z);
		}
			
	}

	fclose(pdbfile);
}
