#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
struct atom_list_t;

struct point {
	float x;
	float y;
	float z;
};

struct atom_list_t {
	unsigned int atom_number;
	float x;
	float y;
	float z;
	struct atom_list_t * next_atom;
};

void list_atoms(struct atom_list_t *atom_list);


int main (int argc, char*argv[])
{

	FILE *pdbfile;
	char line[82];
	int cutoff=0, num_atoms=0, count=0, i;
	float max_xpair[2], max_ypair[2], max_zpair[2], bbx, bby, bbz, temp;
	float shiftvalx, shiftvaly, shiftvalz ,side;
	struct atom_list_t *head, *atom;
	struct atom_list_t **next;
	struct point new_origin, old_origin;

	if (argc!=4) {
		printf("Usage: ./nblist input_file_name num_atoms cutoff\n");
		exit(1);
	}

	num_atoms = atoi(argv[2]);
	cutoff= atoi(argv[3]);
	printf("Input parameters: File name %s  num_atoms %d  cutoff %d\n\n", argv[1], num_atoms, cutoff);

	pdbfile = fopen(argv[1], "r");
	if (!pdbfile) {
		printf("Can't open file\n");
		exit(1);
	}


	head = calloc(1, sizeof(struct atom_list_t));
	if (!head) {
		printf("Insufficient memory\n");
		exit(1);
	}

	next=&(head->next_atom);
	printf("head->next_atom %llx\n", head->next_atom);

	count=0;		
	while (fgets(line, 82, pdbfile)) {
		if (!(strncmp(line, "ATOM", 4)) || !(strncmp(line, "HETATM", 6))) {
			//printf("%s", line);
			atom=calloc(1, sizeof(struct atom_list_t));
			sscanf(&line[31], "%08f", &((atom)->x));
			sscanf(&line[39], "%08f", &((atom)->y));
			sscanf(&line[47], "%08f", &((atom)->z));
			atom->atom_number = count+1;
			atom->next_atom=NULL;
			*next=atom;
			next=&(atom->next_atom);
			count++;
			}
			//printf("%.3f %.3f %.3f\n",(atom)->x, (atom)->y, (atom)->z);
	}

	printf("head->next_atom %llx\n", head->next_atom);
	if (count != num_atoms) {
		printf ("/n/nEXITING: atoms passed as param are different from actual; passed %d count: %d\n\n", num_atoms, count);
		exit(1);
	}

	list_atoms(head->next_atom);
	fclose(pdbfile);
	return 0;
}


void list_atoms(struct atom_list_t *atom_list)
{
	struct atom_list_t *cur_atom;
	cur_atom=atom_list;
	while(cur_atom) {
		printf("%d %f %f %f\n", cur_atom->atom_number, cur_atom->x, cur_atom->y, cur_atom->z);
		cur_atom=cur_atom->next_atom;
	}
}
