#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct atom;
struct nblist {
	struct atom* nb_atom;
	float distance;
	struct nblist *next;
};

struct atom {
	unsigned int atom_number;
	float x;
	float y;
	float z;
	struct nblist * neighbour_list;
};

void build_nblist(struct atom *head, int num_atoms, int cutoff);
float list_all_neighbour_pairs(struct atom *head, int num_atoms);

#define DIST(A,B) sqrt(pow((A->x - B->x), 2) + pow((A->y - B->y), 2) + pow((A->z - B->z),2))
#define VDW(dist) (0.1948)*(pow((3.2/dist), 12) - 2*pow((3.2/dist),6))

int main (int argc, char*argv[])
{


	FILE *pdbfile;
	char line[82];
	int cutoff=0, num_atoms=0, count=0;
	struct atom* head;
	float vdw_energy=0;

	if (argc!=4) {
		printf("Usage: ./nblist input_file_name num_atoms cutoff\n");
		exit(1);
	}

	num_atoms = atoi(argv[2]);
	cutoff= atoi(argv[3]);
	printf("NBLIST: Input parameters: File name %s  num_atoms %d  cutoff %d\n\n", argv[1], num_atoms, cutoff);

	pdbfile = fopen(argv[1], "r");
	if (!pdbfile) {
		printf("Can't open file\n");
		exit(1);
	}


	head = calloc(num_atoms, sizeof(struct atom)); //CHECK THIS
	if (!head) {
		printf("Insufficient memory\n");
		exit(1);
	}

	while (fgets(line, 82, pdbfile)) {
		if (!(strncmp(line, "ATOM", 4)) || !(strncmp(line, "HETATM", 6))) {
			//printf("%s", line);
			sscanf(&line[31], "%08f", &((head+count)->x));
			sscanf(&line[39], "%08f", &((head+count)->y));
			sscanf(&line[47], "%08f", &((head+count)->z));
			(head+count)->atom_number = count+1;
			//printf("%.3f %.3f %.3f\n",(head+count)->x, (head+count)->y, (head+count)->z);
			count++;
		}	
	}


	if (count != num_atoms) {
		printf ("atoms passed as param are different from actual; passed %d count: %d\n", num_atoms, count);
		exit(1);
	}
	build_nblist(head, num_atoms, cutoff);
	vdw_energy=list_all_neighbour_pairs(head, num_atoms);
	printf("\nvdw_energy:  %f\n", vdw_energy);

	fclose(pdbfile);
	return 0;
}

float list_all_neighbour_pairs(struct atom *head, int num_atoms)
{
	int i=0;
	struct nblist * temp;
	float vdw_energy=0 , dist=0;

	for (i=0; i < num_atoms; i++) {
		temp = (head+i)->neighbour_list;
		while (temp) {
			//printf("%u, %u, %f\n", (head+i)->atom_number, temp->nb_atom->atom_number, temp->distance);
			if((head+i)->atom_number < temp->nb_atom->atom_number) { //makes sure vanderwaal's energy is computed only once
				dist=temp->distance;
				printf("%u  %u %f\n", (head+i)->atom_number, temp->nb_atom->atom_number, temp->distance);
				vdw_energy+=VDW(dist);
			}
			temp = temp->next;
		}
	}
	return vdw_energy;
}

void build_nblist(struct atom *head, int num_atoms, int cutoff)
{
	int i=0, j=0;
	float distance;
	struct atom *cur_atom, *neighbour_candidate;
	struct nblist * nblist_head1, *nblist_head2, *temp;
	for (i=0; i < num_atoms; i++) {
		cur_atom = head + i;
		for (j=i+1; j < num_atoms; j++) {
			neighbour_candidate = head+j;
			distance = DIST(cur_atom, neighbour_candidate);
			if (distance <= cutoff) { //Candidate is a neighbour now
				//printf("Neighbours:%d, %d \n", cur_atom->atom_number, neighbour_candidate->atom_number);
				nblist_head1 = calloc (1, sizeof(struct nblist)); //will be the new HEAD of the list
				nblist_head1->distance = distance;
				nblist_head2 = calloc (1, sizeof(struct nblist)); //--ditto---
				nblist_head2->distance = distance;
				//Update both cur_atom's and candidate's lists, both are each other's neighbours
				if (!cur_atom->neighbour_list) {
					nblist_head1->nb_atom = neighbour_candidate;
					nblist_head1->next = NULL;
					cur_atom->neighbour_list = nblist_head1;
				}
				else {
					temp = cur_atom->neighbour_list;
					nblist_head1->nb_atom = neighbour_candidate;
					nblist_head1->next = temp;
					cur_atom->neighbour_list = nblist_head1;
				}
				if (!neighbour_candidate->neighbour_list) {
					nblist_head2->nb_atom = cur_atom;
					nblist_head2->next = NULL;
					neighbour_candidate->neighbour_list = nblist_head2;
				}
				else { 
					temp = neighbour_candidate->neighbour_list;
					nblist_head2->nb_atom = cur_atom;
					nblist_head2->next = temp;
					neighbour_candidate->neighbour_list = nblist_head2;
				}

			}

		}

		#if 0
		temp = cur_atom->neighbour_list;
		printf ("List for %d\n", cur_atom->atom_number);
		while (temp) {
			printf(">> %u, %u, %f\n", cur_atom->atom_number, temp->nb_atom->atom_number, temp->distance);
			temp = temp->next;
		}
		#endif
	}

}
