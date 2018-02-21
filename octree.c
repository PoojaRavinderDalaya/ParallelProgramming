#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "octree.h"
//remember atom number

int main (int argc, char*argv[])
{

	FILE *pdbfile;
	char line[82];
	int cutoff=0, num_atoms=0, count=0, i;
	float max_xpair[2], max_ypair[2], max_zpair[2], bbx, bby, bbz, fswap;
	float shiftvalx, shiftvaly, shiftvalz ,side, vdw_energy;
	struct atom_list_t *head, *atom;
	struct atom_list_t **next;
	struct point_t new_origin, old_origin;
	struct octree_t *ocroot, *ocroot1;

	if (argc!=4) {
		printf("Usage: ./nblist input_file_name num_atoms cutoff\n");
		exit(1);
	}

	global_num_atoms=num_atoms = atoi(argv[2]);
	cutoff= atoi(argv[3]);
	printf("OCTREE:Input parameters: File name %s  num_atoms %d  cutoff %d\n\n", argv[1], num_atoms, cutoff);

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
	max_xpair[0]=0;
	max_ypair[0]=0;
	max_zpair[0]=0;
	max_xpair[1]=0;
	max_ypair[1]=0;
	max_zpair[1]=0;

	count=0;		
	//get the max min values x,y,z boundaries for building the initial cube
	while (fgets(line, 82, pdbfile)) {
		if (!(strncmp(line, "ATOM", 4)) || !(strncmp(line, "HETATM", 6))) {
			//printf("%s", line);
			atom=calloc(1, sizeof(struct atom_list_t));
			sscanf(&line[31], "%08f", &((atom->co_ords).x));
			sscanf(&line[39], "%08f", &((atom->co_ords).y));
			sscanf(&line[47], "%08f", &((atom->co_ords).z));
			atom->atom_number = count+1;
			atom->next_atom=NULL;
			*next=atom;
			next=&(atom->next_atom);

			if (count == 0){
				max_xpair[0]=(atom->co_ords).x;
				max_ypair[0]=(atom->co_ords).y;
				max_zpair[0]=(atom->co_ords).z;
			}
			else if (count ==1){
				if (max_xpair[0] <= (atom->co_ords).x)
					max_xpair[1] = (atom->co_ords).x;
				else {
					fswap = max_xpair[0];
					max_xpair[1]=fswap;
					max_xpair[0]=(atom->co_ords).x;
				}
				
				if (max_ypair[0] <= (atom->co_ords).y)
					max_ypair[1] = (atom->co_ords).y;
				else {
					fswap = max_ypair[0];
					max_ypair[1]=fswap;
					max_ypair[0]=(atom->co_ords).y;
				}

				if (max_zpair[0] <= (atom->co_ords).z)
					max_zpair[1] = (atom->co_ords).z;
				else {
					fswap = max_zpair[0];
					max_zpair[1]=fswap;
					max_zpair[0]=(atom->co_ords).z;
				}
			}
			else {
				if (((atom->co_ords).x) < max_xpair[0])
					max_xpair[0] = (atom->co_ords).x;
				else if (((atom->co_ords).x) > max_xpair[1])
					max_xpair[1] = (atom->co_ords).x;
				
				if (((atom->co_ords).y) < max_ypair[0])
					max_ypair[0] = (atom->co_ords).y;
				else if (((atom->co_ords).y) > max_ypair[1])
					max_ypair[1] = (atom->co_ords).y;
				
				if (((atom->co_ords).z) < max_zpair[0])
					max_zpair[0] = (atom->co_ords).z;
				else if (((atom->co_ords).z) > max_zpair[1])
					max_zpair[1] = (atom->co_ords).z;
			}
			//printf("%.3f %.3f %.3f\n",(atom->co_ords).x, (atom->co_ords).y, (atom->co_ords).z);
			count++;
		}
	}

	if (count != num_atoms) {
		printf ("/n/nEXITING: atoms passed as param are different from actual; passed %d count: %d\n\n", num_atoms, count);
		exit(1);
	}

	shiftvalx = max_xpair[0];
	shiftvaly = max_ypair[0];
	shiftvalz = max_zpair[0];

	get_shift_distance(&shiftvalx, &shiftvaly, &shiftvalz);
	new_origin.x = max_xpair[0] + shiftvalx -1;
	new_origin.y = max_ypair[0] + shiftvaly -1;
	new_origin.z = max_zpair[0] + shiftvalz -1;

	bbx=max_xpair[1]-max_xpair[0];
	bby=max_ypair[1]-max_ypair[0];
	bbz=max_zpair[1]-max_zpair[0];
	side = MAX(bbx, MAX(bby, bbz));

	shift_origin(head->next_atom, num_atoms, shiftvalx, shiftvaly, shiftvalz);

	printf ("\nmax x pairs: %f %f  y pairs: %f %f z pairs: %f %f\n", max_xpair[0], max_xpair[1], max_ypair[0], max_ypair[1], max_zpair[0], max_zpair[1]);

	printf("\nNew origin: %f %f %f \n\nMax bounding dimensions: %f x %f x %f and side %f\n\n", new_origin.x, new_origin.y, new_origin.z, bbx, bby, bbz, side);

	new_origin.x = (int)new_origin.x;
	new_origin.y = (int)new_origin.y;
	new_origin.z = (int)new_origin.z;
	side=(int)(side+1);
	printf("\nNew origin: %f %f %f \n\nMax bounding dimensions: %f x %f x %f and side %f\n\n", new_origin.x, new_origin.y, new_origin.z, bbx, bby, bbz, side);
	//done with building the initial cube
	//list_atoms(head->next_atom);

	ocroot=calloc(1, sizeof(struct octree_t));
	ocroot->atom_list=head->next_atom;
	ocroot->is_leaf=0;
	ocroot->num_atoms=count;
	ocroot->side=side;
	
	ocroot->corigin.x = new_origin.x;
	ocroot->corigin.y = new_origin.y;
	ocroot->corigin.z = new_origin.z;
	//copy node done
	//printf("expanding node\n");
	expand_node(ocroot, 4, 1); //build the octree
	//printf("expansion done\n");
	printLevelOrder(ocroot);
	printf(" built octree\n\n");
	//count_atoms(ocroot);
	printf("\n");
	printf(">>Global count %d\n", global_count);
	printf(">>Global count2 %d\n", global_count2);
	printf(">> Num leaves %d\n", num_leaves);

	ocroot1=ocroot;
	vdw_energy=get_neighbours(ocroot, ocroot1, cutoff); //actual computation
	printf(">> Num leaves2 %d\n", num_leaves2);
	printf("\n\nvdw energy: %f\n", vdw_energy);

	fclose(pdbfile);
	return 0;
}


void shift_origin (struct atom_list_t *head, int num_atoms, float shiftx, float shifty, float shiftz)
{
	struct atom_list_t *ptr;
	ptr=head;
	float shiftvalx, shiftvaly, shiftvalz;
	printf("\nShifting by %f %f %f\n", shiftx, shifty, shiftz);
	while(ptr)
	 {
		//printf ("%d: (%.3f %.3f %.3f)  ", ptr->atom_number, (ptr->co_ords).x, (ptr->co_ords).y, (ptr->co_ords).z);
		(ptr->co_ords).x = (ptr->co_ords).x + shiftx;
		(ptr->co_ords).y = (ptr->co_ords).y + shifty;
		(ptr->co_ords).z = (ptr->co_ords).z + shiftz;
		//printf (" (%.3f %.3f %.3f)\n",(ptr->co_ords).x, (ptr->co_ords).y, (ptr->co_ords).z);
		ptr=ptr->next_atom;
	}
}

void get_shift_distance(float *shiftx, float *shifty, float *shiftz)
{
	if ((((int)fabs(*shiftx) +1) + *shiftx) < 1)
		*shiftx=(int)fabs(*shiftx) + 2;
	else
		*shiftx= 1;


	if ((((int)fabs(*shifty) + 1) + *shifty) < 1)
		*shifty=(int)fabs(*shifty) + 2;
	else
		*shifty= 1;


	if ((((int)fabs(*shiftz) + 1) + *shiftz) < 1)
		*shiftz=(int)fabs(*shiftz) + 2;
	else
		*shiftz= 1;

}
/*
struct octree_t {
	bool is_leaf; //looking at the box only makes sense when its a leaf
	struct point_t corigin;
	float side;
	int num_atoms;
	struct atom_list_t *atom_list;//leaf node will have a non-empty list
	struct octree_t *next_oct;//list of non-empty chidren
};


struct atom_list_t {
	unsigned int atom_number;
	struct point_t co_ords;
	struct atom_list_t * next_atom;
};
*/

void expand_node(struct octree_t* parent, int k , int alpha)
{
	struct atom_list_t *pptr, *cptr1, *cptr2;
	struct octree_t *subcubes[8], *ocptr1, *ocptr2;
	int i=0;

	//printf("atoms: %d \n", parent->num_atoms);	
	if(parent->num_atoms <= alpha*k) {
		//printf("found_leaf with %d atoms\n", parent->num_atoms);
		global_count=global_count+parent->num_atoms;
		parent->is_leaf=1;
		return;
	}
	else {
	//printf("num:atoms %d and side %f\n", parent->num_atoms, parent->side);
	parent->is_leaf=0;
	while(i <8) {
		subcubes[i]=calloc(1, sizeof(struct octree_t));
		compute_subc_origin(&(parent->corigin), &(subcubes[i]->corigin), i, (float)((parent->side)/2)); //get the subcube origin
		subcubes[i]->side=(float)(parent->side/(float)2);
		//printf("i %d: origin (%f %f %f) side %f\n", i, subcubes[i]->corigin.x, subcubes[i]->corigin.y, subcubes[i]->corigin.z, subcubes[i]->side);
		subcubes[i]->is_leaf=1;
		i++;
	}

	pptr = parent->atom_list;

	while(pptr) {
		//printf("list:origin: (%f  %f  %f) side %f atom_no. %d  pt. co-ords: %f %f %f\n", parent->corigin.x, parent->corigin.y, parent->corigin.z, parent->side, pptr->atom_number, pptr->co_ords.x, pptr->co_ords.y ,pptr->co_ords.z);
		i=get_subcube_index(&(parent->corigin), &(pptr->co_ords), (float)((parent->side)/2));
		cptr2=subcubes[i]->atom_list;
		cptr1=calloc(1, sizeof(struct atom_list_t));
		cptr1->atom_number=pptr->atom_number;
		cptr1->co_ords.x=pptr->co_ords.x;
		cptr1->co_ords.y=pptr->co_ords.y;
		cptr1->co_ords.z=pptr->co_ords.z;
		if(!cptr2)
			subcubes[i]->atom_list=cptr1;
		else {
			cptr1->next_atom=cptr2;
			subcubes[i]->atom_list=cptr1;
		}
		subcubes[i]->num_atoms++;
		pptr=pptr->next_atom;
	}

	i=0;
	while(i<8) {
		ocptr2=subcubes[i];
		if(ocptr2->num_atoms > 0) {
			//printf("Child %d has %d atoms\n",i, ocptr2->num_atoms);
			parent->octchild[i]=ocptr2;
			expand_node(ocptr2, k, alpha);
		}
		else {
			free(ocptr2);
			parent->octchild[i]=NULL;
		}
		i++;
	}
	pptr = parent->atom_list;
	parent->atom_list = NULL;
	free(pptr);
	}
}


void compute_subc_origin(struct point_t* base_orig, struct point_t* subc_orig, int oct_num, float side) //octant_number
{

	float base_x,base_y,base_z;
	base_x=base_orig->x;
	base_y=base_orig->y;
	base_z=base_orig->z;

	//printf("Oct_num: %d\n", oct_num);
	switch(oct_num){
		case(0):
			subc_orig->x=base_x;
			subc_orig->y=base_y;
			subc_orig->z=base_z;
			break;
		case(1):
			subc_orig->x=base_x + side;
			subc_orig->y=base_y;
			subc_orig->z=base_z;
			break;
		case(2):
			subc_orig->x=base_x;
			subc_orig->y=base_y + side;
			subc_orig->z=base_z;
			break;
		case(3):
			subc_orig->x=base_x + side;
			subc_orig->y=base_y + side;
			subc_orig->z=base_z;
			break;
		case(4):
			subc_orig->x=base_x;
			subc_orig->y=base_y;
			subc_orig->z=base_z + side;
			break;
		case(5):
			subc_orig->x=base_x + side;
			subc_orig->y=base_y;
			subc_orig->z=base_z + side;
			break;
		case(6):
			subc_orig->x=base_x;
			subc_orig->y=base_y + side;
			subc_orig->z=base_z + side;
			break;
		case(7):
			subc_orig->x=base_x + side;
			subc_orig->y=base_y + side;
			subc_orig->z=base_z + side;
			break;
		default:
			printf("Invalid point and oct_num %d\n", oct_num);
			exit(1);
	}
}

int get_subcube_index(struct point_t* base_orig, struct point_t* pointc, float side)
{
	
	float base_x,base_y,base_z,sub_x,sub_y,sub_z;
	base_x=base_orig->x;
	base_y=base_orig->y;
	base_z=base_orig->z;

	sub_x=pointc->x;
	sub_y=pointc->y;
	sub_z=pointc->z;

	if((sub_x <= (base_x + side)) && (sub_y <= (base_y + side)) && (sub_z <= (base_z + side)))
		return 0;
	else if((sub_x > (base_x + side)) && (sub_y <= (base_y + side)) && (sub_z <= (base_z + side)))
		return 1;
	else if((sub_x <= (base_x + side)) && (sub_y > (base_y + side)) && (sub_z <= (base_z + side)))
		return 2;
	else if((sub_x > (base_x + side)) && (sub_y > (base_y + side)) && (sub_z <= (base_z + side)))
		return 3;
	else if((sub_x <= (base_x + side)) && (sub_y <= (base_y + side)) && (sub_z > (base_z + side)))
		return 4;
	else if((sub_x > (base_x + side)) && (sub_y <= (base_y + side)) && (sub_z > (base_z + side)))
		return 5;
	else if((sub_x <= (base_x + side)) && (sub_y > (base_y + side)) && (sub_z > (base_z + side)))
		return 6;
	else if((sub_x > (base_x + side)) && (sub_y > (base_y + side)) && (sub_z > (base_z + side)))
		return 7;
	else {
		printf("Exiting\n");
		exit(1);
	}
}

void list_atoms(struct atom_list_t *atom_list)
{
	struct atom_list_t *cur_atom;
	cur_atom=atom_list;
	while(cur_atom) {
		//printf("%d %.3f %.3f %.3f\n", cur_atom->atom_number, cur_atom->x, cur_atom->y, cur_atom->z);
		printf("%.3f  %.3f  %.3f\n", (cur_atom->co_ords).x, (cur_atom->co_ords).y, (cur_atom->co_ords).z);
		cur_atom=cur_atom->next_atom;
	}
}

/*
struct octree_t {
	//struct cube_t *box;//is it's a leaf the box can't contain more than k@ atoms in k-addmissible octree
	bool is_leaf; //looking at the box only makes sense when its a leaf
	struct point_t corigin;
	float side;
	int num_atoms;
	struct atom_list_t *atom_list;//leaf node will have a non-empty list
	struct octree_t *next_oct;//list of non-empty chidren
};
*/
void count_atoms(struct octree_t *parent)
{
	struct octree_t *ptr=parent, *ptr2;
	struct atom_list_t *atom_listu;
	int j=0,i=0;
	//printf("num_atoms %d is_leaf %d\n", ptr->num_atoms, ptr->is_leaf);


	if(parent->is_leaf) {
		//printf("leaf: origin: (%f %f %f) : %f  num_atoms: %d\n", ptr->corigin.x, ptr->corigin.y, ptr->corigin.z, ptr->side, ptr->num_atoms);
		printf("leaf oct_num: %d num_atoms %d\n", parent->oct_num, ptr->num_atoms);
		atom_listu=ptr->atom_list;
		while(atom_listu){
			printf("<%d> ", atom_listu->atom_number);	
			atom_listu=atom_listu->next_atom;
		}
		printf("\n\n");
		num_leaves++;
		global_count2=global_count2+ ptr->num_atoms;
		return;
	}
	while(i<8){
		//printf("on going\n");
		if(ptr->octchild[i]) {
			ptr2=ptr->octchild[i];
			//printf("i %d  origin: (%f %f %f) : %f  num_atoms: %d, isleaf %d\n", i,ptr2->corigin.x, ptr2->corigin.y, ptr2->corigin.z, ptr2->side, ptr2->num_atoms, ptr2->is_leaf);
			count_atoms(ptr->octchild[i]);
		}
		else{
			i++;
			continue;
		}
		
		i++;
		j++;
	}
	//printf("j %d\n",j);
	return;
}

float get_neighbours(struct octree_t* u, struct octree_t* v, int cutoff)
{
	int i,j;
	struct atom_list_t *atom_listu, *atom_listv;
	float dist, energy=0;
	if(!(is_within_cutoff(u, v, cutoff))){
		return 0;
	}
	else if ((u->is_leaf) && (v->is_leaf) && (u->oct_num <= v->oct_num)) {
		num_leaves2++;
		atom_listu=u->atom_list;
		atom_listv=v->atom_list;
		while (atom_listu) {
			atom_listv=v->atom_list;
			while(atom_listv) {
				//printf("potential: %d (%f %f %f)  %d (%f %f %f)\n", atom_listu->atom_number, atom_listu->co_ords.x, atom_listu->co_ords.y, atom_listu->co_ords.z, atom_listv->atom_number, atom_listv->co_ords.x, atom_listv->co_ords.y, atom_listv->co_ords.z);
				//printf("cneighbours oct_num %d  %d  oct: %d %d\n", u->oct_num, atom_listu->atom_number, v->oct_num, atom_listv->atom_number);
				if (u->oct_num == v->oct_num) {
					dist=(DIST((atom_listu->co_ords), (atom_listv->co_ords)));
					if ((dist <= cutoff)  && (atom_listu->atom_number < atom_listv->atom_number)) {
						energy+=VDW(dist);
						printf("neighbours: %d  %d %f\n", atom_listu->atom_number, atom_listv->atom_number, dist);
					}
				}
				else  {
					dist=(DIST((atom_listu->co_ords), (atom_listv->co_ords)));
					if ((dist <= cutoff) && (atom_listu->atom_number != atom_listv->atom_number)) {
						energy+=VDW(dist);
						printf("neighbours: %d  %d %f\n", atom_listu->atom_number, atom_listv->atom_number, dist);
					}
				}

				atom_listv=atom_listv->next_atom;
			}
			atom_listu=atom_listu->next_atom;
		}
		return energy;
	}
	else if (u->is_leaf) {
		i=0;
		while(i<8) {
			if (v->octchild[i])
				energy+= get_neighbours(u, v->octchild[i], cutoff);
			i++;
		}
		return energy;
	}
	else if (v->is_leaf) {
		i=0;
		while(i<8) {
			if (u->octchild[i])
				energy += get_neighbours(u->octchild[i], v, cutoff);
			i++;
		}
		return energy;
	}
	else {
		i=0;j=0;
		while(i<8) {
			if (u->octchild[i]){
				j=0;
				while(j<8) {
					if(v->octchild[j]){
						energy += get_neighbours(u->octchild[i], v->octchild[j], cutoff);
					}
					j++;
				}
			}
			i++;
		}
		return energy;
	}
}

bool is_within_cutoff (struct octree_t* u, struct octree_t* v, int cutoff)
{
	float xd, yd, zd, dist;
	float p1, p2, p3;

	if (u->corigin.x <= v->corigin.x) {
		p1=u->corigin.x;
		p2=p1+ u->side;
		p3=v->corigin.x;
	}
	else {
		p1=v->corigin.x;
		p2=p1+ v->side;
		p3=u->corigin.x;
	}
	xd=p3-p2;

	if (u->corigin.y <= v->corigin.y) {
		p1=u->corigin.y;
		p2=p1+ u->side;
		p3=v->corigin.y;
	}
	else {
		p1=v->corigin.y;
		p2=p1+ v->side;
		p3=u->corigin.y;
	}
	yd=p3-p2;
	if (u->corigin.z <= v->corigin.z) {
		p1=u->corigin.z;
		p2=p1+ u->side;
		p3=v->corigin.z;
	}
	else {
		p1=v->corigin.z;
		p2=p1+ v->side;
		p3=u->corigin.z;
	}
	zd=p3-p2;

	dist=fmaxf(xd, fmaxf(yd,zd));

	//printf("Cpoints oct: %d %.3f %.3f %f side %.3f -- oct: %d %.3f %.3f %.3f side %.3f  ==> %f\n", u->oct_num, u->corigin.x, u->corigin.y, u->corigin.z, u->side, v->oct_num, v->corigin.x,v->corigin.y,v->corigin.z, v->side, dist);

	//printf("Cpoints_oct %d %d ==> %f\n", u->oct_num, v->oct_num, dist);
	if(dist > cutoff)
		return 0;
	else 
		return 1;

}
