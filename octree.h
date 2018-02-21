
#define MAX(X,Y) ((X)>(Y) ? X:Y)
int global_num_atoms;
struct atom_list_t;

struct point_t {
	float x;
	float y;
	float z;
};

struct cube_t {
	struct point_t corigin;
	float side;
	int num_atoms; //number of atoms inside this cube
	struct atom_list_t *cube_atom; //list of atoms this cube contains
};

struct atom_list_t {
	unsigned int atom_number;
	struct point_t co_ords;
	struct atom_list_t * next_atom;
};

struct octree_t {
	//struct cube_t *box;//is it's a leaf the box can't contain more than k@ atoms in k-addmissible octree
	bool is_leaf; //looking at the box only makes sense when its a leaf
	struct point_t corigin;
	float side;
	int num_atoms;
	struct atom_list_t *atom_list;//leaf node will have a non-empty list
	int oct_num;
	struct octree_t *octchild[8];//list of non-empty chidren
};

void shift_origin (struct atom_list_t *head, int num_atoms, float shiftx, float shifty, float shiftz);
void get_shift_distance(float *shiftx, float *shifty, float *shiftz);
void list_atoms(struct atom_list_t *atom_list);

void build_subcubes(struct octree_t* parent);//untill it is subdivided, the node is considered a leaf
void compute_subc_origin(struct point_t* base_orig, struct point_t* subc_orig, int oct_num, float side);
int get_subcube_index(struct point_t* base_orig, struct point_t* pointc, float side);
struct octree_t* construct_admissible_octree(struct atom_list_t *atom_list, int k, int alpha);
void expand_node(struct octree_t* node, int k, int alpha);
//void accuminter(struct octree_t* u, struct octree_t* v, int cutoff);
float get_neighbours(struct octree_t* u, struct octree_t* v, int cutoff);
bool is_within_cutoff (struct octree_t* u, struct octree_t* v, int cutoff);

#define DIST(A,B) (sqrt(pow((A.x - B.x), 2) + pow((A.y - B.y), 2) + pow((A.z - B.z),2)))
#define VDW(dist) (0.1948)*(pow((3.2/dist), 12) - 2*pow((3.2/dist),6))

int global_count=0;
int global_count2=0;
int num_leaves=0;
int num_leaves2=0;

void count_atoms(struct octree_t *parent);

///////////////////


 
 
/* frunction prototypes */
struct octree_t** createQueue(int *, int *);
void enQueue(struct octree_t **, int *, struct octree_t *);
struct octree_t *deQueue(struct octree_t **, int *);
 
/* Given a binary tree, print its nodes in level order
   using array for implementing queue */
void printLevelOrder(struct octree_t* root)
{
    int rear, front;
    int i=0, count=-1;
    struct octree_t **queue = createQueue(&front, &rear);
    struct octree_t *temp_node = root;

 
    while (temp_node)
    {

	temp_node->oct_num=++count;
        //printf("bfs:oct %d origin (%f %f %f) side %f num_atoms %d  leaf %d\n", temp_node->oct_num, temp_node->corigin.x, temp_node->corigin.y, temp_node->corigin.z, temp_node->side, temp_node->num_atoms, temp_node->is_leaf);
	i=0; 
	while(i < 8) {
		if (temp_node->octchild[i])
			enQueue(queue, &rear, temp_node->octchild[i]);
		i++;
	}
 
        /*Dequeue node and make it temp_node*/
        temp_node = deQueue(queue, &front);
    }
}
 
/*UTILITY FUNCTIONS*/
struct octree_t** createQueue(int *front, int *rear)
{
    struct octree_t **queue =
        (struct octree_t **)malloc(sizeof(struct octree_t*)*global_num_atoms);
 
    *front = *rear = 0;
    return queue;
}
 
void enQueue(struct octree_t **queue, int *rear, struct octree_t *new_node)
{
    queue[*rear] = new_node;
    (*rear)++;
}
 
struct octree_t *deQueue(struct octree_t **queue, int *front)
{
    (*front)++;
    return queue[*front - 1];
}
 
