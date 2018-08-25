/*
	Parallel code of:
	Barnes-Hut algorithm to calculate forces.	
	This algorithm approximates distant particles as a single point particle
	with position = center of mass of the constinuent particles
	and mass = total mass of these particles.
	Thus it will compute the approximate force on each particle in O(log n) time
	and thus for n particles it will take O(n*log n) time.

	Authors:
		Nisarg S Patel	: 201501134	
		Aagam P Mehta	: 201501133
*/
#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>	//To use bool data type.
#include<math.h>
#include "omp.h"	//OpenMP Directive

struct Point		//Defines a point
{
    float x;
    float y;
};

struct Node		//Defines a node(particle) which will have position as a point and mass as data. 
{
    struct Point pos;
    float data;
    double force_x;
    double force_y;
};

/*
	Barnes Hut algorithm will have 3 steps:
	1) Build a QuadTree.
	2) Traverse the QuadTree from the leaves to the root, computing centre mass and total mass for each parent Tree.
	3) For each particle, traverse the tree from the root, calculate the force.
*/

struct Quad
{
	// Hold details of the boundary of this Tree
	struct Point topLeft;
    struct Point botRight;
    
	/*Data for this tree
		Center of mass in both x and y directions
		Total mass inside this tree section
		Number of particles inside this tree section	
	*/
    float center_of_mass_x;
    float center_of_mass_y;
    float total_mass;
    float number_of_particles;
 
    // Contains details of node(particle) if there is just one particle
    struct Node *n;
 
    // Children of this tree
    struct Quad *topLeftTree;
    struct Quad *topRightTree;
    struct Quad *botLeftTree;
    struct Quad *botRightTree;
};

void initialize_Quad(struct Quad*,struct Point,struct Point);	//Funtion to initialize the QuadTree
void insert(struct Quad*, struct Node*);						//Function to insert a node(particle) in the quad tree
void initialize_Point(struct Point, float, float);				//Function to initialize a point
void initialize_Node(struct Node*, struct Point, float);		//Function to initialize a node
struct Node* search(struct Quad*, struct Point);				//Search function of a point in the Quadtree
bool inBoundary(struct Quad*, struct Point);					//To check if the point lies inside the Tree boundries
void calc_center_of_mass(struct Quad*);							//To calculate center of mass of each Quadtree and its respective subtrees
void calc_force(struct Node**, struct Quad*, int);				//To calculate force on a system of particles	
void calc_node_force(struct Node*, struct Quad*);				//To calculate total force on each particle
void force_of_2_body(struct Node*, struct Quad*);				//To calculate the force on a particle due to another particle

void initialize_Quad(struct Quad* Quadtree, struct Point tL, struct Point bR)	//Funtion to initialize the QuadTree
{
	/*
		tL is the topLeft point that this tree can represent.
		bL is the bottomLeft point that this tree can represent.
		Thus we will have a square quadtree.
		The side of this quadtree should be a power of 2. 
	*/
	Quadtree->topLeft=tL;
    Quadtree->botRight=bR;
    Quadtree->n=NULL;
	Quadtree->topLeftTree  = NULL;
    Quadtree->topRightTree = NULL;
    Quadtree->botLeftTree  = NULL;
    Quadtree->botRightTree = NULL;
    Quadtree->center_of_mass_x = 0;
    Quadtree->center_of_mass_y = 0;
    Quadtree->total_mass = 0;
    Quadtree->number_of_particles = 0;
}

void insert(struct Quad* Quadtree, struct Node* node)	//Function to insert a node(particle) in the quad tree
{
	//If the node(particle) is  NULL 
	if (node == NULL)
        return;
 
    // Current quadtree cannot contain it
    if (!inBoundary(Quadtree, node->pos))
        return;
 
    // We are at a quadtree of unit area
    // We cannot subdivide this quad further
    if (abs(Quadtree->topLeft.x - Quadtree->botRight.x) <= 1 &&
        abs(Quadtree->topLeft.y - Quadtree->botRight.y) <= 1)
    {
        if (Quadtree->n == NULL)
            Quadtree->n = node;
            //printf("Point added with value %f\n", Quadtree->n->data);
        return;
    }
 
    if ((Quadtree->topLeft.x + Quadtree->botRight.x) / 2 >= node->pos.x)
    {
        // Indicates topLeftTree
        if ((Quadtree->topLeft.y + Quadtree->botRight.y) / 2 >= node->pos.y)
        {
            if (Quadtree->topLeftTree == NULL)
            {
				Quadtree->topLeftTree = (struct Quad*)malloc(sizeof(struct Quad));
                struct Point tL = {Quadtree->topLeft.x, Quadtree->topLeft.y};
                struct Point bR = {(Quadtree->topLeft.x+Quadtree->botRight.x)/2,(Quadtree->topLeft.y+Quadtree->botRight.y)/2};
				initialize_Quad(Quadtree->topLeftTree,tL,bR);
        	}
			//printf("In topLeftTree\n");    
            insert(Quadtree->topLeftTree,node);
        }
 
        // Indicates botLeftTree
        else
        {
            if (Quadtree->botLeftTree == NULL)
            {
            	Quadtree->botLeftTree = (struct Quad*)malloc(sizeof(struct Quad));
                struct Point tL = {Quadtree->topLeft.x, (Quadtree->topLeft.y+Quadtree->botRight.y)/2};
                struct Point bR = {(Quadtree->topLeft.x+Quadtree->botRight.x)/2,Quadtree->botRight.y};
                initialize_Quad(Quadtree->botLeftTree,tL,bR);
			}
			//printf("In botLeftTree\n");
            insert(Quadtree->botLeftTree,node);
        }
    }
    else
    {
        // Indicates topRightTree
        if ((Quadtree->topLeft.y + Quadtree->botRight.y) / 2 >= node->pos.y)
        {
            if (Quadtree->topRightTree == NULL)
            {
            	Quadtree->topRightTree = (struct Quad*)malloc(sizeof(struct Quad));
                struct Point tL = {(Quadtree->topLeft.x+Quadtree->botRight.x)/2, Quadtree->topLeft.y};
                struct Point bR = {Quadtree->botRight.x,(Quadtree->topLeft.y+Quadtree->botRight.y)/2};
                initialize_Quad(Quadtree->topRightTree,tL,bR);            	
			}
			//printf("In topRightTree\n");
            insert(Quadtree->topRightTree,node);
        }
 
        // Indicates botRightTree
        else
        {
            if (Quadtree->botRightTree == NULL)
            {
            	Quadtree->botRightTree = (struct Quad*)malloc(sizeof(struct Quad));
                struct Point tL = {(Quadtree->topLeft.x+Quadtree->botRight.x)/2, (Quadtree->topLeft.y+Quadtree->botRight.y)/2};
                struct Point bR = {Quadtree->botRight.x,Quadtree->botRight.y};
                initialize_Quad(Quadtree->botRightTree,tL,bR);
			}
			//printf("In botRightTree\n");
			insert(Quadtree->botRightTree,node);
        }
    }
}

void initialize_Point(struct Point P, float x1 , float y1)		//Function to initialize a point
{
	printf("%f%f\n",x1,y1);
	P.x=x1;
	P.y=y1;
}

void initialize_Node(struct Node* node, struct Point P, float value)	//Function to initialize a node
{
	node->pos=P;
	node->data=value;
	node->force_x=0;
	node->force_y=0;
}

struct Node* search(struct Quad* Quadtree, struct Point p)		//Search function for a Node in the Quadtree
{
    // Current quadtree cannot contain it
    if (!inBoundary(Quadtree, p))
        return NULL;
 	
	/*
		If current quadtree contains a node with this point
		The required node is found
	*/
    if (Quadtree->n != NULL)
        return Quadtree->n;
 
	//Else continue traversing in its one of the subtrees
	
    if ((Quadtree->topLeft.x + Quadtree->botRight.x) / 2 >= p.x)
    {
        // Indicates topLeftTree
        if ((Quadtree->topLeft.y + Quadtree->botRight.y) / 2 >= p.y)
        {
            if (Quadtree->topLeftTree == NULL)
                return NULL;
            //printf("In topLeftTree\n");
            return search(Quadtree->topLeftTree, p);
        }
 
        // Indicates botLeftTree
        else
        {
            if (Quadtree->botLeftTree == NULL)
                return NULL;
           // printf("In botLeftTree\n");
            return search(Quadtree->botLeftTree, p);
        }
    }
    else
    {
        // Indicates topRightTree
        if ((Quadtree->topLeft.y + Quadtree->botRight.y) / 2 >= p.y)
        {
            if (Quadtree->topRightTree == NULL)
                return NULL;
           // printf("In topRightTree\n");
            return search(Quadtree->topRightTree, p);
        }
 
        // Indicates botRightTree
        else
        {
            if (Quadtree->botRightTree == NULL)
                return NULL;
            //printf("In botRightTree\n");
            return search(Quadtree->botRightTree, p);
        }
    }
}

bool inBoundary(struct Quad* Quadtree, struct Point p)	//To check if the point lies inside the Tree boundries
{
    return (p.x >= Quadtree->topLeft.x &&
        p.x <= Quadtree->botRight.x &&
        p.y >= Quadtree->topLeft.y &&
        p.y <= Quadtree->botRight.y);
}

void calc_center_of_mass(struct Quad* Quadtree)		//To calculate center of mass of each Quadtree and its respective subtrees
{
	/*
		If the Tree contains 1 partical:
	*/
	if(Quadtree->n!=NULL)
	{
		Quadtree->center_of_mass_x=Quadtree->n->pos.x;
		Quadtree->center_of_mass_y=Quadtree->n->pos.y;
		Quadtree->total_mass=Quadtree->n->data;
		Quadtree->number_of_particles=1;
		//printf("A Node found of value %d\n", Quadtree->n->data);
	}

	/*
		If any of its subtree contains particle:
	*/	
	else
	{
		if(Quadtree->topLeftTree!=NULL)
		{
			calc_center_of_mass(Quadtree->topLeftTree);
			Quadtree->total_mass+=(float)Quadtree->topLeftTree->total_mass;
			Quadtree->center_of_mass_x+=(float)(Quadtree->topLeftTree->center_of_mass_x*Quadtree->topLeftTree->total_mass);
			Quadtree->center_of_mass_y+=(float)(Quadtree->topLeftTree->center_of_mass_y*Quadtree->topLeftTree->total_mass);
			Quadtree->number_of_particles+=(float)Quadtree->topLeftTree->number_of_particles;
		}
		if(Quadtree->botLeftTree!=NULL)
		{
			calc_center_of_mass(Quadtree->botLeftTree);
			Quadtree->total_mass+=(float)Quadtree->botLeftTree->total_mass;
			Quadtree->center_of_mass_x+=(float)(Quadtree->botLeftTree->center_of_mass_x*Quadtree->botLeftTree->total_mass);
			Quadtree->center_of_mass_y+=(float)(Quadtree->botLeftTree->center_of_mass_y*Quadtree->botLeftTree->total_mass);
			Quadtree->number_of_particles+=(float)Quadtree->botLeftTree->number_of_particles;
		}
		if(Quadtree->topRightTree!=NULL)
		{
			calc_center_of_mass(Quadtree->topRightTree);
			Quadtree->total_mass+=(float)Quadtree->topRightTree->total_mass;
			Quadtree->center_of_mass_x+=(float)(Quadtree->topRightTree->center_of_mass_x*Quadtree->topRightTree->total_mass);
			Quadtree->center_of_mass_y+=(float)(Quadtree->topRightTree->center_of_mass_y*Quadtree->topRightTree->total_mass);
			Quadtree->number_of_particles+=(float)Quadtree->topRightTree->number_of_particles;
		}
		if(Quadtree->botRightTree!=NULL)
		{
			calc_center_of_mass(Quadtree->botRightTree);
			Quadtree->total_mass+=(float)Quadtree->botRightTree->total_mass;
			Quadtree->center_of_mass_x+=(float)(Quadtree->botRightTree->center_of_mass_x*Quadtree->botRightTree->total_mass);
			Quadtree->center_of_mass_y+=(float)(Quadtree->botRightTree->center_of_mass_y*Quadtree->botRightTree->total_mass);
			Quadtree->number_of_particles+=(float)Quadtree->botRightTree->number_of_particles;
		}
		Quadtree->center_of_mass_x=Quadtree->center_of_mass_x/Quadtree->total_mass;
		Quadtree->center_of_mass_y=Quadtree->center_of_mass_y/Quadtree->total_mass;
	}
	
	/*printf statements to get the values for each Quadtrees*/
	//printf("This section has \n");
	//printf("Center of Mass x : %f\n",Quadtree->center_of_mass_x);
	//printf("Center of Mass y : %f\n",Quadtree->center_of_mass_y);
	//printf("Total Mass : %f\n",Quadtree->total_mass);
	//printf("Total Particles : %f\n\n",Quadtree->number_of_particles);
}

void calc_force(struct Node** node, struct Quad* Quadtree, int num)		//To calculate force on a system of particles
{
	int i=0;

	double start=omp_get_wtime();	
	#pragma omp parallel for schedule(guided)
	for(i=0;i<num;i++)
	{
		calc_node_force(node[i], Quadtree);
		/*To calculate the accurate time, comment these printf statements*/
		printf("Force x for node(%0.2f,%0.2f) of mass %0.2f : %0.10lf *G N\n",node[i]->pos.x,node[i]->pos.y, node[i]->data, node[i]->force_x);
		printf("Force y for node(%0.2f,%0.2f) of mass %0.2f : %0.10lf *G N\n\n",node[i]->pos.x,node[i]->pos.y, node[i]->data, node[i]->force_y);
	}
	double end2 = omp_get_wtime();
	double time = end2-start;
	printf("\nTime : %0.10lf\n", time);
}

void calc_node_force(struct Node* node, struct Quad* Quadtree)		//To calculate total force on each particle
{	
	/*
		If the Quadtree contains just one particle, then directly calculate the force between these two particles
	*/	
	if(Quadtree->number_of_particles==1)
	{
		force_of_2_body(node, Quadtree);
	}
	
	/*
		If the Quadtree contains more than 1 particles,
		Then check if the Quadtree is at a relatively larger distance from the node
		A measure to check if the Quadtree is far or near is by introducing a parameter theta

			theta = d/r
			
		where d=the size of the respective quadtree and r is the distance between the node and the center of mass of the quadtree
		
		if d/r is smaller than theta then Quadtree is at a relatively larger distance and can be considered as a point mass.
		if d/r is greater than theta then Quadtree is relatively near
			Thus recurse the same procedure of finding the force on node with all f its child-quadtrees.		
	*/
	else
	{
		double x1=(double)node->pos.x;
		double y1=(double)node->pos.y;
		double x2=(double)Quadtree->center_of_mass_x;
		double y2=(double)Quadtree->center_of_mass_y;
		double r=sqrt(pow(x2-x1,2)+pow(y2-y1,2));
		double d=Quadtree->botRight.x-Quadtree->topLeft.x;
		double theta = 0.5;
		if(d/r<theta)
		{
			force_of_2_body(node, Quadtree);
		}
		else
		{
			if(Quadtree->topLeftTree!=NULL)
				calc_node_force(node, Quadtree->topLeftTree);
			if(Quadtree->botLeftTree!=NULL)
				calc_node_force(node, Quadtree->botLeftTree);
			if(Quadtree->topRightTree!=NULL)
				calc_node_force(node, Quadtree->topRightTree);
			if(Quadtree->botRightTree!=NULL)
				calc_node_force(node, Quadtree->botRightTree);
		}
	}
}

void force_of_2_body(struct Node* node, struct Quad* Quadtree)	//To calculate the force on a particle due to another particle
{

	/*
		Considering the Quadtree at relatively larger distance as a particle
		with position: center of mass
		and mass: Total Mass
		or Quadtree with one node(particle)
		For each particle i(Xi,Yi,Mi) with respect to particle j(Xj,Yj,Mj)
		Force is calculated as:

			Fx = G*Mi*Mj*(Xj-Xi)/(R*R*R)
			Fy = G*Mi*Mj*(Yj-Yi)/(R*R*R)
			
		Where, R=sqrt((Xj-Xi)^2+(Yj-Yi)^2)
		
		Total Force is vector sum of all the forces		
		
	*/

	double G = 1;//6.67*pow(10,-6);
	double x1=(double)node->pos.x;
	double y1=(double)node->pos.y;
	double x2=(double)Quadtree->center_of_mass_x;
	double y2=(double)Quadtree->center_of_mass_y;
	double m1=(double)node->data;
	double m2=(double)Quadtree->total_mass;
	double r=sqrt(pow(x2-x1,2)+pow(y2-y1,2));
	if(r==0)
	{
		return;
	}
	node->force_x += (x2-x1)*G*m1*m2/(r*r*r);
	node->force_y += (y2-y1)*G*m1*m2/(r*r*r);
}

int main()
{
	float height=pow(2,15);		//size of the quadtree
    struct Quad* center = (struct Quad*)malloc(sizeof(struct Quad));
    struct Point tL = {0,0};
    struct Point bR = {height,height};
    initialize_Quad(center,tL,bR);		//Initialize the quadtree
    
	
    int num,i,j,count,k;
	printf("Enter Number of Particles(max 10^9) : ");
	scanf("%d",&num);
	
	printf("Enter Number of Cores : ");
	scanf("%d",&k);
	
	omp_set_num_threads(10);
	
	struct Point* P=(struct Point*)malloc(num* sizeof(struct Point)); //Dynamic allocation of memory to struct Point
	count=0;

	if(num<(height*(height+1)/2))	//Initialize random points
	{    
		for(i=1;;i++)
		{
			for(j=1;j<=i;j++)
			{
				float x=(float)i;
				float y=(float)j;
				P[count].x=x;
				P[count].y=y;
				count++;
				if(count==num)
				{
					break;
				}
			}
			if(count==num)
			{
				break;
			} 
		}
	}
	else
	{
		for(i=1;;i++)
		{
			for(j=1;j<=height;j++)
			{
				float x=(float)i;
				float y=(float)j;
				P[count].x=x;
				P[count].y=y;
				count++;
				if(count==num)
				{
					break;
				}
			}
			if(count==num)
			{
				break;
			} 
		}
	}

	fflush(stdout);		//To clear the output buffer
	fflush(stdin);		//To clear the input buffer

	struct Node** a= (struct Node**)malloc(num*sizeof(struct Node));	

	for(i=0;i<num;i++)
	{
		a[i] = (struct Node*)malloc(sizeof(struct Node));
		initialize_Node(a[i],P[i],2);
	}
	for(i=0;i<num;i++)	//Insert the nodes(particles) in the quadtree
    {
    	insert(center, a[i]);
	}
	
	calc_center_of_mass(center);	//Calculate centre of mass for each Quadtrees
	
	/*These printf statements outputs the values of the overall space(Main Quadtree)
	/*	
	printf("\n\n\n");
	printf("This section has \n");
	printf("Center of Mass x : %f\n",center->center_of_mass_x);
	printf("Center of Mass y : %f\n",center->center_of_mass_y);
	printf("Total Mass : %f\n",center->total_mass);
	printf("Total Particles : %f\n\n",center->number_of_particles);
	*/
	
	calc_force(a,center,num);		//To calculate the force on all the particles

    return 0;
}
