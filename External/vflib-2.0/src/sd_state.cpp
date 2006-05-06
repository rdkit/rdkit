/*------------------------------------------------------------
 * sd_state.cc
 * Interface of sd_state.cc
 * Definition of a class representing a state of the matching
 * process between two ARGs using the Schmidt-Druffel isomorphism
 * algorithm described in:
 * - D. C. Schmidt, L. E. Druffel, "A Fast Backtracking Algorithm
 *   to Test Directed Graphs for Isomorphism Using Distance Matrices",
 *   Journal of ACM, Vol. 23, No. 3, pp. 433-445, 1973.
 *
 * See: argraph.h state.h sd_state.h
 *
 * Author: P. Foggia
 *-----------------------------------------------------------------*/

#define DEBUG_ALGO
#undef DEBUG_ALGO


#include <assert.h>
#include <stdlib.h>

#include "argraph.h"
#include "sd_state.h"

/*--------------------------------------------------------
 * Static prototypes
 -------------------------------------------------------*/
static node_id **create_distance_matrix(Graph *g);
static void compute_initial_partition(Graph *g1, Graph *g2, 
              node_id **dist1, node_id **dist2,
	      node_id *wrk1, node_id *wrk2,
              node_id *cls1, node_id *cls2);
static void count_classes(node_id *cls, int n, node_id *cnt);
static void compose_vectors(int n, node_id *a1, node_id *a2, 
           node_id *b1, node_id *b2, node_id *out1, node_id *out2);

/*----------------------------------------------------------
 * Methods of class SDState
 ---------------------------------------------------------*/

/*--------------------------------------------------------
 * This constructor builds the initial state, computing
 * the distance matrices and the initial partition
 -------------------------------------------------------*/

SDState::SDState(Graph *g1, Graph *g2)
  { assert(g1!=NULL);
    assert(g2!=NULL);

	parent=NULL;

	this->g1=g1;
	this->g2=g2;
	n1=g1->NodeCount();
	n2=g2->NodeCount();
	core_len = orig_core_len = 0;

	if (n1!=n2)
	  { dist1=dist2=NULL;
	    cls1=cls2=NULL;
		cnt1=cnt2=NULL;
		wrk1=wrk2=NULL;
		core1=core2=NULL;
		share_count=NULL;
		dead_end=true;
		return;
	  }

	share_count = new long;
	*share_count = 1;

	dist1=create_distance_matrix(g1);
	dist2=create_distance_matrix(g2);

	cls1=new node_id[n1];
	cls2=new node_id[n1];
	cnt1=new node_id[n1];
	cnt2=new node_id[n1];

	wrk1=new node_id[n1];
	wrk2=new node_id[n1];
	

	compute_initial_partition(g1, g2, dist1, dist2, wrk1, wrk2, cls1, cls2);

	count_classes(cls1, n1, cnt1);
	count_classes(cls2, n1, cnt2);

	dead_end=false;

	int i;
	for(i=0; i<n1; i++)
	  if (cnt1[i] != cnt2[i])
	    dead_end=true;

	core1=new node_id[n1];
	core2=new node_id[n1];
	for(i=0; i<n1; i++)
	  core1[i]=core2[i]=NULL_NODE;
  }

/*-----------------------------------------------------
 * A copy constructor, which builds a new state as a
 * child of the current one.
 ----------------------------------------------------*/
SDState::SDState(const SDState &state)
  { parent=(SDState *)&state;
    core_len=orig_core_len=state.core_len;
	g1=state.g1;
	g2=state.g2;
	n1=state.n1;
	n2=state.n2;
	dead_end=state.dead_end;
	share_count=state.share_count;
	++ *share_count;
	dist1=state.dist1;
	dist2=state.dist2;
	wrk1=state.wrk1;
	wrk2=state.wrk2;
	core1=state.core1;
	core2=state.core2;

	cls1=new node_id[n1];
	cls2=new node_id[n1];
	cnt1=new node_id[n1];
	cnt2=new node_id[n1];
    int i;
	for(i=0; i<n1; i++)
	    cls1[i]=state.cls1[i];
	for(i=0; i<n1; i++)
	    cls2[i]=state.cls2[i];
	for(i=0; i<n1; i++)
		cnt1[i]=state.cnt1[i];
	for(i=0; i<n1; i++)
		cnt2[i]=state.cnt2[i];
	   
  }


/*------------------------------------------------------------
 * Destructor. If the state is the topmost one, deallocate
 * also the shared structures.
 -----------------------------------------------------------*/
SDState::~SDState()
  { delete[] cls1;
    delete[] cls2;
    delete[] cnt1;
    delete[] cnt2;

	if (share_count==NULL)
	  return;

	if (-- *share_count == 0)
	  { delete[] wrk1;
	    delete[] wrk2;
		delete[] core1;
		delete[] core2;
		delete share_count;
		int i;
		for(i=0; i<n1; i++)
		  { delete[] dist1[i];
		    delete[] dist2[i];
		  }
		delete[] dist1;
		delete[] dist2;
	  }

  }

/*-----------------------------------------------------------------------
 * Find the next pair of nodes to be checked for addition to the
 * isomorphism.
 * Returns false if no other pair can be found.
 ----------------------------------------------------------------------*/
bool SDState::NextPair(node_id *pn1, node_id *pn2,
                       node_id prev_n1, node_id prev_n2)
  { if (prev_n1 == NULL_NODE)
      { prev_n1 = core_len;
        prev_n2 = 0;
      }
    else if (prev_n2 == NULL_NODE)
      { prev_n2 = 0;
      }
    else
      { prev_n2 ++;
      }

    assert(core1[prev_n1]==NULL_NODE);

    while (prev_n2 < n2 && 
          (cls1[prev_n1] != cls2[prev_n2] || core2[prev_n2] != NULL_NODE))
      prev_n2 ++;

    if (prev_n2 < n2)
      { *pn1 = prev_n1;
        *pn2 = prev_n2;
	return true;
      }
    else
      return false;
  }


/*---------------------------------------------------------
 * Checks if the addition of pair (node1, node2) to the
 * provisional mapping generates a new consistent
 * mapping. 
 * NOTE: For this algorithm this method only checks
 * semantic attribute compatibility.
 ---------------------------------------------------------*/
bool SDState::IsFeasiblePair(node_id node1, node_id node2)
  { if (!g1->CompatibleNode(g1->GetNodeAttr(node1),
                            g2->GetNodeAttr(node2)))
      return false;
    int i;
    node_id o1, o2;
    void *attr1, *attr2;


    for(i=0; i<g1->OutEdgeCount(node1); i++)
      { o1=g1->GetOutEdge(node1, i, &attr1);
        if ((o2=core1[o1]) != NULL_NODE)
	  { attr2=g2->GetEdgeAttr(node2, o2);
	    if (!g1->CompatibleEdge(attr1, attr2))
	      return false;
	  }
      }


    for(i=0; i<g2->OutEdgeCount(node2); i++)
      { o2=g2->GetOutEdge(node2, i, &attr2);
        if ((o1=core2[o2]) != NULL_NODE)
	  { attr1=g1->GetEdgeAttr(node1, o1);
	    if (!g1->CompatibleEdge(attr1, attr2))
	      return false;
	  }
      }


    for(i=0; i<g1->InEdgeCount(node1); i++)
      { o1=g1->GetInEdge(node1, i, &attr1);
        if ((o2=core1[o1]) != NULL_NODE)
	  { attr2=g2->GetEdgeAttr(o2, node2);
	    if (!g1->CompatibleEdge(attr1, attr2))
	      return false;
	  }
      }


    for(i=0; i<g2->InEdgeCount(node2); i++)
      { o2=g2->GetInEdge(node2, i, &attr2);
        if ((o1=core2[o2]) != NULL_NODE)
	  { attr1=g1->GetEdgeAttr(o1, node1);
	    if (!g1->CompatibleEdge(attr1, attr2))
	      return false;
	  }
      }


    return true;
  }

/*--------------------------------------------------------------
 * Adds a pair (node1, node2) to the current state, recomputing
 * the class vectors and the class count vectors.
 -------------------------------------------------------------*/
void SDState::AddPair(node_id node1, node_id node2)
  { assert(core_len == orig_core_len);

    #ifdef DEBUG_ALGO 
    printf("AddPair\n");
    #endif
  

    // At this point it can be safely assumed that cls1 is equal to
    // parent->cls1 and cls2 to parent->cls2; hence cls1 and cls2 are
    // initially used as working storage.
    int i;
    for(i=0; i<n1; i++)
      cls1[i]=dist1[i][node1];
    for(i=0; i<n2; i++)
      cls2[i]=dist2[i][node2];
    compose_vectors(n1, cls1, cls2, parent->cls1, parent->cls2,  wrk1, wrk2);
    compose_vectors(n1, wrk1, wrk2, dist1[node1], dist2[node2], cls1, cls2);
    count_classes(cls1, n1, cnt1);
    count_classes(cls2, n1, cnt2);

    for(i=0; i<n1; i++)
      if (cnt1[i] != cnt2[i])
        { dead_end=true;
	  break;
	}
    
    if (!dead_end)
      { core1[node1]=node2;
        core2[node2]=node1;
        core_len ++;
      }
  }

/*-----------------------------------------------------------
 * Puts in the output vectors c1 and c2 the partial mapping
 * found at the current state.
 -----------------------------------------------------------*/
void SDState::GetCoreSet(node_id c1[], node_id c2[])
  { int i;
    for(i=0; i<core_len; i++)
      { c1[i]=i;
        c2[i]=core1[i];
      }
  }

/*------------------------------------------------------
 * Creates a copy of the current state
 -----------------------------------------------------*/
State *SDState::Clone()
  { return new SDState(*this);
  }
  

/*-----------------------------------------------
 * Reverts the effect of the last AddPair on
 * shared data structures
 ----------------------------------------------*/
void SDState::BackTrack()
  { if (core_len > orig_core_len)
      { assert(core_len == orig_core_len+1);
        core2[core1[orig_core_len]] = NULL_NODE;
	core1[orig_core_len] = NULL_NODE;
      }
  }


/*---------------------------------------------------------
 * Static functions
 --------------------------------------------------------*/

/*--------------------------------------------------------
 * Allocates and computes the distance matrix using
 * Floyd's algorithm
 -------------------------------------------------------*/
static node_id **create_distance_matrix(Graph *g)
  { int i, j, k;
    int n=g->NodeCount();

    assert(n < NULL_NODE);

    node_id **d;
    d=new node_id *[n];
    for(i=0; i<n; i++)
      d[i]=new node_id[n];
    for(i=0; i<n; i++)
      { d[i][i]=0;
        for(j=0; j<n; j++)
          { if (i!=j)
	      { if (g->HasEdge(i,j))
	          d[i][j]=1;
		else
		  d[i][j]=NULL_NODE;
	       }
	   }
       }

    #ifdef DEBUG_ALGO
    printf("Adj. matrix:\n");
    for(i=0; i<n; i++)
      { printf("    ");
        for(j=0; j<n; j++) printf("%u \t", d[i][j]);
        printf("\n");
      }
    #endif

    bool changed;
    do {
    changed=false;
    for(i=0; i<n; i++)
      for(j=0; j<n; j++)
        for(k=0; k<n; k++)
	  { if ((long)d[i][j]+d[j][k] < d[i][k])
	      { d[i][k]=d[i][j]+d[j][k];
	        changed=true;
	      }
	  }
     } while (changed);

    #ifdef DEBUG_ALGO
    printf("Dist. matrix:\n");
    for(i=0; i<n; i++)
      { printf("    ");
        for(j=0; j<n; j++) printf("%u \t", d[i][j]);
        printf("\n");
      }
    #endif
	    
    return d;
  }

/*--------------------------------------------------------------------
 * Computes the initial partition of the nodes of the two graphs.
 * wrk1 and wrk2 are two vectors, with the same dimension
 * of cls1/cls2, used for temporary results.
 -------------------------------------------------------------------*/
static void compute_initial_partition(Graph *g1, Graph *g2, 
              node_id **dist1, node_id **dist2,
	      node_id *wrk1, node_id *wrk2, node_id *cls1, node_id *cls2)
   { int col;
     int n=g1->NodeCount();
     int i, j;
     node_id *tmp1=new node_id[n];
     node_id *tmp2=new node_id[n];

     #ifdef DEBUG_ALGO
     printf("compute initial partition\n");
     #endif

     for(i=0; i<n; i++)
       cls1[i]=cls2[i]=0;
     for(col=1; col<n; col++)
       { // Compute the column of the `row characteristic matrix'
         for(i=0; i<n; i++)
	   wrk1[i]=wrk2[i]=0;
         for(i=0; i<n; i++)
	   for(j=0; j<n; j++)
	     if (dist1[i][j]==col)
	       wrk1[i] ++;
         for(i=0; i<n; i++)
	   for(j=0; j<n; j++)
	     if (dist2[i][j]==col)
	       wrk2[i] ++;
        
	 compose_vectors(n, cls1, cls2, wrk1, wrk2, tmp1, tmp2);
         
         // Compute the column of the `column characteristic matrix'
         for(i=0; i<n; i++)
	   wrk1[i]=wrk2[i]=0;
         for(i=0; i<n; i++)
	   for(j=0; j<n; j++)
	     if (dist1[j][i]==col)
	       wrk1[i] ++;
         for(i=0; i<n; i++)
	   for(j=0; j<n; j++)
	     if (dist2[j][i]==col)
	       wrk2[i] ++;
        
	 compose_vectors(n, tmp1, tmp2, wrk1, wrk2, cls1, cls2);
         
       }
     delete[] tmp1;
     delete[] tmp2;
   }

/*------------------------------------------------------------
 * Count the occurrences of each class in the partition
 -----------------------------------------------------------*/
static void count_classes(node_id *cls, int n, node_id *cnt)
  { int i;
    for(i=0; i<n; i++)
      cnt[i]=0;
    for(i=0; i<n; i++)
      if (cls[i]<n)
        cnt[cls[i]]++;
  }


/*---------------------------------------------------------------
 * Declarations used to implement efficiently compose_vectors
 --------------------------------------------------------------*/
typedef struct
  { node_id id;
    node_id a;
    node_id b;
  }   sort_data;
typedef int (*compare_fn)(const void *, const void *);
static int compare(sort_data *x, sort_data *y);

/*-----------------------------------------------------------------------
 * Compose the vectors a1 with b1 giving out1 and a2 with b2 giving out2
 * The composition operation is defined by the following properties:
 * -  out1[i] >= 0 && out1[i]<n
 * -  out1[i] == out1[j] iff
 *         a1[i] == a1[j] && b1[i]==b1[j]
 * -  (out2[i] >= 0 && out2[i]<n) || out2[i]==NULL_NODE
 * -  out2[i] == out1[j] iff
 *         a2[i] == a1[j] && b2[i]==b1[j]
 ----------------------------------------------------------------------*/
 /*--------------------------------------------------------------
  * IMPLEMENTATION NOTE:
  * The algorithm is implemented by sorting the input vectors,
  * thus requiring O(n*log(n)) time versus the O(n**2) of the
  * simplistic approach. The vectors are NOT sorted in place 
  * (i.e. the arguments a1,b1,a2,b2 are not modified).
  * The function uses a heap allocated area for storing
  * information needed to sort the input vectors. This area
  * is pointed by a static local variable, so it is shared
  * by different invocations of the function.
  -------------------------------------------------------------*/
static void compose_vectors(int n, node_id *a1, node_id *a2, 
           node_id *b1, node_id *b2, node_id *out1, node_id *out2)
  { int i, j, cl;

    static int vec_size=0;
    static sort_data *vec1=NULL;
    static sort_data *vec2=NULL;

    /*
     * Allocates the auxiliary vectors for sorting the 
     * input information
     */
    if (n>vec_size)
      { delete []vec1;
        delete []vec2;
        vec_size=n;
        vec1=new sort_data[vec_size];
        vec2=new sort_data[vec_size];
      }
    

    #ifdef DEBUG_ALGO
    printf("\ncompose_vectors\n");
    printf("a1: ");
    for(i=0; i<n; i++) printf("%u ", a1[i]);
    printf("\n");
    printf("b1: ");
    for(i=0; i<n; i++) printf("%u ", b1[i]);
    printf("\n");
    printf("a2: ");
    for(i=0; i<n; i++) printf("%u ", a2[i]);
    printf("\n");
    printf("b2: ");
    for(i=0; i<n; i++) printf("%u ", b2[i]);
    printf("\n");
    #endif

    for(i=0; i<n; i++)
      { vec1[i].id=i;
        vec1[i].a=a1[i];
	vec1[i].b=b1[i];
        vec2[i].id=i;
	vec2[i].a=a2[i];
	vec2[i].b=b2[i];
      }

    qsort(vec1, n, sizeof(vec1[0]), (compare_fn)compare);
    qsort(vec2, n, sizeof(vec2[0]), (compare_fn)compare);

    for(i=0; i<n; i++)
      out1[i]=out2[i]=NULL_NODE;

    cl=0;
    i=0;
    j=0;
    while (i<n)
      { out1[vec1[i].id]=cl;
        while (i<n-1 && vec1[i+1].a == vec1[i].a && vec1[i+1].b==vec1[i].b)
	  { ++i;
	    out1[vec1[i].id]=cl;
	  }
	while (j<n && (vec2[j].a < vec1[i].a || 
	                (vec2[j].a == vec1[i].a && vec2[j].b < vec1[i].b)))
	  j++;
	while (j<n && vec2[j].a == vec1[i].a && vec2[j].b == vec1[i].b)
	  { out2[vec2[j].id] = cl;
	    j++;
	  }
	i++;
	cl++;
      }
	

    #ifdef DEBUG_ALGO
    printf("out1: ");
    for(i=0; i<n; i++) printf("%u ", out1[i]);
    printf("\n");
    printf("out2: ");
    for(i=0; i<n; i++) printf("%u ", out2[i]);
    printf("\n");
    #endif
  }

/*---------------------------------------------------
 * Lexicographic comparison for sort_data.
 * Used with qsort to sort the input vectors of
 * compose_vectors().
 --------------------------------------------------*/
static int compare(sort_data *x, sort_data *y)
  { if (x->a < y->a)
      return -1;
    else if (x->a > y->a)
      return +1;
    else if (x->b < y->b)
      return -1;
    else if (x->b > y->b)
      return +1;
    else 
      return 0;
  }
