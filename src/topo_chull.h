/*
This code is described in "Computational Geometry in C" (Second Edition),
Chapter 4.  It is not written to be comprehensible without the 
explanation in that book.

Input: 3n integer coordinates for the points.
Output: the 3D convex hull, in postscript with embedded comments
        showing the vertices and faces.

Compile: gcc -o chull chull.c (or simply: make)

Written by Joseph O'Rourke, with contributions by 
  Kristy Anderson, John Kutcher, Catherine Schevon, Susan Weller.
Last modified: May 2000
Questions to orourke@cs.smith.edu.

--------------------------------------------------------------------
This code is Copyright 2000 by Joseph O'Rourke.  It may be freely 
redistributed in its entirety provided that this copyright notice is 
not removed.
--------------------------------------------------------------------

--------------------------------------------------------------------
Adapted to work with apf data structures.

*/
#include <stdio.h>
#include <math.h>
#include <vector>
#include <cassert>

#include <apf.h>



/*Define Boolean type */
//typedef	enum { FALSE, TRUE }	bool;

/*====================================================================
    macros.h
 
 	macros used to access data structures and perform quick tests.

  ====================================================================*/

/* general-purpose macros */
#define SWAP(t,x,y)	{ t = x; x = y; y = t; }

char *malloc();

#define NEW(p,type)	if ((p=(type *) malloc (sizeof(type))) == NULL) {\
				printf ("Out of Memory!\n");\
				exit(0);\
			}

#define FREE(p)		if (p) { free ((char *) p); p = NULL; }


#define ADD( head, p )  if ( head )  { \
				p->next = head; \
				p->prev = head->prev; \
				head->prev = p; \
				p->prev->next = p; \
			} \
			else { \
				head = p; \
				head->next = head->prev = p; \
			}

#define DELETE( head, p ) if ( head )  { \
				if ( head == head->next ) \
					head = NULL;  \
				else if ( p == head ) \
					head = head->next; \
				p->next->prev = p->prev;  \
				p->prev->next = p->next;  \
				FREE( p ); \
			} 


/* Define vertex indices. */
#define X   0
#define Y   1
#define Z   2

/* Define structures for vertices, edges and faces */
typedef struct tVertexStructure tsVertex;
typedef tsVertex *tVertex;

typedef struct tEdgeStructure tsEdge;
typedef tsEdge *tEdge;

typedef struct tFaceStructure tsFace;
typedef tsFace *tFace;

struct tVertexStructure {
   double      v[3];
   int	    vnum;
   tEdge    duplicate;	        /* pointer to incident cone edge (or NULL) */
   bool     onhull;		/* T iff point on hull. */
   bool	    mark;		/* T iff point already processed. */
   tVertex  next, prev;
};

struct tEdgeStructure {
   tFace    adjface[2];
   tVertex  endpts[2];
   tFace    newface;            /* pointer to incident cone face. */
   bool     del_flag;		/* T iff edge should be delete. */
   tEdge    next, prev;
};

struct tFaceStructure {
   tEdge    edge[3];
   tVertex  vertex[3];
   bool	    visible;	        /* T iff face visible from new point. */
   tFace    next, prev;
};

/* Define flags */
#define ONHULL   	true
#define REMOVED  	true
#define VISIBLE  	true
#define PROCESSED	true
#define SAFE		100000000		/* Range of safe coord values. */


// ext_shell used to preserve volume by a naive approach.
class cvx_hull {
  private:
  public:

    /* Global variable definitions */
    tVertex vertices;
    tEdge edges;
    tFace faces;
    bool debug;
    bool check;

    /* Function declarations */
    tVertex MakeNullVertex( void );
    void    ReadVertices( std::vector<apf::Vector3>* vert_in );
    void    Print( void );
    void    SubVec( double a[3], double b[3], double c[3]);
    void    SubVec( double a[3], double b[3], apf::Vector3 &c);
    void    DoubleTriangle( void );
    void    ConstructHull( void );
    bool	AddOne( tVertex p );
    int     VolumeSign(tFace f, tVertex p);
    int 	Volumei( tFace f, tVertex p );
    tFace	MakeConeFace( tEdge e, tVertex p );
    void    MakeCcw( tFace f, tEdge e, tVertex p );
    tEdge   MakeNullEdge( void );
    tFace   MakeNullFace( void );
    tFace   MakeFace( tVertex v0, tVertex v1, tVertex v2, tFace f );
    void    CleanUp( tVertex *pvnext );
    void    CleanEdges( void );
    void    CleanFaces( void );
    void    CleanVertices( tVertex *pvnext );
    bool	Collinear( tVertex a, tVertex b, tVertex c );
    void    CheckEuler(int V, int E, int F );
    void	PrintPoint( tVertex p );
    void    Checks( void );
    void	Consistency( void );
    void	Convexity( void );
    void	PrintOut( tVertex v );
    void	PrintVertices( void );
    void	PrintEdges( void );
    void	PrintFaces( void );
    void	CheckEndpts ( void );
    void	EdgeOrderOnFaces ( void );

    std::vector<int> find_chull_vertices(std::vector<apf::Vector3>* vert_in);
    std::vector<int> get_vert_id();
    void get_face_id(std::vector<std::vector<int> > * f_id);
    void get_edge_id(std::vector<std::vector<int> > * e_id);
    cvx_hull();

    ~cvx_hull() {
      CleanUp(&vertices->prev);
    }
};

