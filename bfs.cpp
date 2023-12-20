#include <stdlib.h>
#include <string.h>
#include "bfs.h"
#include <queue>
#include <limits>
using namespace std;

int get_neighbors(const Grid *grid, Point p, Point neighb[])
{   
    int nrVecini = 0;
    int dir[4][2] = {
           {-1, 0},
           {0, -1},
           {1, 0},
           {0, 1}
    };
 
    for (int i = 0; i < 4; i++) {
        int x, y;
        x = p.row + dir[i][0];
        y = p.col + dir[i][1];

        if (((x < grid->rows && y < grid->cols) && (x >= 0 && y >= 0)) && (grid->mat[x][y] != 1)) {
            neighb[nrVecini].row = x;
            neighb[nrVecini].col = y;
            nrVecini++;
        }
    }
    return nrVecini;
}

int get_cal(const Grid* grid, Point p, Point neighb[])
{
    int nrVecini = 0;
    int dir[8][2] = {
           {-1, -2},
           {-2, -1},
           {-2, 1},
           {-1, 2},
           {1, 2},
           {2, 1},
           {2, -1},
           {1, -2}
    };

    for (int i = 0; i < 8; i++) {
        int x, y;
        x = p.row + dir[i][0];
        y = p.col + dir[i][1];

        if (((x < grid->rows && y < grid->cols) && (x >= 0 && y >= 0)) && (grid->mat[x][y] != 1)) {
            neighb[nrVecini].row = x;
            neighb[nrVecini].col = y;
            nrVecini++;
        }
    }
    return nrVecini;
}

void grid_to_graph(const Grid *grid, Graph *graph)
{
    //we need to keep the nodes in a matrix, so we can easily refer to a position in the grid
    Node *nodes[MAX_ROWS][MAX_COLS];
    int i, j, k;
    Point neighb[4];

    //compute how many nodes we have and allocate each node
    graph->nrNodes = 0;
    for(i=0; i<grid->rows; ++i){
        for(j=0; j<grid->cols; ++j){
            if(grid->mat[i][j] == 0){
                nodes[i][j] = (Node*)malloc(sizeof(Node));
                memset(nodes[i][j], 0, sizeof(Node)); //initialize all fields with 0/NULL
                nodes[i][j]->position.row = i;
                nodes[i][j]->position.col = j;
                ++graph->nrNodes;
            }else{
                nodes[i][j] = NULL;
            }
        }
    }
    graph->v = (Node**)malloc(graph->nrNodes * sizeof(Node*));
    k = 0;
    for(i=0; i<grid->rows; ++i){
        for(j=0; j<grid->cols; ++j){
            if(nodes[i][j] != NULL){
                graph->v[k++] = nodes[i][j];
            }
        }
    }

    //compute the adjacency list for each node
    for(i=0; i<graph->nrNodes; ++i){
        graph->v[i]->adjSize = get_neighbors(grid, graph->v[i]->position, neighb);
        if(graph->v[i]->adjSize != 0){
            graph->v[i]->adj = (Node**)malloc(graph->v[i]->adjSize * sizeof(Node*));
            k = 0;
            for(j=0; j<graph->v[i]->adjSize; ++j){
                if( neighb[j].row >= 0 && neighb[j].row < grid->rows &&
                    neighb[j].col >= 0 && neighb[j].col < grid->cols &&
                    grid->mat[neighb[j].row][neighb[j].col] == 0){
                        graph->v[i]->adj[k++] = nodes[neighb[j].row][neighb[j].col];
                }
            }
            if(k < graph->v[i]->adjSize){
                //get_neighbors returned some invalid neighbors
                graph->v[i]->adjSize = k;
                graph->v[i]->adj = (Node**)realloc(graph->v[i]->adj, k * sizeof(Node*));
            }
        }
    }
}

void free_graph(Graph *graph)
{
    if(graph->v != NULL){
        for(int i=0; i<graph->nrNodes; ++i){
            if(graph->v[i] != NULL){
                if(graph->v[i]->adj != NULL){
                    free(graph->v[i]->adj);
                    graph->v[i]->adj = NULL;
                }
                graph->v[i]->adjSize = 0;
                free(graph->v[i]);
                graph->v[i] = NULL;
            }
        }
        free(graph->v);
        graph->v = NULL;
    }
    graph->nrNodes = 0;
}

void bfs(Graph* graph, Node* s, Operation* op)
{
    queue<Node*> q;
    for (int i = 0; i < graph->nrNodes; i++) {
        if (op != NULL) {
            op->count();
        }
        if (graph->v[i] != s) {
            if (op != NULL) {
                op->count(3);
            }
            graph->v[i]->color = COLOR_WHITE;
            graph->v[i]->dist = INT_MAX;
            graph->v[i]->parent = NULL;
        }
    }
    if (op != NULL) {
        op->count(3);
    }
    s->color = COLOR_GRAY;
    s->dist = 0;
    s->parent = NULL;
    q.push(s);
    while (!(q.empty())) {
        if (op != NULL) {
            op->count();
        }
        Node* t = q.front();
        q.pop();
        for (int i = 0; i < t->adjSize; i++) {
            if (op != NULL) {
                op->count();
            }
            if (t->adj[i]->color == COLOR_WHITE) {
                if (op != NULL) {
                    op->count(3);
                }
                t->adj[i]->color = COLOR_GRAY;
                t->adj[i]->dist = t->dist + 1;
                t->adj[i]->parent = t;
                q.push(t->adj[i]);
            }
        }
        if (op != NULL) {
            op->count();
        }
        t->color = COLOR_BLACK;
    }
    // TOOD: implement the BFS algorithm on the graph, starting from the node s
    // at the end of the algorithm, every node reachable from s should have the color BLACK
    // for all the visited nodes, the minimum distance from s (dist) and the parent in the BFS tree should be set
    // for counting the number of operations, the optional op parameter is received
    // since op can be NULL (when we are calling the bfs for display purposes), you should check it before counting:
    // if(op != NULL) op->count();
}

void prettyPrint(int array[], Point* k, int size, int index, int count) {

    for (int i = 0; i < size; i++) {
        if (array[i] == index) {
            for (int j = 0; j < count; j++) {
                printf("      ");
            }
            printf("(%d, %d)\n", k[i].row, k[i].col);
            prettyPrint(array, k, size, i, count + 1);

        }
    }
}


void print_bfs_tree(Graph *graph)
{
    //first, we will represent the BFS tree as a parent array
    int n = 0; //the number of nodes
    int *p = NULL; //the parent array
    Point *repr = NULL; //the representation for each element in p

    //some of the nodes in graph->v may not have been reached by BFS
    //p and repr will contain only the reachable nodes
    int *transf = (int*)malloc(graph->nrNodes * sizeof(int));
    for(int i=0; i<graph->nrNodes; ++i){
        if(graph->v[i]->color == COLOR_BLACK){
            transf[i] = n;
            ++n;
        }else{
            transf[i] = -1;
        }
    }
    if(n == 0){
        //no BFS tree
        free(transf);
        return;
    }

    int err = 0;
    p = (int*)malloc(n * sizeof(int));
    repr = (Point*)malloc(n * sizeof(Node));
    for(int i=0; i<graph->nrNodes && !err; ++i){
        if(graph->v[i]->color == COLOR_BLACK){
            if(transf[i] < 0 || transf[i] >= n){
                err = 1;
            }else{
                repr[transf[i]] = graph->v[i]->position;
                if(graph->v[i]->parent == NULL){
                    p[transf[i]] = -1;
                }else{
                    err = 1;
                    for(int j=0; j<graph->nrNodes; ++j){
                        if(graph->v[i]->parent == graph->v[j]){
                            if(transf[j] >= 0 && transf[j] < n){
                                p[transf[i]] = transf[j];
                                err = 0;
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
    free(transf);
    transf = NULL;

    if(!err){

        prettyPrint(p, repr, n, -1, 0);
        // TODO: pretty print the BFS tree
        // the parrent array is p (p[k] is the parent for node k or -1 if k is the root)
        // when printing the node k, print repr[k] (it contains the row and column for that point)
        // you can adapt the code for transforming and printing multi-way trees from the previous labs
    }

    if(p != NULL){
        free(p);
        p = NULL;
    }
    if(repr != NULL){
        free(repr);
        repr = NULL;
    }
}

int shortest_path(Graph *graph, Node *start, Node *end, Node *path[])
{
    int nrNoduri = 0;
    bfs(graph, start); 
    Node* temp = end;
    while (temp != NULL) {
        path[nrNoduri] = temp;
        temp = temp->parent;
        nrNoduri++; 
    }

    for (int i = 0; i < nrNoduri / 2; i++) {
        swap(path[i], path[nrNoduri - i - 1]);
    } 
    // TODO: compute the shortest path between the nodes start and end in the given graph
    // the nodes from the path, should be filled, in order, in the array path
    // the number of nodes filled in the path array should be returned
    // if end is not reachable from start, return -1
    // note: the size of the array path is guaranteed to be at least 1000

    return nrNoduri;
}

void generateRandomEdges(Graph* f,int noEdges) {

    int nrMuchii = 0;
    int** matrice = new int* [f->nrNodes];
    for (int i = 0; i < f->nrNodes; i++) {
        matrice[i] = new int[f->nrNodes];
        for (int j = 0; j < f->nrNodes; j++) {
            matrice[i][j] = 0;
        }
    }

    while (nrMuchii != noEdges) {

        int emit = rand() % f->nrNodes;
        int prim = rand() % f->nrNodes;
        if (emit != prim && matrice[emit][prim] == 0) {

            f->v[emit]->adjSize++;
            f->v[emit]->adj = (Node**)realloc(f->v[emit]->adj, f->v[emit]->adjSize * sizeof(Node*));
            f->v[emit]->adj[f->v[emit]->adjSize - 1] = f->v[prim];

            f->v[prim]->adjSize++;
            f->v[prim]->adj = (Node**)realloc(f->v[prim]->adj, f->v[prim]->adjSize * sizeof(Node*));
            f->v[prim]->adj[f->v[prim]->adjSize - 1] = f->v[emit];
            matrice[emit][prim] = 1;
            matrice[prim][emit] = 1;
            nrMuchii++;
        }
    }
}

void performance()
{
    int n, i;
    Profiler p("bfs");

    // vary the number of edges
    for(n=1000; n<=4500; n+=100){
        Operation op = p.createOperation("bfs-edges", n);
        Graph graph;
        graph.nrNodes = 100;
        //initialize the nodes of the graph
        graph.v = (Node**)malloc(graph.nrNodes * sizeof(Node*));
        for(i=0; i<graph.nrNodes; ++i){
            graph.v[i] = (Node*)malloc(sizeof(Node));
            memset(graph.v[i], 0, sizeof(Node));
        }
       generateRandomEdges(&graph, n);
        // TODO: generate n random edges
        // make sure the generated graph is connected

        bfs(&graph, graph.v[0], &op);
        free_graph(&graph);
    }

    // vary the number of vertices
    for(n=100; n<=200; n+=10){
        Operation op = p.createOperation("bfs-vertices", n); 
        Graph graph;
        graph.nrNodes = n;
        //initialize the nodes of the graph
        graph.v = (Node**)malloc(graph.nrNodes * sizeof(Node*));
        for(i=0; i<graph.nrNodes; ++i){
            graph.v[i] = (Node*)malloc(sizeof(Node));
            memset(graph.v[i], 0, sizeof(Node));
        }

        generateRandomEdges(&graph, 4500);
        // TODO: generate 4500 random edges
        // make sure the generated graph is connected

        bfs(&graph, graph.v[0], &op);
        free_graph(&graph);
    }

    p.showReport();
}
