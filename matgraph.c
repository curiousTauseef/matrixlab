#include "matrix.h"
#include <malloc.h>
#define maxV 20
#define unseen 0


MAT_GRAPH mat_graph_creat(void)
{
    MAT_GRAPH g = NULL;

    if((g =(MAT_GRAPH)malloc(sizeof(mat_graph)))==NULL) graph_error(GRAPH_MALLOC);
    if((g->z = (MAT_GNODE) malloc(sizeof(mat_gnode)))==NULL) graph_error(GRAPH_MALLOC);
    g->z->next = g->z;
    g->weighted = 0;
    g->nvertices = 0;
    g->nedges = 0;
    g->adj = NULL;
    g->id = 0;
    g->dad = NULL;
    return g;
}

void mat_graph_adjlist(MAT_GRAPH g, int directed, int weighted, MAT_FILEPOINTER fp)
{
    int j, x, y;
    mtype weight;
    MAT_GNODE t;
    if(fscanf(fp, "%d", &(g->nvertices))!=1) graph_error(GRAPH_READ);
    if(fscanf(fp, "%d", &(g->nedges))!=1) graph_error(GRAPH_READ);
    if(g->adj!=NULL) free(g->adj);
    if((g->adj = (MAT_GNODE *) malloc(sizeof(MAT_GNODE)*g->nvertices+1))==NULL) graph_error(GRAPH_MALLOC);
    if((g->dad = (int *) malloc(sizeof(int)*(g->nvertices+1)))==NULL) graph_error(GRAPH_MALLOC);
    if((g->val = (int *) malloc(sizeof(int)*(g->nvertices+1)))==NULL) graph_error(GRAPH_MALLOC);
    if((g->vseq = (int *) malloc(sizeof(int)*(g->nvertices+1)))==NULL) graph_error(GRAPH_MALLOC);
    for(j=1; j<=g->nvertices; j++)g->adj[j]=g->z;
    for(j=1; j<=g->nedges; j++)
    {
#if mtype_n == 0
        if(weighted)
        {
            if(fscanf(fp, "%d %d %f", &x, &y, &weight)!=3) graph_error(GRAPH_READ);
        }
#else
        if(weighted)
        {
            if(fscanf(fp, "%d %d %lf", &x, &y, &weight)!=3) graph_error(GRAPH_READ);
        }
#endif
        else if(fscanf(fp, "%d %d", &x, &y)!=2) graph_error(GRAPH_READ);
        if(directed==0)
        {
            t = (MAT_GNODE) malloc(sizeof(mat_gnode));
            t->v = x;
            if(weighted)t->weight = weight;
            t->next = g->adj[y];
            g->adj[y] = t;
        }
        t = (MAT_GNODE) malloc(sizeof(mat_gnode));
        t->v = y;
        if(weighted)t->weight = weight;
        t->next = g->adj[x];
        g->adj[x] = t;
    }
}

MAT_GRAPH mat_graph_reverse(MAT_GRAPH g, MAT_GRAPH r)
{
    int j;
    MAT_GNODE s, t;
    if(r==NULL) r = mat_graph_creat();
    r->nvertices = g->nvertices;
    r->nedges = g->nedges;
    if(r->adj!=NULL)
    {
        free(r->adj);
        free(r->dad);
        free(r->val);
        free(r->vseq);
    }
    if((r->adj = (MAT_GNODE *) malloc(sizeof(MAT_GNODE)*r->nvertices+1))==NULL) graph_error(GRAPH_MALLOC); /* SMH to check GNODE pointer */
    if((r->dad = (int *) malloc(sizeof(int)*(r->nvertices+1)))==NULL) graph_error(GRAPH_MALLOC);
    if((r->val = (int *) malloc(sizeof(int)*(r->nvertices+1)))==NULL) graph_error(GRAPH_MALLOC);
    if((r->vseq = (int *) malloc(sizeof(int)*(r->nvertices+1)))==NULL) graph_error(GRAPH_MALLOC);
    for(j=1; j<= r->nvertices; j++)
    {
        r->adj[j] = r->z;
        for(s=g->adj[j]; s!=g->z; s=s->next)
        {
            t = (MAT_GNODE) malloc(sizeof(mat_gnode));
            t->v = j;
            t->weight = s->weight;
            t->next = r->adj[s->v];
            r->adj[s->v] = t;
        }
    }
    return r;
}

void mat_graph_adjm_to_adjl(MAT_GRAPH g, MATRIX a)
{
    int m, n, i, j, x, y, e = 0;
    mtype weight;
    MAT_GNODE t;

    m = MatRow(a);
    n = MatCol(a);
    if(g==NULL) g = mat_graph_creat();
    if(m!=n) graph_error(GRAPH_ELSE);
    g->nvertices = m;

    if(g->adj!=NULL) free(g->adj);
    if((g->adj = (MAT_GNODE *) malloc(sizeof(MAT_GNODE)*g->nvertices+1))==NULL) graph_error(GRAPH_MALLOC);
    if((g->dad = (int *) malloc(sizeof(int)*(g->nvertices+1)))==NULL) graph_error(GRAPH_MALLOC);
    if((g->val = (int *) malloc(sizeof(int)*(g->nvertices+1)))==NULL) graph_error(GRAPH_MALLOC);
    if((g->vseq = (int *) malloc(sizeof(int)*(g->nvertices+1)))==NULL) graph_error(GRAPH_MALLOC);
    for(j=1; j<=g->nvertices; j++)g->adj[j] = g->z;
    for(i=0; i<g->nvertices; i++)
    {
        for(j=0; j<=g->nvertices; j++)
        {
            if(a[i][j]>0)
            {
                x = i+1;
                y = j+1;
                weight = a[i][j];
                t = (MAT_GNODE) malloc(sizeof(mat_gnode));
                t->v = y;
                t->weight = weight;
                t->next = g->adj[x];
                g->adj[x] = t;
                e++;
            }
        }
    }
    g->nedges = e;

}

MAT_INT_QUEUE mat_graph_search(MAT_GRAPH g, int connected, int mst)
{
    int k;
    MAT_INT_QUEUE q = mat_int_queue_creat();
    MAT_INT_PRIORITYQUEUE pq = mat_int_priorityqueue_creat(MAT_PQ_MAX);
    g->id = 0;
    for(k=1; k<=g->nvertices; k++) g->val[k] = unseen;
    for(k=1; k<=g->nvertices; k++)
        if(g->val[k]==unseen)
        {
            if(connected) mat_int_queue_enqueue(q, INT_MIN);
            mat_graph_visit(g, k, connected, mst, pq, q);
        }
    mat_int_priorityqueue_free(pq);
    return q;
}

void mat_graph_visit(MAT_GRAPH g, int k, int connected, int mst, MAT_INT_PRIORITYQUEUE pq, MAT_INT_QUEUE q)
{
    MAT_GNODE t;
    if(mst)
    {
        if(mat_int_priorityqueue_update(pq, k, unseen, PQ_INCREASE)!=0) g->dad[k] = 0;
        while(pq->p)
        {
            k = mat_int_priorityqueue_dequeue(pq).data;
            g->val[k] = 1;
            if(connected) mat_int_queue_enqueue(q, k);
            g->vseq[++g->id] = k;
            for(t=g->adj[k]; t!=g->z; t=t->next)
            {
                if(g->val[t->v]<=unseen)
                {
                    if(mat_int_priorityqueue_update(pq, t->v, -(t->weight), PQ_INCREASE)>0)
                    {
                        g->val[t->v] = -t->weight;
                        g->dad[t->v] = k;
                    }
                }
            }
        }
    }
    else
    {
        ++g->val[k];
        g->vseq[++g->id] = k;
        for(t=g->adj[k]; t!=g->z; t=t->next)
        {
            if(g->val[t->v]==unseen)mat_graph_visit(g, t->v, connected, mst, pq, q);
        }
        if(connected) mat_int_queue_enqueue(q, k);
    }

}

/*
MAT_INT_STACK stk= int_stack_creat();
void mat_graph_visit(int k, int connected)
{
    NODE t;
    int_stack_push(stk, k);
    while(stk->p!=0)
    {
        k =int_stack_pop(stk);
        ++val[k];
        vseq[++id] =k;
        for(t=adj[k]; t!=z; t= t->next)
        {
            if(val[t->v]==unseen)
            {
                int_stack_push(stk, t->v);
                ++val[t->v];
            };
        }
        if(connected) cout<<" "<<k;
    }
}*/

void mat_graph_dumpf(MAT_GRAPH g, int mst, MAT_FILEPOINTER fp)
{
    int k, l;
    if(g->vseq!=NULL)
    {
        for(k=1; k<=g->nvertices; k++)
        {
            if(mst)
            {
                l = g->vseq[k];
                fprintf(fp, "%d:", g->dad[l]);
            }
            fprintf(fp, "%d  ", g->vseq[k]);
        }
        fprintf(fp, "\n");
    }
}

void mat_graph_dump(MAT_GRAPH g, int mst)
{
    mat_graph_dumpf(g, mst, stdout);
}
