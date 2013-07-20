#include <malloc.h>
#include <limits.h>
#include <matrix.h>
#define maxV 20
#define unseen 0

MAT_GRAPH mat_graph_creat(void)
{
    MAT_GRAPH g = NULL;

    if((g =(MAT_GRAPH)malloc(sizeof(mat_graph)))==NULL) graph_error(GRAPH_MALLOC);
    if((g->z = (G_NODE) malloc(sizeof(g_node)))==NULL) graph_error(GRAPH_MALLOC);
    g->z->next = g->z;
    g->weighted = 0;
    g->V = 0;
    g->E = 0;
    g->adj = NULL;
    g->id = 0;
    g->dad = NULL;
    return g;
}

void mat_graph_adjlist(MAT_GRAPH g, int directed, int weighted, FILEPOINTER fp)
{
    int j, x, y;
    mtype weight;
    G_NODE t;
    if(fscanf(fp, "%d", &(g->V))!=1) graph_error(GRAPH_READ);
    if(fscanf(fp, "%d", &(g->E))!=1) graph_error(GRAPH_READ);
    if(g->adj!=NULL) free(g->adj);
    if((g->adj = (G_NODE *) malloc(sizeof(G_NODE)*g->V+1))==NULL) graph_error(GRAPH_MALLOC);
    if((g->dad = (int *) malloc(sizeof(int)*(g->V+1)))==NULL) graph_error(GRAPH_MALLOC);
    if((g->val = (int *) malloc(sizeof(int)*(g->V+1)))==NULL) graph_error(GRAPH_MALLOC);
    if((g->vseq = (int *) malloc(sizeof(int)*(g->V+1)))==NULL) graph_error(GRAPH_MALLOC);
    for(j=1; j<=g->V; j++)g->adj[j]=g->z;
    for(j=1; j<=g->E; j++)
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
            t = (g_node *) malloc(sizeof(g_node));
            t->v = x;
            if(weighted)t->weight = weight;
            t->next = g->adj[y];
            g->adj[y] = t;
        }
        t = (g_node *) malloc(sizeof(g_node));
        t->v = y;
        if(weighted)t->weight = weight;
        t->next = g->adj[x];
        g->adj[x] = t;
    }
}

MAT_GRAPH mat_graph_reverse(MAT_GRAPH g, MAT_GRAPH r)
{
    int j;
    G_NODE s, t;
    if(r==NULL) r = mat_graph_creat();
    r->V = g->V;
    r->E = g->E;
    if(r->adj!=NULL)
    {
        free(r->adj);
        free(r->dad);
        free(r->val);
        free(r->vseq);
    }
    if((r->adj = (G_NODE *) malloc(sizeof(G_NODE)*r->V+1))==NULL) graph_error(GRAPH_MALLOC);
    if((r->dad = (int *) malloc(sizeof(int)*(r->V+1)))==NULL) graph_error(GRAPH_MALLOC);
    if((r->val = (int *) malloc(sizeof(int)*(r->V+1)))==NULL) graph_error(GRAPH_MALLOC);
    if((r->vseq = (int *) malloc(sizeof(int)*(r->V+1)))==NULL) graph_error(GRAPH_MALLOC);
    for(j=1; j<= r->V; j++)
    {
        r->adj[j] = r->z;
        for(s=g->adj[j]; s!=g->z; s=s->next)
        {
            t = (g_node *) malloc(sizeof(g_node));
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
    G_NODE t;

    m = MatRow(a);
    n = MatCol(a);
    if(g==NULL) g = mat_graph_creat();
    if(m!=n) graph_error(GRAPH_ELSE);
    g->V = m;

    if(g->adj!=NULL) free(g->adj);
    if((g->adj = (G_NODE *) malloc(sizeof(G_NODE)*g->V+1))==NULL) graph_error(GRAPH_MALLOC);
    if((g->dad = (int *) malloc(sizeof(int)*(g->V+1)))==NULL) graph_error(GRAPH_MALLOC);
    if((g->val = (int *) malloc(sizeof(int)*(g->V+1)))==NULL) graph_error(GRAPH_MALLOC);
    if((g->vseq = (int *) malloc(sizeof(int)*(g->V+1)))==NULL) graph_error(GRAPH_MALLOC);
    for(j=1; j<=g->V; j++)g->adj[j] = g->z;
    for(i=0; i<g->V; i++)
    {
        for(j=0; j<=g->V; j++)
        {
            if(a[i][j]>0)
            {
                x = i+1;
                y = j+1;
                weight = a[i][j];
                t = (g_node *) malloc(sizeof(g_node));
                t->v = y;
                t->weight = weight;
                t->next = g->adj[x];
                g->adj[x] = t;
                e++;
            }
        }
    }
    g->E = e;

}


INT_QUEUE mat_graph_search(MAT_GRAPH g, int connected, int mst)
{
    int k;
    INT_QUEUE q = int_queue_creat();
    INT_PRIORITYQUEUE pq = int_priorityqueue_creat();
    g->id = 0;
    for(k=1; k<=g->V; k++) g->val[k] = unseen;
    for(k=1; k<=g->V; k++)
        if(g->val[k]==unseen)
        {
            if(connected) int_queue_enqueue(q, INT_MIN);
            mat_graph_visit(g, k, connected, mst, pq, q);
        }
    int_priorityqueue_free(pq);
    return q;
}

void mat_graph_visit(MAT_GRAPH g, int k, int connected, int mst, INT_PRIORITYQUEUE pq, INT_QUEUE q)
{
    G_NODE t;
    if(mst)
    {
        if(int_priorityqueue_update(pq, k, unseen, PQ_INCREASE)!=0) g->dad[k] = 0;
        while(pq->p)
        {
            k = int_priorityqueue_dequeue(pq);
            g->val[k] = 1;
            if(connected) int_queue_enqueue(q, k);
            g->vseq[++g->id] = k;
            for(t=g->adj[k]; t!=g->z; t=t->next)
            {
                if(g->val[t->v]<=unseen)
                {
                    if(int_priorityqueue_update(pq, t->v, -(t->weight), PQ_INCREASE)>0)
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
        if(connected) int_queue_enqueue(q, k);
    }

}

/*
INT_STACK stk= int_stack_creat();
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

void mat_graph_dumpf(MAT_GRAPH g, int mst, FILEPOINTER fp)
{
    int k, l;
    if(g->vseq!=NULL)
    {
        for(k=1; k<=g->V; k++)
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
