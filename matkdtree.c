#include "matrix.h"
#include <string.h>
#include <malloc.h>

__inline mtype __kd_dist(MAT_KDNODE a, MAT_KDNODE b, int dim)
{
    mtype t, d = 0.0;
    while(--dim>=0)
    {
        t = a->x[dim]-b->x[dim];
        d += t*t;
    }
    return d;
}

__inline void __kd_swap(MAT_KDNODE x, MAT_KDNODE y)
{
    mtype tmp[MAT_KDTREE_MAX_DIMS];
    int idx;
    memcpy(tmp, x->x, sizeof(tmp));
    memcpy(x->x, y->x, sizeof(tmp));
    memcpy(y->x, tmp, sizeof(tmp));
    memcpy(&idx, &(x->idx), sizeof(idx));
    memcpy(&(x->idx), &(y->idx), sizeof(idx));
    memcpy(&(y->idx), &idx, sizeof(idx));
}

MAT_KDNODE __kd_find_median(MAT_KDNODE kd_start, MAT_KDNODE kd_end, int idx)
{
    MAT_KDNODE p, store, md = kd_start+(kd_end-kd_start)/2;
    mtype pivot;
    if(kd_end<=kd_start) return NULL;
    if(kd_end==(kd_start+1)) return kd_start;
    while(1)
    {
        pivot = md->x[idx];
        __kd_swap(md, kd_end-1);
        for(store=p=kd_start; p<kd_end; p++)
        {
            if(p->x[idx]<pivot)
            {
                if(p!=store) __kd_swap(p, store);
                store++;
            }
        }
        __kd_swap(store, kd_end-1);
        if(store->x[idx] == md->x[idx]) return md;
        if(store>md) kd_end = store;
        else kd_start = store;
    }
}

MAT_KDTREE mat_kdtree_make_tree(MATRIX a, MAT_KDTREE result)
{
    int m, n, i, j;
    if(result==NULL) if((result = (MAT_KDTREE)malloc(sizeof(mat_kdtree)))==NULL) mat_error(MAT_MALLOC);
    m = MatRow(a);
    n = MatCol(a);
    if(m>MAT_KDTREE_MAX_DIMS) mat_error(MAT_SIZEMISMATCH);
    if(result->_is_allocated)
    {
        if(result->length!=n)
        {
            free(result->__data);
        }
    }
    else
    {
        result->_is_allocated = 0;
        result->__data = calloc(n, sizeof(mat_kdnode));
    }
    result->_is_allocated = 1;

    result->ndims = m;
    result->length = n;
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            result->__data[i].x[j] = a[j][i];
        }
        result->__data[i].idx = i;
    }
    result->kdtree = __mat_kdtree_make_tree(result->__data, n, 0, m);
    return result;
}

int mat_kdtree_free(MAT_KDTREE t)
{
    if(t->_is_allocated) free(t->__data);
    if(t!=NULL) free(t);
    return 0;
}

MAT_KDNODE __mat_kdtree_make_tree(MAT_KDNODE t, int len, int i, int dim)
{
    MAT_KDNODE n;
    if(!len) return 0;
    if((n = __kd_find_median(t, t+len, i)))
    {
        i = (i+1)%dim;
        n->left  = __mat_kdtree_make_tree(t, n-t, i, dim);
        n->right = __mat_kdtree_make_tree(n+1, t+len-(n+1), i, dim);
    }
    return n;
}

MATRIX mat_kdtree_nearest(MAT_KDTREE t, MATRIX a, MATRIX result)
{
    int i, m;
    mat_kdnode nd;
    MAT_KDNODE best = NULL;
    mtype best_dist;
    m = MatRow(a);
    if(t->ndims!=m) mat_error(MAT_SIZEMISMATCH);
    if(result==NULL) if((result = mat_creat(2,1, ZERO_MATRIX))==NULL) mat_error(MAT_MALLOC);
    nd.idx = -1;
    for(i=0; i<m; ++i) nd.x[i] = a[i][0];
    __mat_kdtree_nearest(t->kdtree, &nd, 0, t->ndims, &best, &best_dist);
    if(best!=NULL)
    {
        result[0][0] = best->idx;
        result[1][0] = best_dist;
    }
    return result;
}

void __mat_kdtree_nearest(MAT_KDNODE root, MAT_KDNODE nd, int i, int dim, MAT_KDNODE *best, mtype *best_dist)
{
    mtype d, dx, dx2;
    if(!root) return;
    d = __kd_dist(root, nd, dim);
    dx = root->x[i]-nd->x[i];
    dx2 = dx*dx;
    if(!*best||d<(*best_dist))
    {
        *best_dist = d;
        *best = root;
    }
    if(!*best_dist) return;
    if(++i >= dim) i = 0;
    __mat_kdtree_nearest(dx>0?root->left:root->right, nd, i, dim, best, best_dist);
    if(dx2>=(*best_dist)) return;
    __mat_kdtree_nearest(dx>0?root->right:root->left, nd, i, dim, best, best_dist);
}

