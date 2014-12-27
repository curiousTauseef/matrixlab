#include <string.h>
#include <malloc.h>
#include "matrix.h"

/** \cond HIDDEN_SYMBOLS */

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

/** \endcond */

/** \brief Creates a k-d tree from a data matrix
 *
 * \param[in] A Input data matrix of size \f$ d \times N \f$
 * \param[in] result K-d tree to store the result
 * \return Output k-d tree
 *
 */

MAT_KDTREE mat_kdtree_make_tree(MATRIX A, MAT_KDTREE result)
{
    int m, n, i, j;
    m = MatRow(A);
    n = MatCol(A);
    if(m>MAT_KDTREE_MAX_DIMS) mat_error(MAT_SIZEMISMATCH);
    if(result==NULL)
    {
        if((result = (MAT_KDTREE)malloc(sizeof(mat_kdtree)))==NULL) mat_error(MAT_MALLOC);
    }
    else if(result->_is_allocated) free(result->data);
    result->data = calloc(n, sizeof(mat_kdnode));
    result->_is_allocated = 1;

    result->ndims = m;
    result->length = n;
    for(i=0; i<n; ++i)
    {
        for(j=0; j<m; ++j)
        {
            result->data[i].x[j] = A[j][i];
        }
        result->data[i].idx = i;
    }
    result->kdtree = __mat_kdtree_make_tree(result->data, n, 0, m);
    return result;
}

/** \brief Frees a k-d tree
 *
 * \param[in] t Input k-d tree
 * \return Success
 *
 */

int mat_kdtree_free(MAT_KDTREE t)
{
    if(t->_is_allocated) free(t->data);
    if(t!=NULL) free(t);
    return 1;
}

/** \cond HIDDEN_SYMBOLS */

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

/** \endcond */

/** \brief Computes nearest neighbors
 *
 * \param[in] t Input k-d tree
 * \param[in] A Input data matrix of size \f$ d \times N \f$
 * \param[in] result Matrix to store the result
 * \return Output matrix \f$ B \f$ with index B[0][j] and squared distance B[1][j] for \f$ j=1,2,\cdots, N\f$
 *
 */

MATRIX mat_kdtree_nearest(MAT_KDTREE t, MATRIX A, MATRIX result)
{
    int i, j, m, n;
    mat_kdnode nd;
    MAT_KDNODE best;
    mtype best_dist;
    m = MatRow(A);
    n = MatCol(A);
    if(t->ndims!=m) mat_error(MAT_SIZEMISMATCH);
    if(result==NULL) if((result = mat_creat(2, n, ZERO_MATRIX))==NULL) mat_error(MAT_MALLOC);
    for(j=0; j<n; ++j)
    {
        best = NULL;
        nd.idx = -1;
        for(i=0; i<m; ++i) nd.x[i] = A[i][j];
        __mat_kdtree_nearest(t->kdtree, &nd, 0, t->ndims, &best, &best_dist);
        if(best!=NULL)
        {
            result[0][j] = (mtype)best->idx;
            result[1][j] = best_dist;
        }
        else
        {
            result[0][j] = -1.0;
            result[1][j] = -1.0;
        }
    }
    return result;
}

/** \cond HIDDEN_SYMBOLS */

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

/** \endcond */

/** \brief Computes k nearest neighbors
 *
 * \param[in] t Input k-d tree
 * \param[in] A Input data matrix of size \f$ d \times N \f$
 * \param[in] k Number of neighbors
 * \param[in] result Matrix to store the result
 * \return Output matrix \f$ B \f$ with index B[0][j] and squared distance B[1][j] for \f$ j=1,2,\cdots, N\f$
 *
 */

MATRIX mat_kdtree_k_nearest(MAT_KDTREE t, MATRIX A, int k, MATRIX result)
{
    int i, j, kk, m, n, o, p;
    mat_kdnode nd;
    MATRIX bmin = NULL, bmax = NULL;
    MAT_MTYPE_PRIORITYQUEUE pq = NULL;
    mat_mtypepqnode pqn;
    m = MatRow(A);
    n = MatCol(A);
    if(t->ndims!=m) mat_error(MAT_SIZEMISMATCH);
    if(result==NULL) if((result = mat_creat(2, k*n, ZERO_MATRIX))==NULL) mat_error(MAT_MALLOC);
    #pragma omp parallel for private(i, bmin, bmax, pq, nd, kk, o, p, pqn), shared(t, A, result)
    for(j=0; j<n; ++j)
    {
        bmin = mat_creat(1,t->ndims, UNDEFINED);
        bmax = mat_creat(1, t->ndims, UNDEFINED);

        bmin = mat_fill(bmin, -MTYPE_MAX);
        bmax = mat_fill(bmax, +MTYPE_MAX);

        pq = mat_mtype_priorityqueue_creat(MAT_PQ_MIN);
        nd.idx = k;
        for(i=0; i<m; ++i) nd.x[i] = A[i][j];
        __mat_kdtree_k_nearest(t->kdtree, &nd, 0, t->ndims, pq, bmax, bmin);
        o = k*j;
        p = (pq->p)<k?(pq->p):k;
        for(kk=0; kk<p; ++kk)
        {
            pqn = mat_mtype_priorityqueue_dequeue(pq);
            result[0][o+kk] = pqn.data;
            result[1][o+kk] = pqn.priority;
        }
        for(kk=p; kk<k; ++kk)
        {
            result[0][o+kk] = -1.0;
            result[1][o+kk] = -1.0;
        }
        mat_mtype_priorityqueue_free(pq);
        mat_free(bmin);
        mat_free(bmax);

    }

    return result;
}

/** \cond HIDDEN_SYMBOLS */

void __mat_kdtree_k_nearest(MAT_KDNODE root, MAT_KDNODE nd, int i, int dim, MAT_MTYPE_PRIORITYQUEUE pq, MATRIX bmax, MATRIX bmin)
{
    mtype d, dx, tmp, tbound = 0.0;
    int l, df = 1;
    if(!root) return;
    d = __kd_dist(root, nd, dim);
    dx = root->x[i]-nd->x[i];
    mat_mtype_priorityqueue_enqueue(pq, root->idx, d);
    if(dx>0)
    {
        tmp = bmax[0][i];
        bmax[0][i] = root->x[i];
        __mat_kdtree_k_nearest(root->left, nd, (i+1)%dim, dim, pq, bmax, bmin);
        bmax[0][i] = tmp;
        tmp = bmin[0][i];
        bmin[0][i] = root->x[i];
        if(pq->p>=nd->idx)
        {
            for(l=0; l<dim && df==1; ++l)
            {
                if(nd->x[l]<bmin[0][l])
                {
                    tbound += (nd->x[l]-bmin[0][l])*(nd->x[l]-bmin[0][l]);
                    if(tbound>pq->element[nd->idx].data) df = 0;
                }
                else if(nd->x[l]>bmax[0][l])
                {
                    tbound += (nd->x[l]-bmax[0][l])*(nd->x[l]-bmax[0][l]);
                    if(tbound>pq->element[nd->idx].data) df = 0;
                }
            }
        }
        if(df==1)
        {
            __mat_kdtree_k_nearest(root->right, nd, (i+1)%dim, dim, pq, bmax, bmin);
        }
        bmin[0][i] = tmp;
    }
    else
    {
        tmp = bmin[0][i];
        bmin[0][i] = root->x[i];
        __mat_kdtree_k_nearest(root->right, nd, (i+1)%dim, dim, pq, bmax, bmin);
        bmin[0][i] = tmp;
        tmp = bmax[0][i];
        bmax[0][i] = root->x[i];
        if(pq->p>=nd->idx)
        {
            for(l=0; l<dim && df==1; ++l)
            {
                if(nd->x[l]<bmin[0][l])
                {
                    tbound += (nd->x[l]-bmin[0][l])*(nd->x[l]-bmin[0][l]);
                    if(tbound>pq->element[nd->idx].data) df = 0;
                }
                else if(nd->x[l]>bmax[0][l])
                {
                    tbound += (nd->x[l]-bmax[0][l])*(nd->x[l]-bmax[0][l]);
                    if(tbound>pq->element[nd->idx].data) df = 0;
                }
            }
        }
        if(df==1)
        {
            __mat_kdtree_k_nearest(root->left, nd, (i+1)%dim, dim, pq, bmax, bmin);
        }
        bmax[0][i] = tmp;
    }
}

/** \endcond */



