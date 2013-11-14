#ifndef __MATRIX__
#define __MATRIX__
#define _CRT_SECURE_NO_DEPRECATE

#ifdef __cplusplus
#define __MATRIX__CPP__
extern "C"
{
#endif

#ifndef mtype_n
#define mtype_n 1
#endif

#if mtype_n == 0
#define MTYPE FLOAT
#define mtype float
#else
#define MTYPE DOUBLE
#define mtype double
#endif
#include <stdio.h>
#include <time.h>
#include <omp.h>

#define STACK_MAX 100

#define __6188123NAME2(fun,suffix) fun ## _ ## suffix
#define __1267861NAME1(fun,suffix) __6188123NAME2(fun,suffix)

#define MTYPE_(x) __1267861NAME1(MTYPE, x)
#define mtype_(x) __1267861NAME1(mtype, x)

/******************************************/
#if defined (GCC) || defined (__GNUC__)
typedef FILE *FILEPOINTER;
#else
typedef struct _iobuf *FILEPOINTER;
#endif
struct int_stack
{
    int p;
    int length;
    int *stack;
};
typedef struct int_stack int_stack;
typedef struct int_stack *INT_STACK;

/******************************************/
struct mat_mtype_stack
{
    int p;
    int length;
    mtype *stack;
};
typedef struct mat_mtype_stack mat_mtype_stack;
typedef struct mat_mtype_stack *MAT_MTYPE_STACK;

/******************************************/
struct __intnode
{
    int data;
    struct __intnode *next;
};
typedef struct __intnode qintdata;
typedef struct __intnode *QINTDATA;
struct int_queue
{
    int p;
    QINTDATA head;
    QINTDATA tail;
};
typedef struct int_queue int_queue;
typedef struct int_queue *INT_QUEUE;

/******************************************/
#define __mtype(x) __mtype ## x
struct __mtype_node
{
    mtype data;
    struct __mtype_node *next;
};
typedef struct __mtype_node qmdata;
typedef struct __mtype_node *QMDATA;
struct mat_mtype_queue
{
    int p;
    QMDATA head;
    QMDATA tail;
};
typedef struct mat_mtype_queue mat_mtype_queue;
typedef struct mat_mtype_queue *MAT_MTYPE_QUEUE;

/******************************************/
struct __intpqnode
{
    int data;
    int priority;
};
typedef struct __intpqnode pqintdata;
typedef struct __intpqnode *PQINTDATA;
struct int_priorityqueue
{
    int p;
    int length;
    PQINTDATA element;
};
typedef struct int_priorityqueue int_priorityqueue;
typedef struct int_priorityqueue *INT_PRIORITYQUEUE;

/******************************************/
struct __mtype_pqnode
{
    mtype data;
    int priority;
};
typedef struct __mtype_pqnode pqmdata;
typedef struct __mtype_pqnode *PQMDATA;
struct mat_mtype_priorityqueue
{
    int p;
    int length;
    PQMDATA element;
};
typedef struct mat_mtype_priorityqueue mat_mtype_priorityqueue;
typedef struct mat_mtype_priorityqueue *MAT_MTYPE_PRIORITYQUEUE;

/******************************************/
struct tree_node
{
    mtype element;
    struct tree_node *left, *right;
};
typedef struct tree_node *SEARCH_TREE;

/******************************************/
typedef int *INT_VECTOR;

/******************************************/
typedef struct __mathead
{
    int row;
    int col;
}	MATHEAD;
typedef struct __matbody
{
    MATHEAD head;
    mtype *matrix;
}	MATBODY;
typedef mtype **MATRIX;

/******************************************/
typedef INT_VECTOR *INT_VECSTACK;
typedef MATRIX *MATSTACK;

typedef void **MATVEC_DPOINTER;

/******************************************/
struct bayes_model
{
    int num_of_classes;
    int num_of_features;
    INT_VECTOR class_labels;
    MATRIX class_priors;
    MATSTACK class_means;
    MATSTACK class_covars;
};
typedef struct bayes_model bayes_model;
typedef struct bayes_model *BAYES_MODEL;

/******************************************/
struct perceptron
{
    int num_of_classes;
    int num_of_features;
    INT_VECTOR class_labels;
    MATRIX class_weights;
    int istrained;
    int num_of_iterations;
};
typedef struct perceptron perceptron;
typedef struct perceptron *PERCEPTRON;

struct g_node
{
    int v;
    double weight;
    struct g_node *next;
};
typedef struct g_node g_node;
typedef struct g_node *G_NODE;

struct mat_graph
{
    int V, E;
    int *val, *vseq;
    int id;
    G_NODE *adj, z;
    int *dad;
    int weighted;
    INT_PRIORITYQUEUE pq;
};
typedef struct mat_graph mat_graph;
typedef struct mat_graph *MAT_GRAPH;

#ifndef MAT_KDTREE_MAX_DIMS
#define MAT_KDTREE_MAX_DIMS 3
#endif

struct __kdnode
{
    mtype x[MAT_KDTREE_MAX_DIMS];
    int idx;
    struct __kdnode *left, *right;
};
typedef struct __kdnode mat_kdnode;
typedef struct __kdnode* MAT_KDNODE;

struct __mat_kdtree
{
    int ndims;
    int length;
    int _is_allocated;
    MAT_KDNODE __data;
    MAT_KDNODE kdtree;
};
typedef struct __mat_kdtree mat_kdtree;
typedef struct __mat_kdtree* MAT_KDTREE;


extern clock_t MAT_CLOCK_TIME;
extern unsigned int MAT_SEED;
extern int MAT_SET_SEED;

/******************************************/
#define Int_VecLen(a)	(*(a-1))

#define MatRow(a)	(*((int*)(a-2)))
#define MatCol(a)	(*((int*)(a-1)))
#define MatNumel(a) ((*((int*)(a-2)))*(*((int*)(a-1))))

#define MatStacklength(a) ((int)((a-1)[0][0][0]))
/*  (*(int*)(a-1)) */
#define Int_VecStackLength(a) ((a-1)[0][0])
/* (*(int*)(a-1)) */

#define GEN_NOT_CONVERGED 1
#define GEN_FNOTOPEN 2
#define GEN_FNOTGETMAT 3
#define GEN_SIZEMISMATCH 4
#define GEN_MATH_ERROR 5
#define GEN_MALLOC 6
#define GEN_NOT_FOUND 7
#define GEN_SIZE_ERROR 8
#define GEN_BAD_TYPE 9

#define MAT_MALLOC 1
#define MAT_FNOTOPEN 2
#define MAT_FNOTGETMAT 3
#define MAT_SIZEMISMATCH 4
#define MAT_INVERSE_ILL_COND 5
#define MAT_INVERSE_NOT_SQUARE 6
#define MAT_CHOLESKY_FAILED 7

#define MATSTACK_MALLOC 1
#define MATSTACK_FNOTOPEN 2
#define MATSTACK_FNOTGETMAT 3
#define MATSTACK_SIZEMISMATCH 4
#define MATSTACK_INVERSE_ERROR 5

#define INT_VEC_MALLOC 1
#define INT_VEC_FNOTOPEN 2
#define INT_VEC_FNOTGETINT_VEC 3
#define INT_VEC_SIZEMISMATCH 4

#define INT_VECSTACK_MALLOC 1
#define INT_VECSTACK_FNOTOPEN 2
#define INT_VECSTACK_FNOTGETMAT 3
#define INT_VECSTACK_SIZEMISMATCH 4
#define INT_VECSTACK_INVERSE_ERROR 5

#define STACK_MALLOC 0
#define STACK_EMPTY 1

#define QUEUE_MALLOC 0
#define QUEUE_EMPTY 1

#define PQ_MALLOC 0
#define PQ_EMPTY 1

#define GRAPH_MALLOC 0
#define GRAPH_READ 1
#define GRAPH_ELSE 2

#define UNDEFINED -1
#define ZERO_MATRIX 0
#define UNIT_MATRIX 1
#define ONES_MATRIX 2
#define ZERO_INT_VECTOR 0
#define ONES_INT_VECTOR 2
#define SERIES_INT_VECTOR 3

#define GEN_LESS_THAN 0
#define GEN_GREATER_THAN 1
#define GEN_EQUAL_TO 2
#define GEN_NOT_EQUAL_TO 3
#define GEN_LESS_THAN_EQUAL_TO 4
#define GEN_GREATER_THAN_EQUAL_TO 5

#define COLS 1
#define ROWS 2

#define PQ_NORMAL 0
#define PQ_INCREASE 1
#define PQ_DECREASE 2

#define eps 1E-20
#define Eps 1E-10
#define EPS 1E-5
#define DOUBLE_MAX 1E37

#define MAT_LOSS_HUBER 0
#define MAT_LOSS_BISQUARE 1

#define MAT_MDS_METRIC 0
#define MAT_MDS_NONMETRIC 1

#define MAT_FFT2_FORWARD 1
#define MAT_FFT2_BACKWARD -1


#define MAT_PCA_CORRELATION 0
#define MAT_PCA_COVARIANCE 1
#define MAT_PCA_SUMOFSQUARES 2


/******************************************/
int mats_isnan(mtype x);
int mats_isinf(mtype x);

INT_VECTOR _int_vec_creat(int length);
INT_VECTOR int_vec_creat(int length, int type);
INT_VECTOR int_vec_fill(INT_VECTOR A, int type);
int int_vec_free(INT_VECTOR A);

/******************************************/
INT_VECSTACK _int_vecstack_creat(int length);
INT_VECSTACK int_vecstack_creat(int length);
int int_vecstack_free(INT_VECSTACK A);

MATRIX _mat_creat(int r, int c);
MATRIX mat_creat(int r, int c, int type);
MATRIX mat_creat_diag(MATRIX diag_vals);
MATRIX mat_fill(MATRIX A, int v);
int mat_free(MATRIX A);

/******************************************/
MATSTACK matstack_creat(int length);
MATSTACK _matstack_creat(int length);
int matstack_free(MATSTACK A);
MATSTACK matstack_append(MATSTACK s, MATRIX a);

MATVEC_DPOINTER matvec_creat(void);
int matvec_free(MATVEC_DPOINTER a);

/******************************************/
MATRIX mat_copy(MATRIX A, MATRIX result);
MATRIX mat_xcopy(MATRIX A, int si, int ei, int sj, int ej, MATRIX result);
MATRIX mat_xjoin(MATRIX A11, MATRIX A12, MATRIX A21, MATRIX A22, MATRIX result);
MATRIX mat_colcopy1(MATRIX A, MATRIX B, int colA, int colB);
int fgetmat(MATRIX A, FILEPOINTER fp);

/******************************************/
/* matrix dump */
MATRIX mat_dump(MATRIX A);
MATRIX mat_dumpf(MATRIX A, char *s);
MATRIX mat_fdump(MATRIX A, FILEPOINTER fp);
MATRIX mat_fdumpf(MATRIX A, char *s, FILEPOINTER fp);

INT_VECTOR int_vec_dump(INT_VECTOR a);
INT_VECTOR int_vec_dumpf(INT_VECTOR a, char *s);
INT_VECTOR int_vec_fdump(INT_VECTOR a, FILEPOINTER fp);
INT_VECTOR int_vec_fdumpf(INT_VECTOR a, char *s, FILEPOINTER fp);

/******************************************/
/* vector manipulations */
INT_VECTOR int_vec_copy(INT_VECTOR a);
INT_VECTOR int_vec_unique(INT_VECTOR a);
INT_VECTOR int_vec_append(INT_VECTOR A, int i);
INT_VECTOR int_vec_find(INT_VECTOR a, int rel_type, int n);
INT_VECTOR int_vec_concat(INT_VECTOR A, INT_VECTOR B, INT_VECTOR result);
INT_VECTOR mat_get_sub_vector(INT_VECTOR data, INT_VECTOR indices);

/******************************************/
/* error handling functions */
int gen_error(int err_);
INT_VECTOR int_vec_error(int err_);
INT_VECSTACK int_vecstack_error (int err_);
MATRIX mat_error(int err_);
MATSTACK matstack_error(int err_);
int stack_error( int err_);
int queue_error( int err_);
int pq_error( int err_);
int graph_error(int err_);

/******************************************/
/* basic matrix operations */
mtype mat_mean(MATRIX A);
MATRIX mat_mean_row(MATRIX A);
MATRIX mat_mean_col(MATRIX A);

mtype mat_sum(MATRIX A);
MATRIX mat_sum_row(MATRIX A);
MATRIX mat_sum_col(MATRIX A);

MATRIX mat_abs(MATRIX A);
MATRIX mat_add(MATRIX A, MATRIX B, MATRIX result);
MATRIX mat_adds(MATRIX A, mtype s, MATRIX result);

INT_VECTOR int_vec_add(INT_VECTOR A, INT_VECTOR B, INT_VECTOR result);
INT_VECTOR int_vec_adds(INT_VECTOR A, int s, INT_VECTOR result);

INT_VECTOR int_vec_sub(INT_VECTOR A, INT_VECTOR B, INT_VECTOR result);
INT_VECTOR int_vec_subs(INT_VECTOR A, int s, INT_VECTOR result);

INT_VECTOR int_vec_mul(INT_VECTOR A, INT_VECTOR B, INT_VECTOR result);
INT_VECTOR int_vec_muls(INT_VECTOR A, int s, INT_VECTOR result);

INT_VECTOR int_vec_div(INT_VECTOR A, INT_VECTOR B, INT_VECTOR result);
INT_VECTOR int_vec_divs(INT_VECTOR A, int s, INT_VECTOR result);


MATRIX mat_sub(MATRIX A, MATRIX B, MATRIX result);
MATRIX mat_subs(MATRIX A, mtype s, MATRIX result);

MATRIX mat_mul(MATRIX A, MATRIX B, MATRIX result);
MATRIX mat_mul_fast(MATRIX A, MATRIX B, MATRIX result);
MATRIX mat_mul_dot(MATRIX A, MATRIX B, MATRIX result);
MATRIX mat_muls(MATRIX A, mtype s, MATRIX result);

MATRIX mat_div_dot(MATRIX A, MATRIX B, MATRIX result);
MATRIX mat_divs(MATRIX A, mtype s, MATRIX result);

/******************************************/
/* norm calculations */
mtype mat_innerprod(MATRIX A, MATRIX B);
mtype mat_norm_inf(MATRIX A);
mtype mat_norm_one(MATRIX A);
mtype mat_norm_p(MATRIX A, mtype p);

/******************************************/
mtype mat_diagmul(MATRIX A);
MATRIX mat_tran(MATRIX A, MATRIX result);
MATRIX mat_inv(MATRIX a, MATRIX result);
MATRIX mat_pinv(MATRIX a, MATRIX result);
MATRIX mat_wpinv(MATRIX a, MATRIX w, MATRIX result);
MATRIX mat_reg_inv(MATRIX a, mtype r_constant, MATRIX result);
MATRIX mat_symtoeplz(MATRIX R);

/******************************************/
/* linear system equation solver */
int mat_lu(MATRIX A, MATRIX P);
void mat_backsubs1(MATRIX A, MATRIX B, MATRIX C, MATRIX P, int xcol);
MATRIX mat_lsolve(MATRIX a, MATRIX b);
MATRIX mat_cholesky(MATRIX a, MATRIX result);

MATRIX mat_submat(MATRIX A, int i, int j);
mtype mat_cofact(MATRIX A, int i, int j);
mtype mat_det(MATRIX a);
mtype mat_minor(MATRIX A, int i, int j);

MATSTACK mat_qr(MATRIX A, MATSTACK qr);
MATRIX mat_durbin(MATRIX R);
MATRIX mat_lsolve_durbin(MATRIX A, MATRIX B);

/******************************************/
/* sorting algorithms */
mtype mat_median(MATRIX A);
mtype mat_order_statistic(MATRIX A, int k);
void __quicksort(MATRIX a, int l, int r, int offset, MATRIX ind);
MATSTACK mat_qsort(MATRIX A, int dim, MATSTACK result);
MATVEC_DPOINTER mat_max(MATRIX A, int dim);
MATVEC_DPOINTER mat_min(MATRIX A, int dim);

/******************************************/
/* random number generation */
MATRIX mat_rand(int r, int c, MATRIX result);
MATRIX mat_randn(int r, int c, MATRIX result);
MATRIX mat_randexp(int r, int c, mtype mu, MATRIX result);
INT_VECTOR int_vec_permute_vect(int n, int k, INT_VECTOR result);
MATRIX mat_randfun(int r, int c, mtype (*fun)(mtype), mtype xmin, mtype xmax, MATRIX result);
void mat_set_seed(int seed);
mtype mat_randfun_(mtype (*fun)(mtype), mtype xmin, mtype xmax);
mtype mat_rand_(void);
mtype mat_randn_(void);
mtype mat_randexp_(mtype mu);
MATRIX mat_randperm(int m, int n, MATRIX result);
MATRIX mat_randperm_(int n, MATRIX result);
INT_VECTOR int_vec_randperm(int n, INT_VECTOR result);

/******************************************/
/* linear regressions */
MATRIX mat_least_squares(MATRIX A, MATRIX Y, MATRIX result);
MATRIX mat_w_least_squares(MATRIX A, MATRIX Y, MATRIX w, MATRIX result);
MATRIX mat_rob_least_squares(MATRIX A, MATRIX Y, int lossfunc);

MATRIX mat_linear_ls_fit(MATRIX data, MATRIX Y, int degree);
MATRIX mat_robust_fit(MATRIX data, MATRIX Y, int degree, int lossfunc);

/******************************************/
/* matrix manipulations */
MATRIX mat_concat(MATRIX A, MATRIX B , int dim);
MATRIX mat_get_sub_matrix_from_rows(MATRIX data, INT_VECTOR indices, MATRIX submat);
MATRIX mat_get_sub_matrix_from_cols(MATRIX data, INT_VECTOR indices, MATRIX submat);
MATRIX mat_pick_row(MATRIX data, int r);
MATRIX mat_pick_col(MATRIX data, int c);
INT_VECSTACK mat_find(MATRIX a, int rel_type, mtype x);

MATRIX mat_fliplr(MATRIX A, MATRIX result);
MATRIX mat_flipud(MATRIX A, MATRIX result);

/******************************************/
/* distance tools */
MATRIX mat_calc_dist3_sq(MATRIX data, MATRIX curr_data);
INT_VECTOR mat_find_within_dist(MATRIX data, MATRIX curr_data, mtype range);

void _cart2pol(mtype x, mtype y, mtype *rho, mtype *th);
MATRIX mat_cart2pol(MATRIX A, int dim);

mtype _dist3_sq(mtype ux, mtype uy, mtype uz, mtype vx, mtype vy, mtype vz);

/******************************************/
/* function tools */
__inline mtype huber_wt(mtype x, mtype k);
__inline mtype bisquare_wt(mtype x, mtype k);
__inline mtype logplusone(mtype x);
mtype arcsinh(mtype x);
mtype arccosh(mtype x);
mtype arctanh(mtype x);

MATRIX mat_bisquare_wt(MATRIX A, mtype k, mtype sigma);
MATRIX mat_huber_wt(MATRIX A, mtype k, mtype sigma);
MATRIX mat_gfunc(MATRIX A, mtype (*pt2func)(mtype), MATRIX result);
MATRIX mat_bsxfun(MATRIX a, MATRIX b, MATRIX result, mtype (*pt2func)(mtype, mtype));

/******************************************/
/* statistical analysis tools */
MATSTACK mat_corcol(MATRIX data);
MATSTACK mat_covcol(MATRIX data);
MATRIX mat_scpcol(MATRIX data);
void tred2(MATRIX a, MATRIX d, MATRIX e);
void tqli(MATRIX d, MATRIX e, MATRIX z);
MATSTACK mat_pca(MATRIX data, int pca_type);
MATSTACK mat_eig_sym(MATRIX symmat, MATSTACK result);

/******************************************/
/* print helpers */
void mat_nextline(void);
void mat_fnextline(FILEPOINTER fp);

/******************************************/
/* signal domain transforms */
int __powerof2(int width, int *m, int *twopm);

MATSTACK mat_fft2(MATSTACK c, int dir);
int __fft(int dir, int m, mtype *x, mtype *y);

/* filtering functions */
MATRIX mat_conv2(MATRIX a, MATRIX mask, MATRIX scratch, MATRIX result);

/******************************************/
/* matrix vector conversions */
INT_VECTOR mat_2int_vec(MATRIX a);
MATRIX int_vec2_mat(INT_VECTOR a, int direction);
MATRIX mat_vectorize(MATRIX a, MATRIX result);
MATRIX mat_vectorize_tr(MATRIX a, MATRIX result);

/******************************************/
/* integration functions */
mtype mat_int_simpson(mtype (*func)(mtype), int n, mtype lower, mtype upper);
mtype _lint(mtype *x, mtype (*fx) (mtype), mtype x0, mtype xn, mtype f0, mtype f2, mtype f3, mtype f5, mtype f6, mtype f7, mtype f9, mtype fl4, mtype hmin, mtype hmax, mtype re, mtype ae);
mtype mat_int_qadrat(mtype(*fx)(mtype), mtype lower, mtype upper);

/******************************************/
/* polynomials */
MATRIX mat_evalpoly(MATRIX a, mtype x, int direction);
MATRIX mat_dpoly(MATRIX a, int direction);
MATRIX mat_devalpoly(MATRIX a, mtype x, int direction);

/******************************************/
/* pattern recognition */
BAYES_MODEL pat_bayes_model_creat(void);
int pat_bayes_model_free(BAYES_MODEL a);
PERCEPTRON pat_perceptron_creat(void);
int pat_perceptron_free(PERCEPTRON a);

BAYES_MODEL pat_bayes_classifier_train(MATRIX data, INT_VECTOR labels);
INT_VECTOR pat_bayes_classifier_test(MATRIX data, BAYES_MODEL b_model);
PERCEPTRON pat_perceptron_train(MATRIX data, INT_VECTOR labels, int num_of_iterations);
INT_VECTOR pat_perceptron_test(MATRIX data, PERCEPTRON p_model);
PERCEPTRON pat_perceptron_train_(MATRIX data1, MATRIX data2, PERCEPTRON p_model, int class_num);

/******************************************/
/* basic search algorithms */
SEARCH_TREE mat_bs_make_null(void);
SEARCH_TREE mat_bs_free(SEARCH_TREE T);
SEARCH_TREE mat_bs_find(mtype x, SEARCH_TREE T);
SEARCH_TREE mat_bs_find_min(SEARCH_TREE T);
SEARCH_TREE mat_bs_find_max(SEARCH_TREE T);
SEARCH_TREE mat_bs_insert(mtype x, SEARCH_TREE T);
SEARCH_TREE mat_bs_delete(mtype x, SEARCH_TREE T);
int mat_bs_inorder(SEARCH_TREE T, int index, mtype **ordered);


/******************************************/
/* relational functions not reqd */
int gen_gt(mtype a);
int gen_lt(mtype a);
int gen_eq(mtype a);
mtype gen_abs_ceil(mtype a);

/******************************************/
/* textfunctions */
int mat_isnumeric(FILEPOINTER fp);
int mat_go_next_word(FILEPOINTER fp);
int mat_count_words_in_line(FILEPOINTER fp, int *count);
int mat_read_word(FILEPOINTER fp, char *c_word);
MATRIX mat_dlmread(const char *fname);
void mat_dlmwrite(const char *fname, MATRIX a);

/******************************************/
/* mat_timer_functions */
void mat_tic(void);
double mat_toc(void);
void mat_toc_print(void);

/******************************************/

/* data structures */
INT_STACK int_stack_creat(void);
int int_stack_free(INT_STACK s);
void int_stack_push(INT_STACK s, int value);
int int_stack_pop(INT_STACK s);
int int_stack_is_empty(INT_STACK s);

MAT_MTYPE_STACK mat_mtype_stack_creat(void);
int mat_mtype_stack_free(MAT_MTYPE_STACK s);
void mat_mtype_stack_push(MAT_MTYPE_STACK s, mtype value);
mtype mat_mtype_stack_pop(MAT_MTYPE_STACK s);
int mat_mtype_stack_is_empty(MAT_MTYPE_STACK s);

INT_QUEUE int_queue_creat(void);
int int_queue_free(INT_QUEUE s);
void int_queue_enqueue(INT_QUEUE s, int value);
int int_queue_dequeue(INT_QUEUE s);
int int_queue_is_empty(INT_QUEUE s);

MAT_MTYPE_QUEUE mat_mtype_queue_creat(void);
int mat_mtype_queue_free(MAT_MTYPE_QUEUE s);
void mat_mtype_queue_enqueue(MAT_MTYPE_QUEUE s, mtype value);
mtype mat_mtype_queue_dequeue(MAT_MTYPE_QUEUE s);
int mat_mtype_queue_is_empty(MAT_MTYPE_QUEUE s);

INT_PRIORITYQUEUE int_priorityqueue_creat(void);
void int_priorityqueue_enqueue(INT_PRIORITYQUEUE H, int data, int priority);
int int_priorityqueue_dequeue(INT_PRIORITYQUEUE H);
int int_priorityqueue_free(INT_PRIORITYQUEUE H);
int int_priorityqueue_update(INT_PRIORITYQUEUE H, int data, int priority, int type);
int int_priorityqueue_is_empty(INT_PRIORITYQUEUE H);

/******************************************/
/* mds */
MATRIX mat_mds(MATRIX d, int dims, int type, MATRIX result);
MATRIX _mat_mds_metric(MATRIX d, int dims, MATRIX result);
MATRIX _mat_mds_nonmetric(MATRIX d, int dims, MATRIX result);
mtype _mat_sqr(mtype x);
mtype _mat_sqrt(mtype x);

/******************************************/
/* graph */
MAT_GRAPH mat_graph_creat(void);
void mat_graph_adjlist(MAT_GRAPH g, int directed, int weighted, FILEPOINTER fp);
INT_QUEUE mat_graph_search(MAT_GRAPH g, int connected, int mst);
void mat_graph_visit(MAT_GRAPH g, int k, int connected, int mst, INT_PRIORITYQUEUE pq, INT_QUEUE q);
void mat_graph_dumpf(MAT_GRAPH g, int mst, FILEPOINTER fp);
void mat_graph_dump(MAT_GRAPH g, int mst);
void mat_graph_adjm_to_adjl(MAT_GRAPH g, MATRIX a);
MAT_GRAPH mat_graph_reverse(MAT_GRAPH g, MAT_GRAPH r);

/******************************************/
/* kdtree */
MAT_KDTREE mat_kdtree_make_tree(MATRIX a, MAT_KDTREE result);
int mat_kdtree_free(MAT_KDTREE t);
MATRIX mat_kdtree_nearest(MAT_KDTREE t, MATRIX a, MATRIX result);
MAT_KDNODE __mat_kdtree_make_tree(MAT_KDNODE t, int len, int i, int dim);
MAT_KDNODE __kd_find_median(MAT_KDNODE kd_start, MAT_KDNODE kd_end, int idx);
void __mat_kdtree_nearest(MAT_KDNODE root, MAT_KDNODE nd, int i, int dim, MAT_KDNODE *best, mtype *best_dist);


#ifdef __cplusplus
}
#endif


#endif

