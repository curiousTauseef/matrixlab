/* ---------------------------------------------------------------------------
** This software is furnished "as is", without technical support,
** and with no warranty, express or implied, as to its usefulness for
** any purpose.
**
** matrix.h
** C Matrix Library
** Originally adapted from Small Matrix Toolbox for C programmers, ver. 0.4
** by Patrick Ko Shu-pui
**
** Author: Sk. Mohammadul Haque
** Copyright (c) 2013 Sk. Mohammadul Haque
** For more details and updates, visit http://mohammadulhaque.alotspace.com
** -------------------------------------------------------------------------*/
/*! \mainpage Matrixlab
 *
 * \section intro_sec Introduction
 * Matrixlab is a generic C library for matrix routines. It contains over 250 functions for matrix operations. Many of the functions are multi-threaded.
 *
 * The functions are categorized as
 * \section install_sec Installation
 *
 * \subsection step1 Step 1: Opening the box
 *
 * etc...
 */

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
#define MTYPE FLOAT /**< \def MTYPE is FLOAT. */
#define mtype float /**< \def mtype is float. */
#else
#define MTYPE DOUBLE /**< \def MTYPE is DOUBLE. */
#define mtype double /**< \def mtype is double. */
#endif
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#define STACK_MAX 100
#define MAT_PI (3.14159265359) /**< def MAT_PI pi */
#define __6188123NAME2(fun,suffix) fun ## _ ## suffix
#define __1267861NAME1(fun,suffix) __6188123NAME2(fun,suffix)

#define MTYPE_(x) __1267861NAME1(MTYPE, x)
#define mtype_(x) __1267861NAME1(mtype, x)

/******************************************/
#if defined (GCC) || defined (__GNUC__)
typedef FILE *MAT_FILEPOINTER; /**< \typedef FILEPOINTER is pointer to FILE */
#else
typedef struct _iobuf *MAT_FILEPOINTER; /**< \typedef FILEPOINTER is pointer to _iobuf */
#endif

/** \brief Integer Stack Structure
 *
 */

struct mat_int_stack
{
    int p; /**< Current stack position */
    int length; /**< Total allocated stack length */
    int *stack; /**< Pointer to stack data */
};
typedef struct mat_int_stack mat_int_stack; /**< Integer Stack */
typedef struct mat_int_stack *MAT_INT_STACK; /**< Integer Stack Pointer */

/******************************************/
/** \brief Mtype Stack Structure
 *
 */

struct mat_mtype_stack
{
    int p; /**< Current stack position */
    int length; /**< Total allocated stack length */
    mtype *stack; /**< Pointer to stack data */
};
typedef struct mat_mtype_stack mat_mtype_stack; /**< Mtype Stack */
typedef struct mat_mtype_stack *MAT_MTYPE_STACK; /**< Mtype Stack Pointer */

/******************************************/

/** \brief Integer Queue Node Structure
 *
 */

struct mat_qintnode
{
    int data; /**< Integer node data */
    struct mat_qintnode *next; /**< Pointer to next node */
};
typedef struct mat_qintnode mat_qintnode; /**< Integer Queue Node */
typedef struct mat_qintnode *MAT_QINTNODE; /**< Integer Queue Node Pointer */

/** \brief Integer Queue Structure
 *
 */

struct mat_int_queue
{
    int p; /**< Current queue position */
    MAT_QINTNODE head; /**< Queue head node */
    MAT_QINTNODE tail; /**< Queue tail node */
};
typedef struct mat_int_queue mat_int_queue; /**< Integer Queue */
typedef struct mat_int_queue *MAT_INT_QUEUE; /**< Integer Queue Pointer */

/******************************************/
#define __mtype(x) __mtype ## x
/** \brief Mtype Queue Node Structure
 *
 */

struct mat_qmtypenode
{
    mtype data; /**< Mtype node data */
    struct mat_qmtypenode *next; /**< Pointer to next node */
};
typedef struct mat_qmtypenode mat_qmtypenode; /**< Mtype Queue Node */
typedef struct mat_qmtypenode *MAT_QMTYPENODE; /**< Mtype Queue Node Pointer */
/** \brief Mtype Queue Structure
 *
 */

struct mat_mtype_queue
{
    int p; /**< Current queue position */
    MAT_QMTYPENODE head; /**< Queue head node */
    MAT_QMTYPENODE tail; /**< Queue tail node */
};
typedef struct mat_mtype_queue mat_mtype_queue; /**< Mtype Queue */
typedef struct mat_mtype_queue *MAT_MTYPE_QUEUE; /**< Mtype Queue Pointer */

/******************************************/
/** \brief Integer Priority Queue Node Structure
 *
 */

struct mat_intpqnode
{
    int data; /**< Integer node data */
    int priority; /**< Node priority */
};
typedef struct mat_intpqnode mat_intpqnode; /**< Integer Priority Queue Node */
typedef struct mat_intpqnode *MAT_INTPQNODE; /**< Integer Priority Queue Node Pointer */
/** \brief Integer Priority Queue Structure
 *
 */

struct mat_int_priorityqueue
{
    int p; /**< Current priority queue position */
    int length; /**< Total allocated priority queue length*/
    MAT_INTPQNODE element; /**< Pointer to priority queue data */
};
typedef struct mat_int_priorityqueue mat_int_priorityqueue; /**< Integer Priority Queue */
typedef struct mat_int_priorityqueue *MAT_INT_PRIORITYQUEUE; /**< Integer Priority Queue Pointer */

/******************************************/
/** \brief Mtype Priority Queue Node Structure
 *
 */

struct mat_mtypepqnode
{
    mtype data; /**< Mtype node data */
    int priority; /**< Node priority */
};
typedef struct mat_mtypepqnode mat_mtypepqnode; /**< Mtype Priority Queue Node */
typedef struct mat_mtypepqnode *MAT_MTYPEPQNODE; /**< Mtype Priority Queue Node Pointer */
/** \brief Mtype Priority Queue Structure
 *
 */

struct mat_mtype_priorityqueue
{
    int p; /**< Current priority queue position */
    int length; /**< Total allocated priority queue length*/
    MAT_MTYPEPQNODE element; /**< Pointer to priority queue data */
};
typedef struct mat_mtype_priorityqueue mat_mtype_priorityqueue; /**< Mtype Priority Queue */
typedef struct mat_mtype_priorityqueue *MAT_MTYPE_PRIORITYQUEUE; /**< Mtype Priority Queue Pointer */

/******************************************/
/** \brief Search Tree Node Structure
 *
 */

struct mat_tree_node
{
    mtype element; /**< Search tree node data */
    struct mat_tree_node *left; /**< Pointer to left child node */
    struct mat_tree_node *right;/**< Pointer to right child node */
};
typedef struct mat_tree_node mat_tree_node; /**< Search Tree Node */
typedef struct mat_tree_node *MAT_TREE_NODE; /**< Search Tree Node Pointer */
typedef struct mat_tree_node *MAT_TREE; /**< Search Tree Pointer */

/******************************************/
typedef int *INT_VECTOR; /**< Integer Vector */

/******************************************/
/** \cond HIDDEN_SYMBOLS */
typedef struct __mathead
{
    int row;
    int col;
} MATHEAD;

typedef struct __matbody
{
    MATHEAD head;
    mtype *matrix;
} MATBODY;
/** \endcond */
typedef mtype **MATRIX; /**< Mtype Matrix */

/******************************************/
typedef INT_VECTOR *INT_VECSTACK; /**< Integer Vector Stack */
typedef MATRIX *MATSTACK; /**< Mtype Matrix Stack */

typedef void **MATVEC_DPOINTER; /**< Mtype Matrix - Integer Vector Pair */

/******************************************/
/** \brief Bayes Classifier Model Structure
 *
 */

struct mat_bayes_model
{
    int num_of_classes; /**< Number of training class */
    int num_of_features; /**< Number of training features */
    INT_VECTOR class_labels; /**< Training data class label vector */
    MATRIX class_priors; /**< Training data prior information */
    MATSTACK class_means; /**< Training data class means */
    MATSTACK class_covars; /**< Training data class covariances */
};
typedef struct mat_bayes_model mat_bayes_model; /**< Bayes Classifier Model */
typedef struct mat_bayes_model *MAT_BAYES_MODEL; /**< Bayes Classifier Model Pointer */

/******************************************/
/** \brief Perceptron Classifier Model Structure
 *
 */

struct mat_perceptron
{
    int num_of_classes; /**< Number of training classes */
    int num_of_features; /**< Number of training features */
    INT_VECTOR class_labels; /**< Training data class label vector */
    MATRIX class_weights; /**< Trained Classifier Weights */
    int istrained; /**< Is trained */
    int num_of_iterations; /**< Number of training iterations */
};
typedef struct mat_perceptron mat_perceptron; /**< Perceptron Classifier Model */
typedef struct mat_perceptron *MAT_PERCEPTRON; /**< Perceptron Classifier Model Pointer */

/******************************************/
/** \brief Graph Node Structure
 *
 */

struct mat_gnode
{
    int v; /**< Value */
    double weight; /**< Node weight */
    struct mat_gnode *next; /**< Ponter to next node */
};
typedef struct mat_gnode mat_gnode; /**< Graph Node */
typedef struct mat_gnode *MAT_GNODE; /**< Graph Node Pointer */

/******************************************/
/** \brief Graph Structure
 *
 */

struct mat_graph
{
    int nvertices; /**< Number of vertices */
    int nedges; /**< Number of edges */
    int *val; /**<  */
    int *vseq;
    int id;
    MAT_GNODE *adj;
    MAT_GNODE z;
    int *dad;
    int weighted;
    MAT_INT_PRIORITYQUEUE pq;
};
typedef struct mat_graph mat_graph;
typedef struct mat_graph *MAT_GRAPH;

#ifndef MAT_KDTREE_MAX_DIMS
#define MAT_KDTREE_MAX_DIMS 3
#endif

struct mat_kdnode
{
    mtype x[MAT_KDTREE_MAX_DIMS];
    int idx;
    struct mat_kdnode *left, *right;
};
typedef struct mat_kdnode mat_kdnode;
typedef struct mat_kdnode* MAT_KDNODE;

struct mat_kdtree
{
    int ndims;
    int length;
    int _is_allocated;
    MAT_KDNODE data;
    MAT_KDNODE kdtree;
};
typedef struct mat_kdtree mat_kdtree;
typedef struct mat_kdtree* MAT_KDTREE;


extern clock_t MAT_CLOCK_TIME;
extern unsigned int MAT_SEED;
extern int MAT_SET_SEED;
extern MATSTACK mat_cheby_series_table, mat_legendre_series_table, mat_binom_series_table;

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
#define INT_VECSTACK_FNOTGETINT_VEC 3
#define INT_VECSTACK_SIZEMISMATCH 4

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

#define ROWS 0
#define COLS 1

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

INT_VECTOR __int_vec_creat(int len);
INT_VECTOR int_vec_creat(int len, int type);
INT_VECTOR int_vec_fill(INT_VECTOR A, int val);
INT_VECTOR int_vec_fill_type(INT_VECTOR A, int type);
int int_vec_free(INT_VECTOR A);

/******************************************/
INT_VECSTACK __int_vecstack_creat(int len);
INT_VECSTACK int_vecstack_creat(int len);
int int_vecstack_free(INT_VECSTACK A);

MATRIX __mat_creat(int r, int c);
MATRIX mat_creat(int r, int c, int type);
MATRIX mat_creat_diag(MATRIX diag_vals, MATRIX result);
MATRIX mat_fill(MATRIX A, mtype val);
MATRIX mat_fill_type(MATRIX A, int type);
int mat_free(MATRIX A);

/******************************************/
MATSTACK matstack_creat(int len);
MATSTACK __matstack_creat(int len);
int matstack_free(MATSTACK A);
MATSTACK matstack_append(MATSTACK s, MATRIX a);

MATVEC_DPOINTER matvec_creat(void);
int matvec_free(MATVEC_DPOINTER a);

/******************************************/
MATRIX mat_copy(MATRIX A, MATRIX result);
MATRIX mat_xcopy(MATRIX A, int si, int ei, int sj, int ej, MATRIX result);
MATRIX mat_xjoin(MATRIX A11, MATRIX A12, MATRIX A21, MATRIX A22, MATRIX result);
MATRIX mat_rowcopy(MATRIX A, int rowa, int rowb, MATRIX result);
MATRIX mat_colcopy(MATRIX A, int cola, int colb, MATRIX result);
int mat_fgetmat(MATRIX A, MAT_FILEPOINTER fp);

/******************************************/
/* matrix dump */
void mat_dump(MATRIX A);
void mat_dumpf(MATRIX A, const char *s);
void mat_fdump(MATRIX A, MAT_FILEPOINTER fp);
void mat_fdumpf(MATRIX A, const char *s, MAT_FILEPOINTER fp);

void int_vec_dump(INT_VECTOR a);
void int_vec_dumpf(INT_VECTOR a, const char *s);
void int_vec_fdump(INT_VECTOR a, MAT_FILEPOINTER fp);
void int_vec_fdumpf(INT_VECTOR a, const char *s, MAT_FILEPOINTER fp);

/******************************************/
/* vector manipulations */
INT_VECTOR int_vec_copy(INT_VECTOR a, INT_VECTOR result);
INT_VECTOR int_vec_unique(INT_VECTOR a);
INT_VECTOR int_vec_append(INT_VECTOR a, int i);
INT_VECTOR int_vec_find(INT_VECTOR a, int rel_type, int n);
INT_VECTOR int_vec_concat(INT_VECTOR a, INT_VECTOR b, INT_VECTOR result);
INT_VECTOR mat_get_sub_vector(INT_VECTOR a, INT_VECTOR indices);

/******************************************/
/* error handling functions */
int gen_error(int err_);
INT_VECTOR int_vec_error(int err_);
INT_VECSTACK int_vecstack_error (int err_);
MATRIX mat_error(int err_);
MATSTACK matstack_error(int err_);
int stack_error(int err_);
int queue_error(int err_);
int pq_error(int err_);
int graph_error(int err_);

/******************************************/
/* basic matrix operations */
mtype mat_mean(MATRIX A);
MATRIX mat_mean_row(MATRIX A, MATRIX result);
MATRIX mat_mean_col(MATRIX A, MATRIX result);

mtype mat_sum(MATRIX A);
MATRIX mat_sum_row(MATRIX A, MATRIX result);
MATRIX mat_sum_col(MATRIX A, MATRIX result);

MATRIX mat_abs(MATRIX A, MATRIX result);

INT_VECTOR int_vec_add(INT_VECTOR A, INT_VECTOR B, INT_VECTOR result);
INT_VECTOR int_vec_adds(INT_VECTOR A, int s, INT_VECTOR result);
INT_VECTOR int_vec_sub(INT_VECTOR A, INT_VECTOR B, INT_VECTOR result);
INT_VECTOR int_vec_subs(INT_VECTOR A, int s, INT_VECTOR result);
INT_VECTOR int_vec_mul(INT_VECTOR A, INT_VECTOR B, INT_VECTOR result);
INT_VECTOR int_vec_muls(INT_VECTOR A, int s, INT_VECTOR result);
INT_VECTOR int_vec_div(INT_VECTOR A, INT_VECTOR B, INT_VECTOR result);
INT_VECTOR int_vec_divs(INT_VECTOR A, int s, INT_VECTOR result);

MATRIX mat_add(MATRIX A, MATRIX B, MATRIX result);
MATRIX mat_adds(MATRIX A, mtype s, MATRIX result);

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
MATRIX mat_inv(MATRIX A, MATRIX result);
MATRIX mat_pinv(MATRIX A, MATRIX result);
MATRIX mat_wpinv(MATRIX A, MATRIX w, MATRIX result);
MATRIX mat_reg_inv(MATRIX A, mtype r, MATRIX result);
MATRIX mat_symtoeplz(MATRIX R, MATRIX result);

/******************************************/
/* linear system equation solver */
int mat_lu(MATRIX A, MATRIX P);
void mat_backsubs1(MATRIX A, MATRIX B, MATRIX C, MATRIX P, int xcol);
MATRIX mat_lsolve(MATRIX A, MATRIX b, MATRIX result);
MATRIX mat_cholesky(MATRIX A, MATRIX result);
MATRIX mat_conjgrad(MATRIX A, MATRIX b, MATRIX x0, mtype tol, int miters, MATRIX result);

MATRIX mat_submat(MATRIX A, int i, int j, MATRIX result);
mtype mat_cofact(MATRIX A, int i, int j);
mtype mat_det(MATRIX A);
mtype mat_minor(MATRIX A, int i, int j);

MATSTACK mat_qr(MATRIX A, MATSTACK qr);
MATRIX mat_durbin(MATRIX R, MATRIX result);
MATRIX mat_lsolve_durbin(MATRIX A, MATRIX B, MATRIX result);

/******************************************/
/* sorting algorithms */
mtype mat_median(MATRIX A);
mtype mat_order_statistic(MATRIX A, int k);
void __mat_quicksort(MATRIX A, int l, int r, int offset, MATRIX ind);
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
mtype __mat_randfun(mtype (*fun)(mtype), mtype xmin, mtype xmax);
mtype __mat_rand(void);
mtype __mat_randn(void);
mtype __mat_randexp(mtype mu);
MATRIX mat_randperm(int m, int n, MATRIX result);
MATRIX mat_randperm_n(int n, MATRIX result);
INT_VECTOR int_vec_randperm(int n, INT_VECTOR result);

/******************************************/
/* linear regressions */
MATRIX mat_least_squares(MATRIX A, MATRIX Y, MATRIX result);
MATRIX mat_w_least_squares(MATRIX A, MATRIX Y, MATRIX w, MATRIX result);
MATRIX mat_rob_least_squares(MATRIX A, MATRIX Y, int lossfunc, MATRIX result);

MATRIX mat_linear_ls_fit(MATRIX A, MATRIX Y, int deg, MATRIX result);
MATRIX mat_robust_fit(MATRIX A, MATRIX Y, int deg, int lossfunc, MATRIX result);

/******************************************/
/* matrix manipulations */
MATRIX mat_concat(MATRIX A, MATRIX B , int dim);
MATRIX mat_get_sub_matrix_from_rows(MATRIX A, INT_VECTOR indices, MATRIX result);
MATRIX mat_get_sub_matrix_from_cols(MATRIX A, INT_VECTOR indices, MATRIX result);
MATRIX mat_pick_row(MATRIX A, int r, MATRIX result);
MATRIX mat_pick_col(MATRIX A, int c, MATRIX result);
INT_VECSTACK mat_find(MATRIX A, int rel_type, mtype x);

MATRIX mat_fliplr(MATRIX A, MATRIX result);
MATRIX mat_flipud(MATRIX A, MATRIX result);

/******************************************/
/* distance tools */
MATRIX mat_calc_dist_sq(MATRIX A, MATRIX d, MATRIX result);
INT_VECTOR mat_find_within_dist(MATRIX A, MATRIX d, mtype range);

void __mat_cart2pol(mtype x, mtype y, mtype *rho, mtype *th);
void __mat_pol2cart(mtype rho, mtype th, mtype *x, mtype *y);
MATRIX mat_cart2pol(MATRIX A, int dim, MATRIX result);
MATRIX mat_pol2cart(MATRIX A, int dim, MATRIX result);

/******************************************/
/* function tools */
mtype __mat_addfunc(mtype x, mtype y);
mtype __mat_subfunc(mtype x, mtype y);
mtype __mat_mulfunc(mtype x, mtype y);
mtype __mat_divfunc(mtype x, mtype y);
mtype __mat_sqrfunc(mtype x);
mtype __mat_sqrtfunc(mtype x);


mtype __mat_huber_wt(mtype x, mtype k);
mtype __mat_bisquare_wt(mtype x, mtype k);
mtype __mat_logplusone(mtype x);
mtype __mat_arcsinh(mtype x);
mtype __mat_arccosh(mtype x);
mtype __mat_arctanh(mtype x);

MATRIX mat_bisquare_wt(MATRIX A, mtype k, mtype sigma, MATRIX result);
MATRIX mat_huber_wt(MATRIX A, mtype k, mtype sigma, MATRIX result);
MATRIX mat_gfunc(MATRIX A, mtype (*pt2func)(mtype), MATRIX result);
MATRIX mat_bsxfun(MATRIX A, MATRIX B, mtype (*func)(mtype, mtype), MATRIX result);

/******************************************/
/* statistical analysis tools */
MATSTACK mat_corcol(MATRIX data);
MATSTACK mat_covcol(MATRIX data);
MATRIX mat_scpcol(MATRIX data);
void mat_tred2(MATRIX a, MATRIX d, MATRIX e);
void mat_tqli(MATRIX d, MATRIX e, MATRIX z);
MATSTACK mat_pca(MATRIX data, int pca_type);
MATSTACK mat_eig_sym(MATRIX symmat, MATSTACK result);

/******************************************/
/* print helpers */
void mat_nextline(void);
void mat_fnextline(MAT_FILEPOINTER fp);

/******************************************/
/* signal domain transforms */
int __mat_powerof2(int width, int *m, int *twopm);

MATSTACK mat_fft2(MATSTACK c, int dir, MATSTACK result);
int __mat_fft(int dir, int m, mtype *x, mtype *y);

/* filtering functions */
MATRIX mat_conv2(MATRIX A, MATRIX mask, MATRIX scratch, MATRIX result);

/******************************************/
/* matrix vector conversions */
INT_VECTOR mat_2int_vec(MATRIX a);
MATRIX int_vec2_mat(INT_VECTOR a, int dir);
MATRIX mat_vectorize(MATRIX a, MATRIX result);
MATRIX mat_vectorize_tr(MATRIX a, MATRIX result);

/******************************************/
/* integration functions */
mtype mat_int_trapezoid(mtype (*func)(mtype), int n, mtype lower, mtype upper);
mtype mat_int_simpson(mtype (*func)(mtype), int n, mtype lower, mtype upper);
mtype __mat_lint(mtype *x, mtype (*func) (mtype), mtype x0, mtype xn, mtype f0, mtype f2, mtype f3, mtype f5, mtype f6, mtype f7, mtype f9, mtype fl4, mtype hmin, mtype hmax, mtype re, mtype ae);
mtype mat_int_qadrat(mtype(*func)(mtype), mtype lower, mtype upper);

/******************************************/
/* polynomials */
MATRIX mat_poly_eval(MATRIX A, mtype x, int dir, MATRIX result);
MATRIX mat_poly_diff(MATRIX A, int dir, MATRIX result);
MATRIX mat_poly_diff_eval(MATRIX A, mtype x, int dir, MATRIX result);
MATRIX mat_poly_add(MATRIX A, MATRIX B, MATRIX result);
MATRIX mat_poly_mul(MATRIX A, MATRIX B, MATRIX result);
MATSTACK mat_poly_div(MATRIX A, MATRIX B, MATSTACK result);
MATRIX mat_poly_scale(MATRIX A, mtype s, MATRIX result);
MATRIX mat_poly_shift(MATRIX A, int s, MATRIX result);
void mat_cheby_init();
void mat_legendre_init();
void mat_binom_init();
MATRIX mat_cheby(int n);
MATRIX mat_legendre(int n);
mtype mat_binom(int n, int k);
MATRIX mat_cheby_coeffs_to_poly(MATRIX coeffs, MATRIX result);
MATRIX mat_cheby_approx(mtype (*f)(mtype), mtype a, mtype b, int n, MATRIX result);

/******************************************/
/* pattern recognition */
MAT_BAYES_MODEL mat_bayes_model_creat(void);
int mat_bayes_model_free(MAT_BAYES_MODEL a);
MAT_PERCEPTRON mat_perceptron_creat(void);
int mat_perceptron_free(MAT_PERCEPTRON a);

MAT_BAYES_MODEL mat_bayes_classifier_train(MATRIX data, INT_VECTOR labels);
INT_VECTOR mat_bayes_classifier_test(MATRIX data, MAT_BAYES_MODEL b_model);
MAT_PERCEPTRON mat_perceptron_train(MATRIX data, INT_VECTOR labels, int num_of_iterations);
INT_VECTOR mat_perceptron_test(MATRIX data, MAT_PERCEPTRON p_model);
MAT_PERCEPTRON mat_perceptron_train_(MATRIX data1, MATRIX data2, MAT_PERCEPTRON p_model, int class_num);
MATVEC_DPOINTER mat_kmeans(MATRIX data, int k, int iters, MATVEC_DPOINTER result);

/******************************************/
/* basic search algorithms */
MAT_TREE mat_bs_make_null(void);
MAT_TREE mat_bs_free(MAT_TREE T);
MAT_TREE mat_bs_find(mtype x, MAT_TREE T);
MAT_TREE mat_bs_find_min(MAT_TREE T);
MAT_TREE mat_bs_find_max(MAT_TREE T);
MAT_TREE mat_bs_insert(mtype x, MAT_TREE T);
MAT_TREE mat_bs_delete(mtype x, MAT_TREE T);
int mat_bs_inorder(MAT_TREE T, int index, mtype **ordered);


/******************************************/
/* relational functions not reqd */
int gen_gt(mtype a);
int gen_lt(mtype a);
int gen_eq(mtype a);
mtype gen_abs_ceil(mtype a);

/******************************************/
/* textfunctions */
int mat_isnumeric(MAT_FILEPOINTER fp);
int mat_go_next_word(MAT_FILEPOINTER fp);
int mat_count_words_in_line(MAT_FILEPOINTER fp, int *count);
int mat_read_word(MAT_FILEPOINTER fp, char *c_word);
MATRIX mat_dlmread(const char *fname);
void mat_dlmwrite(const char *fname, MATRIX A);

/******************************************/
/* mat_timer_functions */
void mat_tic(void);
double mat_toc(void);
void mat_toc_print(void);

/******************************************/

/* data structures */
MAT_INT_STACK mat_int_stack_creat(void);
int mat_int_stack_free(MAT_INT_STACK s);
void mat_int_stack_push(MAT_INT_STACK s, int value);
int mat_int_stack_pop(MAT_INT_STACK s);
int mat_int_stack_is_empty(MAT_INT_STACK s);

MAT_MTYPE_STACK mat_mtype_stack_creat(void);
int mat_mtype_stack_free(MAT_MTYPE_STACK s);
void mat_mtype_stack_push(MAT_MTYPE_STACK s, mtype value);
mtype mat_mtype_stack_pop(MAT_MTYPE_STACK s);
int mat_mtype_stack_is_empty(MAT_MTYPE_STACK s);

MAT_INT_QUEUE mat_int_queue_creat(void);
int mat_int_queue_free(MAT_INT_QUEUE s);
void mat_int_queue_enqueue(MAT_INT_QUEUE s, int value);
int mat_int_queue_dequeue(MAT_INT_QUEUE s);
int mat_int_queue_is_empty(MAT_INT_QUEUE s);

MAT_MTYPE_QUEUE mat_mtype_queue_creat(void);
int mat_mtype_queue_free(MAT_MTYPE_QUEUE s);
void mat_mtype_queue_enqueue(MAT_MTYPE_QUEUE s, mtype value);
mtype mat_mtype_queue_dequeue(MAT_MTYPE_QUEUE s);
int mat_mtype_queue_is_empty(MAT_MTYPE_QUEUE s);

MAT_INT_PRIORITYQUEUE mat_int_priorityqueue_creat(void);
void mat_int_priorityqueue_enqueue(MAT_INT_PRIORITYQUEUE H, int data, int priority);
int mat_int_priorityqueue_dequeue(MAT_INT_PRIORITYQUEUE H);
int mat_int_priorityqueue_free(MAT_INT_PRIORITYQUEUE H);
int mat_int_priorityqueue_update(MAT_INT_PRIORITYQUEUE H, int data, int priority, int type);
int mat_int_priorityqueue_is_empty(MAT_INT_PRIORITYQUEUE H);

/******************************************/
/* mds */
MATRIX mat_mds(MATRIX d, int dims, int type, MATRIX result);
MATRIX __mat_mds_metric(MATRIX d, int dims, MATRIX result);
MATRIX __mat_mds_nonmetric(MATRIX d, int dims, MATRIX result);

/******************************************/
/* graph */
MAT_GRAPH mat_graph_creat(void);
void mat_graph_adjlist(MAT_GRAPH g, int directed, int weighted, MAT_FILEPOINTER fp);
MAT_INT_QUEUE mat_graph_search(MAT_GRAPH g, int connected, int mst);
void mat_graph_visit(MAT_GRAPH g, int k, int connected, int mst, MAT_INT_PRIORITYQUEUE pq, MAT_INT_QUEUE q);
void mat_graph_dumpf(MAT_GRAPH g, int mst, MAT_FILEPOINTER fp);
void mat_graph_dump(MAT_GRAPH g, int mst);
void mat_graph_adjm_to_adjl(MAT_GRAPH g, MATRIX a);
MAT_GRAPH mat_graph_reverse(MAT_GRAPH g, MAT_GRAPH r);

/******************************************/
/* kdtree */
MAT_KDTREE mat_kdtree_make_tree(MATRIX A, MAT_KDTREE result);
int mat_kdtree_free(MAT_KDTREE t);
MATRIX mat_kdtree_nearest(MAT_KDTREE t, MATRIX A, MATRIX result);
MAT_KDNODE __mat_kdtree_make_tree(MAT_KDNODE t, int len, int i, int dim);
MAT_KDNODE __mat_kd_find_median(MAT_KDNODE kd_start, MAT_KDNODE kd_end, int idx);
void __mat_kdtree_nearest(MAT_KDNODE root, MAT_KDNODE nd, int i, int dim, MAT_KDNODE *best, mtype *best_dist);

/******************************************/
/* sparse */
MATSTACK mat_omp(MATRIX A, MATRIX b, int k, mtype tol, MATSTACK result);

#ifdef __cplusplus
}
#endif


#endif

