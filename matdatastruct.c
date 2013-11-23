#include <limits.h>
#include <malloc.h>
#include "matrix.h"

SEARCH_TREE mat_bs_make_null(void)
{
    return NULL;
}

SEARCH_TREE mat_bs_free(SEARCH_TREE T)
{
    if(T!=NULL)
    {
        mat_bs_free(T->left);
        mat_bs_free(T->right);
        free(T);
    }
    return NULL;
}

SEARCH_TREE mat_bs_find( mtype x, SEARCH_TREE T)
{
    if(T==NULL) return NULL;
    if(x<T->element) return(mat_bs_find(x, T->left));
    else if(x>T->element) return(mat_bs_find(x, T->right));
    else return T;
}

SEARCH_TREE mat_bs_find_min(SEARCH_TREE T)
{
    if(T!=NULL)
        while(T->left!= NULL) T = T->left;
    return T;
}

SEARCH_TREE mat_bs_find_max(SEARCH_TREE T)
{
    if(T!=NULL)
        while(T->right!=NULL) T = T->right;
    return T;
}

SEARCH_TREE mat_bs_insert(mtype x, SEARCH_TREE T)
{
    if(T==NULL)
    {
        T = (SEARCH_TREE) malloc (sizeof(struct tree_node));
        if(T==NULL) gen_error(GEN_MALLOC);
        else
        {
            T->element = x;
            T->left = T->right = NULL;
        }
    }
    else if(x<T->element) T->left = mat_bs_insert(x, T->left);
    else if(x>T->element) T->right = mat_bs_insert(x, T->right);
    return T;
}

SEARCH_TREE mat_bs_delete(mtype x, SEARCH_TREE T)
{
    SEARCH_TREE tmp_cell, child = 0;
    if(T==NULL) gen_error(GEN_NOT_FOUND);
    else if(x<T->element) T->left = mat_bs_delete(x, T->left);
    else if(x>T->element) T->right = mat_bs_delete(x, T->right);
    else if(T->left && T->right)
    {
        tmp_cell = mat_bs_find_min(T->right);
        T->element = tmp_cell->element;
        T->right = mat_bs_delete(T->element, T->right);
    }
    else
    {
        tmp_cell = T;
        if(T->left==NULL) child = T->right;
        if(T->right==NULL) child = T->left;
        free(tmp_cell);
        return child;
    }
    return T;
}

int mat_bs_inorder(SEARCH_TREE T, int index, mtype **p_ordered)
{
    if(T!=NULL)
    {
        index = mat_bs_inorder(T->left, index, p_ordered);
        if((*p_ordered = (mtype *)realloc(*p_ordered, sizeof(mtype)*(index+ 1)))==NULL) gen_error(GEN_MALLOC);
        (*p_ordered)[index++] = T->element;
        index = mat_bs_inorder(T->right, index, p_ordered);
    }
    return index;
}

INT_STACK int_stack_creat(void)
{
    INT_STACK s;
    if((s = (INT_STACK)malloc(sizeof(struct int_stack)))==NULL) stack_error(STACK_MALLOC);
    s->p = 0;
    s->length = STACK_MAX;
    if((s->stack = (int*)malloc(sizeof(int)*(STACK_MAX)))==NULL) stack_error(STACK_MALLOC);
    return s;
}

int int_stack_free(INT_STACK s)
{
    if (s==NULL) return 0;
    free(s->stack);
    return 1;
}

void int_stack_push(INT_STACK s, int value)
{
    if(s->p>=s->length)
    {
        if((s->stack = (int *)realloc(s->stack, sizeof(int)*(s->length+STACK_MAX)))==NULL) stack_error(STACK_MALLOC);
        s->length += STACK_MAX;
    }
    s->stack[s->p++] = value;
}

int int_stack_pop(INT_STACK s)
{
    if(s->p>0)
    {
        return s->stack[(--s->p)];
    }
    else return stack_error(STACK_EMPTY);
}

int int_stack_is_empty(INT_STACK s)
{
    return ((int)(s->p==0));
}

MAT_MTYPE_STACK mat_mtype_stack_creat(void)
{
    MAT_MTYPE_STACK s;
    if((s= (MAT_MTYPE_STACK)malloc(sizeof(struct mat_mtype_stack)))==NULL) stack_error( STACK_MALLOC);
    s->p = 0;
    s->length = STACK_MAX;
    if((s->stack=(mtype*)malloc(sizeof(mtype)*(STACK_MAX)))==NULL) stack_error( STACK_MALLOC);
    return s;
}

int mat_mtype_stack_free(MAT_MTYPE_STACK s)
{
    if(s==NULL) return 0;
    free(s->stack);
    return 1;
}

void mat_mtype_stack_push(MAT_MTYPE_STACK s, mtype value)
{
    if(s->p>=s->length)
    {
        if((s->stack = (mtype *)realloc(s->stack, sizeof(mtype)*(s->length+STACK_MAX)))==NULL) stack_error(STACK_MALLOC);
        s->length += STACK_MAX;
    }
    s->stack[s->p++] =  value;
}

mtype mat_mtype_stack_pop(MAT_MTYPE_STACK s)
{
    if(s->p>0)
    {
        return s->stack[(--s->p)];
    }
    else stack_error(STACK_EMPTY);
    return 0;
}

int mat_mtype_stack_is_empty(MAT_MTYPE_STACK s)
{
    return ((int)(s->p==0));
}

INT_QUEUE int_queue_creat(void)
{
    INT_QUEUE s;
    if((s = (INT_QUEUE)malloc(sizeof(struct int_queue)))==NULL) queue_error(QUEUE_MALLOC);
    s->p = 0;
    s->head = NULL;
    s->tail = NULL;
    return s;
}

int int_queue_free(INT_QUEUE s)
{
    if(s==NULL) return 0;
    if(s->head!=NULL)
    {
        QINTDATA t;
        while(s->head!=NULL)
        {
            t = s->head;
            s->head = s->head->next;
            free(t);
        }
        s->tail = NULL;
    }
    return 1;
}

void int_queue_enqueue(INT_QUEUE s, int value)
{
    QINTDATA t;
    if((t = (QINTDATA)malloc(sizeof(qintdata)))==NULL) queue_error(QUEUE_MALLOC);
    t->next=NULL;
    t->data = value;
    if(s->head==NULL)
    {
        s->tail = t;
        s->head = s->tail;
    }
    else
    {
        s->tail->next = t;
        s->tail =t;
    }
    ++(s->p);
}

int int_queue_dequeue(INT_QUEUE s)
{
    int value;
    QINTDATA t;
    if(s->head==NULL) queue_error(QUEUE_EMPTY);
    t = s->head->next;
    value = s->head->data;
    free(s->head);
    --(s->p);

    if(t==NULL)
    {
        s->head = NULL;
        s->tail = NULL;
    }
    else s->head = t;
    return value;
}

int int_queue_is_empty(INT_QUEUE s)
{
    return ((int)(s->head==NULL));
}

MAT_MTYPE_QUEUE mat_mtype_queue_creat(void)
{
    MAT_MTYPE_QUEUE s;
    if((s = (MAT_MTYPE_QUEUE)malloc(sizeof(struct mat_mtype_queue)))==NULL) queue_error(QUEUE_MALLOC);
    s->p = 0;
    s->head = NULL;
    s->tail = NULL;
    return s;
}

int mat_mtype_queue_free(MAT_MTYPE_QUEUE s)
{
    if(s == NULL) return 0;
    if(s->head!= NULL)
    {
        QMDATA t;
        while(s->head!=NULL)
        {
            t = s->head;
            s->head = s->head->next;
            free(t);
        }
        s->tail = NULL;
    }
    return 1;
}

void mat_mtype_queue_enqueue(MAT_MTYPE_QUEUE s, mtype value)
{
    QMDATA t;
    if((t = (QMDATA)malloc(sizeof(qmdata)))==NULL) queue_error(QUEUE_MALLOC);
    t->next = NULL;
    t->data = value;
    if(s->head==NULL)
    {
        s->tail = t;
        s->head = s->tail;
    }
    else
    {
        s->tail->next = t;
        s->tail = t;
    }
    ++(s->p);
}

mtype mat_mtype_queue_dequeue(MAT_MTYPE_QUEUE s)
{
    mtype value;
    QMDATA t;
    if(s->head == NULL) queue_error(QUEUE_EMPTY);
    t = s->head->next;
    value = s->head->data;
    free(s->head);
    --(s->p);
    if(t==NULL)
    {
        s->head = NULL;
        s->tail = NULL;
    }
    else s->head = t;
    return value;
}

int mat_mtype_queue_is_empty(MAT_MTYPE_QUEUE s)
{
    return ((int)(s->head==NULL));
}

#ifndef INT_MAX
#define INT_MAX 2147483647
#define INT_MIN (-INT_MAX-1)
#endif

INT_PRIORITYQUEUE int_priorityqueue_creat()
{
    INT_PRIORITYQUEUE H;
    if((H = (INT_PRIORITYQUEUE) malloc(sizeof(int_priorityqueue)))==NULL) pq_error(PQ_MALLOC);
    if((H->element = (PQINTDATA) malloc((STACK_MAX+1)*sizeof(pqintdata)))==NULL) pq_error(PQ_MALLOC);
    H->length = STACK_MAX;
    H->p = 0;
    H->element[0].priority = INT_MAX;
    H->element[0].data = 0;
    return H;
}

void int_priorityqueue_enqueue(INT_PRIORITYQUEUE H, int data, int priority)
{
    int i;
    if(H->length==H->p)
    {
        if((H->element = (PQINTDATA)realloc(H->element, sizeof(pqintdata)*(H->length+STACK_MAX+1)))==NULL) pq_error(PQ_MALLOC);
        H->length += STACK_MAX;
    }
    i = ++(H->p);
    while(H->element[i/2].priority<priority)
    {
        H->element[i] = H->element[i/2];
        i /= 2;
    }
    H->element[i].priority = priority;
    H->element[i].data = data;
}

int int_priorityqueue_dequeue(INT_PRIORITYQUEUE H)
{
    int i, child;
    pqintdata min_element, last_element;
    if(H->p==0)
    {
        pq_error(PQ_EMPTY);
        return H->element[0].data;
    }
    min_element = H->element[1];
    last_element = H->element[H->p--];
    for(i=1; i*2<=H->p; i=child)
    {
        child = i*2;
        if((child!=H->p) && (H->element[child+1].priority>H->element[child].priority)) ++child;
        if(last_element.priority<H->element[child].priority) H->element[i] = H->element[child];
        else break;
    }
    H->element[i] = last_element;
    return min_element.data;
}

int int_priorityqueue_free(INT_PRIORITYQUEUE H)
{
    if(H==NULL) return 0;
    free(H->element);
    H->element = NULL;
    free(H);
    return 1;
}

int int_priorityqueue_update(INT_PRIORITYQUEUE H, int data, int priority, int type)
{
    int i, child, p = 0;
    for(i=1; i<=H->p; ++i) if(H->element[i].data==data) p = i;
    if(p!=0)
    {
        if((type==0 || type==2) && H->element[p].priority>priority)
        {
            for(i=p; (i*2)<=H->p; i=child)
            {
                child = i*2;
                if((child!=H->p)&&(H->element[child+1].priority<H->element[child].priority)) ++child;
                if(priority>H->element[child].priority) H->element[i] = H->element[child];
                else break;
            }
            H->element[i].data = data;
            H->element[i].priority = priority;
            return -1;
        }
        else if((type<2) && H->element[p].priority<priority)
        {
            while(H->element[p/2].priority <priority)
            {
                H->element[p] = H->element[p/2];
                p/= 2;
            }
            H->element[p].priority = priority;
            H->element[p].data = data;
            return 1;
        }
        else return 0;
    }
    else
    {
        int_priorityqueue_enqueue(H, data, priority);
        return 2;
    }
}

int int_priorityqueue_is_empty(INT_PRIORITYQUEUE H)
{
    return ((int)(H->p==0));
}

