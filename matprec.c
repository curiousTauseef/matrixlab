#include <stdlib.h>
#include "matrix.h"


MAT_BAYES_MODEL mat_bayes_classifier_train(MATRIX data, INT_VECTOR labels)
{
    int i;
    MAT_BAYES_MODEL b_model = mat_bayes_model_creat();
    MATRIX subdata = NULL;
    INT_VECTOR tmp = NULL;
    MATSTACK covstack = NULL;
    if(data==NULL) gen_error(GEN_NOT_FOUND);
    b_model->class_labels = int_vec_unique(labels);
    b_model->num_of_classes = Int_VecLen(b_model->class_labels);
    b_model->num_of_features = MatCol(data);

    b_model->class_means = matstack_creat(b_model->num_of_classes);
    b_model->class_covars = matstack_creat(b_model->num_of_classes);

    for(i=0; i<b_model->num_of_classes; ++i)
    {
        tmp = int_vec_find(labels, GEN_EQUAL_TO, b_model->class_labels[i]);
        subdata = mat_get_sub_matrix_from_rows(data, tmp, NULL);
        covstack = mat_covcol(subdata);
        b_model->class_means[i] = covstack[0];
        b_model->class_covars[i] = covstack[1];
        mat_free(subdata);
        int_vec_free(tmp);
    }
    return b_model;
}

INT_VECTOR mat_bayes_classifier_test(MATRIX data, MAT_BAYES_MODEL b_model)
{
    int n, i, j;
    mtype *w0 = NULL, tmp0;
    INT_VECTOR results = NULL;
    MATSTACK cov_invs = NULL;
    MATSTACK ws = NULL;
    MATRIX tmp1 = NULL, tmp2 = NULL, tmp3 = NULL, results_raw = NULL;
    MATVEC_DPOINTER p = NULL;
    n = MatRow(data);
    results_raw = mat_creat(b_model->num_of_classes, n, UNDEFINED);
    results = int_vec_creat( n, UNDEFINED);

    cov_invs = matstack_creat(b_model->num_of_classes);
    ws = matstack_creat(b_model->num_of_classes);
    w0 = (mtype *)malloc(sizeof(mtype)*b_model->num_of_classes);

    for(i=0; i<b_model->num_of_classes; ++i)
    {
        tmp0 = mat_det(b_model->class_covars[i]);
        if(tmp0==0)
        {
            cov_invs[i] = mat_reg_inv(b_model->class_covars[i], (mtype)EPS, NULL);
            tmp0 = (mtype)eps;
        }
        else cov_invs[i] = mat_inv(b_model->class_covars[i], NULL);
        tmp1 = mat_tran(b_model->class_means[i], tmp1);
        ws[i] = mat_mul(cov_invs[i], tmp1, NULL);
        tmp2 = mat_mul(b_model->class_means[i], ws[i], NULL);
        w0[i] = - tmp2[0][0] - (mtype)log(fabs(tmp0));
        if(b_model->class_priors!=NULL) w0[i] += (mtype)(2*log(b_model->class_priors[0][i]));
    }
    mat_free(tmp1);


    /* code below and above need to be optimized */
    tmp1 = NULL;
    tmp2 = NULL;
    for(i=0; i<b_model->num_of_classes; ++i)
    {
        tmp1 = mat_mul(data, ws[i], tmp1);
        tmp1 = mat_muls(tmp1, 2, tmp1);
        tmp2 = mat_mul(data, cov_invs[i], tmp2);
        tmp2 = mat_mul_dot(tmp2, data, tmp2);
        tmp3 = mat_sum_row(tmp2, tmp3);
        tmp1 = mat_sub(tmp1, tmp3, tmp1);
        tmp1 = mat_adds(tmp1, w0[i], tmp1);
        for(j=0; j<n; ++j) results_raw[i][j] = tmp1[j][0];
    }
    mat_free(tmp1);
    mat_free(tmp2);
    mat_free(tmp3);
    p = mat_max(results_raw, 1);
    for (i=0; i<n; ++i) results[i] = b_model->class_labels[((INT_VECTOR)p[1])[i]];
    mat_free(results_raw);
    matstack_free(cov_invs);
    matstack_free(ws);
    matvec_free(p);
    free(w0);
    return results;
}

MAT_PERCEPTRON mat_perceptron_train(MATRIX data, INT_VECTOR labels, int num_of_iterations)
{
    int n, i, tmp;
    MATRIX data_ = NULL;
    MAT_PERCEPTRON p_model = mat_perceptron_creat();
    MATRIX subdata1, subdata2;
    INT_VECTOR tmp1 = NULL, tmp2 = NULL;
    if(num_of_iterations<1) p_model->num_of_iterations = num_of_iterations;
    if(data==NULL) gen_error(GEN_NOT_FOUND);
    p_model->class_labels = int_vec_unique(labels);
    p_model->num_of_classes = Int_VecLen(p_model->class_labels);
    if(p_model->num_of_classes<2)
    {
        p_model->istrained = 0;
        return p_model;
    }
    else p_model->istrained = 1;
    p_model->num_of_features = MatCol(data);
    n = MatRow(data);
    subdata1 = mat_creat(n, 1, ONES_MATRIX);
    data_ = mat_concat(data, subdata1, 2);
    mat_free(subdata1);
    if(p_model->num_of_classes==2)
    {
        tmp = 1;
        p_model->class_weights = mat_creat(MatCol(data)+1, 1, UNDEFINED);
    }
    else
    {
        tmp = p_model->num_of_classes;
        p_model->class_weights = mat_creat(MatCol(data)+1, p_model->num_of_classes, UNDEFINED);
    }

    for(i=0; i<tmp; ++i)
    {
        tmp1 = int_vec_find(labels, GEN_EQUAL_TO, p_model->class_labels[i]);
        tmp2 = int_vec_find(labels, GEN_NOT_EQUAL_TO, p_model->class_labels[i]);
        subdata1 = mat_get_sub_matrix_from_rows(data_, tmp1, NULL);
        subdata2 = mat_get_sub_matrix_from_rows(data_, tmp2, NULL);

        p_model = mat_perceptron_train_(subdata1, subdata2, p_model, i);

        mat_free(subdata1);
        int_vec_free(tmp1);
        mat_free(subdata2);
        int_vec_free(tmp2);
    }
    mat_free(data_);
    return p_model;
}

MAT_PERCEPTRON mat_perceptron_train_(MATRIX data1, MATRIX data2, MAT_PERCEPTRON p_model, int class_num)
{
    mtype alpha;
    MATRIX weight = NULL, er_data = NULL, tmp1 = NULL, tmp2 = NULL, tmp3 = NULL, tmp21 = NULL, tmp22 = NULL;
    INT_VECSTACK indices1 = NULL, indices2 = NULL;
    int err_num, iters = 0, tn_data, i;
    err_num = MatRow(data1) + MatRow(data2);
    alpha = 1.0f/((mtype)err_num);
    tn_data = err_num;

    weight = mat_randn(MatCol(data1), 1, NULL);
    while(err_num>0 && iters<(p_model->num_of_iterations*tn_data))
    {
        tmp1 = mat_mul(data1, weight, tmp1);
        indices1 = mat_find(tmp1, GEN_GREATER_THAN, 0);
        tmp2 = mat_mul(data2, weight, tmp2);
        indices2 = mat_find(tmp2, GEN_LESS_THAN_EQUAL_TO, 0);
        err_num = ((indices1[0][0]>-1)?(Int_VecLen(indices1[0])):0)
                  +((indices2[0][0]>-1)?(Int_VecLen(indices2[0])):0);
        if(indices1[0][0]==-1 && indices2[0][0]==-1) break;
        else if(indices2[0][0]==-1)
        {
            tmp21 = mat_get_sub_matrix_from_rows(data1, indices1[0], NULL);
            tmp21 = mat_muls(tmp21, -1, tmp21);
            er_data = mat_sum_col(tmp21, NULL);
            mat_free(tmp21);
        }
        else if(indices1[0][0]==-1)
        {
            tmp21 = mat_get_sub_matrix_from_rows(data2, indices2[0], NULL);
            er_data = mat_sum_col(tmp21, NULL);
            mat_free(tmp21);
        }
        else
        {
            tmp21 = mat_get_sub_matrix_from_rows(data1, indices1[0], NULL);
            tmp22 = mat_get_sub_matrix_from_rows(data2, indices2[0], NULL);
            tmp21 = mat_muls(tmp21, -1,tmp21);
            tmp3 = mat_concat(tmp21, tmp22, COLS);
            er_data = mat_sum_col(tmp3, NULL);
            mat_free(tmp21);
            mat_free(tmp22);
            mat_free(tmp3);
        }
        er_data = mat_muls(er_data, alpha, er_data);
        tmp3 = mat_tran(er_data, NULL);

        weight = mat_add(weight, tmp3, weight);
        weight = mat_divs(weight, (float)(mat_norm_p(weight,2)+Eps), weight);
        mat_free(tmp3);
        alpha = 0.985f * alpha;
        mat_free(er_data);
        int_vecstack_free(indices1);
        int_vecstack_free(indices2);
        iters++;
    }
    mat_free(tmp1);
    mat_free(tmp2);

    for(i=0; i<=(p_model->num_of_features); i++) p_model->class_weights[i][class_num] = weight[i][0];
    mat_free(weight);

    return p_model;
}

INT_VECTOR mat_perceptron_test(MATRIX data, MAT_PERCEPTRON p_model)
{
    int n, i;
    MATRIX tmp0, data00;
    MATVEC_DPOINTER p = NULL;
    INT_VECTOR results = NULL;
    if(data ==NULL||p_model==NULL||p_model->istrained==0) gen_error(GEN_NOT_FOUND);
    n = MatRow(data);
    tmp0 = mat_creat(n, 1, ONES_MATRIX);
    data00 = mat_concat(data, tmp0, 2);
    mat_free(tmp0);
    tmp0 = mat_mul(data00, p_model->class_weights, NULL);
    results = int_vec_creat( n, UNDEFINED);
    if(p_model->num_of_classes==2)
    {
        for(i=0; i<n; ++i)
        {
            if(tmp0[i][0]<0) results[i] = p_model->class_labels[0];
            else results[i] = p_model->class_labels[1];
        }
    }
    else
    {
        p = mat_min(tmp0, 2);
        mat_free(tmp0);
        for(i=0; i<n; ++i) results[i] = p_model->class_labels[((INT_VECTOR)p[1])[i]];
        matvec_free(p);
    }
    mat_free(data00);
    return results;
}

MATVEC_DPOINTER mat_kmeans(MATRIX data, int k, int iters, MATVEC_DPOINTER result)
{
    MATRIX data_centers = NULL, sdata = NULL, curr_distances = NULL, curr_diff = NULL, curr_point = NULL;
    INT_VECTOR data_membership = NULL, curr_cluster_indices = NULL;
    MATVEC_DPOINTER mv = NULL;
    int d, l, i, N;
    d = MatRow(data);
    N = MatCol(data);
    if(result==NULL)
    {
        if((result = matvec_creat())==NULL) mat_error(MAT_MALLOC);
        result[0] = mat_rand(d, k, NULL);
        result[1] = int_vec_creat(N, UNDEFINED);
    }
    data_centers = (MATRIX)result[0];
    data_membership = (INT_VECTOR)result[1];

    for(l=0; l<iters; ++l)
    {
        for(i=0; i<N; ++i)
        {
            curr_point = mat_pick_col(data, i, curr_point);
            curr_diff = mat_bsxfun(data_centers, curr_point, __mat_subfunc, curr_diff);
            curr_diff = mat_gfunc(curr_diff, __mat_sqrfunc, curr_diff);
            curr_distances = mat_sum_col(curr_diff, curr_distances);
            mv = mat_min(curr_distances, COLS);
            data_membership[i]= ((INT_VECTOR)mv[1])[0];
            matvec_free(mv);
        }
        for(i=0; i<k; ++i)
        {
            curr_cluster_indices = int_vec_find(data_membership, GEN_EQUAL_TO, i);
            sdata = mat_get_sub_matrix_from_cols(data, curr_cluster_indices, NULL);

            curr_point = mat_mean_row(sdata, curr_point);
            mat_free(sdata);
            int_vec_free(curr_cluster_indices);
            data_centers = mat_colcopy(curr_point, 0, i, data_centers);
        }
    }
    return result;
}


