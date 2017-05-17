#include "matrix.h"
#include <ctype.h>
#include <stdlib.h>


/** \brief Checks if current word in an opened file is numeric or not
 *
 * \param[in] fp Pointer to an opened file
 * \return Zero/non-zero
 *
 */

int mat_isnumeric(MAT_FILEPOINTER fp)
{
    char flag = 0;
    int ch;
    while((ch = getc(fp))!=EOF)
    {
        if(ch!=' ')
        {
            if((ch>=48) &&(ch<=57))
            {
                ungetc(ch, fp);
                if(flag==1) ungetc('.', fp);
                if(flag==2) ungetc('-', fp);
                if(flag==3)
                {
                    ungetc('.', fp);
                    ungetc('-', fp);
                }
                return 1;
            }
            else
            {
                if(ch=='.')
                {
                    if(flag==2) flag = 3;
                    else flag = 1;
                }
                else if(ch=='-') flag = 2;

                else
                {
                    ungetc(ch, fp);
                    return 0;
                }
            }
        }
        else
        {
            if(flag>0) /* inserted condition for initial space in line */
            {
                ungetc(ch, fp);
                return 0;
            }
        }
    }
    if(ch==EOF)
    {
        ungetc(ch, fp); /*edited just now */

        return -1;
    }
    return 0;
}

/** \brief Moves to next word in an opened file
 *
 * \param[in] fp Pointer to an opened file
 * \return EOF reached
 *
 */

int mat_go_next_word(MAT_FILEPOINTER fp)
{
    int flag = 0, ch;
    while((flag<2)&&((ch = getc(fp))!=EOF))
    {
        if((ch=='\v')||(ch=='\n')||(ch=='\t')||isspace(ch)||(ch==',')||(ch=='!')||(ch=='(')||(ch==')')||(ch=='{')||(ch=='}')||(ch=='[')||(ch==']'))
        {
            if(flag==0) flag = 1;
        }
        else if(flag==1) flag = 2;
    }
    if(ch!=EOF)
    {
        ungetc(ch, fp);
        return 0;
    }
    else return 1;
}

/** \brief Count words in current line in an opened file
 *
 * \param[in] fp Pointer to an opened file
 * \param[in] count Pointer to output count
 * \return EOF reached
 *
 */

int mat_count_words_in_line(MAT_FILEPOINTER fp, int *count)
{
    int flag = -1, ch;
    *count = 0;
    while((flag<3) && ((ch = getc(fp))!=EOF))
    {
        if((ch=='\v')||(ch=='\n'))
        {
            if(flag==0)
            {
                (*count)++;
                flag = 3;
            }
            if(flag==-1) flag = 4;/*  line included to handle empty line */
            else flag = 2;
        }
        else if(isspace(ch)||(ch=='\t')||(ch==',')||(ch=='!')||(ch=='(')||(ch==')')||(ch=='{')||(ch=='}')||(ch=='[')||(ch==']'))
        {
            if(flag==0)/* changed from  <=0 to ==0 to skip initial space */
            {
                flag = 1;
                (*count)++;
            }
        }
        else if(flag==1) flag = 0;
        else if(flag==2) flag = 3;
        else flag = 0;
    }
    if(flag!=-1&&flag!=4) ungetc(ch, fp);
    if(ch==EOF)
    {
        if(flag==0) (*count)++;
        return 1;
    }
    else return 0;
}

/** \brief Reads a matrix from a file
 *
 * \param[in] fname Filename to read from
 * \return Output matrix
 *
 */

MATRIX mat_dlmread(const char *fname)
{
    int m = 0, n = 0, i, j, flag = 0, tmp = 0, k = 0;
    mtype in_value = 0;
    char c_word[100];
    MAT_FILEPOINTER fp = NULL;
    MATRIX data = NULL;

    if((fp = fopen(fname,"rb"))==NULL) gen_error(GEN_FNOTOPEN);
    while(!flag)
    {
        k = mat_isnumeric(fp);
        if(k==1)
        {
            flag = mat_count_words_in_line(fp, &tmp);
            if(tmp!=0) n++; /*improved skip non-data line*/
        }
        else if(k==0)
        {
            flag = mat_count_words_in_line(fp, &tmp);
        }
        else flag = 1;
    }
    m = tmp;
    fclose(fp);
    if(m==0||n==0) return mat_error(MAT_FNOTGETMAT);
    fp = fopen(fname,"rb");
    data = mat_creat(n, m, ZERO_MATRIX);
    for(i=0; i<n; ++i)
    {
        k = mat_isnumeric(fp);
        if(k==0)
        {
            do
            {
                flag = mat_count_words_in_line(fp, &tmp);
                k = mat_isnumeric(fp);
            }
            while(k==0);
        }
        for(j=0; j<m; ++j)
        {
            mat_read_word(fp, c_word);
#if mtype_n == 0
            in_value = (float)strtod(c_word, NULL);
#elif mtype_n == 1
            in_value = strtod(c_word, NULL);
#endif
            /* fscanf(fp, "%lf", &in_value); */
            /* go_next_word(fp); */
            data[i][j] = in_value;
        }
    }
    fclose(fp);
    return data;
}

/** \brief Reads current word from an opened file
 *
 * \param[in] fp Pointer to an opened file
 * \param[in] c_word Pointer to word read
 * \return EOF reached
 *
 */

int mat_read_word(MAT_FILEPOINTER fp, char *c_word)
{
    int flag = 0, t = 0, ch;
    while((flag<3)&&((ch = getc(fp))!=EOF))/*no need for state 3 to be corrected*/
    {
        if((ch=='\v')||(ch=='\n')||(ch=='\t')||isspace(ch)||(ch==',')||(ch=='!')||(ch=='(')||(ch==')')||(ch=='{')||(ch=='}')||(ch=='[')||(ch==']'))
        {
            if(flag!=0) flag = 2;
        }
        else if(flag<2)
        {
            flag = 1;
            c_word[t++] = ch;
        }
        else if(flag==2) flag = 3; /* reached next word */ /*to be corrected for deleting state 3*/
    }
    c_word[t] = '\0';
    if(ch!=EOF)
    {
        ungetc(ch, fp);
        return 1;
    }
    else return 0;
}

/** \brief Writes a matrix to a file
 *
 * \param[in] fname Filename to write into
 * \param[in] A Input matrix
 *
 */

void mat_dlmwrite(const char *fname, MATRIX A)
{
    MAT_FILEPOINTER fp;
    fp = fopen(fname, "wb");
    if(fp==NULL) gen_error(GEN_FNOTOPEN);
    mat_fdump(A, fp);
}

