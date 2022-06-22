

/* type matrice*/

typedef struct
{
    int dimy,dimx;
    double *tab;
} matrice;

void mult_blas(matrice *A, matrice *B, matrice *C, short int trans)
{
    if(trans == 0)
    {/*
     * C->dimx = B->dimx;
     * C->dimy = A->dimy;
     * C->tab = mxMalloc(sizeof(double)*C->dimx*C->dimy);*/
        int m,n,k;
        m = A->dimy;
        n = B->dimx;
        k = A->dimx;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, 1.0, A->tab, m, B->tab, k, 0.0, C->tab, m);
    }else{
        int m,n,k;
        m = A->dimy;
        n = B->dimy;
        k = A->dimx;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                m, n, k, 1.0, A->tab, m, B->tab, n, 0.0, C->tab, m);
    }
}



double count_val(matrice * matr, double val)
{
    int i, j, s = 0;
    for(i=0;i<matr->dimx;i++)
    {
        for(j=0;j<matr->dimy;j++)
        {
            if(matr->tab[i*matr->dimy+j] == val){
                s++;
            }
        }
    }
    return s;
}

void remplir_matrice(matrice *matr, double val)
{
    int i, j;
    for(i=0;i<matr->dimx;i++)
    {
        for(j=0;j<matr->dimy;j++)
        {
            matr->tab[i*matr->dimy+j] = val;
            
            /*mexPrintf("element[%d][%d] = %f\n",j,i,matr[i*dimy+j]);*/
        }
    }
}

void remplir_matrice_rand(matrice *matr,int max_value,double pivot)
{
    int i, j;
    srand (time (NULL));
    for(i=0;i<matr->dimx;i++)
    {
        for(j=0;j<matr->dimy;j++)
        {
            matr->tab[i*matr->dimy+j] = (float)rand()/(float)(RAND_MAX/max_value);
        }
    }
}

/* rend la transpose de M*/
matrice transp(matrice *matr)
{
    matrice m;
    m.dimy = matr->dimx;
    m.dimx = matr->dimy;
    m.tab = mxMalloc(sizeof(double)*m.dimx*m.dimy);
    
    int i,j;
    for(i=0;i<m.dimx;i++)
    {
        for(j=0;j<m.dimy;j++)
        {
            m.tab[i*m.dimy + j] = matr->tab[i+ m.dimx*j];
        }
    }
    return m;
}


matrice mult_mat(matrice * m1, matrice * m2)
{
    int i,j,k;
    matrice m;
    if(m1->dimx == m2->dimy)
    {
        m.dimx = m2->dimx;
        m.dimy = m1->dimy;
        m.tab = mxMalloc(sizeof(m)*m.dimx*m.dimy);
        
        for(i=0; i<m.dimy; i++){
            for(j=0; j<m.dimx; j++){
                m.tab[i*m2->dimx + j] = 0;
                for(k=0; k<m1->dimx; k++)
                {
                    m.tab[i+ m.dimy*j]+=m1->tab[i*m1->dimx + k]*
                            m2->tab[j*m1->dimx + k];
                    
                }
            }
        }
        return m;
    }else{
        mexPrintf("mult: probleme de taille\n");
        return m;
    }
}



double sum_mat(matrice M)
{
    double s=0;
    int i,j;
    for (i=0;i<M.dimy;i++)
    {
        for(j=0;j<M.dimx;j++)
        {
            s = s + M.tab[i*M.dimx + j];
        }
    }
    return s;
}

/* rend la somme des mat1 et mat2 si add_sub >=0 et la diff sinon*/

matrice add_sub_mat(matrice *mat1, matrice *mat2, short int add_sub)
{
    matrice M;
    if (mat1->dimx == mat2->dimx && mat1->dimy == mat2->dimy)
    {
        M.dimx = mat1->dimx;
        M.dimy = mat1->dimy;
        M.tab = mxMalloc(sizeof(double)*M.dimx*M.dimy);
        
        int i,j;
        for (i=0;i<M.dimy;i++)
        {
            for (j=0;j<M.dimx;j++)
            {
                M.tab[i*M.dimx + j]=mat1->tab[i* M.dimx + j] + add_sub*mat2->tab[i* M.dimx + j];
            }
        }
    }else{
        mexPrintf("add: probleme de taille\n");
    }
    return M;
}

/* rend la somme des mat1 et mat2 si add_sub >=0 et la diff sinon*/

void add_mat(matrice *mat1, matrice *mat2, matrice *M, short int add_sub)
{
    int i,j;
    for (i=0;i<M->dimy;i++)
    {
        for (j=0;j<M->dimx;j++)
        {
            M->tab[i*M->dimx + j]=mat1->tab[i* M->dimx + j] +
                    add_sub*mat2->tab[i* M->dimx + j];
        }
    }
    
}

matrice abs_mat(matrice matr)
{
    int i, j;
    matrice M;
    M.dimx = matr.dimx;
    M.dimy = matr.dimy;
    M.tab = mxMalloc(sizeof(double)*M.dimx*M.dimy);
    for(i=0;i<matr.dimx*matr.dimy;i++)
    {
        M.tab[i*matr.dimy+j] = fabs(matr.tab[i*matr.dimy+j]);
        /*mexPrintf("element[%d][%d] = %f\n",j,i,matr[i*dimy+j]);*/
        
    }
    return M;
}

void sign_mat(matrice *matr, matrice * M)
{
    int i, j;
    for(i=0;i<matr->dimx*matr->dimy;i++)
    {
        if(matr->tab[i] < 0)
        {
            M->tab[i] = -1;
            
        }else if(matr->tab[i] > 0){
            M->tab[i] = 1;
        }else{
            M->tab[i] = 0;
        }
        
    }
}