// #include <mex.h>
#include <math.h>
#include <float.h>


union DtB {
double doubleValue;
unsigned char asByte[sizeof(double)];
} ;


static void proj_bucket(double* y, double* x, 
const unsigned int length, const double a){
	union DtB * r1 = malloc(sizeof(union DtB) * length);
    union DtB *r1ptr = r1;
    union DtB * r2 = (union DtB *) x;
    union DtB * auxSwap;
    int     plength; 
    int illength;
    double p;
    union DtB  tau;
    int currentLength;
    int count = 0;
    int t[257];
    double s[257];
    double minS[257];
    double maxS[257];  
    union DtB * tmpswap;
    int * tmp;        
    tmp = &t[0];
    tmp ++;
    int bucketSize;
    int start;

    int i; 
    int j;
    int over = 0;

    illength = length;       
    int depth = 7;

    for (i = 0; i < 257; ++i){
        t[i] = 0;
        s[i] = 0.;
        minS[i] = DBL_MAX;
        maxS[i] = DBL_MIN;
    }

    for (i=0; i<length; i++){
        r1[i].doubleValue = y[i];
        tmp[r1[i].asByte[depth]]++;
        s[r1[i].asByte[depth]] += r1[i].doubleValue;
        minS[r1[i].asByte[depth]] = (minS[r1[i].asByte[depth]] < r1[i].doubleValue )? minS[r1[i].asByte[depth]] : r1[i].doubleValue;
        maxS[r1[i].asByte[depth]] = (maxS[r1[i].asByte[depth]] > r1[i].doubleValue )? maxS[r1[i].asByte[depth]] : r1[i].doubleValue;
    }

    tau.doubleValue = - a;
    illength = length;
    for (depth = 7; depth >= 0; depth --){
        
        for (i = 1; i < 256; ++i){                  // Count sort.
            tmp[i] = tmp[i] + tmp[i-1];
        }
        for (i = 0; i < illength; ++i){
            r2[t[r1[i].asByte[depth]]++] = r1[i];           
        }

        tmpswap = r2;
        r2 = r1;
        r1 = tmpswap;
        currentLength = illength;

        for (i = 255; i >= 0; --i){ // t[i] is the starting point of the i+1 values (because of the ++ )
            start = (i == 0)? 0 : t[i-1];
            bucketSize = currentLength - start; 
            currentLength -= bucketSize;
            if (bucketSize == 0){               
                continue;
            }
            if (tau.doubleValue / count > maxS[i]){ // Best possible remaining value is dominatied: end
                over = 1;
                break;
            }
            if ((tau.doubleValue + s[i]) / (count + bucketSize) < minS[i]){  // try keeping the min of b
                tau.doubleValue += s[i];
                count += bucketSize;
                continue;
            }
            r1 += start;
            r2 += start;
            illength = bucketSize;
            break;
        }  
        depth--;
        if (depth < 0 || over == 1 || i < 0)
        {
            break;
        }
        for (i = 0; i < 257; ++i){
            t[i] = 0;
            s[i] = 0.;
            minS[i] = DBL_MAX;
            maxS[i] = DBL_MIN;
        }
        for (i = 0; i <illength; ++i){   
            tmp[r1[i].asByte[depth]]++;
            s[r1[i].asByte[depth]] += r1[i].doubleValue;
            minS[r1[i].asByte[depth]] = (minS[r1[i].asByte[depth]] < r1[i].doubleValue )? minS[r1[i].asByte[depth]] : r1[i].doubleValue;
            maxS[r1[i].asByte[depth]] = (maxS[r1[i].asByte[depth]] > r1[i].doubleValue )? maxS[r1[i].asByte[depth]] : r1[i].doubleValue;
        }
        depth++;
    }
    tau.doubleValue /= count;
    // printf("%f\n",tau.doubleValue);

    for (i=0; i<length; i++)
        x[i]=(y[i]>tau.doubleValue ? y[i]-tau.doubleValue : 0.0); 
    free(r1ptr);
}

static void proj_bucket_filter(double* y, double* x, 
const unsigned int length, const double a){
	union DtB * r1 = malloc(sizeof(union DtB) * length);
    union DtB *r1ptr = r1;
    union DtB * r2 = (union DtB *) x;
    union DtB * auxSwap;
    int     plength; 
    int illength;
    double p;
    union DtB  tau;
    int currentLength;
    int count = 0;
    int t[257];
    double s[257];
    double minS[257];
    double maxS[257];  
    union DtB * tmpswap;
    int * tmp;        
    tmp = &t[0];
    tmp ++;
    int bucketSize;
    int start;

    int i; 
    int j;
    int over = 0;

    illength = length;       
    int depth = 7;

    for (i = 0; i < 257; ++i){
        t[i] = 0;
        s[i] = 0.;
        minS[i] = DBL_MAX;
        maxS[i] = DBL_MIN;
    }
    
    p = ((*r1).doubleValue=*y)-a;
    i = j = 1;
    start = 0;
    plength=-1; 
    tmp[r1[start].asByte[depth]]++;
    s[r1[start].asByte[depth]] += r1[start].doubleValue;
    minS[r1[start].asByte[depth]] = r1[start].doubleValue;
    maxS[r1[start].asByte[depth]] = r1[start].doubleValue;

    for (; j<length; i++, j++) 
        if (y[j]>p) {
            if ((p += ((r1[i].doubleValue=y[j])-p)/(i-plength)) <= y[j]-a) {
                p=y[j]-a;
                plength = i-1;
            }
            tmp[r1[i].asByte[depth]]++;
            s[r1[i].asByte[depth]] += r1[i].doubleValue;
            minS[r1[i].asByte[depth]] = (minS[r1[i].asByte[depth]] < r1[i].doubleValue )? minS[r1[i].asByte[depth]] : r1[i].doubleValue;
            maxS[r1[i].asByte[depth]] = (maxS[r1[i].asByte[depth]] > r1[i].doubleValue )? maxS[r1[i].asByte[depth]] : r1[i].doubleValue;
        }else{
            i--;
        } 

    tau.doubleValue = - a;
    illength = i;
    for (depth = 7; depth >= 0 && over == 0; depth --){
        
        for (i = 1; i < 256; ++i){                  // Count sort.
            tmp[i] = tmp[i] + tmp[i-1];
        }
        for (i = 0; i < illength; ++i){
            r2[t[r1[i].asByte[depth]]++] = r1[i];           
        }

        tmpswap = r2;
        r2 = r1;
        r1 = tmpswap;
        currentLength = illength;

        for (i = 255; i >= 0; --i){ // t[i] is the starting point of the i+1 values (because of the ++ )
            start = (i == 0)? 0 : t[i-1];
            bucketSize = currentLength - start; 
            currentLength -= bucketSize;
            if (bucketSize == 0){               
                continue;
            }
            if (tau.doubleValue / count > maxS[i]){ // Best possible remaining value is dominatied: end
                over = 1;
                break;
            }
            if ((tau.doubleValue + s[i]) / (count + bucketSize) < minS[i]){  // try keeping the min of b
                tau.doubleValue += s[i];
                count += bucketSize;
                continue;
            }
            r1 += start;
            r2 += start;
            illength = bucketSize;
            break;
        }  
        depth--;
        if (depth < 0 || over == 1|| i < 0)
        {
            break;
        }
        for (i = 0; i < 257; ++i){
            t[i] = 0;
            s[i] = 0.;
            minS[i] = DBL_MAX;
            maxS[i] = DBL_MIN;
        }
        start = illength-1;
        plength=-1; 
        p = (r1[start]).doubleValue - a;
        tmp[r1[start].asByte[depth]]++;
        s[r1[start].asByte[depth]] += r1[start].doubleValue;
        minS[r1[start].asByte[depth]] = r1[start].doubleValue;
        maxS[r1[start].asByte[depth]] = r1[start].doubleValue;
        for (i = illength-2; i >= 0; --i){
            if (r1[i].doubleValue > p && r1[i].doubleValue > tau.doubleValue / count) {
                if ((p += (r1[i].doubleValue - p)/((illength - i) - plength)) <= r1[i].doubleValue - a) {
                    p=r1[i].doubleValue - a;
                    plength = (illength - i) -1;
                }
                tmp[r1[i].asByte[depth]]++;
                s[r1[i].asByte[depth]] += r1[i].doubleValue;
                minS[r1[i].asByte[depth]] = (minS[r1[i].asByte[depth]] < r1[i].doubleValue )? minS[r1[i].asByte[depth]] : r1[i].doubleValue;
                maxS[r1[i].asByte[depth]] = (maxS[r1[i].asByte[depth]] > r1[i].doubleValue )? maxS[r1[i].asByte[depth]] : r1[i].doubleValue;
            }else{
                r1[i] = r1[--illength];
            } 
        }
        depth++;
    }
    tau.doubleValue /= count;
    // printf("%f\n",tau.doubleValue);

    for (i=0; i<length; i++)
        x[i]=(y[i]>tau.doubleValue ? y[i]-tau.doubleValue : 0.0); 
    free(r1ptr);
}


static void proj_pivot(double* y, double* x, 
const unsigned int length, const double a){
	double*	aux = (double*)malloc(length*sizeof(double));
	double*  aux0 = aux;
	int		auxlength=1; 
	int		auxlengthold=-1;	
	double	tau=(*aux=*y)-a;
	int 	i=1;
	for (; i<length; i++) 
		if (y[i]>tau) {
			if ((tau+=((aux[auxlength]=y[i])-tau)/(auxlength-auxlengthold))
			<=y[i]-a) {
				tau=y[i]-a;
				auxlengthold=auxlength-1;
			}
			auxlength++;
		} 
	if (auxlengthold>=0) {
		auxlength-=++auxlengthold;
		aux+=auxlengthold;
		while (--auxlengthold>=0) 
			if (aux0[auxlengthold]>tau) 
				tau+=((*(--aux)=aux0[auxlengthold])-tau)/(++auxlength);
	}
	do {
		auxlengthold=auxlength-1;
		for (i=auxlength=0; i<=auxlengthold; i++)
			if (aux[i]>tau) 
				aux[auxlength++]=aux[i];	
			else 
				tau+=(tau-aux[i])/(auxlengthold-i+auxlength);
	} while (auxlength<=auxlengthold);
// printf("%f\n",tau);

	for (i=0; i<length; i++)
		x[i]=(y[i]>tau ? y[i]-tau : 0.0); 
	if (x==y) free(aux0);
} 




