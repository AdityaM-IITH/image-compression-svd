#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "svd.h"

double** crmat(int m,int n){double **r=(double**)calloc(m,sizeof(double*));for(int i=0;i<m;i++)r[i]=(double*)calloc(n,sizeof(double));return r;}
double *crvec(int n){double *r=(double*)malloc(n*sizeof(double));return r;}
double vnorm(double *v,int n){double s=0.0;for(int i=0;i<n;i++)s+=v[i]*v[i];return sqrt(s);}
double dot(double *a,double *b,int n){double s=0.0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
void Axv(double **a,double *x,double *y,int m,int n){for(int i=0;i<m;i++)y[i]=dot(a[i],x,n);}
void scl(double *v,double s,int n){for(int i=0;i<n;i++)v[i]*=s;}
void ypx(double *x,double *y,double a,int n){for(int i=0;i<n;i++)y[i]+=a*x[i];}
void cpyv(double *x,double *y,int n){for(int i=0;i<n;i++)y[i]=x[i];}

double **trp(double **A,int m,int n){
    double **r=(double**)malloc(n*sizeof(double*));
    for(int i=0;i<n;i++){r[i]=(double*)malloc(m*sizeof(double));
        for(int j=0;j<m;j++)r[i][j]=A[j][i];
    }
    return r;
}

int loop(double **A,double *Q,double *v,int k,int m,int n,double *a,double *b){
    double *v1=crvec(m);
    double **AT=trp(A,m,n);
    int i;
    for(i=0;i<k;i++){
        for(int j=0;j<n;j++)v[j]=0.0;
        for(int j=0;j<m;j++)v1[j]=0.0;
        Axv(A,&Q[i*n],v1,m,n);
        Axv(AT,v1,v,n,m);
        if(i>0)ypx(&Q[(i-1)*n],v,-b[i-1],n);
        a[i]=dot(&Q[i*n],v,n);
        ypx(&Q[i*n],v,-a[i],n);
        for (int s = 0; s <= i; s++) {
            double c = dot(&Q[s * n], v, n);
            ypx(&Q[s * n], v, -c, n);
        }
        if(i<k-1){
            b[i]=vnorm(v,n);
            if(b[i]<1e-12){b[i]=0.0; i++; break;}
            scl(v,1.0/b[i],n);
        }
        cpyv(v,&Q[(i+1)*n],n);
    }
    free(v1);
    for(int j=0;j<n;j++)free(AT[j]);
    free(AT);
    return i;
}

int GRot(double *a,double *b,int k,double *ev,int mit,double lim){
    if(k<=0)return 0;
    b[k-1]=0.0;
    for(int j=0;j<k;++j){
        int it=0;
        while(1){
            int p=j;
            while(p<k-1){
                double tol=lim*(fabs(a[p])+fabs(a[p+1]));
                if(fabs(b[p])<=tol){b[p]=0.0;break;}
                ++p;
            }
            if(p==j)break;
            if(++it>mit){printf("max");return -1;}
            double al=a[p-1],ac=a[p],bl=b[p-1];
            double d=(al-ac)/2.0;
            double den=fabs(d)+sqrt(d*d+bl*bl);
            double sgn=d>=0?1.0:-1.0;
            double u=ac - sgn*(bl*bl)/(den==0.0?1e-300:den);
            double x=a[j]-u;
            double z=b[j];
            for(int i=j;i<=p-1;++i){
                double r=hypot(x,z);
                double c=1.0,s=0.0;
                if(r!=0){c=x/r;s=z/r;}
                if(i>j)b[i-1]=r;
                double ai=a[i],ai1=a[i+1],bi=b[i];
                double t1=c*c*ai + 2.0*c*s*bi + s*s*ai1;
                double t2=s*s*ai - 2.0*c*s*bi + c*c*ai1;
                a[i]=t1; a[i+1]=t2;
                double bin=(c*c - s*s)*bi + c*s*(ai1 - ai);
                b[i]=bin;
                if(i<p-1){
                    double bi1=b[i+1];
                    double xn=b[i];
                    double zn=-s*bi1;
                    b[i+1]=c*bi1;
                    x=xn; z=zn;
                }else{x=b[i]; z=0.0;}
                if(ev){
                    for(int rix=0;rix<k;++rix){
                        double e0=ev[rix*k+i],e1=ev[rix*k+i+1];
                        ev[rix*k+i]=c*e0 - s*e1;
                        ev[rix*k+i+1]=s*e0 + c*e1;
                    }
                }
            }
            b[p-1]=x;
        }
    }
    for(int i=0;i<k-1;++i){
        int l=i; double dp=a[i];
        for(int j=i+1;j<k;++j)if(a[j]>dp){l=j;dp=a[j];}
        if(l!=i){
            double t=a[i]; a[i]=a[l]; a[l]=t;
            if(ev){
                for(int r=0;r<k;++r){
                    double te=ev[r*k+i]; ev[r*k+i]=ev[r*k+l]; ev[r*k+l]=te;
                }
            }
        }
    }
    return 0;
}

double* channel(double *A1,int m,int n,int k_in,int mxg){
    int k=k_in;
    double **A=crmat(m,n);
    for(int i=0;i<m;i++)for(int j=0;j<n;j++)A[i][j]=A1[i*n+j];
    double *Q=crvec(n*(k+1)); for(int i=0;i<n*(k+1);i++)Q[i]=0.0;
    unsigned u=1234567u;
    for(int i=0;i<n;i++){u=1664525u*u+1013904223u; double x=((int)(u>>8)&65535)-32768.0; Q[i]=x;}
    double nq=vnorm(Q,n); if(nq>1e-12) scl(Q,1.0/nq,n);
    double *a=crvec(k+1),*b=crvec(k+1); for(int i=0;i<k+1;i++){a[i]=0.0;b[i]=0.0;}
    double *v=crvec(n);
    double *ev=crvec(k*k); for(int i=0;i<k;i++)for(int j=0;j<k;j++)ev[i*k+j]=(i==j);
    int rnk=loop(A,Q,v,k,m,n,a,b);
    if(rnk<k)k=rnk;
    GRot(a,b,k,ev,10000,1e-12);
    double *S=crvec(k);
    for(int i=0;i<k;i++){double sv=a[i]; if(sv<0 && sv>-1e-10) sv=0.0; S[i]=sqrt(sv>0?sv:0.0);}
    double *V=crvec((size_t)n*k);
    for(int j=0;j<k;j++){double *vj=&V[j*n]; for(int t=0;t<n;t++)vj[t]=0.0; for(int t=0;t<k;t++){double c=ev[t*k+j]; ypx(&Q[t*n],vj,c,n);}}
    double *U=crvec((size_t)m*k),*tmp=crvec(m);
    for(int j=0;j<k;j++){
        double *vj=&V[j*n],*uj=&U[j*m];
        for(int i=0;i<m;i++)tmp[i]=0.0;
        Axv(A,vj,tmp,m,n);
        double sj_actual=vnorm(tmp,m);
        S[j]=sj_actual;
        if(sj_actual>1e-12){for(int i=0;i<m;i++)uj[i]=tmp[i]/sj_actual;}
        else{for(int i=0;i<m;i++)uj[i]=0.0;}
    }
    double *Ah=crvec((size_t)m*n); for(size_t t=0;t<(size_t)m*n;t++)Ah[t]=0.0;
    for(int j=0;j<k;j++){double *vj=&V[j*n],*uj=&U[j*m]; double sj=S[j];
        for(int i=0;i<m;i++){double uis=uj[i]*sj; double *row=&Ah[(size_t)i*n];
            for(int c=0;c<n;c++)row[c]+=uis*vj[c];}}
    double mn=Ah[0],mxA=Ah[0];
    for(size_t t=1;t<(size_t)m*n;t++){if(Ah[t]<mn)mn=Ah[t]; if(Ah[t]>mxA)mxA=Ah[t];}
    double rng=mxA-mn; int outmx=mxg>0?mxg:255;
    if(rng>1e-12){double sc=outmx/rng; for(size_t t=0;t<(size_t)m*n;t++){double y=(Ah[t]-mn)*sc; if(y<0)y=0; if(y>outmx)y=outmx; Ah[t]=y;}}
    else{for(size_t t=0;t<(size_t)m*n;t++)Ah[t]=0.0;}
    for(int i=0;i<m;i++)free(A[i]); free(A);
    free(Q); free(a); free(b); free(v); free(ev); free(S); free(V); free(U); free(tmp);
    return Ah;
}
