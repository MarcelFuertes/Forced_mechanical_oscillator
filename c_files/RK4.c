//code for RUNGA KUTTA 4-th ORDER
#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#define N 300

double g(double v);
double f(double x, double v, double t);
double solExacta(double t);
void rk_4(double h, double t[N], double x[N],double v[N], double tn);
void error(double h, double t[N], double x[N],double v[N], double tn);


int main(){
    FILE *fp;
    FILE *fp2;
    fp=fopen("RK4-points.txt","w");
    fp2=fopen("RK4Errors.txt","w");
    double t[N],x[N],v[N],tn,h;
    y[0]=1;
    t[0]=0;
    printf("Comencem la iteració desde t0=0, fins quin t_n vols arribar?:\n");
    scanf("%lf",&tn);
    printf("Introdueix el pas h:\n");
    scanf("%lf",&h);
    rk_4(h,t,x,v,tn);
    error(h,t,x,v,tn);
    int i;
    for(i=0;i<=(tn/h)+h;i++){
        fprintf(fp,"%.16G  \t %.16G \n",t[i],x[i]);
    }
    for(i=0;i<=(tn/h)+h;i++){
        fprintf(fp2,"%.16G  \t %.16G \n ",t[i], fabs(solExacta(t[i])-x[i]));
    }
    fclose(fp);
    fclose(fp2);
    return 0;
}


double g(double v){
    return v;
}

double f(double x, double v, double t){ //pendent de definir

}

double solExacta(double t){ //pendent de canviar
    return 1000/(1+999*exp(2*t));
}

void rk_4(double h, double t[N], double x[N], double v[N], double tn){
    int i=0;
    double k1,k2,k3,k4;
    double l1,l2,l3,l4;
    l1=h*f(x[0],v[0],t[0]);
    k1=h*g(v[0]);
    l2=h*f(x[0]+k1/2.,v[0]+l1/2.,t[0]+h/2.);
    k2=h*g(v[0]+l1/2.);
    l3=h*f(x[0]+k2/2.,v[0]+l2/2.,t[0]+h/2.);
    k3=h*g(v[0]+l2/2.);
    l4=h*f(x[0]+k3,v[0]+l3,t[0]+h);
    k4=h*g(v[0]+k3);
    while(t[i]<=tn+h){
        i++;
        x[i]=x[i-1]+(k1+2*k2+2*k3+k4)/6.;
        v[i]=v[i-1]+(l1+2*l2+2*l3+l4)/6.;
        t[i]=t[i-1]+h;
        l1=h*f(x[i],v[i],t[i]);
        k1=h*g(v[i]);
        l2=h*f(x[i]+k1/2.,v[i]+l1/2.,t[i]+h/2.);
        k2=h*g(v[i]+l1/2.);
        l3=h*f(x[i]+k2/2.,v[i]+l2/2.,t[i]+h/2.);
        k3=h*g(v[i]+l2/2.);
        l4=h*f(x[i]+k3,v[i]+l3,t[i]+h);
        k4=h*g(v[i]+k3);
    }
}

void error(double h, double t[N], double x[N], double v[N], double tn){
    int i=0;
    double exacta[N];
    printf("Iteration \t t_i \t\t SolExacta \t y_i \t\t error \t \n ");
    while(t[i]<=tn+h/10.){
    exacta[i]=solExacta(t[i]);
    printf("%i \t\t %lf \t %.16G \t %.16G \t %.16G \t \n", i,t[i], exacta[i], x[i], fabs(exacta[i]-x[i]));
    i++;
    t[i]=t[i-1]+h;
    }
    
}







