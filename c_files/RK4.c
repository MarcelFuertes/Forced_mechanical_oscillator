//aquesta és la resolució del primer problema utilitzant RUNGA KUTTA ORDRE 4
#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#define N 300

double f(double y);
double solExacta(double x);
void rk_4(double h, double x[N], double y[N], double xn);
void error(double h, double x[N], double y[N], double xn);


int main(){
    FILE *fp;
    FILE *fp2;
    fp=fopen("PuntsRK4.txt","w");
    fp2=fopen("RK4Errors.txt","w");
    double x[N],y[N],xn,h;
    y[0]=1;
    x[0]=0;
    printf("Comencem la iteració desde x0=0, fins quin x_n vols arribar?:\n");
    scanf("%lf",&xn);
    printf("Introdueix el pas h:\n");
    scanf("%lf",&h);
    rk_4(h,x,y,xn);
    error(h,x,y,xn);
    int i;
    for(i=0;i<=(xn/h)+h;i++){
        fprintf(fp,"%.16G \n",y[i]);
    }
    for(i=0;i<=(xn/h)+h;i++){
        fprintf(fp2,"%.16G  \t %.16G \n ",x[i], fabs(solExacta(x[i])-y[i]));
    }
    fclose(fp);
    fclose(fp2);
    return 0;
}


double f(double y){
    return 0.002*y*(y-1000);
}
double solExacta(double x){
    return 1000/(1+999*exp(2*x));
}

void rk_4(double h, double x[N], double y[N], double xn){
    int i=0;
    double k1,k2,k3,k4;
    k1=h*f(y[0]);
    k2=h*f(y[0]+k1/2.);
    k3=h*f(y[0]+k2/2.);
    k4=h*f(y[0]+k3);
    while(x[i]<=xn+h){
        i++;
        y[i]=y[i-1]+(k1+2*k2+2*k3+k4)/6.;
        k1=h*f(y[i]);
        k2=h*f(y[i]+k1/2.);
        k3=h*f(y[i]+k2/2.);
        k4=h*f(y[i]+k3);
        x[i]=x[i-1]+h;
    }
}

void error(double h, double x[N], double y[N], double xn){
    int i=0;
    double exacta[N];
    printf("Iteració \t x_i \t\t SolExacta \t y_i \t\t error \t \n ");
    while(x[i]<=xn+h/10.){        
    exacta[i]=solExacta(x[i]);
    printf("%i \t\t %lf \t %.16G \t %.16G \t %.16G \t \n", i,x[i], exacta[i], y[i], fabs(exacta[i]-y[i]));
    i++;
    x[i]=x[i-1]+h;
    }
    
}







