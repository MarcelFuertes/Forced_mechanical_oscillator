//code for RUNGA KUTTA 4-th ORDER
#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#define PI 3.141592653589793
//#define N 300
#define B 2.
#define Wo 10.
#define M 0.1
#define Fo 1.
#define cte 0.2

double g(double v);
double f(double x, double v, double t);
double solExacta(double t);
void rk_4(double h, double t[N], double x[N],double v[N], double tn);
void error(double h, double t[N], double x[N],double v[N], double tn);


int main(){
    FILE *fp;
    FILE *fp2;
    fp=fopen("RK4-points.txt","w"); //we create the txt file that contains the points of the RK4 algorithm
    fp2=fopen("RK4Errors.txt","w"); //we create the txt file that contains the exact errors.
    printf("Introduce the value for h:\n");
    scanf("%lf",&h);
	double tn = 10*PI/Wo;
	int N;
	double value=tn/h;
	N = (int)value+1;
    double t[N],x[N],v[N],tn,h;
    x[0]=0.;
	v[0]=1.;
    t[0]=0.;

    rk_4(h,t,x,v,tn);
    error(h,t,x,v,tn);
    int i;
    for(i=0;i<=(tn/h)+h;i++){
        fprintf(fp,"%.16G  \t %.16G \n",t[i],x[i]);
    }
    for(i=0;i<=(tn/h)+h;i++){ //this should be in the error function
        fprintf(fp2,"%.16G  \t %.16G \n ",t[i], fabs(solExacta(t[i])-x[i]));
    }
    fclose(fp);
    fclose(fp2);
    return 0;
}


double g(double v){
    return v;
}

double f(double x, double v, double t){
    double w=cte*Wo;
    return -(B/M)*v-Wo*Wo*x+(Fo/M)*cos(w*t);
}

double solExacta(double t){
	double lambda1,lambda2,alpha,theta;
	double w=cte*Wo;
	lambda1=(B+sqrt(B*B-4.*M*M*Wo*Wo))/(2.*M);
	lambda2=(B-sqrt(B*B-4.*M*M*Wo*Wo))/(2.*M);
	alpha=Fo/sqrt(M*M*(Wo*Wo-w*w)*(Wo*Wo-w*w)+B*B*w*w);
	theta=atan((B*w)/(M*(Wo*Wo-w*w)));
		
	double c1,c2;
	c1=(alpha*lambda2*cos(theta)+alpha*w*sin(theta)-1.)/(lambda1-lambda2);
	c2=(alpha*lambda1*cos(theta)+alpha*w*sin(theta)-1.)/(lambda2-lambda1);

    return c1*exp(-t*lambda1)+c2*exp(-t*lambda2)+alpha*cos(w*t-theta);
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
        t[i]=t[0]+i*h;
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

void error(double h, double t[N], double x[N], double v[N], double tn){//this should save the points in a txt file
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







