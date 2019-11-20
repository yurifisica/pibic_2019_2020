#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#define length(a) sizeof(a)/sizeof(a[0])
const float pi = 3.1415926;
const float mu=4*3.1415926E-7;
const float eps= 8.85E-12;

double tgh(double x)
{
    if (x<=86)
    {
        return (exp(x)-exp(-x))/(exp(x)+exp(-x));
    } else
    {
        return 1;
    }
}

complex double mult_complex(complex double Z1, complex double Z2)
{
complex double Z = (creal(Z1)*creal(Z2)-cimag(Z1)*cimag(Z2))+(creal(Z1)*cimag(Z2)+cimag(Z1)*creal(Z2))*I;
return Z;
}

complex double div_complex(complex double Z1, complex double Z2)
{
    complex double Z = (creal(Z1)*creal(Z2)-cimag(Z1)*cimag(Z2)+(-creal(Z1)*cimag(Z2)+cimag(Z1)*creal(Z2))*I)/(pow(creal(Z2),2)-pow(cimag(Z2),2));
    return Z;
}

//cálculo do campo primário
void campo_prim(int n,double sigma[],int n_h, double h[], double frq, double zp[], int np,complex double Hp[])
{
    FILE *file;
    file=fopen("a.dat","w");
    int cont,ic;
    if (file==NULL)
    {
        printf("caminho inválido");
    }
    double w=2*pi*frq;
    complex double Z=I*w*mu;
    complex double H[n];
    H[0]=1;
    double z[n_h+1];
    z[0]=0;
    for (cont=0;cont<n_h;cont++)
    {
        z[cont+1]=z[cont]+h[cont];
    }
    complex double Y[n+1];
    complex double u[n+1];
    complex double Zint[n+1];
    Y[0]=I*w*eps;
    u[0]=I*sqrt(-Z*Y[0]);
    Zint[0]=u[0]/Y[0];
    printf("%4.12f %4.12f",creal(u[0]),cimag(u[0]));
    for (cont=1;cont<=n;cont++)
    {
        Y[cont]=sigma[cont-1]+I*w*eps;
        u[cont]=I*sqrt(-Z*Y[cont]);
        Zint[cont]=u[cont]/Y[cont];
    }
    complex double Zap[n];
    Zap[n-1]=Zint[n];
    for (cont=n-2;cont>=0;cont--)
    {
        complex double aux=u[cont+1]*h[cont];
        Zap[cont]=Zint[cont+1]*(Zap[cont+1]+Zint[cont+1]*tgh(aux))/(Zint[cont+1]+Zap[cont+1]*tgh(aux));
//                 printf("%f %f\n", creal(Zap[cont]),cimag(Zap[cont]));
    }
    complex double rtm[n-1];
    rtm[0]=(Zint[0]-Zap[0])/(Zint[0]+Zap[0]);
    for (cont=1;cont<n;cont++)
    {
        complex double aux=u[cont]*h[cont-1];
        rtm[cont]=(Zint[cont]-Zap[cont])/(Zint[cont]+Zap[cont]);
        H[cont]=H[cont-1]*(1+rtm[cont-1])*exp(-aux)/(1+rtm[cont]*exp(-2*aux));
    }
    H[n]=H[n-1]*(1+rtm[n-1]);
//    printf("%f %f", creal(H[n]),cimag(H[n]));
    int aa=(np+1)*sizeof(complex double);
    //complex double *Hp = malloc((np+1)*sizeof(complex double));
    int cmd=0;
    for (cont=0;cont<np;cont++)
    {
        int flag=0;
        for (ic=0;ic<n;ic++)
        {
            //printf("%f\n",z[ic]);
            if (zp[cont]<z[ic] && flag==0)
            {
                cmd=ic;
                flag=1;
            }
        }
       //printf("%f\n",zp[0]);
        if (zp[cont]>z[0])
        {
            if (zp[cont]<z[n_h])
            {
                complex double aux=u[cmd]*(zp[cont]-z[cmd]);
                Hp[cont]=H[cmd]*(exp(-aux)+rtm[cmd]*exp(aux));
            }
            else
            {
                complex double aux=u[n]*(zp[cont]-z[cmd]);
                Hp[cont]=H[n]*exp(-aux);
            }
        }
        else
        {
            complex double aux=u[0]*zp[cont];
            Hp[cont]=H[0]*(exp(-aux)+rtm[0]*exp(aux));
        }
       fprintf(file,"%i %f %f\n",cont,creal(Hp[cont]),cimag(Hp[cont]));
    }

    fclose(file);
    return;
}

//programa principal
int main()
{
    printf("%1.14f",eps);
	int n_cam=2,cont,j,i;//numero de camadas
    double frq=5*pow(10,3);//frequência
    double espessura[1]={200};//espessura das camadas intermediárias
    size_t n_h=length(espessura);
    double sigma[2]={pow(10,-3),0.5};
    size_t n_s=length(sigma);
    double z0=0;
    double zn=espessura[0]*2;//profundidade da malha
    int nz=37;//numero de pontos em z
    double dz=abs(zn-z0)/(nz-1);
    double z[nz];
    for (cont=0;cont<nz;cont++)
    {
        z[cont]=cont*dz;
    }
    nz=length(z);
    double x0=0;
    double xn=10000;//largura da malha
    int nx=37;//numero de pontos em x, na malha reduzida
    double x[nx+4],xt[nx];
    x[0]=-xn+xn/4;x[1]=-xn+xn/2;
    double dx=abs(xn-x0)/(nx-1);
    for (cont=0;cont<nx;cont++)
    {
        xt[cont]=cont*dx;
        x[cont+2]=cont*dx;
    }
    x[nx+3]=xn+xn/4;x[nx+4]=xn+xn/2;
    nx=length(x);//numero de pontos em x, na malha completa
    int nx2=length(xt);
    int np=nx*nz;
    int ne=(nx-1)*(nz-1)*2;
    double xv[np],zv[np];
    complex double Y[np],Yj[np],sig[np];
    double w =2*pi*frq;
    Y[0]=I*w*eps;
    for (cont=1;cont<=n_cam;cont++)
    {
        Y[cont]=sigma[cont-1]+I*w*eps;
    }
    int ic=0;
    for (j=0;j<nz;j++)
    {
        for(i=0;i<nx;i++)
        {
            if (z[j]>=espessura[0])
            {
                sig[ic]=sigma[1];
                Yj[ic]=sigma[1]+I*w*eps;
                Y[ic]=sigma[1]+I*w*eps;
            }
            else
            {
                sig[ic]=sigma[0];
                Yj[ic]=sigma[0]+I*w*eps;
                Y[ic]=sigma[0]+I*w*eps;
            }
            xv[ic]=x[i];zv[ic]=z[j];
            ic++;
        }
    }
    double sigma1=0.5;//solo embaixo dos prédios
    double sigma2=4.8;//água salgada
    double sigma3=2.8;//proporção de 1/2 de areia e cimento (prédios), com fibras de aço carbono
    ic=0;
    complex double Hp[np];
    campo_prim(n_s,sigma,n_h,espessura,frq,zv,np,Hp);
    //free(Hp);
    return 0;
}



