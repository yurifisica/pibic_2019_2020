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
    u[0]=I*csqrt(-Z*Y[0]);
    Zint[0]=u[0]/Y[0];
    for (cont=1;cont<=n;cont++)
    {
        Y[cont]=sigma[cont-1]+I*w*eps;
        u[cont]=I*csqrt(-Z*Y[cont]);
        Zint[cont]=u[cont]/Y[cont];
    }
    complex double Zap[n];
    Zap[n-1]=Zint[n];
    for (cont=n-2;cont>=0;cont--)
    {
        complex double aux=u[cont+1]*h[cont];
        Zap[cont]=Zint[cont+1]*(Zap[cont+1]+Zint[cont+1]*tgh(aux))/(Zint[cont+1]+Zap[cont+1]*tgh(aux));
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
    int cmd=0;
    for (cont=0;cont<np;cont++)
    {
        int flag=0;
        for (ic=0;ic<n;ic++)
        {
            if (zp[cont]<z[ic] && flag==0)
            {
                cmd=ic;
                flag=1;
            }
        }
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
	int n_cam=2,cont,j,i;//numero de camadas
    double frq=5*pow(10,3);//frequência
    double espessura[1]={200};//espessura das camadas intermediárias
    size_t n_h=length(espessura);
    double sigma[2]={pow(10,-3),0.5};
    size_t n_s=length(sigma);
    double z0=0;
    double zn=espessura[0]*2;//profundidade da malha
    int nz=37;//numero de pontos em z
    double dz1=abs(zn-z0)/(nz-1);
    double z[nz];
    for (cont=0;cont<nz;cont++)
    {
        z[cont]=cont*dz1;
    }
    nz=length(z);
    double x0=0;
    double xn=10000;//largura da malha
    int nx=37;//numero de pontos em x, na malha reduzida
    double x[nx+4],xt[nx];
    x[0]=-xn+xn/4;x[1]=-xn+xn/2;
    double dx1=abs(xn-x0)/(nx-1);
    for (cont=0;cont<nx;cont++)
    {
        xt[cont]=cont*dx1;
        x[cont+2]=cont*dx1;
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
    complex double Z=I*w*mu;
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
    double xv2[nz*nx2];
    for (j=0;j<nz;j++)
    {
      for (i=0;i<nx2;i++)
      {
        xv2[ic]=x[i];
        ic=ic+1;
      }
    }
    //definição de coeficientes da edo
    complex double a(double x, double z){return 1;}
    complex double b(double x, double z){return 0;}
    complex double c(double x, double z){return 1;}
    complex double d(double x, double z){return 0;}
    complex double e(double x, double z){return 0;}
    complex double f(double Z, double Y){return -Z*Y;}
    complex double g(double Z, complex double Yj, complex double Yp, complex double Hp){return -Z*(Yj-Yp)*Hp;}
    ic=0;
    int cc=0;
    int el[(nz-1)*(nx-1)*2][4];
    cont=1;
    for (j=0;j<(nz-1);j++)
    {
      for (i=0;i<nx;i++)
      {
        if(i<nx)
        {
          el[cc][0]=ic+1;
          el[cc][1]=ic+2;
          el[cc][2]=ic+nx+1;
          el[cc][3]=cc;
          cc=cc+1;
          el[cc][0]=ic+2;
          el[cc][1]=ic+nx+2;
          el[cc][2]=ic+nx+1;
          el[cc][3]=cc;
        } else
        {
          cont++;
        }
        ic++;
      }
    }
    //condições de fronteira
    int vc[np],qt_front;
    qt_front=0;
    complex double Hs[np];
    for (i=0;i<np;i++)
    {
      vc[i]=0;
      if (xv[i]==xv[1] || xv[i]==xv[np] || zv[i]==zv[1]|| zv[i]==zv[np])
      {
        Hs[i]=0; vc[i]=1;
        qt_front++;
      }
    }

    //definição dos meios secundários
    complex double M[np];
    //meio secundário 1 (solo)
    for (i=0;i<np;i++)
    {
      if (xv[i]>=3000 && xv[i]>=10000 && zv[i]>5000 && zv[i]<=6500)
      {
        Yj[i]=sigma1+I*w*eps; sig[i]=sigma1;
      }
      M[i]=0;
    }

    //meio secundário 2 (mar)
    for (i=0;i<np;i++)
    {
      if (xv[i]>=0 && xv[i]<=3000)
      {
        if ((zv[i]>=espessura[0]) && (zv[i]<(5*espessura[0])/4))
        {
          Yj[i]=sigma2+I*w*eps; sig[i]=sigma2;
        }
      }
    }
    ic=0;
    int flag=0;
    for (j=1;j<(nz-1);j++)
    {
      for (i=1;i<(nx-1);i++)
      {
        ic++;
        if (xv[i]>=2000 && x[i]<=5000)
        {
          if (z[j]>=(((z[j]-z[j-1])/(x[i]-x[i-1]))*x[i]) && z[j]<=(3*espessura[0])/2 && z[j]>(13*espessura[0])/12)
          {
            Yj[i]=sigma2+I*w*eps; sig[i]=sigma2;
          }
        }
      }
    }
    //meio secundário 3 (prédios)
    for (i=0;i>np;i++)
    {
        if (xv[i]>=6500 && xv[i]<=7000 && zv[i]>=espessura[0]-50 && zv[i]<espessura[0]) {
            Yj[i]=sigma3+I*w*eps;
            sig[i]=sigma3;
        }
        if (xv[i]>=7500 && xv[i]<=8000 && zv[i]>=espessura[0]-50 && zv[i]<espessura[0]) {
            Yj[i]=sigma3+I*w*eps;
            sig[i]=sigma3;
        }
        if (xv[i]>=8500 && xv[i]<=9500 && zv[i]>=espessura[0]-50 && zv[i]<espessura[0]) {
            Yj[i]=sigma3+I*w*eps;
            sig[i]=sigma3;
        }
    }

    //construção da matriz e vetor da direita
    complex double K[np][np];
    complex double Kr[np-qt_front][np-qt_front];
    complex double k[3][3],m[3];
    int np2=(nx-2)*(nz-2);
    complex double Mr[np-qt_front];
    double dt;
    double dx[3],dz[3];
    int i1,j1;
    int jc=0;
    complex double al[3];
    complex double bl[3];
    complex double cl[3];
    complex double dl[3];
    complex double El[3];
    complex double fl[3];
    complex double gl[3];
    for (i=0;i<ne;i++){
        dx[0]=zv[el[i][1]]-zv[el[i][2]];
        dx[1]=zv[el[i][2]]-zv[el[i][0]];
        dx[2]=zv[el[i][0]]-zv[el[i][1]];
        dz[0]=xv[el[i][2]]-xv[el[i][1]];
        dz[1]=xv[el[i][0]]-xv[el[i][2]];
        dz[2]=xv[el[i][1]]-xv[el[i][0]];
        dt=abs(xv[el[i][1]]*zv[el[i][2]]+xv[el[i][0]]*zv[el[i][1]]+xv[el[i][2]]*zv[el[i][0]]-xv[el[i][1]]*zv[el[i][0]]-xv[el[i][2]]*zv[el[i][1]]-xv[el[i][0]]*zv[el[i][2]])/2;//delta
        for (cont=0;cont<3;cont++){
          al[cont]=a(xv[el[i][cont]],zv[el[i][cont]]);
          bl[cont]=b(xv[el[i][cont]],zv[el[i][cont]]);
          cl[cont]=c(xv[el[i][cont]],zv[el[i][cont]]);
          dl[cont]=d(xv[el[i][cont]],zv[el[i][cont]]);
          El[cont]=e(xv[el[i][cont]],zv[el[i][cont]]);
          fl[cont]=f(Z,Yj[el[i][cont]]); gl[cont]=g(Z,Yj[el[i][cont]],Y[el[i][cont]],Hp[el[i][cont]]);
        }
        //matriz de cada elemento
        for (i=0;i<3;i++)
        {
          for(j=0;j<3;j++)
          {
            k[i][j]+=(-1/(4*dt))*(dx[1]*dx[1]*al[j]+dx[1]*dz[1]*bl[j]+dz[1]*dz[1]*cl[j])+(1/6)*(dx[1]*dl[j]+dz[1]*El[j])+(dt/12)*2*fl[j];
          }
        }
        //vetor da direita em cada elemento
        m[0]=(dt/12)*(2*gl[0]+gl[1]+gl[2]);
        m[1]=(dt/12)*(gl[0]+2*gl[1]+gl[2]);
        m[2]=(dt/12)*(gl[0]+gl[1]+2*gl[2]);

        //Matriz global vetor da direita
        for (j1=0;j1<3;j1++){
            for (i1=0;i1<3;i1++){
                K[el[i][j1]][el[i][i1]]=K[el[i][j1]][el[i][i1]]+k[j1][i1];
            }
            M[el[i][j1]]=M[el[i][j1]]+m[j1];
        }
    }

    //redução do sistema - Inclusão das condições de fronteira
    ic=0;
    jc=0;
    for (i=0;i<np;i++)
    {
      if (vc[i]!=1)
      {
        Mr[ic]=M[i];
        for (j=0;j<np;j++)
        {
          if (vc[i]!=1)
          {
            Kr[ic][jc]=K[i][j];
            jc++;
          } else
          {
            Mr[ic]-=K[i][j]*Hs[j];
          }
        }
        ic++;
      }
    }


    //solução por eliminação gaussiana
    //1-triangularização
    complex double Ur[np-qt_front];
    complex double temp;
    for (i=0;i<(np-qt_front-1);i++)
    {
      for (j=i+1;j<np-qt_front;j++)
      {
        temp=Kr[j][i]/Kr[i][i];
        for (int l=i+1;l<np-qt_front;l++)
        {
          Kr[j][l]-=temp*Kr[i][l];
        }
        Mr[j]-=temp*Mr[i];
      }
    }
    //2-retrossubstituição
    Ur[np-qt_front]=Mr[np-qt_front]/Kr[np-qt_front][np-qt_front];
    for (i=(np-qt_front-1);i>=0;i--)
    {
      Ur[i]=Mr[i];
      for (j=i+1;j<np-qt_front;j++)
      {
        Ur[i]-=Kr[i][j]*Ur[j];
      }
      Ur[i]=Ur[i]/Kr[i][i];
    }
    //Incorporando as condições de fronteira
    complex double U[np];
    ic=0;
    for (i=0;i<np;i++)
    {
      if (vc[i]!=1)
      {
        U[i]=Ur[ic];
        ic++;
      }
      else
      {
        U[i]=Hs[i];
      }
    }
    complex double Htotal[np];
    FILE *file;
    file=fopen("H_sec.dat","w");
    for (i=0;i<np;i++)
    {
      Htotal[i]=U[i]+Hp[i];
      fprintf(file,"%i %f %f\n",cont,creal(Htotal[i]),cimag(Htotal[i]));
    }

    return 0;
}
