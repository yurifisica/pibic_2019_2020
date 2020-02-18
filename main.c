#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#define length(a) sizeof(a)/sizeof(a[0])
const float pi = 3.1415926;
const float mu=4*3.1415926E-7;
const float eps= 8.85E-12;

// double complex ctanh(double complex x)
// {
//     if (creal(x)<=86)
//     {
//         return (cexp(x)-cexp(-x))/(cexp(x)+cexp(-x));
//     } else
//     {
//         return 1;
//     }
// }

double complex mult_complex(double complex Z1, double complex Z2)
{
  double complex Z = (creal(Z1)*creal(Z2)-cimag(Z1)*cimag(Z2))+(creal(Z1)*cimag(Z2)+cimag(Z1)*creal(Z2))*I;
  return Z;
}

double complex div_complex(double complex Z1, double complex Z2)
{
  double complex Z = (creal(Z1)*creal(Z2)-cimag(Z1)*cimag(Z2)+(-creal(Z1)*cimag(Z2)+cimag(Z1)*creal(Z2))*I)/(pow(creal(Z2),2)-pow(cimag(Z2),2));
  return Z;
}

//calculo do campo primario
void campo_prim(int n,double sigma[],int n_h, double h[], double frq, double zp[], int np,double complex Hp[])
{
  FILE *file;
  file=fopen("a.dat","w");
  int cont,ic;
  if (file==NULL)
  {
    printf("caminho invalido");
  }
  double w=2*pi*frq;
  double complex Z=I*w*mu;
  double complex H[n];
  H[0]=1;
  double z[n_h+1];
  z[0]=0;
  for (cont=0;cont<n_h;cont++)
  {
    z[cont+1]=z[cont]+h[cont];
  }
  double complex Y[n+1];
  double complex u[n+1];
  double complex Zint[n+1];
  Y[0]=I*w*eps;
  u[0]=I*csqrt(-Z*Y[0]);
  Zint[0]=u[0]/Y[0];
  for (cont=1;cont<=n;cont++)
  {
    Y[cont]=sigma[cont-1]+I*w*eps;
    u[cont]=I*csqrt(-Z*Y[cont]);
    Zint[cont]=u[cont]/Y[cont];
  }
  double complex Zap[n];
  Zap[n-1]=Zint[n];
  for (cont=n-2;cont>=0;cont--)
  {
    double complex aux=u[cont+1]*h[cont];
    Zap[cont]=Zint[cont+1]*(Zap[cont+1]+Zint[cont+1]*ctanh(aux))/(Zint[cont+1]+Zap[cont+1]*ctanh(aux));
  }
  double complex rtm[n-1];
  rtm[0]=(Zint[0]-Zap[0])/(Zint[0]+Zap[0]);
  for (cont=1;cont<n;cont++)
  {
    double complex aux=u[cont]*h[cont-1];
    rtm[cont]=(Zint[cont]-Zap[cont])/(Zint[cont]+Zap[cont]);
    H[cont]=H[cont-1]*(1+rtm[cont-1])*cexp(-aux)/(1+rtm[cont]*cexp(-2*aux));
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
        double complex aux=u[cmd]*(zp[cont]-z[cmd]);
        Hp[cont]=H[cmd]*(cexp(-aux)+rtm[cmd]*cexp(aux));
      }
      else
      {
        double complex aux=u[n]*(zp[cont]-z[cmd]);
        Hp[cont]=H[n]*cexp(-aux);
      }
    }
    else
    {
      double complex aux=u[0]*zp[cont];
      Hp[cont]=H[0]*(cexp(-aux)+rtm[0]*cexp(aux));
    }
    fprintf(file,"%i %f %f\n",cont,creal(Hp[cont]),cimag(Hp[cont]));
  }
  fclose(file);
  return;
}
double complex a(double x, double z){return 1;}
double complex b(double x, double z){return 0;}
double complex c(double x, double z){return 1;}
double complex d(double x, double z){return 0;}
double complex e(double x, double z){return 0;}
double complex f(double complex Z, double complex Y){return -Z*Y;}
double complex g(double complex Z, double complex Yj, double complex Yp, double complex Hp){return -Z*(Yj-Yp)*Hp;}
//programa principal
int main(int argc, char *argv )
{
  int n_cam=2,cont,j,i;//numero de camadas
  double frq=5*pow(10,3);//frequência
  double espessura[1]={200};//espessura das camadas intermediarias
  size_t n_h=length(espessura);
  double sigma[2];
  sigma[0]=pow(10,-3);sigma[1]=0.5;
  size_t n_s=length(sigma);
  double z0=0;
  double zn=espessura[0]*2;//profundidade da malha
  int nz=20;//numero de pontos em z
  double dz1=abs(zn-z0)/(nz-1);
  double z[nz];
  for (cont=0;cont<nz;cont++)
  {
    z[cont]=cont*dz1;
  }
  nz=length(z);
  double x0=0;
  double xn=10000;//largura da malha
  int nx=20;//numero de pontos em x, na malha reduzida
  double x[nx+4],xt[nx];
  x[0]=-xn+xn/4;x[1]=-xn+xn/2;
  double dx1=abs(xn-x0)/(nx-1);
  for (cont=0;cont<nx;cont++)
  {
    xt[cont]=cont*dx1;
    x[cont+2]=cont*dx1;
  }
  x[nx+2]=xn+xn/4;x[nx+3]=xn+xn/2;
  nx=length(x);//numero de pontos em x, na malha completa
  int nx2=length(xt);
  int np=nx*nz;
  int ne=(nx-1)*(nz-1)*2;
  double xv[np],zv[np];
  double complex Y[np],Yj[np],sig[np];
  double w =2*pi*frq;
  Y[0]=I*w*eps;
  double complex Z=I*w*mu;
  double complex K[np][np];
  double complex k[3][3],m[3];
  int np2=(nx-2)*(nz-2);
  double dt;
  double dx[3],dz[3];
  int i1,j1,l;
  int jc=0;
  double sigma1=0.5;//solo embaixo dos predios
  double sigma2=4.8;//agua salgada
  double sigma3=2.8;//proporcao de 1/2 de areia e cimento (predios), com fibras de aco carbono
  double complex Hp[np];
  double complex al[3];
  double complex bl[3];
  double complex cl[3];
  double complex dl[3];
  double complex El[3];
  double complex fl[3];
  double complex gl[3];
  double xv2[nz*nx2];
  double complex Hs[np];
  int vc[np],qt_front;
  int cc;
  int el[(nz-1)*(nx-1)*2][4];
  double complex M[np];
  int flag;
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


  ic=0;
  campo_prim(n_s,sigma,n_h,espessura,frq,zv,np,Hp);
  for (j=0;j<nz;j++)
  {
    for (i=0;i<nx2;i++)
    {
      xv2[ic]=x[i];
      ic++;
    }
  }
  //definicao de coeficientes da edo
  cc=-1;
  ic=0;
  FILE *fileaux;
  fileaux=fopen("el.dat","w");

  for (j=0;j<(nz-1);j++)
  {
    for (i=0;i<nx;i++)
    {
      if(i<nx-1)
      {
        cc++;
        el[cc][0]=ic+1;
        el[cc][1]=ic+2;
        el[cc][2]=ic+nx+1;
        el[cc][3]=cc;
        fprintf(fileaux, "%i\t",el[cc][0]);
        fprintf(fileaux, "%i\t",el[cc][1]);
        fprintf(fileaux, "%i\t",el[cc][2]);
        fprintf(fileaux, "%i\n",el[cc][3]);
        cc++;
        el[cc][0]=ic+2;
        el[cc][1]=ic+nx+2;
        el[cc][2]=ic+nx+1;
        el[cc][3]=cc;
        fprintf(fileaux, "%i\t",el[cc][0]);
        fprintf(fileaux, "%i\t",el[cc][1]);
        fprintf(fileaux, "%i\t",el[cc][2]);
        fprintf(fileaux, "%i\n",el[cc][3]);
      } else
      {
        cont++;
      }
      ic++;
    }
  }
  printf("%i\n",cc);
  printf("%i\n",(nz-1)*(nx-1)*2);
  printf("%i %i\n",nz,nx);
  //condicoes de fronteira
  qt_front=0;
  for (i=0;i<np;i++)
  {
    vc[i]=0;
    if (xv[i]==xv[1] || xv[i]==xv[np] || zv[i]==zv[1]|| zv[i]==zv[np])
    {
      Hs[i]=0; vc[i]=1;
      qt_front++;
    }
  }
  double complex Kr[np-qt_front][np-qt_front];
  double complex Mr[np-qt_front];
  //definicao dos meios secundarios
  //meio secundario 1 (solo)
  for (i=0;i<np;i++)
  {
    if (xv[i]>=3000 && xv[i]>=10000 && zv[i]>5000 && zv[i]<=6500)
    {
      Yj[i]=sigma1+I*w*eps; sig[i]=sigma1;
      //        printf("\n%i %f",i,creal(Yj[i]-Y[i]));
    }

    M[i]=0;
  }

  //meio secundario 2 (mar)
  for (i=0;i<np;i++)
  {
    if (xv[i]>=0 && xv[i]<=3000)
    {
      if ((zv[i]>=espessura[0]) && (zv[i]<(5*espessura[0])/4))
      {
        Yj[i]=sigma2+I*w*eps; sig[i]=sigma2;
        //          printf("\n%i %f",i,creal(Yj[i]-Y[i]));
      }
    }
  }
  ic=0;
  flag=0;
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
          //            printf("\n%i %f",i,creal(Yj[i]-Y[i]));
        }
      }
    }
  }
  //meio secundario 3 (predios)
  for (i=0;i<np;i++)
  {
    if (xv[i]>=6500 && xv[i]<=7000 && zv[i]>=espessura[0]-50 && zv[i]<espessura[0]) {
      Yj[i]=sigma3+I*w*eps;
      sig[i]=sigma3;
      //            printf("\n%i %f",i,creal(Yj[i]-Y[i]));
    }
    if (xv[i]>=7500 && xv[i]<=8000 && zv[i]>=espessura[0]-50 && zv[i]<espessura[0]) {
      Yj[i]=sigma3+I*w*eps;
      sig[i]=sigma3;
      //            printf("\n%i %f",i,creal(Yj[i]-Y[i]));
    }
    if (xv[i]>=8500 && xv[i]<=9500 && zv[i]>=espessura[0]-50 && zv[i]<espessura[0]) {
      Yj[i]=sigma3+I*w*eps;
      sig[i]=sigma3;
      //            printf("\n%i %f",i,creal(Yj[i]-Y[i]));
    }
  }

  // construcao da matriz e vetor da direita
  for (i=0;i<ne;i++){
    dx[0]=zv[el[i][1]]-zv[el[i][2]];
    dx[1]=zv[el[i][2]]-zv[el[i][0]];
    dx[2]=zv[el[i][0]]-zv[el[i][1]];
    dz[0]=xv[el[i][2]]-xv[el[i][1]];
    dz[1]=xv[el[i][0]]-xv[el[i][2]];
    dz[2]=xv[el[i][1]]-xv[el[i][0]];
    dt=abs(xv[el[i][1]]*zv[el[i][2]]+xv[el[i][0]]*zv[el[i][1]]+xv[el[i][2]]*zv[el[i][0]]-xv[el[i][1]]*zv[el[i][0]]-xv[el[i][2]]*zv[el[i][1]]-xv[el[i][0]]*zv[el[i][2]])/2;//delta
    //printf("%f\n",dt);
    for (cont=0;cont<3;cont++){
      al[cont]=a(xv[el[i][cont]],zv[el[i][cont]]);
      bl[cont]=b(xv[el[i][cont]],zv[el[i][cont]]);
      cl[cont]=c(xv[el[i][cont]],zv[el[i][cont]]);
      dl[cont]=d(xv[el[i][cont]],zv[el[i][cont]]);
      El[cont]=e(xv[el[i][cont]],zv[el[i][cont]]);
      fl[cont]=f(Z,Yj[el[i][cont]]);
      gl[cont]=g(Z,Yj[el[i][cont]],Y[el[i][cont]],Hp[el[i][cont]]);
      //fprintf(fileaux,"%f\n",creal(fl[cont]));
      // if (cabs(gl[cont])!= 0)
      // {
      // printf("\n%i %f %f",el[i][cont],creal(gl[cont]),cimag(gl[cont]));
      // }
      //printf("\n%i %f %f",el[i][cont],creal(gl[cont]),cimag(gl[cont]));
    }
    //matriz de cada elemento
    for (cc=0;cc<3;cc++)
    {
      for(j=0;j<3;j++)
      {
        k[cc][j]+=(-1/(4*dt))*(dx[1]*dx[1]*al[j]+dx[1]*dz[1]*bl[j]+dz[1]*dz[1]*cl[j])+(1/6)*(dx[1]*dl[j]+dz[1]*El[j])+(dt/12)*2*fl[j];
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
        M[el[i][j1]]=M[el[i][j1]]+m[j1];
        //printf("\n%5.15f %5.15f",creal(K[el[i][j1]][el[i][i1]]),cimag(K[el[i][j1]][el[i][i1]]));
        //printf("\n%5.15f %5.15f",creal(M[el[i][j1]]),cimag(M[el[i][j1]]));
      }
    }
  }

  //reducao do sistema - Inclusao das condicoes de fronteira
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


  //solucao por eliminacao gaussiana
  //1-triangularizacao
  double complex Ur[np-qt_front];
  double complex temp;
  for (i=0;i<(np-qt_front-1);i++)
  {
    for (j=i+1;j<np-qt_front;j++)
    {
      temp=Kr[j][i]/Kr[i][i];
      for (l=i;l<np-qt_front;l++)
      {
        Kr[j][l]-=temp*Kr[i][l];
      }
      Mr[j]-=temp*Mr[i];
    }
  }
  //2-retrossubstituicao
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
  //Incorporando as condicoes de fronteira
  double complex U[np];
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
  double complex Htotal[np];
  FILE *file;
  file=fopen("H_sec.dat","w");
  for (i=0;i<np;i++)
  {
    Htotal[i]=U[i]+Hp[i];
    //printf("%i %f %f\n",i,creal(Htotal[i]),cimag(Htotal[i]));
    fprintf(file,"%i %f %f\n",i,creal(Htotal[i]),cimag(Htotal[i]));
  }

  return 0;
}
