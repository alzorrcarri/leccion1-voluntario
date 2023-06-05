/* SIMULACIÓN DE LA DINÁMICA MOLECULAR DE UN GAS CON UN POTENCIAL DE LENNARD-JONES */
/* EJERCICIO VOLUNTARIO 1 LECCIÓN 1 SISTEMA SOLAR */

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>


using namespace std;

/* Definimos unas funciones necesarias en nuestro programa */
double distancia (double x1, double y1, double x2, double y2, double dif_coord[]);
double posicion (double x);
double fuerza (double x1, double y1, double x2, double y2, double dif_coord[]);

int main (void)
{
    int n; n=6;
    int cont;
    double x[n], y[n], vx[n], vy[n], ax[n], ay[n], wx[n], wy[n], dif_coord[2];
    double t,h,energia,T,V,f,r,r6,aux,temp;
    int i,j,k;
    bool insertar;

    ofstream fichero_densidad;

    srand(time(NULL));

    fichero_densidad.open("/Users/alzorrcarri/Documents/cphys/leccion1-voluntario/apartado8/densidad_gas.dat");


    /* Damos valores iniciales a las posiciones de las 16 partículas dentro de la cuadrícula 4x4 */
    /* En el caso del gas, consideramos 6 partículas dispuestas aletoriamente. */

    cout << "X" << '\t' << "Y" << '\t' << "V_X" << '\t' << "V_Y" << endl;
   /*for ( i = 0; i < 4; i++)
   {
      for ( j = 0; j < 4; j++)
      {
         x[4*i+j]=i+0.5;
         y[4*i+j]=j+0.5;
      }
      
   }*/
   for ( i = 0; i <= n-1; i++)
    {
      do
      {
          x[i] = rand() % 100000;
          x[i] = x[i]*4/100000;
          y[i] = rand() % 100000;
          y[i] = y[i]*4/100000;
          insertar=true;
          for ( j = 0; j < i; j++)
          {
            if (distancia(x[i],y[i],x[j],y[j],dif_coord)<=1.1)
            {
               insertar=false;
            }
          }
          
      } while (insertar==false);
    }
        
      /* Consideramos velocidades iniciales nulas. */
   for ( i = 0; i <= n-1; i++)
    {
      vx[i] = 0;
      vy[i] = 0;
      cout << x[i] << '\t' << y[i] << '\t' << vx[i] << '\t' << vy[i] << endl; // Vemos si todo es correcto.
    }

    /* Damos formato al fichero de escritura de datos */
    fichero_densidad << "t";
    for ( i = 1; i <= 40; i++)
    {
      fichero_densidad << '\t' << ((1.0+(i-1)/20.0)+(1.0+i/20.0))/2.0;
    }
    fichero_densidad << '\n';
    
   
    /* Iniciamos el proceso iterativo de cálculo de las aceleraciones, posiciones, vector auxiliar, nuevas aceleraciones
    y velocidades según el algoritmo de Verlet. */
   h=0.001;
   t=0;
   for ( i = 0; i <= n-1; i++)
    {
        ax[i]=0;
        ay[i]=0;
    }

   /* En primer lugar el cálculo de la aceleración a t=0. */
   for ( i = 0; i <= n-1; i++)
   {
      for ( j = 0; j <= n-1; j++)
      {
         if (j!=i)
         {
            f=fuerza(x[i],y[i],x[j],y[j],dif_coord);
            ax[i]=ax[i]+f*dif_coord[0];
            ay[i]=ay[i]+f*dif_coord[1];
         }   
      }
   }

   /* Cálculo de la función de correlación. */
   /* En el caso de la fase sólida, dejamos el sistema evolucionar hasta el equilibrio. */
   /* Para la fase líquida, calentamos el sistema hasta que su tempertura sea mayor a la temperatura crítica caluclada en el apartado 7: 1.383694. */
   /* Para la fase gaseosa, calentamos el gas muchísimo. */

   for ( k = 0; k <= 30000; k++)
   {
      /* Cálculo de la temperatura del sistema. */
      temp=0;
      for ( i = 0; i <= n-1; i++)
      {
         temp=temp+vx[i]*vx[i]+vy[i]*vy[i];
      }
      temp=temp/(2*16.0);


      if (t>=10.0 && temp > 5)
      {
         if (k%10==0)
         {
            fichero_densidad << t;
            for ( i = 1; i <= 40; i++)
            {
               cont=0;
               for ( j = 1; j <= n-1; j++)
               {
                  aux=distancia(x[j],y[j],x[0],y[0],dif_coord);
                  if (1.0+(i-1)/20.0 < aux && aux < 1.0+i/20.0)
                  {
                     cont++;
                  }
               }
               fichero_densidad << '\t' << cont;
            }
            fichero_densidad << '\n';
            
         }
         
      }

      /* Calentamos el sistema en ciertos instantes de tiempo si la temperatura no es la suficiente para estar en fase líquida. */
      if ((k*h==10.0) || (k*h==12.0) || (k*h==14.0) || (k*h==16.0) || (k*h==18.0) || (k*h==20.0) && temp < 5)
      {
         for ( i = 0; i < n-1; i++)
         {
            vx[i]=2*vx[i];
            vy[i]=2*vy[i];
         }
         
      }

      /* Actualizamos el tiempo y calculamos las nuevas posiciones y velocidades de los planetas. */
      t=t+h;
      
      /* ALGORTIMO DE VERLET */

      for (i = 0; i <= n-1; i++)
      {
         x[i]=x[i]+h*vx[i]+h*h*ax[i]/2;
         x[i]=posicion(x[i]);
         y[i]=y[i]+h*vy[i]+h*h*ay[i]/2;
         y[i]=posicion(y[i]);
      }
      for (i = 0; i <= n-1; i++)
      {
         wx[i]=vx[i]+h*ax[i]/2;
         wy[i]=vy[i]+h*ay[i]/2;
         ax[i]=0; ay[i]=0;
      }

      for ( i = 0; i <= n-1; i++)
      {
         for ( j = 0; j <= n-1; j++)
         {
            if (j!=i)
            {
               f=fuerza(x[i],y[i],x[j],y[j],dif_coord);
               ax[i]=ax[i]+f*dif_coord[0];
               ay[i]=ay[i]+f*dif_coord[1];
            }
         }
      }
      for ( i = 0; i <= n-1; i++)
      {
         vx[i]=wx[i]+h*ax[i]/2;
         vy[i]=wy[i]+h*ay[i]/2;
      }
   }

    fichero_densidad.close();
    return 0;
}

/* Definición de las funciones empleadas */
/* La función distancia calcula la distancia entre dos partículas teniendo en cuenta las condiciones de contorno periódicas. */
/* Usamos que los arrays siempre se pasan por referencia en las funciones. */
double distancia (double x1, double y1, double x2, double y2, double dif_coord[])
{
    double difx, dify, difx2, dify2;
    double d1,d2,d3,d4;
    double r;

    difx=x1-x2;
    dify=y1-y2;
    if (x1<=x2)
    {
      difx2=x1-x2+4;
    } else difx2=x1-x2-4;
    if (y1<=y2)
    {
      dify2=y1-y2+4;
    } else dify2=y1-y2-4;

    // Calculamos las cuatro posibles distancias y escogemos la mínima.
    d1=difx*difx+dify*dify;
    d2=difx*difx+dify2*dify2;
    d3=difx2*difx2+dify*dify;
    d4=difx2*difx2+dify2*dify2;

    if (d1<=d2 && d1<=d3 && d1<=d4)
    {
      r=sqrt(d1);
      dif_coord[0]=difx;
      dif_coord[1]=dify;
    }
    if (d2<=d1 && d2<=d3 && d2<=d4)
    {
      r=sqrt(d2);
      dif_coord[0]=difx;
      dif_coord[1]=dify2;
    }
    if (d3<=d1 && d3<=d2 && d3<=d4)
    {
      r=sqrt(d3);
      dif_coord[0]=difx2;
      dif_coord[1]=dify;
    }
    if (d4<=d1 && d4<=d2 && d4<=d3)
    {
      r=sqrt(d4);
      dif_coord[0]=difx2;
      dif_coord[1]=dify2;
    }
    return r;   
    
}

/* La función posición comprubea que todas las partículas queden dentro de la caja, y aquéllas que están fuera, las mete dentro
siguiendo las condiciones de contorno periódicas. */
double posicion (double x)
{
    if (x>4.0 || x<0.0)
    {
        x=x-4*floor(x/4.0);
    }
    return x;
}

/* La función fuerza devuelve el valor del módulo de la fuerza entre dos partículas según la distancia que las separa.*/
double fuerza (double x1, double y1, double x2, double y2, double dif_coord[])
{
   double r,r6,f;
   r=distancia(x1,y1,x2,y2,dif_coord);
   r6=pow(r,6);
   f=24*(2.0/r6-1)/(r6*r*r);
   return f;
}
