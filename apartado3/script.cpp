/* SIMULACIÓN DE LA DINÁMICA MOLECULAR DE UN GAS CON UN POTENCIAL DE LENNARD-JONES */
/* EJERCICIO VOLUNTARIO 1 LECCIÓN 1 SISTEMA SOLAR */

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <time.h>

#define PI 3.14159265


using namespace std;

/* Definimos unas funciones necesarias en nuestro programa */
double distancia (double x1, double y1, double x2, double y2, double dif_coord[]);
double posicion (double x, double v, double& momentosparedes);
double fuerza (double x1, double y1, double x2, double y2, double dif_coord[]);

int main (void)
{
    int n; n=16;
    double x[n], y[n], vx[n], vy[n], ax[n], ay[n], wx[n], wy[n], dif_coord[2];
    double t,h,energia,T,V,f,r,r6,omega,temp,aux;
    double momentosparedes,aux2; // Suma de los momentos de cada partícula cuando atravesiesan las paredes. La variable auxiliar la uso para obtener la presión del gas.
    int i,j,k;
    bool insertar;
    //ofstream fichero_particulas;
    ofstream fichero_energia;
    //ofstream fichero_velocidades;
    ofstream fichero_resultados;
    srand(time(NULL));

    //fichero_particulas.open("/Users/alzorrcarri/Documents/cphys/leccion1-voluntario/particles_positions_data.dat");
    fichero_energia.open("/Users/alzorrcarri/Documents/cphys/leccion1-voluntario/presion_tiempo.dat");
    //fichero_velocidades.open("/Users/alzorrcarri/Documents/cphys/leccion1-voluntario/velocity_time.dat");
    fichero_resultados.open("/Users/alzorrcarri/Documents/cphys/leccion1-voluntario/resultados.dat");


    /* Damos valores iniciales a las posiciones y velocidades de las partículas de forma aleatoria. Las posiciones iniciales
    siguen una distribución uniforme en el cuadrado, y las velocidades tienen todas módulo uno y ángulo aleatorio uniforme en 0-2pi. */
    /* Para que las partículas no salgan disparadas al principio, imponemos una speración mínima entre ellas desde el principio. */

    cout << "X" << '\t' << "Y" << '\t' << "V_X" << '\t' << "V_Y" << endl;
    for ( i = 0; i <= n-1; i++)
    {
      do
      {
          x[i] = rand() % 100000;
          x[i] = x[i]*10/100000;
          y[i] = rand() % 100000;
          y[i] = y[i]*10/100000;
          insertar=true;
          for ( j = 0; j < i; j++)
          {
            if (distancia(x[i],y[i],x[j],y[j],dif_coord)<=1.1)
            {
               insertar=false;
            }
          }
          
      } while (insertar==false);
        
      /* Consideramos velocidades iniciales de módulo 1,2,3,4 y 5 (así se consigue variar la temperatura del sistema). */
      omega = rand() % 100000;
      omega = omega*2*PI/100000;
      vx[i] = 5*cos(omega);
      vy[i] = 5*sin(omega);
      cout << x[i] << '\t' << y[i] << '\t' << vx[i] << '\t' << vy[i] << endl; // Vemos si todo es correcto.
    }
   
    /* Iniciamos el proceso iterativo de cálculo de las aceleraciones, posiciones, vector auxiliar, nuevas aceleraciones
    y velocidades según el algoritmo de Verlet. */
   h=0.002;
   t=0;
   temp=0;
   momentosparedes=0; aux2=0;
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

   /* ALGORTIMO DE VERLET */

   for ( k = 0; k < 30000; k++)
   {
      /* Escribimos las posiciones de las partículas, el módulo de la velocidad y la energía del sistema en los ficheros correspondientes antes de pasar al siguiente tiempo. */
      for ( i = 0; i <= n-1; i++)
      {
         //fichero_particulas << x[i] << ',' << y[i] << endl;
      }
      //fichero_particulas << endl;

      if (k*h>=20 && k*h<=50)
      {  
         aux2=aux2+momentosparedes/(4.0*10*t);

         aux=0;
         for ( i = 0; i <= n-1; i++)
         {
            aux=aux+vx[i]*vx[i]+vy[i]*vy[i];
            //fichero_velocidades << i+1 << '\t' << vx[i] << '\t' << vy[i] << '\t' << sqrt(vx[i]*vx[i]+vy[i]*vy[i]) << '\n';
         }
         aux=aux/20.0;
         temp=temp+aux;
      }
      

      energia=0; T=0; V=0;
      for ( i = 0; i <= n-1; i++)
      {
         T=T+(vx[i]*vx[i]+vy[i]*vy[i])/2.0;
         for ( j = 0; j <= n-1; j++)
         {
            if (j!=i)
            {
               r=distancia(x[i],y[i],x[j],y[j],dif_coord);
               r6=pow(r,6);
               V=V+4.0*(1/r6-1)/r6;
            }
         } 
      }
      V=V/2.0;
      energia=T+V;
      /* Escribimos la presión en las paredes en el archivo de las energías. La presión es momentosparedes entre el tiempo y entre 4*10 */
      fichero_energia << t << '\t' << momentosparedes/(4.0*10*t) << endl;
      

      /* Actualizamos el tiempo y calculamos las nuevas posiciones y velocidades de los planetas. */
      t=t+h;

      for (i = 0; i <= n-1; i++)
      {
         x[i]=x[i]+h*vx[i]+h*h*ax[i]/2;
         x[i]=posicion(x[i],vx[i],momentosparedes);
         y[i]=y[i]+h*vy[i]+h*h*ay[i]/2;
         y[i]=posicion(y[i],vy[i],momentosparedes);
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

    temp=temp/(30.0/h);
    aux2=aux2/(30.0/h);
    fichero_resultados << "Temperatura del sistema estimada por el teorema de equipartición:" << '\t' << temp/2.0 <<'\n';
    fichero_resultados << "La presión media del sistema entre los instantes 20 seg y 50 seg es:" << '\t' << aux2 << '\n';

    //fichero_particulas.close();
    fichero_energia.close();
    //fichero_velocidades.close();
    fichero_resultados.close();
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
      difx2=x1-x2+10;
    } else difx2=x1-x2-10;
    if (y1<=y2)
    {
      dify2=y1-y2+10;
    } else dify2=y1-y2-10;

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
double posicion (double x, double v, double& momentosparedes)
{
    if (x>10.0 || x<0.0)
    {
        momentosparedes=momentosparedes+2*abs(v);
        x=x-10*floor(x/10.0);
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
