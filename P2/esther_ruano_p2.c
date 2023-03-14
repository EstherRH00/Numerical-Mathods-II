#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

#define num 30000 /*nombre maxim de punts calculats a cada corba*/
#define delta 0.005 /*distancia entre dos punts consecutius de la corba*/
#define tol 1e-5 /*no dividim entre 0*/
#define pre 1e-8 /*precisio a la correccio*/
#define iter 3 /*maxim d'iterats a la correccio*/

/* DEFINICIO DE FUNCIONS */
void pred_corr(double x_inicial, double y_inicial, int * d, FILE * f);
double norma_2(double x, double y);
double dist_euc(double x,double  y,double  a,double  b);
double dx_dist_euc(double x,double  y,double  a,double  b);
double dy_dist_euc(double x,double  y,double  a,double  b);
double h(double p, double q,double x0, double y0,  int * d);
double f(double p, double q, int * d);
double dx_h(double p, double q, int * d);
double dy_h(double p, double q, int * d);
double newton_raphson(double x0, double y0, double x_i, double y_i, int * digits, double * x_k);

int main () {
	int i, j;
	
	/* 
	 * Variables:
	 * d = {d[0], d[1], d[2], d[3], d[4], d[5], d[6], c8} - els nombres del niub
	 * bx = 0.25 * (d[0] + d[2] + d[4] + d[6]) - coordenada x del baricentre
	 * by = 0.25 * (d[1] + d[3] + d[5] + d[7]) - coordenada y del baricentre
	 * c = {0.2, 0.1, -0.05, -0.15} - constants per aconseguir els 16 punts inicials
	 */
	
	int * d = (int *)malloc(8*sizeof(int));
	double * c = (double *)malloc(4*sizeof(double));
	
	
	/* Omplir d */
	printf("Introdueix el teu niub, separat per espais (per exemple: 2 0 1 5 0 2 5 5):\n");
    for (i = 0; i < 8; i++) {
        scanf("%d", &d[i]);
    }
	
	/* Calcular baricentre (B) */
	double bx = 0.25 * (d[0] + d[2] + d[4] + d[6]);
	double by = 0.25 * (d[1] + d[3] + d[5] + d[7]);
	
	/* Omplir c */
	c[0] = 0.2;
	c[1] = 0.1;
	c[2] = - 0.05;
	c[3] = - 0.15;
	
	/* (px, py) - Punts inicials a partir del niub: (d[0], d[1]), (d[2], d[3]), (d[4], d[5]), (d[6], d[7]) */
	double px;
	double py;
	
	/* Fitxer resultat */
	FILE * f = fopen("p2_esther_ruano.txt", "w");
	
	for(i = 0; i < 4; i++){
		px = d[2*i];
		py = d[2*i+1];
		for(j = 0; j < 4; j++){
			/* 
			 * Per cada una de les constants, calculo el punt inicial
			 * (px - c(px-bx), py - c(py-by) 
			 */
			printf("Calculant la corba pel punt inicial: (%.4f, %.4f)\n", px + c[j] * (px-bx),py + c[j] * (py - by));
			fprintf(f,"Iniciant plot %d %d\n", i, j);
			pred_corr(px + c[j] * (px-bx),py + c[j] * (py - by), d, f);
			fprintf(f,"\nAcabant plot %d %d\n", i, j);
		}
		
	} 
	
	fclose(f);
	
	/*
	 * Per fer el plot:
	 * plot "p2_esther_ruano.txt" u 2:3 w l
	 */
	 
	return 0;
}

/*
 * Calcula una de les corbes fent prediccio - correccio
 */ 
void pred_corr(double x_inicial, double y_inicial, int * d, FILE * f){
	
	/* 
	 * Variables:
	 * (x0, y0): punt inicial; es manté sempre per tal de poder cridar la funció h
	 * (x_i, y_i): inicialment es el punt inicial, despres va canviant per cada iteracio de pred-corr
	 * z = {z0, z1}: en cada iteracio, es inicialment el punt predit i va canviat amb les diferents correccions
	 * (hx, hy): derivades d'h al punt(x_i, y_i)
	 * norm: norma del vector (hx, hy)
	 * fi: variable booleana que controla que el vector tangent no sigui massa petit o la correccio no convergeixi
	 *      (la seva funcio es no posar un break)
	 */
	 
	double x_0 = x_inicial;
	double y_0 = y_inicial;
	
	double x_i = x_inicial;
	double y_i = y_inicial;
	
	double * z = (double *)malloc(2*sizeof(double));
	
	double hx, hy, norma;
	
	int i = 0;
	
	bool fi = false; /* Valdra true si no s'ha pogut fer la prediccio (valor tangent massa proper al nul) 
						o si no s'ha pogut fer la correccio (NR no convergeix en iter iteracions) */
						
	/* Pararem quan:
	 * 	(a) Haguem fet massa iteracions
	 * 	(b) Siguem massa aprop del punt inicial, ja hem tancat la corba (*) 
	 * 	(c) Hi hagi hagut algun incident (variable fi = True)
	 *
	 * 	(*) Cal ignorar la condicio el pimer cop que entrem al bucle, ja que partim de l'inicial
	 * 		A mes, pot passar que la primera volta no s'allunyi prou del punt inicial, per tant
	 *		podem, o be ignorar-ho tambe a la primera volta, o demanar que la distancia sigui menor
	 *		que delta - pre, es a dir, que en el "pas" entre punts afegim tolerancia a l'error permes
	 *		en la correccio
	 */
	while ( i < num && ( i == 0 || dist_euc(x_i, y_i, x_0, y_0) >= (delta - pre)) && !fi){
		
		/* Predicció */
		/* 1. el gradient d'h*/
		hx = dx_h(x_i, y_i, d);
		hy = dy_h(x_i, y_i, d);
		
		norma = norma_2(hx, hy);
		
		/* 2. Predicció inicial*/
		if(norma > tol){
			
			
			/* Valors predits */
			z[0] = x_i + delta * hy/norma;
			z[1] = y_i - delta * hx/norma;
			
			/* 3. Correcció */
			int j = 0;
			double err = 1;
			while(j < iter && (err >= pre) && !fi){
				err = newton_raphson(x_0, y_0, x_i, y_i, d, z);
				j++;
			}
			
			if(err > pre){
				fi = true;
				printf("No s'ha pogut fer la correccio, despres de %d iteracions l'error era %.6f\n", iter, err);
			}
			
			/* Actualitzar el punt */
			x_i = z[0];
			y_i = z[1];
			
			/* Desar al fitxer */
			fprintf(f,"%d\t%+.6le\t%+.6le\n", i, x_i, y_i);
				
		} else {
			fi = true;
			printf("No s'ha pogut fer la prediccio, el vector tangent es massa proper a zero, te norma %.6f\n",norma);
		}
		
		i ++;
 
	}
	free(z);
}

/*
 * Calcula la norma euclidiana del vector(x, y) 
 */
double norma_2(double x, double y){
	return sqrt(pow(x,2) + pow(y, 2));
}

/*
 * Calcula la distancia euclidea entre els punts (x, y) i (a, b)
 */
double dist_euc(double x,double  y,double  a,double  b){
	return sqrt(pow(x-a,2) + pow(y-b, 2));
}

/*
 * Calcula la derivada respecte x de la distancia euclidea entre els punts (x, y) i (a, b)
 */
double dx_dist_euc(double x,double  y,double  a,double  b){
	return (double) (x-a)/sqrt(pow(x-a,2) + pow(y-b, 2));
}

/*
 * Calcula la derivada respecte y de la distancia euclidea entre els punts (x, y) i (a, b)
 */
double dy_dist_euc(double x,double  y,double  a,double  b){
	return (double) (y-b)/sqrt(pow(x-a,2) + pow(y-b, 2));
}

/*
 * Calcula f al punt (p, q)
 */
double f(double p, double q, int * d){
	return dist_euc(p,q,d[0],d[1])*dist_euc(p,q,d[2],d[3])*dist_euc(p,q,d[4],d[5])*dist_euc(p,q,d[6],d[7]);
}


/*
 * Calcula h(p,q) = f(p,q) - f(x0, y0) al punt (p, q)
 */
double h(double p, double q, double x0, double y0, int * d){
	return f(p, q, d) - f(x0, y0, d);
}

/*
 * Calcula la derivada de h respecte x al punt (p, q) aplicant la regla del producte
 */
double dx_h(double p, double q, int * d){
	return dx_dist_euc(p,q,d[0],d[1])*dist_euc(p,q,d[2],d[3])*dist_euc(p,q,d[4],d[5])*dist_euc(p,q,d[6],d[7]) +
	dist_euc(p,q,d[0],d[1])*dx_dist_euc(p,q,d[2],d[3])*dist_euc(p,q,d[4],d[5])*dist_euc(p,q,d[6],d[7])+
	dist_euc(p,q,d[0],d[1])*dist_euc(p,q,d[2],d[3])*dx_dist_euc(p,q,d[4],d[5])*dist_euc(p,q,d[6],d[7])+
	dist_euc(p,q,d[0],d[1])*dist_euc(p,q,d[2],d[3])*dist_euc(p,q,d[4],d[5])*dx_dist_euc(p,q,d[6],d[7]);
}

/*
 * Calcula la derivada de h respecte y al punt (p, q) aplicant la regla del producte
 */
double dy_h(double p, double q, int * d){
	return dy_dist_euc(p,q,d[0],d[1])*dist_euc(p,q,d[2],d[3])*dist_euc(p,q,d[4],d[5])*dist_euc(p,q,d[6],d[7]) +
	dist_euc(p,q,d[0],d[1])*dy_dist_euc(p,q,d[2],d[3])*dist_euc(p,q,d[4],d[5])*dist_euc(p,q,d[6],d[7])+
	dist_euc(p,q,d[0],d[1])*dist_euc(p,q,d[2],d[3])*dy_dist_euc(p,q,d[4],d[5])*dist_euc(p,q,d[6],d[7])+
	dist_euc(p,q,d[0],d[1])*dist_euc(p,q,d[2],d[3])*dist_euc(p,q,d[4],d[5])*dy_dist_euc(p,q,d[6],d[7]);
}

/*
 * Metode de Newton Raphson:
 * 
 *  x_(k+1) = x_(k) − [DF(x_k))]^−1 * F(x_(k))
 * 
 * Amb:
 * 	 		[  h(x, y)                        ]    [ p ]
 * 		F = [ (x-xi)^2 + (y-yi)^2 - delta^2   ]  = [ q ]
 *
 *       	 [ dx_h(x, y)        dy_h(x, y)   ]    [ a    b ]
 *		DF = [ 2(x-xi) 			2(y-yi)       ]	 = [ c    d ]
 *
 *                              [  d    -b ]
 *      DF^-1 = 1/(a*d - b*c) * [ -c     a ]
 *
 * Paramentres:
 * 	x0, y0: punt original per poder cridar h
 * 	x_i, y_i: punt des del qual hem fet la prediccio
 * 	digits = {d1, d2, d3, d4, d5, d6, d7, d8}
 *  x_k = {x_k, x_k+1} punt al qual volem aplicar la correccio
 * 
 * Rebo x_k a la variable "x_k" i hi desaré x_k+1
 * 
 * Retorno l'error com la diferencia entre x_k+1 i x_k
 */
double newton_raphson(double x0, double y0, double x_i, double y_i, int * digits, double * x_k){
	/* Coeficients DF */
	double a = dx_h(x_k[0], x_k[1], digits);
	double b = dy_h(x_k[0], x_k[1], digits);
	double c = 2*(x_k[0] - x_i);
	double d = 2*(x_k[1] - y_i);
	
	/* Coeficients de F */
	double p = h(x_k[0], x_k[1], x0, y0, digits);
	double q = pow((x_k[0] - x_i), 2) + pow((x_k[1] - y_i), 2) - pow(delta, 2);
	
	/* Calculo resultat x_k+1 */
	double x_k1 = x_k[0] - (double) 1 / (a*d - b*c) * (d * p - b * q);
	double y_k1 = x_k[1] - (double) 1 / (a*d - b*c) * (-c * p + a * q);
	
	/* Calculo error */
	double err = dist_euc(x_k1, y_k1, x_k[0], x_k[1]);
	
	/* Actualitzo x_k amb x_k+1 */
	x_k[0] = x_k1;
	x_k[1] = y_k1;
	
	return err;
}
