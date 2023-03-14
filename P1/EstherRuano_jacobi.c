#include <stdio.h>
#include <stdlib.h>
#include <math.h>

 
#define A 1
#define B 0
#define C 1
#define D 0
#define E 0
#define F 0

/* DEFINICIO DE FUNCIONS */
double itera_jacobi(double h, double ** mat, double coeficients[], int n);

double g(double x, double y);
double p(double x, double y);

void inicialitza(double h, double ** arr, double coeficients[], int n);
void inicialitza_zero(double ** arr, int n, double h);

void print_file(double** arr, int n, double h);

int main () {
	/*Entrades de l'usuari*/
	int n;
	double tolerancia, max_iteracions;
	char zeros;
	printf("Introdueix n: \n");
	scanf("%d", &n);
	
	printf("Introdueix la tolerancia (error): (p.ex 1e-8)\n");
	scanf("%lf", &tolerancia);

	printf("Introdueix el maxim d'iteracions: (p.ex 1e6)\n");
	scanf("%lf", &max_iteracions);
	
	printf("Vols inicialitzar a zeros?(y/n)\n");
	scanf("\n%c", &zeros);


	/* definim el pas */
	double h = 1./(n+1);
	
	/* definim la matriu resultat
	 p00  	p10  p20				pn0
	 p01 	u11, u21,       ,un1	pn1
			...
	 p0n	u1n,            ,unn	pnn
	 p0n+1							pn+1n+1
	 */
	int i, j;
	double err_prev = 1, error=1;
	double** arr = (double**)malloc((n+2) * sizeof(double*));
    for (i = 0; i < (n+2); i++){
        arr[i] = (double*)malloc((n+2) * sizeof(double));
	}
 
	/* valors per multiplicar les u */
	double coeficients[] = {(double)B/(4*h*h), (double)C/(h*h)-(double)E/(2*h), (double)A/(h*h)-(double)D/(2*h), (double)C/(h*h)+(double)E/(2*h), (double)A/(h*h)+(double)D/(2*h), -(double)(2*A +2*C)/(h*h) + F};
	
	if(zeros == 'y'){
		inicialitza_zero(arr, n, h);
	}else{
		inicialitza(h, arr, coeficients, n);
	}
	
	/* Print dels valors inicials
	printf("Inicialitzem: \n");
	for (i = 0; i < (n+2); i++){
        for (j = 0; j < (n+2); j++){
            printf("%.3f ", arr[i][j]);
		}
		printf("\n");
	}*/
 
	/* Iterem Jacobi */
	int k = 0;
	while(k < max_iteracions && error > tolerancia){
		error = itera_jacobi(h, arr, coeficients, n);
		k++;
		if(k%1000 == 0)
			printf("Iteracio: %d, Ratio d'errors: %f\n", k, (double)error/err_prev);
		err_prev = error;		
	}
	
	/* Desar resultat */
	print_file(arr, n, h);
	
	/* per veure si convergeix a p: miro diferencies al quadrat de valors*/
	double err = 0;
	for (i = 0; i < (n+2); i++){
        for (j = 0; j < (n+2); j++){
            err += (arr[i][j]-p(i*h, j*h))*(arr[i][j]-p(i*h, j*h));
		}
	}
	printf("Error real (norma euclidiana): %.10lf\n", err);
	
	/* display */
	if(k == max_iteracions){
		printf("No ha convergit en  %d iteracions\n", k);
	} else {
		printf("Ha convergit en %d iteracions\n", k);
	}
 
    /* alliberem memoria */
    for (i = 0; i < (n+2); i++){
        free(arr[i]);
	}
 
    free(arr);

	return 0;
}

/* inicialitza amb b/diag: b son els valors de g - "falses incognites" */
void inicialitza(double h, double ** arr, double coeficients[], int n){
	int i, j;
	/* Inicialitzem: perimetre: valors de p */
    for (i = 0; i < (n+2); i++){
        arr[i][0] = p(i*h, 0);
        arr[i][n+1] = p(i*h, 1);
	}
	for (j = 1; j < (n+1); j++){
        arr[0][j] = p(0, j*h);
        arr[n+1][j] = p(1, j*h);
	}
	
	
	/* interior: valors de g-falses_incognites/diag (bi/diag) */
	for (i = 1; i < (n+1); i++){
        for (j = 1; j < (n+1); j++){
			double suma = g(i*h, j*h);
			/* primera columna */
			if(i == 1){
				suma -= coeficients[0]*arr[i-1][j-1];
				suma -= coeficients[2]*arr[i-1][j];
				suma += coeficients[0]*arr[i-1][j+1];
			}
			/* ultima columna */
			if(i == n){
				suma += coeficients[0]*arr[i+1][j-1];
				suma -= coeficients[4]*arr[i+1][j];
				suma -= coeficients[0]*arr[i+1][j+1];
			}
			/* primera fila, ojo amb les cantonades */
			if(j == 1){
				if(i != 1){
					suma -= coeficients[0]*arr[i-1][j-1];
				} else if (i != n){
					suma += coeficients[0]*arr[i+1][j-1];
				}
				suma -= coeficients[1]*arr[i][j-1];
			}
			/* ultima fila, ojo amb les cantonades */
			if(j == n){
				if(i != 1){
					suma += coeficients[0]*arr[i-1][j+1];
				} else if (i != n){
					suma -= coeficients[0]*arr[i+1][j+1];
				}
				suma -= coeficients[3]*arr[i][j+1];
			}
			arr[i][j] = suma/coeficients[5];
		}
	}
}

/* inicialitza a zero */
 void inicialitza_zero(double ** arr, int n, double h){
	 int i, j;
	 
	 /* Inicialitzem: perimetre: valors de p */
    for (i = 0; i < (n+2); i++){
        arr[i][0] = p(i*h, 0);
        arr[i][n+1] = p(i*h, 1);
	}
	for (j = 1; j < (n+1); j++){
        arr[0][j] = p(0, j*h);
        arr[n+1][j] = p(1, j*h);
	}
	 
	/* interior a 0 */
	 for (i = 1; i < (n+1); i++){
		for (j = 1; j < (n+1); j++){
			arr[i][j] = 0;
		}
	}
 }
 
/* calcula la matriu */
double itera_jacobi(double h, double ** mat, double coeficients[], int n) {
	int i, j;
	double error, err_tmp;
	
	/* Com que es jacobi, he de copiar el vector resultat */
	double** arr = (double**)malloc((n+2) * sizeof(double*));
    for (i = 0; i < (n+2); i++){
        arr[i] = (double*)malloc((n+2) * sizeof(double));
	}
 
    /* Inicialitzem la copia */
    for (i = 0; i < (n+2); i++){
        for (j = 0; j < (n+2); j++){
            arr[i][j] = mat[i][j];
		}
	}
	
	/* Fare mat = D^-1(L+U) * arr + D^-1*b */
	for (i = 1; i < (n+1); i++){
        for (j = 1; j < (n+1); j++){
			double suma = g(i*h, j*h);
			suma -= coeficients[0]*(arr[i-1][j-1]-arr[i+1][j-1]-arr[i-1][j+1]+arr[i+1][j+1]);
			suma -= coeficients[1]*arr[i][j-1];
			suma -= coeficients[2]*arr[i-1][j];
			suma -= coeficients[3]*arr[i][j+1];
			suma -= coeficients[4]*arr[i+1][j];
			
			mat[i][j] = (double)suma/coeficients[5];
			
			err_tmp = fabs(mat[i][j]-arr[i][j]);
			
			if((i==1&&j==1) || err_tmp > error){
				error = err_tmp;
			}
		}
	}
	
	/* alliberem memoria */
    for (i = 0; i < (n+2); i++){
        free(arr[i]);
	}
 
    free(arr);
	return error;
}

/* funcio perimetre */
double p(double x, double y) {
	return pow(x,5)*pow(y,4);
}

/* funcio resultat */
double g(double x, double y){
	return A*20*pow(x,3)*pow(y,4) + C*12*pow(x,5)*pow(y,2) + B*20*pow(x,4)*pow(y,3) + D*5*pow(x,4)*pow(y,4) + E*4*pow(x,5)*pow(y,3) + F*pow(x,5)*pow(y,4);
}

/* Desa resultats en un arxiu: jacobi_result.txt*/
void print_file(double** arr, int n, double h){
	int i, j;
	FILE * f = fopen("jacobi_result.txt", "w");
	fprintf(f, "x\t\ty\t\tu\n");
	for (i = 0; i < (n+2); i++){
        for (j = 0; j < (n+2); j++){
            fprintf(f,"%+.6le\t%+.6le\t%+.6le\n", i*h, j*h, arr[i][j]);
		}
		fprintf(f,"\n");
	}
	fclose(f);
}