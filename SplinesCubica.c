//****************************************************************//
//**		Autor: Leonardo Marcos Santiago						  //
//**		Programa de segmentacion Cúbica						  //
//****************************************************************//

#include <stdio.h>
#include <stdlib.h>
#define N 10//Numero de puntos
#define Max 100
float X[N]={-0.102,0.915,1.876,3.049,4.067,5.007,6.041,6.864,7.959,9.108};//valores de x
float Y[N]={9.052,8.804,8.953,8.334,8.002,7.414,6.665,5.613,4.093,0.3};//valores de y
void Cancela(float Mat[Max][Max],int Filas, int FilaP, int FilaC);
void Divide(float Mat[Max][Max], int Filas, int Fila);
void DespM(float Mat[Max][Max],int Filas);
void LlenaM(float Mat[Max][Max],int Filas);
void DespM1(float Mat[Max][Max],int Filas);
void Cambiar_Fila(float M[Max][Max],int Filas,int i);
void Gauss();
void Generador_Polinomios_Splines(float M[Max][Max]);
//funcion para la segmentacion cùbica
int main(){
	float M[(N-1)*4][(N-1)*4+1];//matriz generada
	int i,a,j;//a para obtener valores del vector y
	a=0;
	for(i=0;i<(N-1)*4;i++){//llenar matriz con ceros
		for(j=0;j<(N-1)*4+1;j++){
			M[i][j]=0;
		}
	}
	for(i=0;i<(N-1)*2;i++){//guardar valores de Y en la ultima columna de matriz generada
		if(i!=0 && i!=(N-1)*2-1){
			M[i][(N-1)*4]=Y[a];//si el valor de y no es el primero ni el ultimo, entonces repetir dos veces
			i++;
			M[i][(N-1)*4]=Y[a++];
		}else{
			M[i][(N-1)*4]=Y[a++];
		}
	}
	//******************************************************************
	printf("\nValores guardados en la ultima columna\n");
	for(i=0;i<(N-1)*4;i++){
		printf("\n%f",M[i][(N-1)*4]);
	}
	//******************************************************************
		//primera derivada total(N-1)*2// obteniendo primera derivada
	a=0;
	j=-4;
	for(i=0;i<(N-1)*2;i++){
		if(i%2==0){
			j=j+4;
		}
		else{
			j=j;
			a++;
		}
		M[i][j]=(X[a])*(X[a]*X[a]);
		M[i][j+1]=(X[a]*X[a]);
		M[i][j+2]=X[a];
		M[i][j+3]=1;
		
	}
	//******************************************************************
	//imprime la matriz completa
	printf("Primera derivada");
	for(i=0;i<(N-1)*4;i++){
		printf("\n");
		for(j=0;j<(N-1)*4+1;j++)
			printf("[%f]",M[i][j]);
	}//*****************************************************************
	//igualando primeras derivadas para discontinuidad,
	j=-4;
	a=1;
	int b=0;
	for(i=(N-1)*2;b<N-2;i++,b++){
		j=j+4;
		M[i][j]=3*X[a]*X[a];
		M[i][j+1]=2*X[a];
		M[i][j+2]=1;
		M[i][j+3]=0;
		M[i][j+4]=-3*X[a]*X[a];
		M[i][j+5]=-2*X[a];
		M[i][j+6]=-1;
		M[i][j+7]=0;
		a++;
	}
	int derivada1=i;//guardar valor de i para segunda derivada
	//******************************************************************
	//imprime la matriz completa
	printf("\n");
	for(i=0;i<(N-1)*4;i++){
		printf("\n");
		for(j=0;j<(N-1)*4+1;j++)
			printf("[%f]",M[i][j]);
	}//*****************************************************************
	//igualando Segundas derivadas para discontinuidad,
	j=-4;
	a=1;
	b=0;
	for(i=derivada1;b<N-2;i++,b++){
		j=j+4;
		M[i][j]=6*X[a];
		M[i][j+1]=2;
		M[i][j+2]=0;
		M[i][j+3]=0;
		M[i][j+4]=-6*X[a];
		M[i][j+5]=-2;
		M[i][j+6]=0;
		M[i][j+7]=0;
		a++;
	}
	//imprime la matriz completa
	printf("\nMatriz con segunda derivada");
	for(i=0;i<(N-1)*4;i++){
		printf("\n");
		for(j=0;j<(N-1)*4+1;j++)
			printf("[%f]",M[i][j]);
	}//*****************************************************************
	//Agregando condiciones S´´(X0)=0 and S´´(Xn)=0
	//para punto inicial donde S´´(X0)=0
	M[(N-1)*4-2][0]=3*X[0];
	M[(N-1)*4-2][1]=1;
	//para punto final donde S´´(Xn)=0
	M[(N-1)*4-1][(N-1)*4-4]=3*X[N-1];
	M[(N-1)*4-1][(N-1)*4-3]=1;
	
	printf("\nMatriz generado para resolver");
	for(i=0;i<(N-1)*4;i++){
		printf("\n");
		for(j=0;j<(N-1)*4+1;j++)
			printf("[%f]",M[i][j]);
	}//*****************************************************************
	FILE *out;
	out=fopen("Resolver.dat","w");
	for(i=0;i<(N-1)*4;i++){
		for(j=0;j<=(N-1)*4;j++){
			//printf("[%f]",M[i][j]);
			fprintf(out,"%f\t",M[i][j]);
		}
		fprintf(out,"\n");
	}
	fclose(out);
	printf("\nResolviendo con GaussJordan\n");
	//El archivo Resolver.dat 
	//reslverlo con gauss de manera automática
	Gauss();//método para resolver la matriz	
	return 0;
}
//********************************************************************************************
////******************************************************************************************
////******************************************************************************************
////******************************************************************************************
//********************************************************************************************
void Gauss(){
  int Filas=(N-1)*4;
  int i,j;
  float Matriz[Max][Max];
  LlenaM(Matriz,Filas);
  DespM(Matriz,Filas);

  for(i=0;i<=Filas;i++)
  {
	if(Matriz[i][i]==0)//No se puede dividir entre 0, hacer intercambio de filas
		Cambiar_Fila(Matriz,Filas,i);
    Divide(Matriz,Filas,i+1);
    for(j=i+1;j<Filas;j++){
      //printf("En ciclo de cancelación\n");
      Cancela(Matriz,Filas,i+1,j+1);
    }
  }
  for(i=Filas-1;i>=0;i--){
    for(j=i-1;j>=0;j--){
      Cancela(Matriz,Filas,i+1,j+1);
    }
  }
	printf("\nMatriz Resuelta\n");
	DespM1(Matriz,Filas);//imrime en pantalla
	DespM(Matriz,Filas);//imprime en archivo, con nombre Resuelto.dat
	Generador_Polinomios_Splines(Matriz);
}
void LlenaM(float Mat[Max][Max],int Filas){
  int i,j;
 	FILE *fp;
 	fp = fopen ( "Resolver.dat", "r" );
	for (i=0;i<=Filas;i++){
		for(j=0;j<=Filas;j++){
			fscanf(fp, "%f" ,&Mat[i][j]);
			printf("%f\t",Mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
 	fclose ( fp );
 	printf("\nMatriz leida del método de Segmentacion cùbica\n");
 	DespM1(Mat,Filas);
} // fin de LlenaM
void DespM(float Mat[Max][Max],int Filas){
	////************procedimiento para imprimir en *****///// archivo** /////en resultado
	FILE *out;
	out=fopen("Resuelto.dat","w");
	int i,j;
	for(i=0;i<Filas;i++){
		for(j=0;j<=Filas;j++){
			//printf("%f ",Mat[i][j]);
			fprintf(out,"%f\t",Mat[i][j]);
		}
		fprintf(out,"\n");
	}
	fclose(out);
} // fin de DespM
void Cambiar_Fila(float M[Max][Max],int Filas,int Fila){
	int i,j;
	float aux;
	for(i=Fila+1;i<Filas;i++){
		if(M[i][Fila]!=0){
			for(j=0;j<=Filas;j++){
				aux=M[Fila][j];
				M[Fila][j]=M[i][j];
				M[i][j]=aux;
				
			}
		}
	}
		//printf("\nFila intercambiada\n");
}
void DespM1(float Mat[Max][Max],int Filas){
	////************procedimiento para imprimir en pantalla en resultado
	int i,j;
	for(i=0;i<Filas;i++){
		for(j=0;j<=Filas;j++){
			printf("%f ",Mat[i][j]);
		}
		printf("\n");
	}
} // fin de DespM
void Divide(float Mat[Max][Max], int Filas, int Fila){
  int i;
  float Aux;

  Aux=Mat[Fila-1][Fila-1];
  for(i=Fila-1;i<=Filas;i++){
    Mat[Fila-1][i]=Mat[Fila-1][i]/Aux;
  }
}
void Cancela(float Mat[Max][Max],int Filas, int FilaP, int FilaC){
  int i;
  float Aux;
  Aux=Mat[FilaC-1][FilaP-1];
	#ifdef TEST
	  printf("Valor de coeficiente a anular: %f\n",Aux);
	  printf("Función Cancela: Aux = %f\n",Aux);
	  DespM(Mat,Filas);
	#endif
  for(i=FilaP-1;i<=Filas;i++){
    Mat[FilaC-1][i]=Mat[FilaC-1][i]-Aux*Mat[FilaP-1][i];
  }
}

void Generador_Polinomios_Splines(float M[Max][Max]){
	int i,a;
	FILE *Salida;
	Salida =fopen("PolSegCubico.dat","w");
	fprintf(Salida,"set parametric\n");
	fprintf(Salida,"set trange[0:1]\n");  
	float delta[N];
	for(i=0;i<N-1;i++){
		delta[i]=X[i+1]-X[i];
	}
	//Generando primer polinomio
	printf("\nPolinomios Generados\n");
	fprintf(Salida,"plot t*%f+(%f), %f*(t*%f+%f)*(t*%f+%f)*(t*%f+%f)+%f*(t*%f+%f)*(t*%f+%f)+%f*(t*%f+%f)+(%f)\n",
					delta[0],X[0],M[0][(N-1)*4],delta[0],X[0],delta[0],X[0],delta[0],X[0],M[1][(N-1)*4],delta[0],X[0],delta[0],X[0],
					M[2][(N-1)*4],delta[0],X[0],M[3][(N-1)*4]);//imprimir en archivo para graficar en gnuplot
	printf("%f*x*x*x+(%f*x*x)+(%f*x)+(%f)\n",M[0][(N-1)*4],M[1][(N-1)*4],M[2][(N-1)*4],M[3][(N-1)*4]);//Para imprimir en pantalla
	a=1;
	for(i=4;i<(N-1)*4;i+=4){
		fprintf(Salida,"replot t*%f+(%f), %f*(t*%f+%f)*(t*%f+%f)*(t*%f+%f)+%f*(t*%f+%f)*(t*%f+%f)+%f*(t*%f+%f)+(%f)\n",
					delta[a],X[a],M[i][(N-1)*4],delta[a],X[a],delta[a],X[a],delta[a],X[a],M[i+1][(N-1)*4],delta[a],X[a],delta[a],X[a],
					M[i+2][(N-1)*4],delta[a],X[a],M[i+3][(N-1)*4]);//imprimir en archivo para graficar en gnuplot
		printf("%f*x*x*x+(%f*x*x)+(%f*x)+(%f)\n",M[i][(N-1)*4],M[i+1][(N-1)*4],M[i+2][(N-1)*4],M[i+3][(N-1)*4]);//Para imprimir en pantalla
		a++;
	}
	fprintf(Salida,"replot 'puntos.dat'\n");
	
	//printf("t*%f+(%f), %f *(t*%f+(%f)) + %f\n",delta[a],X[a],);
	fclose(Salida);
	Salida =fopen("puntos.dat","w");
	for(i=0;i<N;i++){
		fprintf(Salida,"%f\t%f\n",X[i],Y[i]);
	}
	fclose(Salida);
}
