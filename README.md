# interpolacion_splines
Programa para generar polinomios, por segmentación o por splines cuadrática y cúbica
listo para graficar en gnuplot
Pasos
      tener instalado gnuplot en linux
           1. gcc -o cubica SplinesCubica.c --compilando el programa
           2. gnuplot  --abriendo gnuplot
                    3. load "PolSegCubico" --graficando archivo generado
son los mismos pasos para ambos programas, cambiando nombre de archivo generado por PolSegCuad
y SplinesCubica.c por SplinesCuadratica.c
