# Introducción al Método de Cuadratura Gaussiana

El propósito es implementar el **método de cuadratura de Gauss-Legendre** para aproximar integrales definidas de funciones. Este método de integración numérica evalúa la integral utilizando puntos estratégicamente seleccionados, conocidos como puntos de Gauss, y unos coeficientes ponderados, que se llaman pesos. 

El método es altamente eficiente al evaluar integrales, incluso para funciones complejas que no poseen una solución analítica sencilla. Además, logra esta precisión con un número relativamente pequeño de iteraciones (\(N\)). En esta tarea, se implementa el método de cuadratura de Gauss-Legendre en el módulo `cuadrature.py`.

### Propósito del código

- **Objetivo**: Aproximar la integral de funciones complejas mediante el método de cuadratura de Gauss-Legendre.
- **Método**: 
  - Implementación de funciones que calculan los puntos y pesos de Gauss-Legendre para integrales definidas en un intervalo \([a, b]\).
  - Transformación de los puntos y pesos desde el intervalo estándar \([-1, 1]\) al intervalo de integración deseado mediante una función de escalado.
  - Uso de recursividad para determinar el número mínimo de iteraciones (\(N\)) necesarias para alcanzar un resultado fiable, respecto a una tolerancia epsilon, así como una comparación con un valor de referencia o solución exacta.
