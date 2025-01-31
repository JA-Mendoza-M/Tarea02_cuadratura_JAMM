# Descripción del Método de Cuadratura Gaussiana

El método de cuadratura de Gauss-Legendre es una técnica avanzada de integración numérica que utiliza un conjunto de puntos y pesos para aproximar el valor de una integral definida de una función. A diferencia de otros métodos de integración, como los de Riemann o trapezoidal, la cuadratura de Gauss-Legendre no utiliza particiones uniformes del intervalo de integración, lo que lo hace más eficiente y preciso.

## Fórmula de Cuadratura de Gauss-Legendre

La fórmula de cuadratura de Gauss-Legendre para un intervalo $[a, b]$ es la siguiente:

$$
I = \sum_{i=1}^{N} w_i f(x_i)
$$


- Donde $x_i$ son los puntos de integración (ceros de los polinomios de Legendre).
- Los pesos corresppondientes a los puntos se representan con $w_i$ .
- Además $f(x_i)$ es el valor de la función que estamos integrando en los puntos de integración.

## Diferencias con los Métodos de Newton-Cotes

### Métodos de Newton-Cotes:
- Los puntos de muestreo son **equidistantes**.
- Una ecuación de Newton-Cotes de orden $N$ es *exacta* para un polinomio de grado $N$.
- Un polinomio de orden $N$ aproxima una función bien comportada mejor que un polinomio de orden $N-1$.

### Cuadratura Gaussiana:
- Los puntos de muestreo **no son equidistantes**, lo que introduce más grados de libertad.
- Es exacta para un polinomio de orden $(2N - 1)$.
- Utiliza los ceros de los polinomios de Legendre $P_N(x)$ como puntos de muestreo:
  - Los pesos se calculan como:
    $$
    w_k = \left[\frac{2}{1 - x_k^2} \left(\frac{dP_N}{dx}\right)^{-2}\right]_{x=x_k}, \quad \text{donde } P_N(x_k) = 0.
    $$

## Ventajas del Método

- Alta precisión con un número relativamente bajo de puntos de integración.
- Convergencia rápida: el error decrece como $cte/N^2$ al aumentar los puntos de muestreo.
- Es ideal para funciones bien comportadas en el intervalo de integración.

## Limitaciones del Método

- No es adecuado para funciones con singularidades o comportamientos complicados; se requieren más puntos de muestreo en esas regiones.
- Evaluar el error de forma precisa puede ser complejo.

## Polinomios de Legendre

Los polinomios de Legendre son un sistema de polinomios ortogonales definidos en el intervalo $[-1, 1]$. Satisfacen la relación de ortogonalidad:

$$
\int_{-1}^1 P_N(x) P_M(x) \, dx = \frac{2\delta_{MN}}{2N+1}.
$$

Se definen de forma recursiva:

- $P_0(x) = 1$, $P_1(x) = x$

- La relación de recurrencia es:
  $$
  (N+1)P_{N+1}(x) = (2N+1)xP_N(x) - NP_{N-1}(x).
  $$

Alternativamente, pueden definirse mediante la fórmula de Rodrigues:
$$
P_N(x) = \frac{1}{2^N N!} \frac{d^N}{dx^N}\left[(x^2 - 1)^N\right].
$$

## Uso de la Cuadratura Gaussiana

Una vez que se tienen los puntos de muestreo $x_k$ y los pesos $w_k$, la integral se evalúa como:

$$
\int_a^b f(x) \, dx \approx \sum_{k=1}^{N+1} w_k f(x_k).
$$

