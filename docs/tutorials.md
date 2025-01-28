# Integración Numérica con Cuadratura Gaussiana

Este tutorial explica cómo implementar y usar un método de integración numérica basado en cuadratura Gaussiana en Python, y proporciona un ejemplo práctico paso a paso.

---

## Objetivo del Ejemplo

Calcular la integral de la función:

\[
f(x) = x^6 - x^2 \sin(2x)
\]

en el intervalo \([1, 3]\), con una tolerancia de error de \(10^{-6}\), utilizando la cuadratura de Gauss-Legendre.

---

## Explicación general del código para realizar ediciones:

### Paso 1: Definir la Función a Integrar

Primero, se define la función que deseas integrar. En este caso:

```python
import numpy as np

def func_grado6(xvar):
    """
    Se define la función matemática x^6 - x^2 * sin(2x) a integrar.

    Argumentos:
        xvar (tipo float): Valores en los que se evalúa la función.

    Devuelve:
        Un valor tipo float que es el resultado de la evaluación de la función.

    Ejemplo:
        >>> x = np.array([1, 2, 3])
        >>> func_grado6(x)
        array([-0.90929743, 62.72789228, 717.97165513])
    """
    return xvar**6 - xvar**2 * np.sin(2 * xvar)


```

### Paso 2: Implementar el Cálculo de Puntos y Pesos

Se define la función para calcular los puntos y pesos de la cuadratura Gaussiana en el intervalo estándar \([-1, 1]\):

```python
def gaussxw(N):
    """
    Calcula los puntos y pesos de cuadratura Gauss-Legendre para el intervalo [-1, 1].

    Argumentos:
        N (tipo int): Número de puntos de cuadratura.

    Devuelve:
        Una tupla (tuple):
            - np.ndarray: Puntos en el intervalo estándar [-1, 1].
            - np.ndarray: Pesos en el intervalo estándar [-1, 1].

    Ejemplo:
        >>> puntos, pesos = gaussxw(3)
        >>> print(puntos)
        [ 0.77459667  0.          -0.77459667]
        >>> print(pesos)
        [0.55555556 0.88888889 0.55555556]
    """
    # Aproximación inicial
    a = np.linspace(3, 4 * (N - 1), N) / ((4 * N) + 2)
    x = np.cos(np.pi * a + 1 / (8 * N * N * np.tan(a)))

    # Creación de los puntos utilizando el método de Newton
    epsilon = 1e-15
    delta = 1.0
    while delta > epsilon:
        p0 = np.ones(N, dtype=float)
        p1 = np.copy(x)
        for k in range(1, N):
            p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
        dp = (N + 1) * (p0 - x * p1) / (1 - x * x)
        dx = p1 / dp
        x -= dx
        delta = np.max(np.abs(dx))

    # Cálculo de los pesos
    w = 2 * (N + 1) * (N + 1) / (N * N * (1 - x * x) * dp * dp)
    return x, w
```

### Paso 3: Escalar los puntos y pesos

Se escalan los puntos y pesos calculados al intervalo deseado \([a, b]\):

```python
def escalado(lim_inf, lim_sup, puntos, pesos):
    """
    Escala los puntos y pesos de [-1, 1] al intervalo [lim_inf, lim_sup].

    Argumentos:
        lim_inf (tipo float): Límite inferior de integración.
        lim_sup (tipo float): Límite superior de integración.
        puntos (np.ndarray): Puntos en el intervalo [-1, 1].
        pesos (np.ndarray): Pesos en el intervalo [-1, 1].

    Devuelve:
        Tupla (tuple):
            - np.ndarray: Puntos escalados al intervalo [lim_inf, lim_sup].
            - np.ndarray: Pesos escalados al intervalo [lim_inf, lim_sup].

    Ejemplo:
        >>> puntos = np.array([-0.77459667, 0.0, 0.77459667])
        >>> pesos = np.array([0.55555556, 0.88888889, 0.55555556])
        >>> escalado(1, 3, puntos, pesos)
        (array([1.11270167, 2.0, 2.88729833]), array([0.55555556, 0.88888889, 0.55555556]))
    """
    return 0.5 * (lim_sup - lim_inf) * puntos + 0.5 * (lim_sup + lim_inf), 0.5 * (lim_sup - lim_inf) * pesos
```

### Paso 4: Calcular la integral aproximada

Se define una función que use los puntos, pesos y la función dada para calcular la integral aproximada:

```python
def Aprox_integ(func, pesos, puntos):
    """
    Calcula la integral aproximada utilizando cuadratura Gaussiana.

    Argumentos:
        func (function): Función a integrar.
        pesos (np.ndarray): Pesos de cuadratura.
        puntos (np.ndarray): Puntos de colocación.

    Devuelve:
         Un float: Resultado de la integral aproximada.

    Ejemplo:
        >>> puntos = np.array([1.11270167, 2.0, 2.88729833])
        >>> pesos = np.array([0.55555556, 0.88888889, 0.55555556])
        >>> Aprox_integ(func_grado6, pesos, puntos)
        309.45365734
    """
    return np.sum(pesos * func(puntos))
```

## Tutorial de ejecución del algoritmo completo:

Configurar el intervalo, la tolerancia y el valor correcto a encontrar, además de llamar la función `encontrar_N_aproximado`:

```python
# Configuración inicial
inf, sup = 1, 3  # Límites inferior y superior de la integral
resultado_correcto = 317.3442467  # Resultado correcto de la integral (calculado numéricamente)

print("Buscando el valor de N necesario para alcanzar la tolerancia...\n")
N_optimo = encontrar_N_aproximado(func_grado6, inf, sup, resultado_correcto, 1e-6)
```

---

### Salida esperada

Al ejecutar el código, el programa iterará sobre diferentes valores de \(N\) y mostrará resultados similares a los siguientes:

```
Buscando el valor de N necesario para alcanzar la tolerancia...

N = 1, Resultado = 134.05441996, Error = 1.83289827e+02
N = 2, Resultado = 598.24047746, Error = 2.80896231e+02
N = 3, Resultado = 375.49436118, Error = 5.81501145e+01
N = 4, Resultado = 317.34539033, Error = 1.14363416e-03
N = 5, Resultado = 448.11612481, Error = 1.30771878e+02
N = 6, Resultado = 414.07734548, Error = 9.67330988e+01
N = 7, Resultado = 317.34424667, Error = 2.77739218e-08

Con N = 7, se alcanza el resultado correcto: 317.34424667
```

---




