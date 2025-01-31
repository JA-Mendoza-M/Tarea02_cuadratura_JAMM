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

def func_grado6(xvar: np.ndarray) -> np.ndarray:
    """Calcula la función matemática x^6 - x^2 * sin(2x).

    Parameters:
        xvar (np.ndarray): Valores en los que se evalúa la función.

    Returns:
        np.ndarray: Resultado de la evaluación de la función.

    Example:
        >>> x = np.array([1, 2, 3])
        >>> func_grado6(x)
        array([-0.90929743, 62.72789228, 717.97165513])
        
    """
    return xvar**6 - xvar**2 * np.sin(2 * xvar)

```

### Paso 2: Implementar el Cálculo de Puntos y Pesos

Se define la función para calcular los puntos y pesos de la cuadratura Gaussiana en el intervalo estándar \([-1, 1]\):

```python
def gaussxw(N: int) -> tuple[np.ndarray, np.ndarray]:
    """Calcula los puntos y pesos de cuadratura Gauss-Legendre en el intervalo [-1, 1].

    Parameters:
        N (int): Número de puntos de cuadratura.

    Returns:
        tuple[np.ndarray, np.ndarray]: 
            - Puntos en el intervalo [-1, 1].
            - Pesos en el intervalo [-1, 1].

    Example:
        >>> puntos, pesos = gaussxw(3)
        >>> puntos
        array([ 0.77459667,  0., -0.77459667])
        >>> pesos
        array([0.55555556, 0.88888889, 0.55555556])
        
    """
    a = np.linspace(3, 4 * (N - 1), N) / ((4 * N) + 2)
    x = np.cos(np.pi * a + 1 / (8 * N * N * np.tan(a)))

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

    w = 2 * (N + 1) ** 2 / (N**2 * (1 - x**2) * dp**2)
    return x, w
```

### Paso 3: Escalar los puntos y pesos

Se escalan los puntos y pesos calculados al intervalo deseado \([a, b]\):

```python
def escalado(lim_inf: float, lim_sup: float, puntos: np.ndarray, pesos: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Escala los puntos y pesos de [-1, 1] al intervalo [lim_inf, lim_sup].

    Parameters:
        lim_inf (float): Límite inferior de integración.
        lim_sup (float): Límite superior de integración.
        puntos (np.ndarray): Puntos en el intervalo [-1, 1].
        pesos (np.ndarray): Pesos en el intervalo [-1, 1].

    Returns:
        tuple[np.ndarray, np.ndarray]: 
            - Puntos escalados al intervalo [lim_inf, lim_sup].
            - Pesos escalados al intervalo [lim_inf, lim_sup].

    Example:
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
def Aprox_integ(func: callable, pesos: np.ndarray, puntos: np.ndarray) -> float:
    """Calcula la integral aproximada utilizando cuadratura Gaussiana.

    Parameters:
        func (callable): Función a integrar.
        pesos (np.ndarray): Pesos de cuadratura.
        puntos (np.ndarray): Puntos de colocación.

    Returns:
        float: Resultado de la integral aproximada.

    Example:
        >>> puntos = np.array([1.11270167, 2.0, 2.88729833])
        >>> pesos = np.array([0.55555556, 0.88888889, 0.55555556])
        >>> Aprox_integ(func_grado6, pesos, puntos)
        309.45365734
        
    """
    return np.sum(pesos * func(puntos))
```

### Paso 5: Encontrar en N necesario

Se define una función que itera hasta encontrar un valor cuyo error sea igual o menor al valor de tolerancia:

```python
def encontrar_N_aproximado(func: callable, lim_inf: float, lim_sup: float, resultado_correcto: float, tolerancia: float) -> int:
    """Encuentra el valor de N necesario para aproximar la integral con una tolerancia específica.
    
    Parameters:
        func (callable): Función a integrar.
        lim_inf (float): Límite inferior de integración.
        lim_sup (float): Límite superior de integración.
        resultado_correcto (float): Valor de referencia para la integral.
        tolerancia (float): Tolerancia permitida para el error relativo.

    Returns:
        int: Número de puntos (N) necesario para alcanzar la tolerancia.

    Example:
        >>> func = lambda x: x**2
        >>> a, b = 0, 1
        >>> resultado_correcto = 1/3
        >>> tolerancia = 1e-6
        >>> encontrar_N_aproximado(func, a, b, resultado_correcto, tolerancia)
        3
        
    """
    N = 1
    error = float("inf")
    while error > tolerancia:
        puntos, pesos = gaussxw(N)
        puntos_esc, pesos_esc = escalado(lim_inf, lim_sup, puntos, pesos)
        resultado = Aprox_integ(func, pesos_esc, puntos_esc)
        error = abs(resultado - resultado_correcto) / resultado_correcto
        print(f"N = {N}, Resultado = {resultado:.8f}, Error = {error:.8e}")
        if error <= tolerancia:
            print(f"\nCon N = {N}, se alcanza el resultado correcto: {resultado:.8f}")
            return N
        N += 1

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




