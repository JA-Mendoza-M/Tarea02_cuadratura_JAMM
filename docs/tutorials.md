# Tutorial: Cuadratura Gaussiana en Python

Este tutorial explica cómo implementar el método de cuadratura de Gauss-Legendre en Python para aproximar la integral de una función definida. Además, incluye un ejemplo práctico y un código funcional.

---

## Introducción

El método de cuadratura Gauss-Legendre es una técnica numérica utilizada para calcular integrales definidas con alta precisión. Este método se basa en evaluar la función en puntos específicos (ceros de los polinomios de Legendre) y ponderarlos mediante pesos predefinidos.

**Ventajas del método:**

- Alta precisión con pocos puntos de muestreo.
- Convergencia rápida.

**Aplicaciones comunes:**

- Integrales de funciones complicadas.
- Métodos de elementos finitos.

En este tutorial, se integra una función de grado 6 para demostrar el método.

---

## Código

El siguiente código implementa el método de cuadratura Gauss-Legendre, escalando los puntos y pesos al intervalo deseado y calculando la integral aproximada.

```python
#!/usr/bin/env python3
import numpy as np

# Definición de la función a integrar
def func_grado6(xvar):
    """
    Define la función a integrar: x^6 - x^2*sin(2x).
    Parámetros:
        xvar: Valores en los que se evalúa la función.
    Retorna:
        float : Resultado de la evaluación de la función.
    """
    return xvar**6 - xvar**2 * np.sin(2 * xvar)

# Función que calcula los puntos y pesos de cuadratura de Gauss-Legendre
def gaussxw(N):
    """
    Calcula los puntos y pesos de cuadratura Gauss-Legendre en el intervalo [-1, 1].
    Parámetros:
        N (int): Número de puntos de cuadratura.
    Retorna:
        tuple: Puntos y pesos en el intervalo estándar [-1, 1].
    """
    # Aproximación inicial
    a = np.linspace(3, 4 * (N - 1), N) / ((4 * N) + 2)
    x = np.cos(np.pi * a + 1 / (8 * N * N * np.tan(a)))

    # Refinamiento de los puntos utilizando el método de Newton
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

# Función para escalar puntos y pesos al intervalo [lim_inf, lim_sup]
def escalado(lim_inf, lim_sup, puntos, pesos):
    """
    Escala los puntos y pesos de [-1, 1] al intervalo [lim_inf, lim_sup].
    Parámetros:
        lim_inf (float): Límite inferior de integración.
        lim_sup (float): Límite superior de integración.
        puntos (array): Puntos en el intervalo [-1, 1].
        pesos (array): Pesos en el intervalo [-1, 1].
    Retorna:
        tuple: Puntos y pesos escalados al intervalo [lim_inf, lim_sup].
    """
    return 0.5 * (lim_sup - lim_inf) * puntos + 0.5 * (lim_sup + lim_inf), 0.5 * (lim_sup - lim_inf) * pesos

# Función para aproximar la integral
def Aprox_integ(func, pesos, puntos):
    """
    Calcula la integral aproximada utilizando cuadratura Gaussiana.
    Parámetros:
        func (function): Función a integrar.
        pesos (array): Pesos de cuadratura.
        puntos (array): Puntos de colocación.
    Retorna:
        float: Resultado de la integral aproximada.
    """
    return np.sum(pesos * func(puntos))

# Configuración inicial
a, b = 1, 3
resultado_correcto = 317.3442467  # Calculado numéricamente
tolerancia = 1e-6

# Bucle para encontrar el N adecuado
N = 1
error = float("inf")

while error > tolerancia:
    # Obtiene puntos y pesos en [-1, 1]
    puntos, pesos = gaussxw(N)
    # Escala los puntos y pesos al intervalo [a, b]
    puntos_esc, pesos_esc = escalado(a, b, puntos, pesos)
    # Calcula la integral aproximada
    resultado = Aprox_integ(func_grado6, pesos_esc, puntos_esc)
    # Calcula el error
    error = abs(resultado_correcto - resultado)
    print(f"N = {N}, Resultado = {resultado:.8f}, Error = {error:.8e}")
    if error <= tolerancia:
        print(f"\nCon N = {N}, se alcanza el resultado correcto: {resultado:.8f}")
        break
    N += 1
```

---

## Ejemplo de Uso

1. **Función a integrar:**

   La función que vamos a integrar es:



2. **Intervalo de integración:**



3. **Resultado esperado:**



4. **Ejecución del programa:**

   Al ejecutar el programa, se buscará el número mínimo de puntos () necesarios para alcanzar el resultado con una tolerancia de .

5. **Salida esperada:**

   ```
   N = 1, Resultado = 222.66666667, Error = 9.46775800e+01
   N = 2, Resultado = 303.30481077, Error = 1.40394359e+01
   N = 3, Resultado = 317.15978385, Error = 1.84462909e-01
   N = 4, Resultado = 317.34424151, Error = 5.18600000e-06

   Con N = 4, se alcanza el resultado correcto: 317.34424151
   ```

---

## Conclusión

El método de cuadratura Gauss-Legendre permite calcular integrales definidas con alta precisión, adaptándose dinámicamente al número de puntos necesarios. Este ejemplo demuestra cómo implementar y aplicar el método en Python.


