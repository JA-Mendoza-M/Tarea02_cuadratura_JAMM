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



