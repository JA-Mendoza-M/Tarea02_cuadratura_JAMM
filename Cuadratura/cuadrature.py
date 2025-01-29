import numpy as np

# Definición de la función a integrar
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

# Función que calcula los puntos y pesos de cuadratura de Gauss-Legendre
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

# Función para escalar puntos y pesos al intervalo [lim_inf, lim_sup]
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

# Función para aproximar la integral
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

# Función para encontra la cantidad de N necesarias para aproximar el valor.
def encontrar_N_aproximado(func, lim_inf,lim_sup,resultado_correcto, tolerancia):
    """
    Encuentra el valor de N necesario para aproximar la integral con una tolerancia específica.
    
    Argumentos:
        func (function): Función a integrar.
        lim_inf (float): Límite inferior de integración.
        lim_sup (float): Límite superior de integración.
        resultado_correcto (float): Valor correcto de la integral para comparar.
        tolerancia (float): Tolerancia permitida para el error.
    
    Devuelve:
        int: Número de puntos (N) necesario para alcanzar la tolerancia.
    
    Ejemplo:
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
        error = abs(resultado-resultado_correcto)/resultado_correcto
        print(f"N = {N}, Resultado = {resultado:.8f}, Error = {error:.8e}")
        if error <= tolerancia:
            print(f"\nCon N = {N}, se alcanza el resultado correcto: {resultado:.8f}")
            return N
        N += 1
        
# Configuración inicial
inf, sup = 1, 3  # Límites inferior y superior de la integral
resultado_correcto = 317.3442467  # Resultado correcto de la integral (calculado numéricamente)

print("Buscando el valor de N necesario para alcanzar la tolerancia...\n")
N_optimo = encontrar_N_aproximado(func_grado6, inf, sup, resultado_correcto, 1e-6)
