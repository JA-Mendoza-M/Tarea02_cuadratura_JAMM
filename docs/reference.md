# Reference Documentation

Este documento contiene la referencia de las funciones definidas en el módulo, con descripciones detalladas, parámetros, valores de retorno y ejemplos de uso.

## Funciones

### `func_grado6(xvar)`

Define la función matemática a integrar.

#### Parámetros:
- `xvar` (float o array): Valores en los que se evalúa la función.

#### Retorno:
- `float`: Resultado de la evaluación de la función.

#### Ejemplo de uso:
```python
import numpy as np
from modulo import func_grado6

# Evaluar la función en un punto
resultado = func_grado6(2.0)
print(resultado)  # Output: 62.234248557783406

# Evaluar la función en un array
x = np.array([1.0, 2.0, 3.0])
resultados = func_grado6(x)
print(resultados)  # Output: [ 0.          62.23424856 729.40019833 ]
```

---

### `gaussxw(N)`

Calcula los puntos y pesos de cuadratura Gauss-Legendre en el intervalo estándar [-1, 1].

#### Parámetros:
- `N` (int): Número de puntos de cuadratura.

#### Retorno:
- `tuple`: Puntos (`array`) y pesos (`array`) en el intervalo estándar [-1, 1].

#### Ejemplo de uso:
```python
from modulo import gaussxw

# Calcular puntos y pesos para N = 4
N = 4
puntos, pesos = gaussxw(N)
print("Puntos:", puntos)
print("Pesos:", pesos)
```

---

### `escalado(lim_inf, lim_sup, puntos, pesos)`

Escala los puntos y pesos de [-1, 1] al intervalo [lim_inf, lim_sup].

#### Parámetros:
- `lim_inf` (float): Límite inferior de integración.
- `lim_sup` (float): Límite superior de integración.
- `puntos` (array): Puntos en el intervalo [-1, 1].
- `pesos` (array): Pesos en el intervalo [-1, 1].

#### Retorno:
- `tuple`: Puntos (`array`) y pesos (`array`) escalados al intervalo [lim_inf, lim_sup].

#### Ejemplo de uso:
```python
from modulo import escalado, gaussxw

# Calcular puntos y pesos en [-1, 1]
puntos, pesos = gaussxw(4)

# Escalar al intervalo [0, 2]
puntos_esc, pesos_esc = escalado(0, 2, puntos, pesos)
print("Puntos escalados:", puntos_esc)
print("Pesos escalados:", pesos_esc)
```

---

### `Aprox_integ(func, pesos, puntos)`

Calcula la integral aproximada utilizando cuadratura Gaussiana.

#### Parámetros:
- `func` (function): Función a integrar.
- `pesos` (array): Pesos de cuadratura.
- `puntos` (array): Puntos de colocación.

#### Retorno:
- `float`: Resultado de la integral aproximada.

#### Ejemplo de uso:
```python
import numpy as np
from modulo import func_grado6, Aprox_integ, gaussxw, escalado

# Configuración
lim_inf, lim_sup = 1, 3
N = 4

# Calcular puntos y pesos en [-1, 1]
puntos, pesos = gaussxw(N)

# Escalar al intervalo [lim_inf, lim_sup]
puntos_esc, pesos_esc = escalado(lim_inf, lim_sup, puntos, pesos)

# Aproximar la integral
resultado = Aprox_integ(func_grado6, pesos_esc, puntos_esc)
print("Resultado de la integral:", resultado)
```

---

### Código de Ejemplo Completo

A continuación, un ejemplo completo que usa todas las funciones descritas para aproximar la integral de la función `func_grado6` en el intervalo [1, 3]:

```python
import numpy as np
from modulo import func_grado6, gaussxw, escalado, Aprox_integ

# Configuración inicial
a, b = 1, 3
resultado_correcto = 317.3442467  # Valor calculado numéricamente
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

Este archivo proporciona una referencia exhaustiva y fácil de usar para las funciones del módulo. Cada función está completamente documentada y respaldada por ejemplos prácticos para facilitar su comprensión e implementación.


