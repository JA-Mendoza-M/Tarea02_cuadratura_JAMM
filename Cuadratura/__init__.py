"""
Este paquete proporciona el método de integración numérica utilizando la cuadratura de Gauss-Legendre.

Módulos exportables de este paquete:

- `func_grado6`: Define la función a integrar: x^6 - x^2 * sin(2x).
- `gaussxw`: Calcula los puntos y los pesos de la cuadratura de Gauss-Legendre para un intervalo dado.
- `escalado`: Escala el intervalo [a, b] para la integración numérica.
- `Aprox_integ`: Aproxima la integral sumando áreas bajo la curva utilizando el método de cuadratura de Gauss-Legendre.
"""

