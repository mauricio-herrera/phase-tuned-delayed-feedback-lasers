# Resultados esperados y criterios de aceptación

Este documento dice qué números debería obtener, cuánto tarda cada corrida, cómo decidir si
un resultado es aceptable y dónde va cada cifra en el manuscrito.

Todos los valores de las secciones 2 y 3 provienen de una corrida piloto con los presets por
defecto en 1 núcleo, numba 0.68 y Δt = 0.02. Con las mismas semillas, su máquina debería
reproducirlos (posiblemente con diferencias en el último dígito). Si aumenta `--batches` o
`--T`, los valores centrales cambiarán dentro de la incertidumbre y los errores se reducirán.

## 1. Valores objetivo (teoría lineal, exactos)

Son las razones respecto al láser solitario. La simulación debe acercarse a ellas a medida
que crece la estadística.

| preset | caso | κ | θ/π | ancho de línea | RIN | portadores |
|---|---|---|---|---|---|---|
| tau20 | trade-off | 0.008 | −0.15584 | 0.471 | 1.565 | 1.438 |
| tau20 | win-win | 0.009 | 0.01795 | 0.750 | 0.480 | 0.518 |
| long200 | trade-off | 0.005 | −0.15584 | 0.070 | 1.369 | 1.352 |
| long200 | en fase | 0.005 | 0 | 0.246 | 0.761 | 0.764 |
| long1000 | trade-off | 0.001 | −0.15584 | 0.070 | 1.069 | 1.069 |
| long1000 | en fase | 0.001 | 0 | 0.249 | 0.935 | 0.934 |

## 2. Resultados de la corrida piloto

Formato: valor simulado ± 2 errores estándar (predicción lineal). Las dos últimas columnas
comparan la varianza absoluta simulada con la teoría lineal, sin dividir por el solitario.

**tau20**: τ = 20, T = 40000, burn-in = 400, 128 realizaciones, ventana de ajuste de fase 7.5τ–40τ.

| caso | κ | θ/π | ancho de línea | RIN | portadores | ancho de línea abs. sim/lin | RIN abs. sim/lin |
|---|---|---|---|---|---|---|---|
| solitary | 0 | 0.0000 | 1.000 ± 0.043 (1.000) | 1.000 ± 0.016 (1.000) | 1.000 ± 0.016 (1.000) | 0.962 ± 0.029 | 1.015 ± 0.011 |
| tradeoff | 0.008 | -0.1558 | 0.474 ± 0.021 (0.471) | 1.568 ± 0.029 (1.565) | 1.442 ± 0.027 (1.438) | 0.968 ± 0.031 | 1.017 ± 0.015 |
| win_win | 0.009 | 0.0180 | 0.791 ± 0.034 (0.750) | 0.477 ± 0.006 (0.480) | 0.515 ± 0.007 (0.518) | 1.014 ± 0.031 | 1.008 ± 0.007 |

**long200**: τ = 200, T = 100000, burn-in = 2000, 128 realizaciones, ventana de ajuste de fase 3τ–12τ.

| caso | κ | θ/π | ancho de línea | RIN | portadores | ancho de línea abs. sim/lin | RIN abs. sim/lin |
|---|---|---|---|---|---|---|---|
| solitary | 0 | 0.0000 | 1.000 ± 0.055 (1.000) | 1.000 ± 0.009 (1.000) | 1.000 ± 0.009 (1.000) | 1.003 ± 0.039 | 1.016 ± 0.006 |
| tradeoff | 0.005 | -0.1558 | 0.070 ± 0.004 (0.070) | 1.371 ± 0.013 (1.369) | 1.353 ± 0.013 (1.352) | 1.010 ± 0.042 | 1.017 ± 0.008 |
| win_win | 0.005 | 0.0000 | 0.238 ± 0.014 (0.246) | 0.756 ± 0.006 (0.761) | 0.760 ± 0.007 (0.764) | 0.969 ± 0.045 | 1.009 ± 0.006 |

**long1000**: τ = 1000, T = 400000, burn-in = 10000, 128 realizaciones, ventana de ajuste de fase 3τ–12τ.

| caso | κ | θ/π | ancho de línea | RIN | portadores | ancho de línea abs. sim/lin | RIN abs. sim/lin |
|---|---|---|---|---|---|---|---|
| solitary | 0 | 0.0000 | 1.000 ± 0.065 (1.000) | 1.000 ± 0.004 (1.000) | 1.000 ± 0.004 (1.000) | 1.031 ± 0.047 | 1.014 ± 0.003 |
| tradeoff | 0.001 | -0.1558 | 0.066 ± 0.004 (0.070) | 1.069 ± 0.005 (1.069) | 1.068 ± 0.005 (1.069) | 0.972 ± 0.046 | 1.014 ± 0.003 |
| win_win | 0.001 | 0.0000 | 0.248 ± 0.017 (0.249) | 0.933 ± 0.004 (0.935) | 0.933 ± 0.004 (0.934) | 1.024 ± 0.051 | 1.012 ± 0.003 |

**Control de paso temporal** (`long1000 --dt 0.01 --tag dt001`, solo el caso en fase; las
razones no se calculan porque no se simuló el solitario):

**long1000_dt001**: τ = 1000, T = 400000, burn-in = 10000, 64 realizaciones, ventana de ajuste de fase 3τ–12τ.

| caso | κ | θ/π | ancho de línea | RIN | portadores | ancho de línea abs. sim/lin | RIN abs. sim/lin |
|---|---|---|---|---|---|---|---|
| win_win | 0.001 | 0.0000 | — | — | — | 1.041 ± 0.067 | 1.005 ± 0.005 |

**Registro cuatro veces más largo** (`long1000 --T 1.6e6 --tag T16`, solo el caso en fase; las razones no se calculan porque no se simuló el solitario):

**long1000_T16**: τ = 1000, T = 1.6e+06, burn-in = 10000, 128 realizaciones, ventana de ajuste de fase 3τ–12τ.

| caso | κ | θ/π | ancho de línea | RIN | portadores | ancho de línea abs. sim/lin | RIN abs. sim/lin |
|---|---|---|---|---|---|---|---|
| win_win | 0.001 | 0.0000 | — | — | — | 1.029 ± 0.026 | 1.012 ± 0.002 |

**Sensibilidad a la ventana de ajuste de la fase:**

| preset | ventana (τ) | caso | razón de ancho de línea | 2SE | abs. sim/lin |
|---|---|---|---|---|---|
| long200 | 2-8 | solitary | 1.000 | 0.040 | 0.996 |
| long200 | 2-8 | tradeoff | 0.070 | 0.003 | 0.999 |
| long200 | 2-8 | win_win | 0.244 | 0.011 | 0.987 |
| long200 | 3-12 | solitary | 1.000 | 0.055 | 1.003 |
| long200 | 3-12 | tradeoff | 0.070 | 0.004 | 1.010 |
| long200 | 3-12 | win_win | 0.238 | 0.014 | 0.969 |
| long200 | 5-20 | solitary | 1.000 | 0.071 | 1.006 |
| long200 | 5-20 | tradeoff | 0.070 | 0.005 | 1.015 |
| long200 | 5-20 | win_win | 0.235 | 0.018 | 0.958 |
| long1000 | 2-8 | solitary | 1.000 | 0.053 | 1.022 |
| long1000 | 2-8 | tradeoff | 0.068 | 0.004 | 0.981 |
| long1000 | 2-8 | win_win | 0.247 | 0.013 | 1.012 |
| long1000 | 3-12 | solitary | 1.000 | 0.065 | 1.031 |
| long1000 | 3-12 | tradeoff | 0.066 | 0.004 | 0.972 |
| long1000 | 3-12 | win_win | 0.248 | 0.017 | 1.024 |
| long1000 | 5-20 | solitary | 1.000 | 0.085 | 1.039 |
| long1000 | 5-20 | tradeoff | 0.067 | 0.006 | 0.993 |
| long1000 | 5-20 | win_win | 0.245 | 0.023 | 1.021 |

## 3. Observaciones de la corrida piloto

1. **Las razones concuerdan con la teoría lineal.** Para long200 y long1000, todas las
   razones (ancho de línea, RIN y portadores) están dentro de 2 errores estándar (|z| ≤ 2).
   En tau20, el caso win-win da z ≈ +2.4 en la razón de ancho de línea. Esto se debe a que el
   ancho de línea absoluto del solitario simulado salió 3.8 % bajo (z ≈ −2.6), mientras que
   el win-win, en absoluto, coincide (1.014 ± 0.031). Es una fluctuación del denominador y no
   afecta lo que se inserta, que es solo long200 y long1000.
2. **Sesgo de Euler en las varianzas absolutas.** Las varianzas absolutas de RIN y de
   portadores salen alrededor de 1.5 % por encima de la teoría en todos los casos, incluido el
   solitario, con Δt = 0.02. Con Δt = 0.01 el exceso se reduce aproximadamente a la mitad.
   Es el sesgo de primer orden del método de Euler y se cancela en las razones.
3. **Deriva residual de frecuencia (corregida).** En la primera versión del análisis, el
   ancho de línea absoluto del modo en fase con τ = 1000 salía alrededor de 9 % alto y crecía
   con la ventana de ajuste. La causa era un pequeño desplazamiento de la frecuencia media
   respecto al marco rotante del ECM (≈ 2.5×10⁻⁶ en unidades de 1/τ_p), común a todas las
   realizaciones. Ese desplazamiento añade un término (δω·t)² al desplazamiento cuadrático
   medio de la fase.

   `aggregate.py` ahora estima δω a partir del incremento medio del ensamble y lo resta. No
   lo hace por realización, para no restar también parte de la difusión de fase. La columna
   `mean_freq_offset` de `summary_*.csv` informa δω. Con esta corrección, el ancho de línea
   absoluto pasa de 1.089 ± 0.051 a 1.024 ± 0.051 (1.041 ± 0.067 con Δt = 0.01) y concuerda
   con la teoría. El texto del suplemento que inserta `update_manuscript.py` describe el
   procedimiento. Si quiere ver el efecto, use
   `python aggregate.py --no-detrend`.

4. **Comprobado: la difusión de fase es estable a lags largos.** Con T = 4×10⁵, el modo en
   fase de τ = 1000 parecía aumentar su tasa de difusión para lags mayores que unos 20τ
   (1.17 ± 0.20 en 30–60τ y 1.30 ± 0.34 en 60–100τ). Para descartar un efecto físico, repetí
   ese caso con registros cuatro veces más largos (`--T 1.6e6 --tag T16`, 128 realizaciones,
   unos 23 min en 1 núcleo). La tasa de difusión de fase, en valor absoluto respecto a la
   teoría lineal, queda en:

   | ventana | 3–12τ | 12–30τ | 30–60τ | 60–100τ | 100–200τ |
   |---|---|---|---|---|---|
   | T = 4×10⁵ | 1.02 ± 0.05 | 1.01 ± 0.11 | 1.17 ± 0.20 | 1.30 ± 0.34 | — |
   | T = 1.6×10⁶ | 1.03 ± 0.03 | 1.03 ± 0.05 | 1.01 ± 0.08 | 1.00 ± 0.13 | 0.96 ± 0.22 |

   El aumento era fluctuación estadística. No hace falta añadir nada al suplemento. Si quiere
   ver el resultado con más precisión en su máquina, use `--T 1.6e6 --tag T16`; el resumen
   está en `lag_windows_long1000_T16.csv`.

## 4. Criterios de aceptación antes de insertar resultados

| comprobación | dónde mirarla | aceptar si |
|---|---|---|
| razones frente a teoría lineal | columnas z de `aggregate.py` | \|z\| ≤ 2 en casi todos los casos y ninguno \|z\| > 3 (`update_manuscript.py` avisa) |
| sin saltos de modo | columna «deriva<I>» y `max_meanI_drift` | < 10⁻³ (típico: 10⁻⁵) |
| ventana de fase | `lag_sensitivity.csv` | las razones con 2–8τ, 3–12τ y 5–20τ difieren menos que su 2SE |
| difusión de fase bien definida | `diag_<preset>.png` y `lag_windows_<preset>.csv` | fluctúa alrededor de 1 dentro de la franja gris y en las ventanas sucesivas |
| paso temporal | variante `--dt 0.01` | razones iguales dentro de 2SE; absolutos más cerca de 1 |
| varianzas absolutas | columnas «abs. sim/lin» | RIN y portadores entre 1.00 y 1.02 con Δt = 0.02 (sesgo de Euler) |

## 5. Precisión frente a esfuerzo

El error estándar de las varianzas de RIN y de portadores escala como 1/√(N·T). El de la
difusión de fase escala aproximadamente como √(L/(N·T)), donde L es la ventana de lags. Con
los presets por defecto, las incertidumbres (2 errores estándar) de las razones son:

| preset | ancho de línea | RIN |
|---|---|---|
| long200 | ≈ ±6 % | ≈ ±1 % |
| long1000 | ≈ ±6–7 % | ≈ ±0.5 % |

Para llevar el ancho de línea a alrededor de ±3 %, multiplique N·T por 4. Dos formas de
hacerlo:
- **Más realizaciones:** `--batches 32`, que añade lotes a los ya existentes.
- **Registros más largos:** `--T` cuatro veces mayor, como variante `--tag` en una carpeta
  nueva. Use esta opción si prefiere menos lotes largos.

**Tiempos medidos** en 1 núcleo, para 16 realizaciones por lote y Δt = 0.02. Numba reparte
las 16 realizaciones de cada lote entre los núcleos, así que con hasta 16 núcleos el tiempo se
divide aproximadamente por el número de núcleos. Con más núcleos, use una variante con
`--batch-size 32` o `64` (con `--tag`, porque no se pueden mezclar tamaños de lote).

| preset | por lote y caso | total (8 lotes × casos) |
|---|---|---|
| smoke | < 1 s | < 1 min |
| tau20 | 4 s | ≈ 2 min |
| long200 | 10 s | ≈ 4 min |
| long1000 | 37–41 s | ≈ 16 min |
| long1000 con Δt = 0.01 | 82 s | — |
| fixedkt (opcional) | ≈ 20 s | ≈ 11 min |

## 6. Dónde va cada número

| resultado | archivo | ubicación | cómo se inserta |
|---|---|---|---|
| razones long200 y long1000, ±2SE, frente a lineal | `main.tex` | Sec. 7, segundo párrafo, frases que empiezan «Nonlinear stochastic simulations for τ=200…» | `update_manuscript.py` |
| N, T, burn-in, ventana de ajuste, corrección de deriva | `supplement.tex` | Sec. S10, párrafo «Nonlinear simulations» | `update_manuscript.py` |
| filas de la Tabla S14 | `supplement.tex` | Tabla S14 (`tab:sdelong`) | `update_manuscript.py` (desde `tableS14_rows.tex`) |
| ítem «Nonlinear validation» | carta | Reviewer 3 | `update_manuscript.py` |
| tau20 (control) | — | no se inserta; puede citarse en la carta como verificación adicional de la Tabla S5 | manual, opcional |

Después de insertar, compruebe además:
- La carta, en la respuesta al Revisor 3, dice «Supplement S10 (methods, Table S13, Fig. S2,
  and simulation Table S14)». Esto no cambia.
- Si cambió T o N, la frase de la carta menciona el T del preset long1000 automáticamente.
