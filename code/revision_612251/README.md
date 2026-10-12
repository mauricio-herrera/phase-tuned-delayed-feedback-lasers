# Simulaciones estocásticas largas — revisión JOSA B 612251

Este paquete repite, con registros mucho más largos, las simulaciones no lineales del modelo
de Lang–Kobayashi estocástico que respaldan la Sec. 7 del manuscrito revisado, la Tabla S14
del suplemento y la respuesta al Revisor 3. Al final, inserta los números nuevos en los tres
documentos sin tocar los originales.

## 1. Instalación

```bash
python -m venv venv && source venv/bin/activate      # opcional
python -m pip install --only-binary=:all: -r requirements.txt
```

**Compatibilidad macOS/Python 3.12:** los paquetes ahora estan fijados a NumPy 2.2.6, Numba 0.61.2 y llvmlite 0.44.0; lea `../MACOS_NUMBA_FIX.md`.

`numba` es muy recomendable. Con él, el integrador avanza unos 6×10⁶ pasos·realización por
segundo y por núcleo, y usa todos los núcleos disponibles. Sin `numba` el código funciona
igual (backend `numpy`), pero entre 30 y 100 veces más lento.

## 2. Estructura

```
run_all.sh               ejecuta todo en orden (pasos 3–5)
run_sde.py               simulaciones; un archivo por lote; reanudable
aggregate.py             estadística, comparación con teoría lineal, tablas y diagnósticos
update_manuscript.py     inserta los resultados en main.tex, supplement.tex y la carta
config/cases.json        presets: parámetros, casos, longitud de registro, nº de realizaciones
lk/                      modelo: rev_core.py (ECM por θ, resolvente), lk_model.py, sde_kernel.py
manuscript_current/      copia de los .tex revisados tal como se entregaron
EXPECTED_RESULTS.md      valores esperados, tiempos y criterios de aceptación
results/                 (se crea) raw/<preset>/*.npz, resúmenes CSV, JSON, tablas, figuras
pilot_results/           corrida piloto de referencia (mismo formato que results/); no la borre
```

## 3. Presets

| preset | τ | κ | casos | registro T | realizaciones | uso en el manuscrito |
|---|---|---|---|---|---|---|
| `smoke` | 20 | 0.008 / 0.009 | solitario, trade-off, win-win | 4 000 | 16 | prueba de instalación (≈1 min) |
| `tau20` | 20 | 0.008 / 0.009 | ídem | 40 000 | 128 | control de los casos de la Tabla 1 (no se inserta automáticamente) |
| `long200` | 200 | 0.005 | solitario, θ = −0.156π, θ = 0 | 100 000 | 128 | Sec. 7, Tabla S14, carta R3 |
| `long1000` | 1000 | 0.001 | ídem | 400 000 | 128 | Sec. 7, Tabla S14, carta R3 |
| `fixedkt` | 200 | 0.0009 | solitario, θ = 0.018π | 200 000 | 256 | opcional (efecto de RIN pequeño, 0.933) |

Cada caso se define por (κ, θ/π). θ es la fase de retorno total, de modo que se sigue el ECM
elegido aunque coexistan otros. Los demás parámetros son α = 4, ε = 0.005, P = 1,
σ = 5×10⁻⁴ y Δt = 0.02, iguales a los de la Sec. S5.

## 4. Ejecución

```bash
python run_sde.py --preset smoke          # 1) verificar
bash run_all.sh                           # 2) todo: smoke, tau20, long200, long1000 y agregación
```

Los presets también se pueden lanzar por separado. Cada lote se guarda en
`results/raw/<preset>/<caso>_bNNN.npz` y nunca se recalcula. Por eso:

- **Interrumpir y reanudar:** basta volver a ejecutar el mismo comando.
- **Más realizaciones:** `python run_sde.py --preset long1000 --batches 32` añade lotes a los
  ya existentes.
- **Registros más largos:** como no se pueden mezclar lotes con distinta longitud T, use una
  variante, `python run_sde.py --preset long1000 --T 1e6 --tag T1e6`, y agregue con
  `python aggregate.py --presets long200 long1000_T1e6`.

  Para insertar esa variante en el manuscrito, renombre o copie la carpeta
  `results/raw/long1000_T1e6` como `results/raw/long1000`. La forma más simple es definir un
  preset nuevo en `config/cases.json` con el T deseado.
- **Amplitud de ruido:** `--sigma 2.5e-4 --tag s25`. Sirve para comprobar que una
  discrepancia es o no un efecto no lineal: las correcciones no lineales escalan como σ².
- **Control de paso temporal** (recomendado, una vez):
  `python run_sde.py --preset long200 --dt 0.01 --tag dt001` y luego
  `python aggregate.py --presets long200 long200_dt001`. Las razones deben coincidir dentro
  de la incertidumbre.

Las semillas son deterministas (hash de preset, variante, caso, lote y realización). El mismo
comando reproduce los mismos números en la misma máquina; entre máquinas distintas pueden
aparecer diferencias en los últimos dígitos.

## 5. Agregación y diagnóstico

```bash
python aggregate.py --presets tau20 long200 long1000
```

Compare lo que obtenga con `EXPECTED_RESULTS.md`, que incluye los valores de la corrida
piloto, los criterios de aceptación y los tiempos.

Para cada caso, el script imprime la razón respecto al láser solitario simulado ±2 errores estándar,
la predicción lineal entre paréntesis, el estadístico z y la deriva máxima de ⟨|E|²⟩. Además
escribe en `results/`:

- `summary_<preset>.csv`: todas las cifras, incluidas las razones absolutas simulación/teoría.
- `lag_sensitivity.csv`: el ancho de línea con tres ventanas de ajuste de fase.
- `lag_windows_<preset>.csv`: la tasa de difusión de fase en ventanas sucesivas de lags
  (3–12τ, 12–30τ, …), en valor absoluto respecto a la teoría. Detecta un posible vagabundeo
  lento de la fase.
- `diag_<preset>.png`: tasa local de difusión de fase, d(MSD)/dt / 2D_lin. Debe fluctuar
  alrededor de 1 dentro de la ventana sombreada.

El ajuste de difusión de fase resta un pequeño desplazamiento residual de la frecuencia media,
común a todas las realizaciones (ver `EXPECTED_RESULTS.md`, observación 3). Para ver el
resultado sin esa corrección, use `--no-detrend`.
- `manuscript_numbers.json` y `tableS14_rows.tex`: lo que usa el paso 6.

## 6. Inserción en el manuscrito

```bash
python update_manuscript.py \
   --main   manuscript_current/main.tex \
   --supp   manuscript_current/supplement.tex \
   --letter manuscript_current/response_letter.tex
```

Los archivos de `manuscript_current/` ya contienen los resultados piloto; el script los
reemplaza por los suyos. Si después editó a mano estos archivos, apunte a sus versiones.

El script reemplaza exactamente tres pasajes:

1. En `main.tex`, Sec. 7: el bloque de frases sobre simulaciones estocásticas.
2. En `supplement.tex`, Sec. S10: el párrafo «Nonlinear simulations», que recoge la nueva
   longitud de registro, el número de realizaciones y la ventana de ajuste, y las filas de la
   Tabla S14.
3. En la carta, respuesta al Revisor 3: el ítem «Nonlinear validation».

Escribe `*_updated.tex`, muestra el diff y **no sobrescribe nada**. Si algún |z| > 3 o se
detecta deriva de intensidad (posible salto de modo), imprime un AVISO y marca el texto con
«CHECK». En ese caso no inserte los resultados sin revisarlos.

Después:

1. Revise el diff.
2. Renombre `*_updated.tex` a sus nombres originales.
3. Compile la versión limpia y la marcada (`touch marked.flag`).
4. Suba los scripts y `results/*.csv` al repositorio de GitHub y Zenodo citado en *Data
   availability*.
