# Corriendo GHOST

Para correr todas las simulaciones,
`make` (o `make all`)
compila el solver,
genera los parámetros a partir de `parameters.py`,
y corre todas las simulaciones,
que quedan guardadas en `output/<variante>/`.

## Instalando dependencias

Para instalar las dependencias para GHOST,
`make install`
crea un *conda environment* a partir del `environment.yaml`
que está en el *root* del repositorio.

## Compilando el solver

Para compilar el solver,
`make SOLVER`
copia el directorio `src` al de GHOST,
compila y mueve el binario a este directorio.

## Generando parámetros

Para generar los parámetros de las distintas corridas,
`make parameters`
corre el archivo `parameters.py`,
que lee `parameters.in` y reemplaza ciertos parámetros.
Para cada condición,
guarda un archivo de parámetros en una carpeta dentro de `output/`.

## Corriendo una simulación para una dada condición

Para correr una dada condición `<cond>`,
`make output/<cond>`
hardlinkea el binario SOLVER en el directorio correspoindiente
y corre la simulación.
