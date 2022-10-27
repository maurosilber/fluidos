# Corriendo GHOST

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

## Correr la simulación

Para correr la simulación,
hay que crear la carpeta `output/`
y poner el archivo `parameter.inp` ahí.
`make output/vx.0001.out` corre la simulación.
