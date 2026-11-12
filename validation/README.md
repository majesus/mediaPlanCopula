# Validación (excluida del paquete)

* `validate.R`: comprobaciones numéricas de exactitud (marginales, r1/r2, duplicaciones, conservación de GRPs).
* `simulation/`: estudio de simulación pareado frente a los modelos de **mediaPlanR** (8 procesos generadores, 150 escenarios cada uno; régimen poblacional y régimen de panel de 3000). Resultados en `tab_contrastes.csv` y `tab_descriptivo.csv`. G1 y H1 son el mundo propio de la cópula (circulares); los demás no.
* `real_data/`: validación con 440 programas observados del panel comScore (Cheong, 2007, Apéndices E–G). Los datos de entrada **no se distribuyen** porque proceden de una tesis con derechos de terceros; `real_cheong.R` documenta el formato esperado (`vehicles.csv`, `pairs.csv`, `schedules.csv`). Es una validación en muestra: no mide capacidad predictiva fuera de muestra.

## Alcance e interpretación

* Los veredictos de `simulation/` (por ejemplo «copula mejor (relevante)») son condicionales a los ocho procesos generadores y a los regímenes evaluados, y no son una afirmación general de precisión. Los procesos G1 y H1 son el mundo propio de la cópula (circulares). En H1 la cópula resulta peor que CANEX, CBD, CSD y MSAD en el régimen de panel de 3000 (diferencias de TV negativas en la mayoría de los contrastes), y en H2 los contrastes son en gran parte no concluyentes.
* `real_data/` es una validación en muestra con el panel comScore; la ventaja de TV frente a CANEX, CBD y CSD es de 0,003 a 0,005.
* Una comparación posterior con microdatos AMPS 2011 (22 diarios, n = 25 160), con los mismos modelos de **mediaPlanR** y las mismas entradas (r1, r2, duplicación), dio una ventaja de la cópula de 0,0015 a 0,0021 de TV con una verdad hipergeométrica construida a partir de la frecuencia declarada, diferencias no distinguibles del ruido muestral fuera de muestra (calibración en una mitad y verdad en la otra) y una ventaja de CANEX, CBD, CSD y MSAD con exposición observada de una inserción («ayer»). Esa comparación y sus cotas de identificación no forman parte de este directorio.
* Las entradas de los programas de `real_data/` no se distribuyen (ver arriba).

