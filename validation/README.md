# Validación (excluida del paquete)

* `validate.R`: comprobaciones numéricas de exactitud (marginales, r1/r2, duplicaciones, conservación de GRPs).
* `simulation/`: estudio de simulación pareado frente a los modelos de **mediaPlanR** (8 procesos generadores, 150 escenarios cada uno; régimen poblacional y régimen de panel de 3000). Resultados en `tab_contrastes.csv` y `tab_descriptivo.csv`. G1 y H1 son el mundo propio de la cópula (circulares); los demás no.
* `real_data/`: validación con 440 programas observados del panel comScore (Cheong, 2007, Apéndices E–G). Los datos de entrada **no se distribuyen** porque proceden de una tesis con derechos de terceros; `real_cheong.R` documenta el formato esperado (`vehicles.csv`, `pairs.csv`, `schedules.csv`). Es una validación en muestra: no mide capacidad predictiva fuera de muestra.
