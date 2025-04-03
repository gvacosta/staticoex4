# staticoex4 (development version) (lo que está en progreso y será parte de la versión 1.0.0)

* En DESCRIPTION hay que consultar si en "autores" hay que poner el nombre de alguien más que sea el "R package maintainer"
* Se empezó el release checklist: Hay que finiquitar DESCRIPTION, Readme y otros archivos más.
* Hay que actualizar los datos respecto a la traducción de títulos, y hay datos de meses que podrían estar en inglés también.
* Hay que realizar los testthat()
* Initial CRAN submission (usethis::use_release_issue)


## Primeras acciones realizadas

* Toda la creación del paquete se realizó siguiendo un ejemplo subido en youtube: https://youtube.com/playlist?list=PLfX5C7cc6LRIuAYefr8tJw9JABx656lmc&si=FddQDZU3g7M_n1U6
* Primero se conectó R con Github. Luego se creó el paquete con su estructura básica: Licencia, Documentaciones, Readme, etc, y finalmente se cargó el paquete a github para crear el repositorio.
* Se crearon uno por uno los scripts. Primero las funciones: PrepTables() que sirve para configurar el set de datos de la manera en que la función HJSTATICO los necesitará. Luego se creó la función HJSTATICO. Cada función con sus respectivos argumentos según la configuración necesaria para R packages.
* Luego se cargaron los datasets al paquete. Cada función y dataset con su script y archivo .Rd (documentación) correspondiente.
* Luego de subir todo el repositorio, se crearon las GitHub actions para colocar etiquetas de funcionamiento/rendimiento en los diferentes sistemas operativos.
* Se continuó a la subida del paquete, ahora con la guía de Jim Hester en: https://www.youtube.com/watch?v=-zID-rVDEHQ&t=22s
* Initial CRAN submission (usethis::use_release_issue)
* Se empezó el archivo NEWS.md para documentar los cambios, pasos, actualizaciones en la creación del paquete.
