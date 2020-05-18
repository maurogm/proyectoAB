#' ---
#' title: "Poisson test"
#' author: ""
#' date: "`r Sys.Date()`"
#' output: 
#'   github_document:
#'     toc: true
#' ---

#'<style>
#'  body {
#'    text-align: justify}
#'</style>

#+ r knitr-global-options, include = FALSE
knitr::opts_chunk$set(warning = FALSE, message=FALSE, verbose=TRUE)
knitr::opts_knit$set(root.dir=normalizePath('../'))

# Para compilar el reporte correr: rmarkdown::render('scripts/poisson_tests.R')

#+ r Resto

#' Busco jugar un poco con 2 distintos tests de Poisson para diferencias entre medias, explicados en el paper
#' [A more powerful test for comparing two Poisson means](https://userweb.ucs.louisiana.edu/~kxk4695/JSPI-04.pdf), 
#' de K. Krishnamoorthy y Jessica Thomson.

#' 
#' # Reproducción manual de las posibles formas de realizar el test

#' Carga de dependencias:
library("data.table")
library('tidyverse')

#' Antes que nada, seteo algunos parámetros de ejemplo:
set.seed(999)
alpha <- 0.05
n1 <- 10
n2 <- 15
lambda1 <- 6
lambda2 <- 4

#' Simulo un par de muestras:
muestra1 <- rpois(n1, lambda1)
muestra2 <- rpois(n2, lambda2)


#' ## Usando el paquete 'stats' de R

k1 <- sum(muestra1) # Valores observados de una Poisson(n1 * lambda1)
k2 <- sum(muestra2) # Valores observados de una Poisson(n2 * lambda2)
k <- k1 + k2

#' Resultado de la función estandar de R:
poisson.test(c(k1, k2), c(n1, n2), alternative = 'greater')
#poisson.test(c(850, 780), c(1, 1), alternative = 'greater')

#' ## Implementación manual del C-test

#' Planteo hipótesis nulas de la forma:
#' 
#' H0: lambda1 / lambda2 <= c, o bien
#' 
#' H0: lambda1 - lambda2 <= d

#' 
#' (X1|X1+X2=k) ~ Binomial(n = k, p = lambda1 / (lambda1 + lambda 2)),
#' 
#' O si n1 != n2, p = (n1/n2)\*(lambda1/lambda2) / (1 + (n1/n2)\*(lambda1/lambda2))

p_cociente <- function(n1, n2, cociente)  (n1/n2)*(cociente) / (1 + (n1/n2)*(cociente))

#' Bajo la Hipótesis nula de que el cociente sea igual a c (por ejemplo c = 1, para lambda1=lambda2):
c <- 1
d <- 0
p <- p_cociente(n1, n2, c)

#' Dado que bajo la hipótesis nula y los k casos totales observados, X1 sigue una distribución binomial que conocemos,
#' podemos calcular el p-valor (la probabilidad de que X1 tome un valor tanto o más extremo que el valor observado k1)
#' como la suma de las probabilidades puntuales para todos los valores entre k1 y k:
comb <- function(n, k) factorial(n) / factorial(n-k) / factorial(k) #Función combinatoria
sumando <- function(p, k, i) comb(k, i) * p^i * (1 - p)^(k - i) #Términos de la sumatoria
pvalor <- function(k1, k) {
  seq(k1, k, 1) %>%
  lapply(function(i) sumando(p, k, i)) %>% 
  unlist() %>% 
  sum()
}
pvalor(k1, k)

#' Observamos entonces que el p-valor así calculado es idéntico al que había arrojado la función poisson.test.
#' 
#' **Alternativamente**, en vez de calcular el p-valor, puedo ir calculando la suma cumulativa de probabilidades puntuales,
#' de modo que al llegar a _1 - alpha_ habré encontrado el valor crítico del test.
x1_posibles <- seq(0, k)
probs_acumuladas <- cumsum(sumando(p, k, xs))
crit_val <- x1_posibles[probs_acumuladas > (1 - alpha)][1]
# Rechazo H0 sii k1 > crit_val
crit_val

#' De este modo se obtiene el valor crítico ```r crit_val```, que tal vez podría ser utilizado para un eventual análisis de potencia
#' (bajo el supuesto de que el total de casos sea igual a k). Dado que k1 tiene un valor de ```r k1```, decido no rechazar;
#'  lo que está en consonancia con el p-valor obtenido.
#'
#' ##  Implementación manual del E-test
#' 
#' ### Cálculos preliminares
#' 
#' #### Validación empírica de la varianza de la diferencia
#' 
#' A través de simulaciones veo que la varianza de X1 - X2 es efectivamente lambda1/n1 + lambda2/n2.

nsims <- 10000
sims <- rep(0, nsims)
for (i in 1:nsims) {
  n1_sim <- 1000
  n2_sim <- 100
  X1_sim <- rpois(n1, lambda1) %>% sum()
  X2_sim <- rpois(n2, lambda2) %>% sum()
  
  diff <- X1_sim / n1_sim - X2_sim / n2_sim
  sims[i] <- diff
}
var(sims)
lambda1/n1 + lambda2/n2

#' ### Implementación del E-test

#' Estimador de la varianza de la diferencia entre X1 y X2:
var_hat <- function(x1, x2, n1, n2) x1 / n1^2 + x2 / n2^2
#' Se usa la diferencia estandarizada como estadístico pivot, donde _d_ es la diferencia pautada por la Hipótesis Nula (comunmente 0):
T_x1_x2 <- function(x1, x2, n1, n2, d) (x1 / n1 - x2 / n2 - d) / sqrt(var_hat(x1, x2, n1, n2))

#' Dadas las muestras (k1, k2, n1, n2) obtenidas, el valor observado para el estadístico T es:
T_k1_k2 <- T_x1_x2(k1, k2, n1, n2, d)

#' Bajo la hipótesis nula, lambda1 = lambda2 + d; lo que lleva al siguiente estimador de lambda2:
lambda2k_hat <- (k1 + k2) / (n1 + n2) - d * n1 / (n1 + n2)

#' Si lambda2 da negativo, no se rechaza, dado que mi mejor estimación de lambda2 está dando más grande que lambda1 + d:
if (lambda2k_hat <= 0) {
  reject <- FALSE
  break
}

#' Para una Poisson(lambda), la función de probabilidad puntual es:
p_puntual_poisson <- function(x, lambda) exp(-lambda) * lambda^x / factorial(x)

#' Luego, para el vector (X1, X2) bajo la H0 y con lambda2 estimado por lambda2k_hat, la función de probabilidad puntual es:
p_puntual_x1_x2 <- function(x1, x2) p_puntual_poisson(x2, n2 * lambda2k_hat) * p_puntual_poisson(x1, n1 * (lambda2k_hat + d)) 

#' El p-valor de nuestro test es la posibilidad de haber observado un valor ```T_x1_x2``` tan o más extremo que ```T_k1_k2``` bajo la H0.
#' Y esto puede ser calculado exactamente como la suma de las probabilidades puntuales ```p_puntual_x1_x2``` en los casos en que ```T_k1_k2 >= T_k1_k2```.

#' El cálculo exacto del p-valor implicaría una suma infinita, dado que teóricamente x1 y x2 pueden tomar el valor de cualquier entero positivo.
#' No obstante dado que lo único que nos interesa es saber si el p-valor está por debajo del nivel de significancia _alpha_, podemos dejar de sumar
#' una vez que tengamos certeza sobre si esto será así; ya sea porque la suma del p-valor excedió _alpha_, 
#' o porque la suma sobre ```T_k1_k2 < T_k1_k2``` excedió _1-alpha_.
#' 
#' Dado que esta no pretende ser la implementación final del E-test, no se tendrán consideraciones como buscar dinámicamente los x1, x2 en los que dejar de sumar,
#' ni aprovechar la recursividad de P(X = j+1) = P(X = j) * lambda / (j+1) para hacer más eficientes las cuentas.
#' 
#' En vez de ello, tomamos "a ojo" un rango suficientemente grande y calculamos las probabilidades puntuales correspondientes.

expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
x1_max <- 2 * ceiling(n2 * lambda2k_hat)
x2_max <- 2 * ceiling(n1 * (lambda2k_hat + d))
grid_x1_x2 <- expand.grid.df(1:x1_max, 1:x2_max) %>% setnames(c('x1', 'x2'))
probs <- map2_dbl(grid_x1_x2$x1, grid_x1_x2$x2, p_puntual_x1_x2)
estadisticoT <- map2_dbl(grid_x1_x2$x1, grid_x1_x2$x2, function(x1, x2) T_x1_x2(x1, x2, n1, n2, d))
indicadora <- estadisticoT >= T_k1_k2

sum(probs)
#' Observemos que considerar estos valores de (x1, x2) barre con casi la totalidad de la probabilidad acumulable, 
#' dejando por fuera una cola de tan solo ```r 1 - sum(probs)```, lo cual sería relevante sólo en casos límite.
 
sumandos_consolidados %>% 
  ggplot() +
  geom_contour(aes(x1, x2, z = probs, colour = ..level.., linetype = indicadora)) +
  theme_minimal() +
  ggtitle("Probabilidades puntuales bajo H0 para distintos valores de x1 y x2")

sumandos_consolidados <- data.table(grid_x1_x2, probs, estadisticoT, indicadora)

sumandos_consolidados[, .(suma_probabilidades = sum(probs)), indicadora]
#' En este ejemplo, observamos que la suma de las probabilidades puntuales supera el umbral de 0.95.
#' Por lo tanto, podemos estar seguros de que el p-valor estará por debajo de 0.05, y por lo tanto rechazamos la hipótesis nula.
#' 
#' Es interesante notar que el E-test logra rechazar la Hipótesis Nula, mientras que el C-test no lo había hecho, 
#' al haber llegado a un p-valor de ```r pvalor(k1, k) ```.
#' Esto está en consonancia con la afirmación respecto a mayor potencia por parte de Krishnamoorthy y Thomson.

