Comparación power analysis
================
2020-04-17

Carga de dependencias:

``` r
library("data.table")
library('tidyverse')
library("pwr")
```

Función para cálculo manual de la potencia de un test de diferencia de
proporciones entre muestras de distinto
tamaño:

``` r
power_calculation_manual <- function(thetaA, effect_size, nA, nB, alpha = 0.05, tail = c('two.tailed', 'right', 'left')) {
  thetaB_null <- thetaA
  thetaB_alt <- thetaA + effect_size
  
  sigmaA <- sqrt(thetaA * (1 - thetaA) / nA)
  sigmaB_null <- sqrt(thetaB_null * (1 - thetaB_null) / nB)
  sigmaB_alt  <- sqrt(thetaB_alt * (1 - thetaB_alt) / nB)
  
  muDiff_null <- thetaA - thetaB_null
  sigmaDiff_null <- sqrt(sigmaA^2 + sigmaB_null^2)
  crit_value <- case_when(tail == 'right' ~ qnorm(1 - alpha, muDiff_null, sigmaDiff_null),
                          tail == 'left' ~ qnorm(alpha, muDiff_null, sigmaDiff_null),
                          tail == 'two.tailed' ~ qnorm(1 - alpha / 2, muDiff_null, sigmaDiff_null))
  
  muDiff_alt <- thetaA - thetaB_alt
  sigmaDiff_alt <- sqrt(sigmaA^2 + sigmaB_alt^2)
  
  beta <- pnorm(crit_value, muDiff_alt, sigmaDiff_alt)
  beta <- case_when(tail == 'right' ~ pnorm(crit_value, muDiff_alt, sigmaDiff_alt),
                    tail == 'left' ~ 1 - pnorm(crit_value, muDiff_alt, sigmaDiff_alt),
                    tail == 'two.tailed' ~ pnorm(crit_value, muDiff_alt, sigmaDiff_alt) - pnorm(-crit_value, muDiff_alt, sigmaDiff_alt) )
  power <- 1 - beta
  list(thetaA = thetaA, effect_size = effect_size, nA = nA, nB = nB, alpha = alpha, tail = tail, beta = beta, power = power)
}
```

Datos de ejemplo, tomados del [blog post de Despegar Research sobre AB
testing](https://researchdespegar.wordpress.com/2015/06/19/what-is-hypothesis-te)

``` r
alpha <- 0.05
tail <- 'left'
nA <- 2000000
nB <- 1000000
thetaA <- 0.02
effect_size <- thetaA * 1.01 - thetaA
```

Output del cálculo manual:

``` r
power_calculation_manual(thetaA, effect_size, nA, nB, alpha, tail)
```

    ## $thetaA
    ## [1] 0.02
    ## 
    ## $effect_size
    ## [1] 2e-04
    ## 
    ## $nA
    ## [1] 2e+06
    ## 
    ## $nB
    ## [1] 1e+06
    ## 
    ## $alpha
    ## [1] 0.05
    ## 
    ## $tail
    ## [1] "left"
    ## 
    ## $beta
    ## [1] 0.6832747
    ## 
    ## $power
    ## [1] 0.3167253

Es consistente con el post previamente mencionado.

El output de usar el [paquete de R
pwr](https://cran.r-project.org/web/packages/pwr/index.html) arroja el
mismo
resultado:

``` r
pwr.2p2n.test(h = ES.h(thetaA, thetaA + effect_size), n1 = nA, n2 = nB, sig.level = alpha, alternative = 'less')
```

    ## 
    ##      difference of proportion power calculation for binomial distribution (arcsine transformation) 
    ## 
    ##               h = -0.00142509
    ##              n1 = 2e+06
    ##              n2 = 1e+06
    ##       sig.level = 0.05
    ##           power = 0.3151615
    ##     alternative = less
    ## 
    ## NOTE: different sample sizes

Ahora bien, para correrlo correctamente se debió usar la función ES.h
para calcular el effect size apropiado, tal como indican el manual del
paquete y [la viñeta del
mismo](https://cran.r-project.org/web/packages/pwr/vignettes/pwr-vignette.html).
Si se pasa el effect size sin ser transformado, el resultado es
otro:

``` r
pwr.2p2n.test(h = effect_size, n1 = nA, n2 = nB, sig.level = alpha, alternative = 'less')
```

    ## 
    ##      difference of proportion power calculation for binomial distribution (arcsine transformation) 
    ## 
    ##               h = 2e-04
    ##              n1 = 2e+06
    ##              n2 = 1e+06
    ##       sig.level = 0.05
    ##           power = 0.03529135
    ##     alternative = less
    ## 
    ## NOTE: different sample sizes

Pareciera entonces que la transformación que es necesario realizar para
pasar al parámetro h tiene que ver con la implementación específica
elegida por este paquete, más que con un procedimiento necesario en
general. Este paquete es la implementación del libro [Cohen -
Statistical Power Analysis for the Behavioral
Sciences](https://drive.google.com/open?id=1_fZ6zUU9Zv9-gnBZVlz2Uw80aQNCghK1).
Allí dice (sección 6.2, p.180):

If we were to define j = P 1 - P 2 , and try to use j as our ES, we
would soon discover that the delectability of some given value of j,
under given fixed conditions of a and n, would not be constant, but
would vary depending upon where along the scale of P between zero and
one the value j occurred.

Luego procede a sugerir usar como h la diferencia entre las
probabilidades transformadas por 2\*arcsin(sqrt(p)). No obstante, en la
implementación power\_calculation\_manual no se hizo nada por el estilo,
y sin embargo llegamos al mismo resultado.

En efecto, en ese mismo capítulo un poco más arriba menciona:

It is possible to approach the testing of hypotheses about proportions
by different statistical techniques, including **the classical normal
curve test using a “critical ratio” applied directly to the
proportions** (Edwards, 1972, pp. 42-44; Guilford\&Fruchter, 1978,
pp. 159-161; Blalock, 1972, pp. 228-232), by a chi-square contingency
test (see Chapter 7 and references), by a special case of the
hypergeometric probability distribution ("Fisher’s Exact Method’’) for 2
x 2 tables (Hays, 1981, pp. 552-554; Owen, 1962, pp. 479-496), **or by
means of a normal curve test applied to the arcsine transformation of
the proportions. Despite its unfamiliarity, it is the last of these
alternatives that provides the basis for the approach of this chapter
because of certain advantages it has, particularly from the viewpoint of
power analysis**. However, the results from using any of these
procedures will be the same to a close approximation, particularly when
samples are not small (Cohen, 1970).

Lo bueno es que refuerza y reconoce que su método es análogo al clásico.
Lo malo es que por ahora no me quedan claros los pormenores de su
método, ni cuáles son esas *ciertas ventajas* que menciona.

En la introducción (sección 1.4) habla de la necesidad de normalizar
escalas para desprenderse de unidades, en particular teniendo en cuenta
que después presenta tablas de referencia para consultar potencias y
tamaños de muestra. Así que puede ser que una normalización permita una
ventaja en ese sentido. Pero la verdad que por ahora no me queda claro;
y tampoco sé si vale la pena en esta instancia.
