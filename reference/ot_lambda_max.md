# Maximum finite rate of population increase (Cole's equation)

Solves Cole's (1954) characteristic equation for the maximum intrinsic
rate of natural increase, returning the finite (per-year) multiplication
rate \\\lambda\_{max} = e^{r\_{max}}\\. This is the rate required by the
production model
[`ot_pro()`](https://stangandaho.github.io/offtake/reference/ot_pro.md)
when it is not supplied directly.

## Usage

``` r
ot_lambda_max(b, a, w, interval = c(1e-06, 5))
```

## Arguments

- b:

  Mean annual number of female offspring per female (i.e. litter size
  \\\times\\ litters per year \\\times\\ proportion female).

- a:

  Age at first reproduction (years).

- w:

  Age at last reproduction (years); often approximated by lifespan.

- interval:

  Search interval for \\r\_{max}\\ passed to
  [`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html).

## Value

A numeric vector, recycled to the length of the longest of `b`, `a` and
`w`, giving \\\lambda\_{max}\\ – the maximum finite (per-year) rate of
increase. It is dimensionless and greater than 1 for a growing
population (e.g. `1.30` means the population can grow by at most 30% per
year). `NA` is returned for elements with missing inputs or with no
positive root inside `interval`.

## Details

Cole's equation is \$\$1 = e^{-r} + b\\e^{-r a} - b\\e^{-r (w + 1)}\$\$
where `b` is the mean annual number of *female* offspring per female,
`a` is the age at first reproduction and `w` is the age at last
reproduction (all in years). The equation is solved numerically for
\\r\\ and \\\lambda\_{max}=e^{r}\\ is returned.

## References

Cole, L. C. (1954) The population consequences of life history
phenomena. *The Quarterly Review of Biology* 29, 103-137.

Robinson, J. G. & Redford, K. H. (1986) Intrinsic rate of natural
increase in Neotropical forest mammals. *Oecologia* 68, 516-520.

## Examples

``` r
# Inputs are life-history traits (not columns of a data frame here):
#   b = female young per female per year, a = age at first reproduction,
#   w = age at last reproduction. A duiker breeding from age 1 to 8:
ot_lambda_max(b = 0.6, a = 1, w = 8)
#> [1] 1.584931

# Vectorised over several species at once:
ot_lambda_max(b = c(0.6, 0.5), a = c(1, 2), w = c(8, 12))
#> [1] 1.584931 1.355822
```
