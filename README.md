Medbetareg is a R package containing the software for the method proposed in the article
__Generalized Beta Regression to Elicit Conditional Distributions of Medical Variables__
by A. Magrini, D. Luciani and F. M. Stefanini, to be appeared on _Austrian Journal of Statistics_.

R (The R Project for Statistical Computing) need to be installed on your system for
Medbetareg to work. R can be downloaded from https://www.r-project.org/.

To install Medbetareg, open the console of R and type:
```
install.packages(path_to_file, repos = NULL, type="source")
```
where 'path_to_file' is the directory containing the file 'Medbetareg_1.0.tar.gz'.

If you find any bug, please write to <alessandro.magrini@outlook.com> (Alessandro Magrini)

Below, you find the code to replicate the results in the article.
_________________________________________________________________

Code for expert assessments (Table 1):
```
assess.test <- 'RESP RespRate (0,15,25,40);
  CEV intraShunt (0,2,5,100) 0.5 hp 0.5 hp 5;
  CEV deadSpace (0,0,30,100) 0.5 hp 0.5 hp 5;
  CEV extraShunt (0,0,5,100) 0.5 hp 0.5 hp 5;
  CEV redAlvSpace (0,0,5,100) 0.5 hp 0.5 hp 5;
  BEV Panic 0.25 hp 25;
  BEV Neuromusc 0.6 lp 100;
  INTER intraShunt deadSpace 0.9 hp 5;
  TAU 0.3 n'
```
Computation of the prior:
```
set.seed(10)
prior.test <- newPrior(assess.test, nrep=5000)
```
Marginal quantile summaries of the prior (Table 2):
```
summary(prior.test)
```
Probability density of the predictive distributions inspected during the revision of the prior (Figure 5):
```
set.seed(10)
predictive(prior.test, xcfg1, nrep=50000, title="Configuration 1")
set.seed(10)
predictive(prior.test, xcfg2, nrep=50000, title="Configuration 2")
set.seed(10)
predictive(prior.test, xcfg3, nrep=50000, title="Configuration 3")
set.seed(10)
predictive(prior.test, xcfg4, nrep=50000, title="Configuration 4")
```
