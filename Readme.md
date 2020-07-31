# gggpScatteringIHQCD
C++ code that implements a holographic description of the following physical observables: photon structure function <a href="https://www.codecogs.com/eqnedit.php?latex=F^{\gamma}_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F^{\gamma}_2" title="F^{\gamma}_2" /></a>, the hadronic photon-photon total cross-section <a href="https://www.codecogs.com/eqnedit.php?latex=\sigma\left(\gamma&space;\gamma&space;\rightarrow&space;X\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma\left(\gamma&space;\gamma&space;\rightarrow&space;X\right)" title="\sigma\left(\gamma \gamma \rightarrow X\right)" /></a>,
the proton structure functions <a href="https://www.codecogs.com/eqnedit.php?latex=F^{p}_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F^{p}_2" title="F^{p}_2" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=F^{p}_L" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F^{p}_L" title="F^{p}_L" /></a> and the hadronic photon-proton total cross section <a href="https://www.codecogs.com/eqnedit.php?latex=\sigma\left(\gamma&space;p&space;\rightarrow&space;X&space;\right&space;)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma\left(\gamma&space;p&space;\rightarrow&space;X&space;\right&space;)" title="\sigma\left(\gamma p \rightarrow X \right )" /></a>.

# Dependencies
This packages requires the [Boost](https://www.boost.org/) libraries on your computer to solve numerically a set of differential equations.
 The minimum required version is 1.53.0. 
Another important intermediary step of these computations is the solution of the eigenvalue problem of a square matrix. To do this end we use the [Armadillo](http://arma.sourceforge.net) C++ library.
 The minimum required version is 8.600.1
 
# How to use?
First you will need to create the executable files. These files will be contained inside a folder called execs.
To build just type in your shell
```r
cd path_of_the_directory
cmake .
make
```
After these steps a folder named execs will appear inside path_of_the_directory.

As an example of usage let's say you want to determine the best fit parameters all the available data of <a href="https://www.codecogs.com/eqnedit.php?latex=F^{\gamma}_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F^{\gamma}_2" title="F^{\gamma}_2" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\sigma\left(\gamma&space;\gamma&space;\rightarrow&space;X\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma\left(\gamma&space;\gamma&space;\rightarrow&space;X\right)" title="\sigma\left(\gamma \gamma \rightarrow X\right)" /></a>.
Then type in the shell
```r
./execs/fitGammaGammaAllPars.exe f2photon_data_path sigma_gg_data_path invls a b c d k1 k2 k3 k4
```
The first two parameters are the location of the data we want to use while the rest are the initial guess values used to start the fit.
In case you don't pass input to the program it will run with default values and display a message explaining the usage.

# Acknowledgements
This research received funding from the Simons Foundation grants 488637 (Simons collaboration on the Non-perturbative bootstrap) and from the grant CERN/FIS- PAR/0019/2017. I was also funded by Fundação para a Ciência e a Tecnologia (FCT) under the IDPASC doctorate programme with the fellowship PD/BD/114158/2016.