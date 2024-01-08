# LLPgammaAnalyzer
Long Lived Particle to Photon Analyzer

## Bayesian Photon Shower Shapes
### Install [BayesianClustering](https://github.com/mlazarovits/BayesianClustering)
First, install the BayesianClustering package. Be sure this is done in the CMSSW `src/` directory. This package does not rely on CMSSW, but the paths in this skimmer's `makefile` rely on the directory structure of CMSSW.
```
git clone https://github.com/mlazarovits/BayesianClustering master
```
### Make the shared library
Go to the BayesianClustering package and make the shared library.
``` make lpclib
```
### Update paths in environment variables
The CGAL package needs to be found during compilation.
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/cgal/4.2/lib
```
In order for the skimmer to be able to find the library during compilation, the proper environment variables must be updated to include the package path.
```
export $LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/BayesianClustering/lib
```
This line can also be added at the end of your `~/.bash_profile` to have this set up upon login. The CGAL export line can also be added to this bash script, but be sure it is *before* the BayesianClustering line.

### Include the header
Then, in your executable be sure to include the header `#include BayesianClustering/BayesianClustering.hh`

### Troubleshooting
For more information, see the BayesianClustering README.
