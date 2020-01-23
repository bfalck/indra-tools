# indra-tools
Code repository for the Indra simulations, hosted on the [SciServer](http://sciserver.org).

This branch, `pre_runid`, preserves the old API in which the Indra run had to be specified by 3 `X_Y_Z` integers in the function calls. Each goes from 0 to 7, corresponding to the raveled index of the 512 Indra simulations in an 8x8x8 cube. In the new API, the run can either be specified by the unraveled index (from 0 to 511) or a tuple of the raveled index, (`X`,`Y`,`Z`).


Written by Bridget Falck, 2018-2020