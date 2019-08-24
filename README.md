# CharibdeOptim

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://yashcodes.github.io/CharibdeOptim.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://yashcodes.github.io/CharibdeOptim.jl/dev)
[![Build Status](https://travis-ci.com/yashcodes/CharibdeOptim.jl.svg?branch=master)](https://travis-ci.com/yashcodes/CharibdeOptim.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/yashcodes/CharibdeOptim.jl?svg=true)](https://ci.appveyor.com/project/yashcodes/CharibdeOptim-jl)
[![Codecov](https://codecov.io/gh/yashcodes/CharibdeOptim.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/yashcodes/CharibdeOptim.jl)
[![Coveralls](https://coveralls.io/repos/github/yashcodes/CharibdeOptim.jl/badge.svg?branch=master)](https://coveralls.io/github/yashcodes/CharibdeOptim.jl?branch=master)
[![Build Status](https://api.cirrus-ci.com/github/yashcodes/CharibdeOptim.jl.svg)](https://cirrus-ci.com/github/yashcodes/CharibdeOptim.jl)

This Julia package is for to solve mathematical optimisation problems of both unconstrained as well as constrained kind. The package uses the approach named *Charibde* in which it uses two parallel running algorithms *Differential Evolution* and *Interval Branch & Contract algorithm* to achieve solution of any difficult problem.

The package makes these two algorithms run in parallel either on same worker (processor) or on two different workers while maintaining contact with one another through `channels` or `remotechannels`.

The package also allows us to use just *Interval Branch & Contract* algorithm to solve the problems..

The package can also be used through JuMP syntax.


## Documentation
Documentation for the package is available [here](https://yashcodes.github.io/CharibdeOptim.jl/docs/src/index.html).
 
Example of some difficult problems are present in the Example folder.

## Author

- Yashvardhan Sharma, IIT Kharagpur, West Bengal, India

## Mentors

- [David P. Sanders](http://sistemas.fciencias.unam.mx/~dsanders),
Departamento de Física, Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM)

- [Charlie Vanaret](https://cvanaret.wordpress.com/)
Researcher at Fraunhofer ITWM

## References:
- Hybridization of Interval CP and Evolutionary
Algorithms for Optimizing Difficult Problems
by Charlie Vanaret, Jean-Baptiste Gotteland, Nicolas Durand, Jean-Marc Alliot [Charibde - Archive ouverte HAL](https://hal.archives-ouvertes.fr/hal-01168096/document)
- Hybridisataion of evolutionary algorithms and methods of intervals for the optimization of difficult problems [DOCTORATE OF THE UNIVERSITY OF TOULOUSE](http://ethesis.inp-toulouse.fr/archive/00002966/01/vanaret.pdf) by Charlie Vanaret
