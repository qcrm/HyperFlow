# HyperFlow

![Type3 shock-expansion interaction density gradient](https://github.com/qcrm/HyperFlow/blob/master/media/readme-screenshot.png?raw=true)

HyperFlow is a C++11-based open source modular research framework for the computation of two-dimensional numerical solutions to [hyperbolic conservation laws](https://en.wikipedia.org/wiki/Hyperbolic_partial_differential_equation#Hyperbolic_system_and_conservation_laws).

It implements various [high-resolution schemes](https://en.wikipedia.org/wiki/High-resolution_scheme) via the [finite volume method](https://en.wikipedia.org/wiki/Finite_volume_method). The [method of lines](https://en.wikipedia.org/wiki/Method_of_lines) is used to separate spatial and temporal derivative discretisation.

HyperFlow uses a block-structured axis-aligned Cartesian grid for spatial domain discretisation, allowing reasonably complex domains to be studied.

Its primary application is conducting experiments for supersonic aerodynamics via numerical investigation of the compressible inviscid Euler equations.

# Development Status

| Development   | Details       |
| ------------- | ------------- |
| License       | ![GitHub](https://img.shields.io/github/license/qcrm/HyperFlow?style=flat-square&label=License) |      

HyperFlow is developed by QCRM Ltd.

# Schemes Implemented

- *Gounov First-Order Upwind:* The well-known first-order scheme of Godunov for non-linear hyperbolic systems. The implementation currently makes use of Toro's HLLC approximate Riemann solver with forward Euler time integration in a dimensionally unsplit manner.

- *MUSCL-Hancock:* A second-order 'predictor/corrector' [total-variational diminishing](https://en.wikipedia.org/wiki/Total_variation_diminishing) (TVD) [MUSCL](https://en.wikipedia.org/wiki/MUSCL_scheme)-based scheme of Van Leer and Hancock, using MINBEE or the Van Leer slope limiter to eliminate spurious oscillations. The implementation uses Toro's HLLC approximate Riemann solver with forward Euler time integration in a dimensionally unsplit manner.

# Media

Various video examples of HyperFlow applied to the compressible Euler equations (in both 2D and 3D) can be found at the [QCRM YouTube channel](https://www.youtube.com/channel/UCJrlFNmxBaFf6Ul78bEp4xg/videos).

# Copyright

Copyright (c) 2019-2020 QCRM Ltd
