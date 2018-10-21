CHANGES in v0.7.4

* add Schroder method
* close issue #143 by allowing fns to Newton, Halley to compute f, f/fp, fp/fpp
* add `newton` function to simple.jl
* change find_zeros to identify zeros on [a,b], not (a,b). Closes #141.
* bug fix: issue with quad step after a truncated M-step in find_zero(M,N,...)
* bug fix: verbose argument for Bisection method
* bug fix: unintentional widening of types in intial secant step

CHANGES in v0.7.3

* fix bug with find_zeros and Float32
* speeds up bisection function

CHANGES in v0.7.2

* speed up bisection

CHANGES in v0.7.1

* refactor Bisection method to reduce conditional checks inside loop

* took algorithm from Order0, and made it an alternative for find_zero allowing other non-bracketing methods to be more robust

* In FalsePosition there is a parameter to adjust when a bisection step should be used. This was changed in v0.7.0, the old value is restored. (This method is too sensitive to this parameter. It is recommended that either A42 or AlefeldPotraShi be used as alternatives.
