CHANGES in v0.7.4

* set find_zero(s) to specialize on the function, fzero(s) to not. (#148)
* adjust Steffensen method logic to take secant step or steffensen
  step, rather than modified steffensen step. Seems to improve
  robustness. (#147)
* add Schroder method (order 2 for multiplicity with derivative), King (1B)
  (superlinear for multiplicity, no derivative), Esser (2B) (order 2
  for multipicity, no derivative) (#143, #147)
* close issue #143 by allowing fns to Newton, Halley to compute f, f/fp, fp/fpp
* add `newton` function to simple.jl
* change find_zeros to identify zeros on [a,b], not (a,b). Closes #141.
* bug fix: issue with quad step after a truncated M-step in find_zero(M,N,...)
* bug fix: verbose argument for Bisection method (#139)
* bug fix: unintentional widening of types in initial secant step (#139)

CHANGES in v0.7.3

* fix bug with find_zeros and Float32
* speeds up bisection function

CHANGES in v0.7.2

* speed up bisection

CHANGES in v0.7.1

* refactor Bisection method to reduce conditional checks inside loop

* took algorithm from Order0, and made it an alternative for find_zero allowing other non-bracketing methods to be more robust

* In FalsePosition there is a parameter to adjust when a bisection step should be used. This was changed in v0.7.0, the old value is restored. (This method is too sensitive to this parameter. It is recommended that either A42 or AlefeldPotraShi be used as alternatives.
