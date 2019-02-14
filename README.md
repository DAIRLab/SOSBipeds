# SOSBipeds
Crude port of some old code from the paper "Michael Posa, Twan Koolen, Russ Tedrake. _Balancing and Step Recovery Capturability via Sums-of-Squares Optimization_. Robotics: Science and Systems, 2017." 

Runfiles are:
- LIPM
  - lipm2DNstepCapturability.m (single SOS program for outer approximation)
  - lipm2DNStepAlternations.m (sequence of SOS programs for inner approximation)
- Variable height model
  - lipmHeightVariationCapturability.m (single SOS program for outer approximation)
  - lipmHeightVariationAlternations.m (sequence of SOS programs for inner approximation)
