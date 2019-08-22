# SOSBipeds

## Runfiles
### Inner Approximation
- Without scaling:
    - `lipmSwingLeg1StepAlternations.m`
    - Volume decreases initially and then starts increasing rapidly before becoming infeasible.
- With scaling:
    - `lsl1StepAlternations.m`
    - The problem remains feasible but there is no improvement in the volume of the ROA.

### Outer Approximation
- N-step capturability (outer)
    - `lipmSwingLegNStepCapturability.m`
    - Takes as input the step number `N`  

### Other useful files
- Plot trajectory of
    - `misc/plot_traj.m`
    - Useful for quick testing (but the code is in a crude state)

## Note
The function `W` in inner approximation denotes the 1-step goal region (stepping here will reset the states to 0-step region)
