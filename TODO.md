# Things to do

1. [ ] Add handling for missing data in stream histories. Just do dropmissing() (preferred), or replace with average, or interpolate between points on either side?
2. [ ] Add proper mass balance recon algorithm, or at least use ridge regression. Step detection is already coded, but need ramp detection too.
3. [ ] Should stream histories detect the steps or should this be done externally and the periods between steps/ramps etc fed to a stream history? And then reconciled seperately? Best place to do is probably in the recon function. Need to split stream histories into separate variables for each steady state period, as each period could have different corrections.
4. [ ] We need to determine uncertainty in flowrates to assign weights to the errors. Streams oscilating due to a control input will look like there is high uncertainty, when there is not. How do we detect this? FFT and a high pass filter?
5. [ ] Create Plots/Makie recipes for plotting stuff
6. [X] Add KPIs for histories
7. [ ] Check for missing pretty printing (like on lists...)
8. [ ] Print stream tables for a boundary (PrettyTables?)
9. [ ] Make order of parameters to macros consistent (how did they get so mixed-up anyway?)
10. [ ] Macro to define boundaries
11. [ ] Now that we have active UnitOps, we need a Flowsheet type to manage the whole systems calculations
12. [ ] Define unit tests for CI at some future point
