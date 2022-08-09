# Things to do

1. Add handling for missing data in stream histories. Just do dropmissing() (preferred), or replace with average, or interpolate between points on either side?
2. Add proper mass balance recon algorithm, or at least use ridge regression. Step detection is already coded, but need ramp detection too.
3. Should stream histories detect the steps or should this be done externally and the periods between steps/ramps etc fed to a stream history? And then reconciled seperately? Best place to do is probably in the recon function. Need to split stream histories into separate variables for each steady state period.
4. Wrap the code into a package
5. Create Active and Passive UnitOps - Passive is the current one. Active takes a transformation function to generate the outlet streams for the inlets. Not aiming to make Aspen.jl just yet, but calling simple reactors etc would be useful.
6. We need to determine uncertainty in flowrates to assign weights to the errors. Streams oscilating due to a control input will look like there is high uncertainty, when there is not. How do we detect this? FFT and a high pass filter?
7. Create Plots/Makie recipes for plotting stuff
8. Add KPIs for histories
9. Check for missing pretty printing (like on lists...)
10. NB! NB! NB! Drop CSV and DataFrames as dependencies NB! NB! NB!
11. Print stream tables for a boundary (PrettyTables?)
