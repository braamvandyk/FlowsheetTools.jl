# Things to do

1. Add handling for missing data in stream histories
2. Add proper mass balance recon algorithm, or at least use ridge regression
3. Should stream histories detect the steps or should this be done externally and the periods between steps/ramps etc fed to a stream history? And them reconciled seperately?
4. Wrap the code into a package
5. Create Active and Passive UnitOps - Passive is the current one. Active takes a transformation function to generate the outlet streams for the inlets. Not aiming to make Aspen.jl just yet, but calling simple reactors etc would be useful.
6. We need to determine uncertainty in flowrates to assign weights to the errors. Streams oscilating due to a control input will look like there is high uncertainty, when there is not. How do we detect this? FFT and a high pass filter?
7. Create Plots/Makie recipes for plotting stuff
