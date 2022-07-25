# Things to do

1. We need to be able to keep track of streams to apply corrections etc. Proposal: rather than create streams as variables, force them to be entries in an array? Then ID them by index?
2. Add handling for missing data in stream histories
3. Pull the date/time field into StreamHistories for plots etc
4. Add proper mass balance recon algorithm, or at least use ridge regression
5. Should stream histories detect the steps or should this be done externally and the periods between steps/ramps etc fed to a stream history? And them reconciled seperately?
6. Wrap the code into a package
7. Create Active and Passive UnitOps - Passive is the current one. Active takes a transformation function to generate the outlet streams for the inlets. Not aiming to make Aspen.jl just yet, but calling simple reactors etc would be useful.
8. Switch to StaticArrays?
