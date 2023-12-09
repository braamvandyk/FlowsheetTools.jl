# Things to do

 ## YES
 - Add proper mass balance recon algorithm
 - Add change point detection and stationarity checks to isolate (pseudo)steady states
 - We need to determine uncertainty in flowrates to assign weights to the errors. Streams oscilating due to a control input will look like there is high uncertainty, when there is not. How do we detect this? FFT and a Should stream histories detect the steps or should this be done externally and the periods between steps/ramps etc fed to a stream history? And then reconciled seperately? Best place to do is probably in the recon function. Need to split stream histories into separate variables for eac pass filter?
 - Define unit tests for CI at some future point
 - Add docstrings to functions without them
 - Add some form of feedback when a Flowsheet executes
 - Add a precompile workload
 - Add a stoichiometric reactor to unitops
 - Overload Base.==
 - Overload Base.≈

 ## MAYBE
 - Create Plots/Makie recipes for plotting stuff - mostly handled through using TimeSeries
 - Add callbacks to mermaid diagram to display streams when clicking on them
 - Reconsider using --> vs ~ in macros
 - Replace @assert and @error with throw?
