# Things to do

 ## YES
 - Add proper mass balance recon algorithm
 - Add change point detection and stationarity checks to isolate (pseudo)steady states
 - We need to determine uncertainty in flowrates to assign weights to the errors. Streams oscilating due to a control input will look like there is high uncertainty, when there is not. How do we detect this? FFT and a Should stream histories detect the steps or should this be done externally and the periods between steps/ramps etc fed to a stream history? And then reconciled seperately? Best place to do is probably in the recon function. Need to split stream histories into separate variables for eac pass filter?
 - Define unit tests for CI at some future point
 - Add docstrings to functions without them
 - Add some form of feedback when a Flowsheet executes -> write a Base.show for Flowsheet
 - Add a precompile workload
 - Add a stoichiometric reactor to unitops
 - Overload Base.== or Base.≈

 ### Data clean-up
  - Should be a separate step, before reconciliation
  
 ### Reconciliation
  - In closemb_simple, include variances for flows in loss function
    - Add a calculated variance for streams with history, or a specified one for point values
        - Add a list of default variances for various flow measurement types and allow this to be specified as an option (vararg)
  - In main algorithm:
    - Streams need a status - has flow, is fixed, or unknown, so add status variable
    - Construct the incidence matrix in the boundary
    - Build x, y, z, and f
    - Build A
    - rref(A)
    - discoverability analysis, with appropriate output for undiscoverable streams
    - reconcilability analysis, with appropriate output is there is no redundancy
    - reconcile, if possible
    - We are NOT including non-linear constraints in this round

 ## MAYBE
 - Create Plots/Makie recipes for plotting stuff - mostly handled through using TimeSeries
