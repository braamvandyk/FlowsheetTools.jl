# Things to do

 ## YES
 - Add stroichiometry check to reaction constructor
 - Define unit tests for CI at some future point
 - Don't limit the length of pretty table outputs.

 ### Data clean-up
 - Should be a separate step, before reconciliation, hence separate package
 - Add change point detection and stationarity checks to isolate (pseudo)steady states
  
 ### Reconciliation alternative algorithm
  - In closemb_simple, include variances for flows in loss function
    - Add a calculated variance for streams with history, or a specified one for point values
        - Add a list of default variances for various flow measurement types and allow this to be specified as an option (vararg)
  - In main algorithm:
    - Construct the incidence matrix in the boundary
    - Build x, y, z, and f
    - Build A
    - rref(A)
    - discoverability analysis, with appropriate output for undiscoverable streams
    - reconcilability analysis, with appropriate output is there is no redundancy
    - reconcile, if possible
    - We are NOT including non-linear constraints in this round

