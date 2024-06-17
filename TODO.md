# Things to do

 ## YES
 - Change the setup so you start by defining a Flowsheet, which automatically creates a ComponentList, StreamList, UnitOpList and BoundaryList. Then we just pass the Flowsheet object everywhere.
 
 - Define unit tests for CI at some future point
 - Don't limit the length of pretty table outputs.

 ### Data clean-up
 - Should be a separate step, before reconciliation, hence separate package inside the main package
 - Add change point detection and stationarity checks to isolate (pseudo)steady states
