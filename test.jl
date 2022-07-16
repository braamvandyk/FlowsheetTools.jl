#----Testing-----------------

# Define components via constructor
ethylene = Component("C2H4", ["C", "H"], [2, 4])

# or via macro
hydrogen = @comp begin
    H --> 2
end "H2"

ethane = @comp begin
    C --> 2
    H --> 6
end "Ethane"


# Define stream via constructor
feedcomps = [ethylene, hydrogen]
feedflows = [ethylene.Mr, 2*hydrogen.Mr]
feed = Stream("Feed", feedcomps, feedflows)

# or via macro
prod = @stream begin
    ethylene --> ethylene.Mr*0.1
    ethane --> ethane.Mr*0.9
    hydrogen --> hydrogen.Mr*1.1
end "Product"

# You get the same result
C2c = Stream("C2", [ethylene, ethane], [ethylene.Mr*0.1, ethane.Mr*0.9])
C2m = @stream begin
    ethylene --> ethylene.Mr*0.1
    ethane --> ethane.Mr*0.9
end "C2"
C2c.comps ==C2m.comps
C2c.massflows == C2m.massflows
C2c.moleflows == C2m.moleflows

# Copy a stream to change its name
C2cc = copystream(C2c, "mystream")
# Different name, so not identical
!(C2cc == C2c)
# But all the flows are the same
C2c.atomflows == C2cc.atomflows

# Multiply a stream to change all its flows by a scalar factor
moreC2 = 2.0*C2c
moreC2.totalmassflow ≈ 2.0*C2c.totalmassflow

# Test unit ops and mass balance boundaries
# Define some streams
H2 = Stream("H2", [hydrogen], [hydrogen.Mr*1.1])
C2 = Stream("C2", [ethylene, ethane], [ethylene.Mr*0.1, ethane.Mr*0.9])

# Define some unit ops
reactor = UnitOp("RX101", [feed], [prod])
membrane = UnitOp("MX101", [prod], [C2, H2])

# Define a mass balance boundary
b = BalanceBoundary([reactor, membrane])

# Check the closures
b.atomclosures
b.closure
b.total_in.totalmassflow
b.total_out.totalmassflow

# And some defined KPIs
conversion(b, ethane)
conversion(b, ethylene)
selectivity(b, ethylene, ethane)

# Now test the reconciliations
# Change some flows to simulate measurement errors
feed2 = copystream(0.95*feed, "Feed2")
prod2 = copystream(1.01*prod, "Prod2")

# Define the unit ops and boundary
reactor2 = UnitOp("RX101", [feed2], [prod2])
membrane2 = UnitOp("MX101", [prod2], [C2, H2])
b2 = BalanceBoundary([reactor2, membrane2])

# Get the correction factors on the inlets and outlets
corrections = closemb(b2)

# Test the pretty printing
print(ethane)
print(feed)
print(membrane)
print(b)

# Read in a list of components
comps, count = readcomponentlist(["Ethylene", "Ethane", "hydrogen"])
count == 3

# Read in some stream histories
feedhist = readstreamhistory("FeedStream", comps)
prodhist = readstreamhistory("ProdStream", comps)

# And mix them
combhist = feedhist + prodhist

# And scale the result
combhist = 2.0 * combhist

# Test UnitOpHistory
RX101History = UnitOpHistory("RX101", [feedhist], [prodhist])

# Test BalanceBoundaryHistory
bh = BalanceBoundaryHistory([RX101History])
corrections = closemb(bh)