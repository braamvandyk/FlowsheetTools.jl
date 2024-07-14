# Flowsheet Tools

Basic Julia infrastructure for flowsheets.

1. Components
2. Streams
3. UnitOps
4. Mass balance boundaries
5. KPIs
6. Mass and element balance reconciliation

This is VERY MUCH a work in progress. If you want to collaborate, please let me know, so
we can discuss ideas and do this together.

## Current status

The library is finally in package form and complete to a usable point. The [demos.jl / demos.pynb](https://github.com/braamvandyk/FlowsheetTools.jl/blob/main/demos.ipynb)
file shows the current capabilities. TODO.md lists the things I currently want to do.

## The backstory

The project started off as an exercise in metaprogramming, which resulted in the macros to define
components, and then later streams. And then the thing grew a life of its own.

The vision was to build a basic toolkit to define flowsheets in order to calculate KPIs
and reconcile mass balances. This is mostly there, expect that the mass balance reconciliation
algorithms are a bit basic. More sophisticated methods are available and can be added.

In the process I thought it may also be useful to be able to define a simple function to
calculate the product streams from the feed streams for each unit operation. The idea is not
to build a fully fledged process simulator, but often I need to e.g. fit kinetics on piloting
data for non-trivial reactor configurations and this type of library could come in handy.

Of course, nothing stops us from later adding physical properties, energy balances etc, but
is there REALLY a need for that? After all, we already have Aspen, Pro/II etc.

After the idea came to add "active" unit ops, it occured to me that it would be a really good
idea to have all of the components and streams in some form of global structure, rather than
a scattering of loose variables. This required a significant rewrite, followed by another and then
yet another, as I refined the user interface. We are now at the point where you just define a flowsheet
container, which holds all of the components, streams, unit operations, mass balance boundaries etc.
There is always room for improvement, but at least I now feel that I can show this to people who may
want to use the tool.

Now, you at least have the option to iterate through all components and streams, should you
need to do that, and everything in the flowsheet is nice and consistent - one list of components, 
used in one list of streams, one lost of unit operations etc.

I am sure that I'll keep on tinkering with this for some time still. At the very least to add more
reconciliation algorithms, and maybe a few more unit operations. I am VERY open to ideas and collaborators.

### A note on the reconciliation algorithm

There are many algorithms that will modify the stream compositions to get a perfect closure. I am NOT a fan...
In my work, we deal with complex streams from Fischer-Tropsch Syntheses. This means potentially hundreds of 
components, but practically at least 60+. These stream compositions have several constraints - they must obey
an Anderson-Schultz-Flory distribution, certain ratios of paraffin : olefin : carbonyl : alcohol : acid : ... etc.
Opening up the stream compositions to modification not only make the problem extremely high-dimensional, but, is
extremely likely to break all of these constrains. Better, thus, to only modify the stream flows and do the best
you can with the available analyses.

This is why the existing algorithm looks the way it does. Adding any other algorithm is of course not a problem,
but also not a priority. If anyone feels they really need some other methods, PRs would be more than welcome.
