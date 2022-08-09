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

The story started off as an exercise in meta programming, which ended in the macros to define
components, and then later streams. And then the thing grew a life of its own.

The vision was to build a basic toolkit to define flowsheets in order to calculate KPIs
and reconcile mass balances.

In the process I thought it may also be useful to be able to define a simple function to
calculate the product streams from the feed streams for each unit operation. The idea is not
to build a fully fledged process simulator, but often I need to e.g. fit kinetics on piloting
data for non-trivial reactor configurations and this type of library could come in handy.

Of course, nothing stops us from later adding physical properties, energy balances etc, but
is there REALLY a need for that? After all, we already have Aspen, Pro/II etc.

After the idea came to add "active" unit ops, it occured to me that it would be a really good
idea to have all of the components and streams in some form of global structure, rather than
a scattering of loose variables. This required a significant rewrite and passing the system
component list to all streams and the system stream list to all unit ops.

Now, you at least have the option to iterate through all components and streams, should you
need to do that. I sincerely hope this proves to be useful, since it required a lot of rework,
including modifying the macros, months after I forgot how I wrote them in the first place :-(

The streams hold an array of the names of the components actually present, so you don't
need a bunch of zero flows for the components not present. The system component list is a Dict,
indexed on the component names, which is fairly flexible. I figured we don't need to squeeze
every last nanosecond out of the speed in this specific area - that will more be the case in the
"active" unit ops.
