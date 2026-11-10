"""Parse a graph with CoRAL 3.1.1 and write it back with CoRAL's own writer."""
import sys
from coral.breakpoint.parse_graph import parse_breakpoint_graph
from coral.breakpoint.breakpoint_utilities import output_breakpoint_graph_lr
from coral.datatypes import OutputPCOptions
src, dst = sys.argv[1], sys.argv[2]
with open(src) as fh:
    g = parse_breakpoint_graph(fh)
g.infer_amplicon_intervals(0)
output_breakpoint_graph_lr(g, dst, OutputPCOptions.ALL)
print("intervals", [str(i) for i in g.amplicon_intervals], "seq", len(g.sequence_edges), "conc", len(g.concordant_edges), "disc", len(g.discordant_edges), "pc", len(g.path_constraints))
