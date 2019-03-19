#!/usr/bin/env python3
###############################################################################
# ampliviz: Visualizes and tries to simplify AmpliconArchitect graphs.
# Copyright (C) 2019 Marcus Fedarko
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
###############################################################################
import click
import math
import os
import numpy


# These don't really correspond to anything in particular, they're just
# aesthetically pleasing options that I thought would look nice
# Use of the first three colors was definitely based on my prior experience
# working with Graphviz colors on MetagenomeScope
COLORS = ["salmon", "cornflowerblue", "limegreen", "burlywood", "orchid",
          "navajowhite1", "yellow", "honeydew2"]
CHROMOSOME2COLOR = {}
NUM_CHROMOSOMES_SEEN = 0

EDGE_TYPE2COLOR = {"concordant": "blue", "discordant": "red", "source":
                   "black"}


class GenomicPosition(object):

    def __init__(self, chromosome, pos, strand):

        global COLORS, CHROMOSOME2COLOR, NUM_CHROMOSOMES_SEEN

        self.chromosome = chromosome
        # small hack to assign a unique color to each distinct chromosome
        # note that this is dependent on COLORS having enough values to account
        # for all chromosomes in the input file
        if self.chromosome not in CHROMOSOME2COLOR:
            CHROMOSOME2COLOR[self.chromosome] = COLORS[NUM_CHROMOSOMES_SEEN]
            NUM_CHROMOSOMES_SEEN += 1

        self.pos = pos
        self.strand = strand

    def __repr__(self):
        return "{}:{}{}".format(self.chromosome, self.pos, self.strand)

    def __eq__(self, other):
        c1 = self.chromosome == other.chromosome
        c2 = self.pos == other.pos
        c3 = self.strand == other.strand
        return (c1 and c2 and c3)

    def __hash__(self):
        # NOTE that this only takes into account chromosome and position -- it
        # doesn't account for the strand
        # Based on https://stackoverflow.com/a/2909119
        return hash((self.chromosome, self.pos, self.strand))

    @staticmethod
    def get_pos_details(chromosome_position):
        """Returns the chromosme, position, and strand associated with a
           position as listed in an AmpliconArchitect breakpoint graph file.

           >>> gp1 = GenomicPosition.get_pos_details('chr1:1234-')
           >>> gp1.chromosome
           'chr1'
           >>> gp1.pos
           1234
           >>> gp1.strand
           '-'
           >>> gp2 = GenomicPosition.get_pos_details('super_cool_chr:333+')
           >>> gp2.chromosome
           'super_cool_chr'
           >>> gp2.pos
           333
           >>> gp2.strand
           '+'
        """
        # chromosome should be of the format "chr4" or something similar
        chromosome, pos = chromosome_position.split(":")
        # strand is either + or -.
        strand = pos[-1]
        # pos is just the integer position on the chromosome (with the trailing
        # strand character removed).
        pos = int(pos[:-1])
        return GenomicPosition(chromosome, pos, strand)


class SequenceEdge(object):

    def __init__(self, start_pos, end_pos, pcc, cov, size, num_mapped_reads):
        """Initializes a node based on a sequence line in the breakpoint graph.

           Note that all of the arguments of this function are assumed to be
           strings, since they've been read from a file.
        """

        # Get position details
        self.start_pos = GenomicPosition.get_pos_details(start_pos)
        self.end_pos = GenomicPosition.get_pos_details(end_pos)

        # A node shouldn't span multiple chromosomes. If it does, this is a
        # problem.
        assert self.start_pos.chromosome == self.end_pos.chromosome

        # Get size; assert that it matches with the positional data
        self.size = int(size)
        assert self.size == self.end_pos.pos - self.start_pos.pos

        self.pcc = float(pcc)
        self.coverage = float(cov)
        self.num_mapped_reads = int(num_mapped_reads)

    def size_to_height(self):
        """This should ostensibly generate a "height" value for this node that
           can be used in a DOT file for this graph. However, in practice this
           ended up producing ugly drawings and didn't really help a ton with
           analysis. So this function is currently unused -- it's left here for
           reference.
        """

        # log(1) is 0, and log(x) for x < 1 will be either negative or
        # undefined. So we just use 0.5 as a minimum height.
        if self.size <= 1:
            return 0.5
        initial_log = math.sqrt(math.log(self.size))
        return max(initial_log, 0.5)

    def dot_repr(self, with_attrs=False, color=False):
        """Generates a representation of this node that can be used in a DOT
           file. Set with_attrs to True if this is for a node's own line in the
           DOT file; keep with_attrs as False this is for just mentioning the
           node in an edge line.

           If color is True, this will color the node based on its chromosome.
           Passing color as True if with_attrs is False won't do anything.

           >>> se = SequenceEdge("chr4:2660-", "chr4:9030+", "14.07", \
                                 "200.54321", "6370", "1234")
           >>> se.dot_repr()
           '"chr4:2660--9030"'
           >>> se.dot_repr(with_attrs=True)
           '"chr4:2660--9030" [shape="rect"]'
           >>> # See https://stackoverflow.com/a/13395612
           >>> se.dot_repr(with_attrs=True, color=True)
           '"chr4:2660--9030" [shape="rect", style="filled", \
fillcolor="limegreen"]'
           >>> se.dot_repr(with_attrs=False, color=True)
           '"chr4:2660--9030"'
        """
        attrs = ""
        if with_attrs:
            attrs = []
            attrs.append("shape=\"rect\"")
            if color:
                attrs.append("style=\"filled\"")
                attrs.append("fillcolor=\"{}\"".format(
                    CHROMOSOME2COLOR[self.start_pos.chromosome]
                ))
            attrs = " " + str(attrs).replace("'", "")
        return '"{}:{}--{}"{}'.format(
            self.start_pos.chromosome,
            self.start_pos.pos,
            self.end_pos.pos,
            attrs
        )

    def __repr__(self):
        return "SequenceEdge({}:{}{} to {}:{}{})".format(
            self.start_pos.chromosome, self.start_pos.pos,
            self.start_pos.strand, self.start_pos.chromosome, self.end_pos.pos,
            self.end_pos.strand
        )


class BreakpointEdge(object):

    def __init__(self, edge_type, positions, pcc, num_read_pairs,
                 hom_size, hom_seq):

        self.edge_type = edge_type

        start_pos, end_pos = positions.split("->")
        self.start_pos = GenomicPosition.get_pos_details(start_pos)
        self.end_pos = GenomicPosition.get_pos_details(end_pos)

        # should be updated to point to a corresponding SequenceEdge object
        self.start_node = None
        self.end_node = None
        # should be updated based on what sides of nodes this edge connects to
        self.headport = None
        self.tailport = None

        # will be updated based on analysis of PacBio reads
        self.support_value = 0

        self.pcc = float(pcc)
        # I've seen some graphs where this is a floating-point number, so for
        # safety's sake we just treat this as a float. (We don't use this for
        # anything yet so this isn't a huge deal)
        self.num_read_pairs = float(num_read_pairs)

        # going to be strings
        self.hom_size = hom_size
        self.hom_seq = hom_seq

    def has_matched_nodes(self):
        """Returns True if this edge has its start and end nodes defined."""
        return self.start_node is not None and self.end_node is not None

    def dot_repr(self, color=False, pacbio=False):
        """Requires self.start_node and self.end_node to be defined."""
        assert self.has_matched_nodes()
        attrs = []
        attrs.append("headport=\"{}\"".format(self.headport))
        attrs.append("tailport=\"{}\"".format(self.tailport))
        if color:
            attrs.append("color=\"{}\"".format(
                EDGE_TYPE2COLOR[self.edge_type]
            ))
        if pacbio:
            attrs.append("penwidth={}".format(self.support_value))
        attrs = str(attrs).replace("'", "")
        return "{}->{} {}".format(self.start_node.dot_repr(),
                                  self.end_node.dot_repr(), attrs)

    def __repr__(self):
        return "{} -- {}:{}{} -> {}:{}{}".format(
            self.edge_type, self.start_pos.chromosome, self.start_pos.pos,
            self.start_pos.strand, self.end_pos.chromosome, self.end_pos.pos,
            self.end_pos.strand
        )


@click.command()
@click.option("-i", "--input-graph", required=True, help="AmpliconArchitect "
              "breakpoint graph. The filename for this should be of the "
              "format {out}_amplicon{id}_graph.txt, per the AmpliconArchitect "
              "README.")
@click.option("-o", "--output-prefix", required=True,
              help="Output prefix for the various files generated.")
@click.option("-d", "--output-directory", required=False, help="Directory to "
              "write output files to. If not passed, files will be written to "
              "the current working directory.")
@click.option("-sm", "--simulated-mean", required=False, default=5, type=float,
              help="Mean of the normal distribution to use when simulating "
              "edge support values.")
@click.option("-ss", "--simulated-std-dev", required=False, default=5,
              type=float, help="Standard deviation of the normal distribution "
              "to use when simulating edge support values.")
@click.option("-p", "--pacbio-bam", required=False, help="BAM file describing "
              "the alignment between PacBio long reads and the sequences in "
              "the breakpoint graph. If this is not passed, simulated edge "
              "support values will be generated and used. NOTE THAT THIS "
              "OPTION IS NOT SUPPORTED YET -- only simulated edge support "
              "values can be generated at present.")
def convert_graph(input_graph: str, output_prefix: str, output_directory: str,
                  simulated_mean: float, simulated_std_dev: float,
                  pacbio_bam: str) -> None:
    with open(input_graph, 'r') as input_graph_file:
        on_first_line = True
        # Collections of relevant objects
        seq_edges = []
        bp_edges = []
        start_pos_to_node = {}
        end_pos_to_node = {}
        for line in input_graph_file:
            if on_first_line:
                on_first_line = False
                continue
            else:
                line_parts = line.split("\t")

                if line_parts[0] == "sequence":
                    n = SequenceEdge(*line_parts[1:])
                    seq_edges.append(n)
                    start_pos_to_node[(n.start_pos)] = n
                    end_pos_to_node[(n.end_pos)] = n

                elif line_parts[0] in ["concordant", "discordant", "source"]:
                    e = BreakpointEdge(*line_parts)
                    bp_edges.append(e)
    # end "with"
    #
    # Now, we have seq_edges, bp_edges, start_pos_to_node, and end_pos_to_node:
    # so we have the entire graph structure ready. We can convert this graph
    # to DOT format.
    #
    # 0. Create a new seq edge for the source. (we can do this later -- for
    # now ignoring source edges is totally ok imo)
    #
    # 1. Go through each bp edge. Using seq_edges, convert the bp edge to a
    # connection between two nodes (based on the start and end
    # coordinate, which -- based on the strand -- should match up to either
    # the start or end of a node).
    #
    # 2. Create a DOT file. for a node, dot_repr() should be good BUT we'd
    # want it to use the "size" (just take the log10 of size or something).
    # for an edge, just get the node dot IDs and use those. Optionally can
    # set weight/etc based on attrs but non essential.
    #
    # 3. at this point we can start messing with pacbio read mapping

    # Assert that no GenomicPosition objects are used as both the start
    # node and end node of a BreakpointEdge
    n_shared = set(start_pos_to_node.keys()) & set(end_pos_to_node.keys())
    assert len(n_shared) == 0
    # Assert that there's one start position and one end position for each
    # node
    n_union = set(start_pos_to_node.keys()) | set(end_pos_to_node.keys())
    assert len(n_union) == 2 * len(seq_edges)

    for e in bp_edges:
        # NOTE: so the "tail" of an edge is actually its beginning. I uh
        # definitely thought is was the other way around, which is why i
        # just spent a non-insignificant amount of time messing with dot
        # to get it to work
        if e.start_pos in start_pos_to_node:
            e.start_node = start_pos_to_node[e.start_pos]
            e.tailport = "w"
        elif e.start_pos in end_pos_to_node:
            e.start_node = end_pos_to_node[e.start_pos]
            e.tailport = "e"

        if e.end_pos in end_pos_to_node:
            e.end_node = end_pos_to_node[e.end_pos]
            e.headport = "e"
        elif e.end_pos in start_pos_to_node:
            e.end_node = start_pos_to_node[e.end_pos]
            e.headport = "w"

    # Filter out edges to the "source" vertex from bp_edges for now
    bp_edges_no_src = filter(lambda e: e.has_matched_nodes(), bp_edges)
    bp_edges_list = list(bp_edges_no_src)

    def make_dot_graph(color_n=False, color_e=False, pacbio_e=False,
                       esv=None, split_edges=False):
        """Returns a str of DOT output for the breakpoint graph, with differing
           levels of detail based on the options passed.

        Parameters
        ----------
        color_n: If True, this will color nodes by their reference chromosome.

        color_e: If True, this will color edges by their "type" (concordant vs.
                 discordant vs. source).

        pacbio_e: If True, this will modify edges' thickness by their
                  "support_value" attribute (either the existing support_value
                  attribute or by a support value from the esv list if
                  specified). Edges with a support_value of <= 0 will not be
                  drawn, since they are "unsupported" by the long read data.

        esv: If passed, this should be a list of simulated "support" values
             such that len(esv) >= len(bp_edges_list). This function will
             overwrite the support_value attribute of each edge in
             bp_edges_list with the corresponding value in esv. So only use
             this argument if you want to simulate the effect of PacBio reads
             on the graph visualization. If this is passed, pacbio_e should
             also be passed. (Passing esv but not pacbio_e does nothing.)

        split_edges: If True, this will remove edges with a support value of
                     <= 0 from the graph. (As with esv, this only does
                     something if pacbio_e is passed.)
        """
        dot_output = "digraph bp_graph {\n\trankdir=\"LR\"\n"
        for n in seq_edges:
            node_dot_str = n.dot_repr(with_attrs=True, color=color_n)
            dot_output += "\t{}\n".format(node_dot_str)
        i = 0
        for e in bp_edges_list:
            if pacbio_e:
                # If the user passed in simulated support values, use them
                if esv is not None:
                    e.support_value = esv[i]
                # Don't draw unsupported edges if split_edges is True.
                # (If split_edges is False, these edges will still be
                # incorporated in the layout process.)
                if split_edges and e.support_value <= 0:
                    i += 1
                    continue
            edge_dot_str = e.dot_repr(color=color_e, pacbio=pacbio_e)
            dot_output += "\t{}\n".format(edge_dot_str)
            i += 1
        dot_output += "}"
        return dot_output

    # Get ready to start creating output files
    if output_directory is None:
        output_directory = os.getcwd()
        print('Creating files in current working directory (-d not passed)...')
    else:
        os.makedirs(output_directory, exist_ok=True)
        print('Creating files in output directory {}...'.format(
            output_directory
        ))

    p0 = "{}_0.gv".format(output_prefix)
    p1 = "{}_1.gv".format(output_prefix)
    p2 = "{}_2.gv".format(output_prefix)
    p3 = "{}_3.gv".format(output_prefix)
    p4 = "{}_4.gv".format(output_prefix)
    fn0 = os.path.join(output_directory, p0)
    fn1 = os.path.join(output_directory, p1)
    fn2 = os.path.join(output_directory, p2)
    fn3 = os.path.join(output_directory, p3)
    fn4 = os.path.join(output_directory, p4)

    # Create basic visualizations
    with open(fn0, 'w') as out_file:
        out_file.write(make_dot_graph())
    print('Basic graph ({}) created.'.format(p0))

    with open(fn1, 'w') as out_file:
        out_file.write(make_dot_graph(color_n=True))
    print('Graph with colorized nodes ({}) created.'.format(p1))

    with open(fn2, 'w') as out_file:
        out_file.write(make_dot_graph(color_n=True, color_e=True))
    print('Graph with colorized nodes and edges ({}) created.'.format(p2))

    if pacbio_bam is not None:
        # Below code would be used for working with real long read data
        # edge_support_values = read_pacbio_bam(pacbio_bam, bp_edges_list)
        # for e in bp_edges_list:
        #     e.support_value = edge_support_values[e]
        #     print(e, e.support_value)
        raise NotImplementedError("The -p option isn't supported yet, sorry.")
    else:
        edge_support_values = []
        for i in range(len(bp_edges_list)):
            # The max(..., 0) is just a way of saying: if this edge gets
            # assigned a support value of <= 0, just set it to 0.
            v = max(round(numpy.random.normal(
                loc=simulated_mean, scale=simulated_std_dev
            )), 0)
            edge_support_values.append(v)

    with open(fn3, 'w') as out_file:
        out_file.write(make_dot_graph(color_n=True, color_e=True,
                                      pacbio_e=True, esv=edge_support_values))
    print('Graph with colorization and edge support ({}) created.'.format(p3))

    with open(fn4, 'w') as out_file:
        out_file.write(make_dot_graph(color_n=True, color_e=True,
                                      pacbio_e=True, esv=edge_support_values,
                                      split_edges=True))
    print('Graph with colorization, edge support, and splitting ({}) '
          'created.'.format(p4))


if __name__ == '__main__':
    convert_graph()
