# ampliviz
[![Build Status](https://travis-ci.org/fedarko/ampliviz.svg?branch=master)](https://travis-ci.org/fedarko/ampliviz)

Visualizes breakpoint graphs generated by
[AmpliconArchitect](https://github.com/virajbdeshpande/AmpliconArchitect).

Supports generating simulated edge support values (which would hypothetically
have arisen from long-read data) and using these to adjust the thickness of
edges in the graph visualization. "Unsupported" edges are removed from the graph
visualization, which can simplify it drastically.

This is my class project for [CSE 280A](http://proteomics.ucsd.edu/vbafna/teaching-2/cse280a-algorithms-for-genetics/) (Winter 2019, UC San Diego) with Prof. Vineet Bafna.

## Installation

```bash
git clone https://github.com/fedarko/ampliviz
pip install click numpy
```

## Usage

This assumes that the ampliviz directory is on your `PATH`. You can add the
directory to your `PATH` if you like (then you can just run `ampliviz.py` from
anywhere), or you can just `cd` to the ampliviz directory and then run it using
`./ampliviz.py`. Either way is fine.

```
Usage: ampliviz.py [OPTIONS]

Options:
  -i, --input-graph TEXT       AmpliconArchitect breakpoint graph. The
                               filename for this should be of the format
                               {out}_amplicon{id}_graph.txt, per the
                               AmpliconArchitect README.  [required]
  -o, --output-prefix TEXT     Output prefix for the various files generated.
                               [required]
  -d, --output-directory TEXT  Directory to write output files to. If not
                               passed, files will be written to the current
                               working directory.
  -p, --pacbio-bam TEXT        BAM file describing the alignment between
                               PacBio long reads and the sequences in the
                               breakpoint graph. If this is not passed,
                               simulated edge support values will be generated
                               and used. NOTE THAT THIS OPTION IS NOT
                               SUPPORTED YET -- only simulated edge support
                               values can be generated at present.
  --help                       Show this message and exit.
```

### Using Graphviz to draw generated DOT files

After running ampliviz, it should produce four `.gv` files. These all describe
the breakpoint graph at varying levels of detail. In particular:
  - `{output_prefix}_0.gv` describes the nodes and edges in the graph with
    minimal styling.
  - `{output_prefix}_1.gv` is like the previous graph, but nodes are colored
    by their assigned chromosome.
  - `{output_prefix}_2.gv` is like the previous graph, but edges are colored
    by their "type" (either concordant or discordant).
  - `{output_prefix}_3.gv` is like the previous graph, but the thicknesses of
    edges are altered based on randomly simulated "edge support" values. Edges
    assigned a simulated support value of less than or equal to zero do not
    have their lines drawn, but their arrowhead is still drawn and the edge still
    impacts the layout of the graph.
  - `{output_prefix}_4.gv` is like the previous graph, but edges assigned a
    simulated support value of less than or equal to zero are removed
    from the graph, which can simplify the graph topology and resulting
    visualization significantly.

All of these files can be visualized in [Graphviz](https://www.graphviz.org/)
using the "dot" layout program. If you use ampliviz' `Makefile` to generate
these files, it will automatically try to lay them out as PNG images using dot
(and then open the resulting images).

### Determining the number of connected components in `{output_prefix}_4.gv`

(This requires that Graphviz be installed, in order for its `ccomps` tool to be
available.)

```bash
ccomps {output_prefix}_4.gv | grep -r "subgraph" | wc -l
```

### Caveats

0. This isn't really a "caveat," per se, but the `Makefile` relies on you
   having Graphviz (or at least dot) installed on your system, as well as some
   sort of program called `open` that knows how to open images. If either or
   both of these things are not installed, the Makefile will crash. (As
   mentioned above, you can totally just run ampliviz by itself -- however, the
   bulk of what it does is produce visualizations, so if you don't have dot
   installed then it won't be very useful.)

1. If your AA graph spans more than eight reference chromosomes, ampliviz will run
   out of colors to use for colorizing nodes (for `{output_prefix}_1.gv`, as
   mentioned above). You can fix this very easily by editing the `COLORS` list in
   `ampliviz.py` to include more colors (see [here](http://www.graphviz.org/doc/info/colors.html) for a list of all usable color names). Feel free to change these colors around even if you don't need more -- I'd advise you to pick relatively contrasting colors, but it's completely up to you.
   In the future, I might write some code to autogenerate random colors so that an arbitrary amount of chromosomes can be used.

2. This currently filters out `source` edges in the breakpoint graph, so only
   edges classified as `concordant` and `discordant` will be shown.

## Screenshots

Coming soon (need to find a publicly available AA dataset).

## License

This code is released under the [GNU General Public License (GPL), version 3](https://www.gnu.org/licenses/gpl-3.0.en.html).

## Acknowledgements

Thanks to Jens Luebeck, Nam-Phuong Nguyen, Mehrdad Baktiari, and Prof. Bafna
for advice and help with getting started.
