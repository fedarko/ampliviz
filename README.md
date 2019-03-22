# ampliviz
[![Build Status](https://travis-ci.org/fedarko/ampliviz.svg?branch=master)](https://travis-ci.org/fedarko/ampliviz)

Visualizes breakpoint graphs generated by
[AmpliconArchitect](https://github.com/virajbdeshpande/AmpliconArchitect).

Via the `-p` option, supports generating "edge support values" for each edge
based on a sorted and indexed BAM file mapping long reads to the reference
chromosomes of the graph. These support values are used to adjust the thickness
of edges in the graph visualization, and "unsupported" edges are removed from
the visualization entirely (which can simplify it drastically).

Also supports generating simulated edge support values (which would hypothetically
have arisen from long-read data) and using these to impact the visualization.
These simulated values are randomly drawn from a normal distribution (the
parameters for this distribution are configurable via the `-sm` and `-ss`
options).

This is my class project for [CSE 280A](http://proteomics.ucsd.edu/vbafna/teaching-2/cse280a-algorithms-for-genetics/) (Winter 2019, UC San Diego) with Prof. Vineet Bafna.

## Installation

```bash
git clone https://github.com/fedarko/ampliviz
pip install click numpy pysam
```

## Usage

This assumes that the ampliviz directory is on your `PATH`. You can add the
directory to your `PATH` if you like (then you can just run `ampliviz.py` from
anywhere), or you can just `cd` to the ampliviz directory and then run it using
`./ampliviz.py`. Either way is fine.

```
Usage: ampliviz.py [OPTIONS]

Options:
  -i, --input-graph TEXT          AmpliconArchitect breakpoint graph. The
                                  filename for this should be of the format
                                  {out}_amplicon{id}_graph.txt, per the
                                  AmpliconArchitect README.  [required]
  -o, --output-prefix TEXT        Output prefix for the various files
                                  generated.  [required]
  -d, --output-directory TEXT     Directory to write output files to. If not
                                  passed, files will be written to the current
                                  working directory.
  -sm, --simulated-mean FLOAT     Mean of the normal distribution to use when
                                  simulating edge support values.
  -ss, --simulated-std-dev FLOAT  Standard deviation of the normal
                                  distribution to use when simulating edge
                                  support values.
  -p, --pacbio-bam TEXT           BAM file describing the alignment between
                                  PacBio long reads and the sequences in the
                                  breakpoint graph. If this is not passed,
                                  simulated edge support values will be
                                  generated and used.
  -msv, --max-support-value INTEGER
                                  Maximum accepted edge support value. If any
                                  edge support values (either real or
                                  simulated) exceed this value, then every
                                  edge support value will be scaled to have
                                  this value as a maximum. (This is done so
                                  that edge support values can be used
                                  directly as values for Graphviz' "penwidth"
                                  attribute.)
  -dot, --draw-with-dot           If passed, will try to use dot to convert
                                  each DOT file to a PNG drawing. This will
                                  crash if Graphviz isn't installed.
  --help                          Show this message and exit.
```

### Interpreting the output of ampliviz

After running ampliviz, it should produce five `.gv` files. These all describe
the breakpoint graph at varying levels of detail. In particular:

  - `{output_prefix}_0.gv` describes the nodes and edges in the graph with
    minimal styling.
  - `{output_prefix}_1.gv` is like the previous graph, but nodes are colored
    by their assigned chromosome.
  - `{output_prefix}_2.gv` is like the previous graph, but edges are colored
    by their "type" (either concordant or discordant).
  - `{output_prefix}_3.gv` is like the previous graph, but the thicknesses of
    edges are altered based on the "edge support" values (these are randomly
    simulated if `-p` is not passed, and computed from a specified BAM file if
    `-p` is passed). Edges with a support value of less than or equal to zero
    do not have their lines drawn, but their arrowhead is still drawn and the
    edge still impacts the layout of the graph here.
  - `{output_prefix}_4.gv` is like the previous graph, but edges assigned a
    support value of less than or equal to zero are removed entirely
    from the graph, which can simplify the graph topology and resulting
    visualization significantly.

All of these files can be visualized in [Graphviz](https://www.graphviz.org/)
using the "dot" layout program. If you pass the `-dot` flag when running
ampliviz, then ampliviz will automatically call dot to lay out these graphs.

(If you use the `-dot` flag and Graphviz isn't installed on your system, this
will result in an error being raised. You can check if Graphviz is installed by
running `dot -V` -- if it gives you a version number then it's installed, and
if it gives you an error then it's not installed.)

### Caveats

1. If your AA graph spans more than eight reference chromosomes, ampliviz will run
   out of colors to use for colorizing nodes (for `{output_prefix}_1.gv`, as
   mentioned above). You can fix this very easily by editing the `COLORS` list in
   `ampliviz.py` to include more colors (see [here](http://www.graphviz.org/doc/info/colors.html) for a list of all usable color names). Feel free to change these colors around even if you don't need more -- I'd advise you to pick relatively contrasting colors, but it's completely up to you.
   In the future, I might write some code to autogenerate random colors so that an arbitrary amount of chromosomes can be used.

2. This currently filters out `source` edges in the breakpoint graph, so only
   edges classified as `concordant` and `discordant` will be shown.

3. There isn't a lot of test code for this (yet) -- so the standard GPL
   disclaimers about no warranty, etc. apply.

## Screenshots

Coming soon (need to find a publicly available AA dataset).

## License

This code is released under the [GNU General Public License (GPL), version 3](https://www.gnu.org/licenses/gpl-3.0.en.html).

## Acknowledgements

Thanks to Jens Luebeck, Nam-Phuong Nguyen, Mehrdad Bakhtiari, and Prof. Bafna
for advice and help with getting started.

This project relies on a number of Python libraries, including:

  - [click](https://click.palletsprojects.com/en/7.x/)
  - [flake8](http://flake8.pycqa.org/en/latest/) (used when testing with Travis-CI)
  - [numpy](http://www.numpy.org/)
  - [pysam](https://github.com/pysam-developers/pysam/)
