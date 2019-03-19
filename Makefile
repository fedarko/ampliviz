# You can specify these three variables from the command-line to use arbitrary
# AmpliconArchitect files as inputs: e.g.
# make INGRAPH=your_amplicon1_graph.txt OUTDIR=your_outdir OUTPREF=your_prefix
#
# Currently INGRAPH doesn't point to an actual AA graph. Once I get clearance
# to use a certain AA graph in public, I'll upload it to the repo and use it in
# the tests/Makefile/etc.
INGRAPH = data/test_amplicon1_graph.txt
OUTPREF = test
OUTDIR = data/test
all:
	./ampliviz.py -i $(INGRAPH) -o $(OUTPREF) -d $(OUTDIR)
	dot -Tpng $(OUTDIR)/$(OUTPREF)_0.gv > $(OUTDIR)/$(OUTPREF)_0.png
	dot -Tpng $(OUTDIR)/$(OUTPREF)_1.gv > $(OUTDIR)/$(OUTPREF)_1.png
	dot -Tpng $(OUTDIR)/$(OUTPREF)_2.gv > $(OUTDIR)/$(OUTPREF)_2.png
	dot -Tpng $(OUTDIR)/$(OUTPREF)_3.gv > $(OUTDIR)/$(OUTPREF)_3.png
	dot -Tpng $(OUTDIR)/$(OUTPREF)_4.gv > $(OUTDIR)/$(OUTPREF)_4.png
	open $(OUTDIR)/$(OUTPREF)*.png
