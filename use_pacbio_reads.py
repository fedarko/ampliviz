import pysam


def read_pacbio_bam(pacbio_bam, bp_edges_list):
    """Returns a dict mapping each edge to its support value based on the long
       read data. We say a long read supports an edge if it maps to both the
       start and end position of the edge, regardless of chromosome (this
       method allows for long reads that hit multiple chromosomes).
    """

    pf = pysam.AlignmentFile(pacbio_bam, 'rb')
    # Determine sequence IDs (not just 'chr1', etc -- from what I've seen
    # these can include their positions as well, which is why we look through
    # the SQ header to get the IDs
    seq_ids = [seq_header['SN'] for seq_header in pf.header["SQ"]]
    short_id2seq_id = {}
    long_id2seq_bounds = {}
    for long_id in seq_ids:
        assert ':' in long_id
        short_id2seq_id[long_id.split(':')[0].lower()] = long_id
        # Assumes the reference chromosome ID is of the form
        # "chr1:abc-def", where abc is the observed start position of the
        # chromosome and def is the observed end position of the chromosome.
        bounds = [int(c) for c in long_id.split(':')[1].split('-')]
        long_id2seq_bounds[long_id] = bounds

    edge2support = {}
    for e in bp_edges_list:
        # Need to define a way to get the chromosome ID in the BAM file
        # from the chromosome ID in the AmpliconArchitect graph. This works
        # based on some of the data I've seen so far, but it may break and need
        # to be adjusted to work with other naming conventions.
        start_chr_format = e.start_pos.chromosome.lower().replace('_', '.')
        end_chr_format = e.end_pos.chromosome.lower().replace('_', '.')
        start_chromosome = short_id2seq_id[start_chr_format]
        end_chromosome = short_id2seq_id[end_chr_format]
        # This doesn't really guarantee support *between* the nodes, but it's
        # a reasonable approximation
        chr1bounds = long_id2seq_bounds[start_chromosome]
        chr2bounds = long_id2seq_bounds[end_chromosome]
        assert chr1bounds[1] - chr1bounds[0] > 1
        assert chr2bounds[1] - chr2bounds[0] > 1
        edge_start = e.start_pos.pos - chr1bounds[0]
        edge_end = e.end_pos.pos - chr2bounds[0]
        # print(e.start_pos.pos, chr1bounds, e.end_pos.pos, chr2bounds)
        # From testing, I don't *think* this will cause problems if edge_start
        # is located at the last possible position on a reference chromosome.
        f1 = pf.fetch(start_chromosome, edge_start, edge_start + 1)
        # ...However, if edge_end is 0, that will be a problem -- since -1 is
        # an invalid coordinate. We compensate for this by using (0, 1) as the
        # f2 region in that case (since we've asserted that these chromosomes
        # have at least two base pairs, we should be fine).
        if edge_end == 0:
            edge_end += 1
        f2 = pf.fetch(end_chromosome, edge_end - 1, edge_end)
        r1 = [r.query_name for r in f1]
        r2 = [r.query_name for r in f2]
        supporting_reads = set(r1) & set(r2)
        # Use of set() means that split reads (i.e. with same query name) will
        # only be counted as one, which is perfect.
        # Also, since we only look at intersection of r1 and r2, we don't care
        # about order -- so if theoretically this read hits the end pos then
        # start pos, this will catch that also. All we care about is presence
        # at the extreme positions.
        # print(e, len(supporting_reads))
        edge2support[e] = len(supporting_reads)
    pf.close()
    return edge2support
