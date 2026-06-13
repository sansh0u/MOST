import snapatac2 as snap

adata = snap.pp.import_fragments(
    fragment_file= "/mnt/d/program/program/HYS_sample1_sorted.bed.gz",
    chrom_sizes=snap.genome.hg38,
    sorted_by_barcode=False,
)

peaks = snap.tl.macs3(
    adata,
    groupby="leiden"
)

snap.pp.make_peak_matrix(
    adata,
    use_rep=peaks
)


snap.metrics.frip(
    adata,
    peaks
)