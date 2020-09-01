import marsnpdiff

snps = marsnpdiff.find_snps(["sample_a.bam"], ["sample_b.bam"], {"1": int(1e6)}, ll_threshold=0, chunk_size=100,)
snps.to_csv("marsnpdiff_sample_a_vs_sample_b.tsv", sep="\t")
