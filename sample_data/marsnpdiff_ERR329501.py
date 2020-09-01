import marsnpdiff

snps = marsnpdiff.find_snps(["ERR329501_chr4.bam"], ["GSM1553106_chr4.bam"], {"4": int(100e6)},)
snps.to_csv("marsnpdiff_ERR329501_chr4_vs_GSM1553106_chr4.tsv", sep="\t")
