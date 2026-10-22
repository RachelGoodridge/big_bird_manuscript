from common import count_snps, load_founder_ibis, load_master, negative_controls, positive_controls, read_bim

# load the IBIS segments shared with the founder
seg, coef = load_founder_ibis("tuning/finches")

# fraction of all SNPs that each bird shares IBD with the founder
bim = read_bim("tuning/dataset")
seg["SNP_count"] = count_snps(seg, bim)
fraction_snps_ibd = seg.groupby("ind")["SNP_count"].sum() / len(bim)

# get the birds used for tuning (birds with no segments share 0% IBD)
master = load_master()
master = master[master["meaningful.unique"].isin(coef["ind"])].copy()
master["fraction_snps_ibd"] = master["meaningful.unique"].map(fraction_snps_ibd).fillna(0)

# offspring should share IBD across their whole genome with the founder, so anything less is an error
pos_mean_error = (1 - positive_controls(master)["fraction_snps_ibd"]).mean()

# birds that predate the founder's descendants should share no IBD with him
neg_con = negative_controls(master)
neg_mean_error = neg_con["fraction_snps_ibd"].mean()
neg_max_error = neg_con["fraction_snps_ibd"].max()

print(f"{pos_mean_error},{neg_mean_error},{neg_max_error}")
