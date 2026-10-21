# identify the false positive segments and create bed files
python plot_false_pos_loc.py

# trim the original segments to exclude the false positive regions
bedtools subtract -a trimming/features.bed -b trimming/exclude.bed -sorted > trimming/trimmed.bed

# plot the locations of the trimmed segments
python plot_IBD_loc.py