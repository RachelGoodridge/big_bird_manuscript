#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p trimming trimming_plots output

# create bed files of the negative control segments, the segments to trim, and the chromosome lengths
python scripts/make_trim_beds.py

# count how many negative controls share IBD with the founder at each position
bedtools genomecov -bg -i trimming/neg_control.bed -g trimming/genome.txt > trimming/neg_control_depth.bedgraph

# exclude regions where more than 5% of negative controls share IBD with the founder (likely false positives)
NUM_NEG=$(wc -l < trimming/neg_control_ids.txt)
echo "excluding regions shared with more than $(awk -v n=$NUM_NEG 'BEGIN{print 0.05 * n}') of $NUM_NEG negative controls"
awk -v n=$NUM_NEG '$4 > 0.05 * n' trimming/neg_control_depth.bedgraph | bedtools merge -i - > trimming/exclude.bed

# trim the original segments to exclude the false positive regions
bedtools subtract -a trimming/features.bed -b trimming/exclude.bed > trimming/trimmed.bed

# plot the false positive regions and the locations of the trimmed segments
python scripts/plot_false_pos_loc.py
python scripts/plot_IBD_loc.py
