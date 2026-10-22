import os
import sys
from common import FOUNDER, VCF_SAMPLES, load_master, negative_controls, positive_controls

# every Daphne bird in the master database that is also in the full vcf
master = load_master()
with open(VCF_SAMPLES) as f:
    vcf_birds = {line.strip() for line in f}
master = master[master["meaningful.unique"].isin(vcf_birds)]
all_birds = set(master["meaningful.unique"])

# birds for tuning: the founder, his offspring (positive controls), and birds that predate his descendants (negative controls)
tuning_birds = {FOUNDER} | set(positive_controls(master)["meaningful.unique"]) | set(negative_controls(master)["meaningful.unique"])

# write each list, but never replace an existing file with a different set of birds
def write_birds(path, birds):
    if os.path.exists(path):
        with open(path) as f:
            existing = {line.strip() for line in f if line.strip()}
        if existing == birds:
            print(f"{path} already lists the same {len(birds)} birds, keeping it unchanged")
            return
        sys.exit(f"ERROR: {path} differs from the new selection ({len(existing - birds)} only in file, {len(birds - existing)} only in new); not overwriting")
    with open(path, "w") as f:
        f.writelines(f"{bird}\n" for bird in sorted(birds))
    print(f"wrote {len(birds)} birds to {path}")

write_birds("dataset/Daphne_birds_for_ibis.txt", all_birds)
write_birds("dataset/birds_for_tuning.txt", tuning_birds)
