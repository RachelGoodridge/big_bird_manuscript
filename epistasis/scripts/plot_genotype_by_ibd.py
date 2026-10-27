from common import TRAITS, cell_summary, group_title, load_birds, plot_mean_heatmap, set_plot_theme, split_by_ibd

set_plot_theme()

# mean traits of each two-locus genotype combination, with genotypes in IBD with the founder at the locus shown separately
ANALYSES = [
    ("fortis", [("G29", "G27"), ("G03", "G27")], ["PC(Beak Size)", "Beak:Body Ratio"]),
    ("Big Birds", [("G29", "G27")], ["PC(Body Size)", "PC(Beak Size)", "Beak:Body Ratio"]),
]

birds = split_by_ibd(load_birds(), ["G03", "G27", "G29"])
for group, pairs, traits in ANALYSES:
    group_birds = birds[birds["Species"] == group]
    for locus1, locus2 in pairs:
        for trait in traits:
            means, _, ns = cell_summary(group_birds, locus1, locus2, trait, by_ibd=True)
            print(f"{group} {locus1} x {locus2} {trait}: {int(ns.sum().sum())} birds")
            plot_mean_heatmap(means, ns, locus1, locus2, trait, f"{group_title(group)}  |  {locus1} × {locus2}",
                              f"{group.replace(' ', '_')}_{locus1}_x_{locus2}_{TRAITS[trait]}_mean_by_ibd.png")
