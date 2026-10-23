from common import IBD_LABEL, control_yearly_summary, load_master, plot_controls, set_plot_theme

set_plot_theme()

# mean fraction of SNPs in IBD among fortis / scandens over time for the founder vs. each negative control treated as a founder
summary = control_yearly_summary(load_master())
plot_controls(summary, "mean_ibd", f"Mean {IBD_LABEL}", "supplement/year_vs_mean_IBD_control.png")
