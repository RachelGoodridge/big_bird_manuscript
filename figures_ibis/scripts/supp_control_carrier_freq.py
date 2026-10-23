from common import control_yearly_summary, load_master, plot_controls, set_plot_theme

set_plot_theme()

# IBD carrier frequency among fortis / scandens over time for the founder vs. each negative control treated as a founder
summary = control_yearly_summary(load_master())
plot_controls(summary, "carrier_freq", "IBD Segment Carrier Frequency", "supplement/year_vs_carrier_control.png")
