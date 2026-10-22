import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# plot of segment length
df = pd.read_csv("results.csv")
y_vars = ["NEG_MEAN_ERR", "NEG_MAX_ERR", "POS_MEAN_ERR"]
y_labels = ["Negative Control Mean Error", "Negative Control Maximum Error", "Positive Control Mean Error"]

# Create subplots: 1 row, 3 columns
fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True)

# Loop over each subplot and plot
for i, (ax, y_var, y_label) in enumerate(zip(axes, y_vars, y_labels)):
    sns.lineplot(data=df, x="LENGTH", y=y_var, hue="MARKERS", style="ERROR", dashes=False, markers=True, palette="viridis", ax=ax, legend=(i==0))
    ax.set_xlabel("Minimum Segment Length")
    ax.set_ylabel(y_label)

plt.tight_layout()
plt.savefig("tuning_plots/length.png", dpi=450, bbox_inches="tight", pad_inches=0)
plt.close()


# plot of number of markers
# Create subplots: 1 row, 3 columns
fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True)

# Loop over each subplot and plot
for i, (ax, y_var, y_label) in enumerate(zip(axes, y_vars, y_labels)):
    sns.lineplot(data=df, x="MARKERS", y=y_var, hue="LENGTH", style="ERROR", dashes=False, markers=True, palette="viridis", ax=ax, legend=(i==0))
    ax.set_xlabel("Minimum Number of Markers")
    ax.set_ylabel(y_label)

plt.tight_layout()
plt.savefig("tuning_plots/markers.png", dpi=450, bbox_inches="tight", pad_inches=0)
plt.close()


# plot of number of markers (for specific parameters)
# Create subplots: 1 row, 3 columns
fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True)

# Loop over each subplot and plot
subset = df[(df["LENGTH"]==0.5)]
for i, (ax, y_var, y_label) in enumerate(zip(axes, y_vars, y_labels)):
    sns.lineplot(data=subset, x="MARKERS", y=y_var, hue="LENGTH", style="ERROR", dashes=False, markers=True, palette="viridis", ax=ax, legend=(i==0))
    ax.set_xlabel("Minimum Number of Markers")
    ax.set_ylabel(y_label)

plt.tight_layout()
plt.savefig("tuning_plots/subset_markers.png", dpi=450, bbox_inches="tight", pad_inches=0)
plt.close()


# plot of error rate
# Create subplots: 1 row, 3 columns
fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True)

# Loop over each subplot and plot
for i, (ax, y_var, y_label) in enumerate(zip(axes, y_vars, y_labels)):
    sns.lineplot(data=df, x="ERROR", y=y_var, hue="LENGTH", style="MARKERS", dashes=False, markers=True, palette="viridis", ax=ax, legend=(i==0))
    ax.set_xlabel("Acceptable Error Rate")
    ax.set_ylabel(y_label)

plt.tight_layout()
plt.savefig("tuning_plots/error.png", dpi=450, bbox_inches="tight", pad_inches=0)
plt.close()


# plot of error rate (for specific parameters)
# Create subplots: 1 row, 3 columns
fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True)

# Loop over each subplot and plot
subset = df[(df["LENGTH"]==0.5) & ((df["MARKERS"]==100) | (df["MARKERS"]==200))]
for i, (ax, y_var, y_label) in enumerate(zip(axes, y_vars, y_labels)):
    sns.lineplot(data=subset, x="ERROR", y=y_var, hue="LENGTH", style="MARKERS", dashes=False, markers=True, palette="viridis", ax=ax, legend=(i==0))
    ax.set_xlabel("Acceptable Error Rate")
    ax.set_ylabel(y_label)

plt.tight_layout()
plt.savefig("tuning_plots/subset_error.png", dpi=450, bbox_inches="tight", pad_inches=0)
plt.close()