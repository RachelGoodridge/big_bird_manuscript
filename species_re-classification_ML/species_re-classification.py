import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
import random
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.mixture import GaussianMixture
from sklearn.ensemble import RandomForestClassifier
from scipy.stats import norm
from scipy.optimize import fsolve
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score

# make plot theme more like R
custom_params = {"font.family": "DejaVu Sans", "axes.facecolor": "white", "axes.edgecolor": "#DCDCDC", "axes.linewidth": 1.2, "axes.grid": True,
    "grid.color": "#EAEAEA", "grid.linestyle": "-", "grid.linewidth": 0.5, "xtick.color": "black", "ytick.color": "black", "xtick.major.size": 0,
    "ytick.major.size": 0}
sns.set_style("whitegrid", rc=custom_params)
plt.rcParams.update(custom_params)
color_map_hex = {"fortis": "#1C86EE", "scandens": "#43CD80", "fuliginosa": "#C6E2FF", "hybrid": "#AB82FF"}
# species epithets get italicized in figure text/legends
SPECIES_NAMES = {"fortis", "scandens", "fuliginosa"}

# set random seed
seed = 1
random.seed(seed)
np.random.seed(seed)

####################################################################################
# PCA

# first get the master database and filter
master = pd.read_csv("/local/storage/all_birds_highcov_lowcov/Finch_Database_V2.csv")
master.drop_duplicates(subset=["meaningful.unique"], inplace=True)
master = master[master["Island"] == "Daphne"].copy()

# fill in missing years with the minimum year, year born, or last year
master["First.year"] = master["First.year"].fillna(master["First.year.min"])
master["First.year"] = master["First.year"].fillna(master["Year.born"])
master["First.year"] = master["First.year"].fillna(master["Last.year"])
master["Last.year"] = master["Last.year"].fillna(master["First.year"])

# update the newest lineage classifications
bbs_lines = pd.read_csv("/local/storage/BigBirds/BigBirds.csv")
bbs_lines.drop_duplicates(subset=["meaningful.unique"], inplace=True)
master["lineage"] = master["meaningful.unique"].map(bbs_lines.set_index("meaningful.unique")["lineage"])
master["lineage"] = master["lineage"].fillna(master["Species"])
master = master.drop(columns=["Species"]).rename(columns={"lineage": "Species"})

# filter to only the species of interest
master = master[master["Species"].isin(["scandens", "fortis", "hybrid", "fuliginosa"])]
master.set_index("meaningful.unique", inplace=True)
master = master[["Species", "First.year", "bill.length", "bill.depth", "bill.width", "weight"]].dropna().copy()

# PC1/PC2 of bill dimensions -> beak size / beak pointedness
bill_scaled = StandardScaler().fit_transform(master[["bill.length", "bill.depth", "bill.width"]])
pca = PCA(n_components=2)
bill_pcs = pca.fit_transform(bill_scaled)
exp_var1, exp_var2 = pca.fit(bill_scaled).explained_variance_ratio_
master["PC(Beak Size)"] = bill_pcs[:, 0]
master["PC(Beak Pointedness)"] = bill_pcs[:, 1]

# compute beak dimensions
master["pointy"] = master["bill.length"] / master["bill.depth"]
master["sum"] = master["bill.length"] + master["bill.depth"] + master["bill.width"]

# plot the PCA results
plt.figure(figsize=(6.4, 4.8))
sns.scatterplot(data=master, x="PC(Beak Size)", y="PC(Beak Pointedness)", hue="Species", palette=color_map_hex, s=40, alpha=0.7)
plt.xlabel(f"PC1, Beak Size ({round(exp_var1*100, 2)}%)")
plt.ylabel(f"PC2, Beak Pointedness ({round(exp_var2*100, 2)}%)")
leg = plt.legend(title=None)
for text in leg.get_texts():
    if text.get_text() in SPECIES_NAMES:
        text.set_fontstyle("italic")
sns.despine(top=True, right=True, left=True, bottom=True)
plt.tight_layout()
plt.savefig("beak_traits_pca.png", dpi=450, bbox_inches="tight", pad_inches=0.05)
plt.close()

# plot PC1 vs. sum of beak traits
plt.figure(figsize=(6.4, 4.8))
sns.scatterplot(data=master, x="PC(Beak Size)", y="sum", hue="Species", palette=color_map_hex, s=40, alpha=0.7)
plt.xlabel("PC1, Beak Size")
plt.ylabel("Sum of Beak Traits")
leg = plt.legend(title=None)
for text in leg.get_texts():
    if text.get_text() in SPECIES_NAMES:
        text.set_fontstyle("italic")
sns.despine(top=True, right=True, left=True, bottom=True)
plt.tight_layout()
plt.savefig("pc1_vs_sum.png", dpi=450, bbox_inches="tight", pad_inches=0.05)
plt.close()

# plot PC2 vs. pointedness
plt.figure(figsize=(6.4, 4.8))
sns.scatterplot(data=master, x="PC(Beak Pointedness)", y="pointy", hue="Species", palette=color_map_hex, s=40, alpha=0.7)
plt.xlabel("PC2, Beak Pointedness")
plt.ylabel("Beak Pointedness")
leg = plt.legend(title=None)
for text in leg.get_texts():
    if text.get_text() in SPECIES_NAMES:
        text.set_fontstyle("italic")
sns.despine(top=True, right=True, left=True, bottom=True)
plt.tight_layout()
plt.savefig("pc2_vs_pointy.png", dpi=450, bbox_inches="tight", pad_inches=0.05)
plt.close()

####################################################################################
# Random Forest

# remove hybrids
master = master[master["Species"] != "hybrid"].copy()

# get the admixture info
data = pd.read_csv("/SSD/reg259/big_bird_manuscript/ancestry_estimates/whole_data_new/admixture/finches.9.Q", sep="\\s+", header=None)
fam = pd.read_csv("/SSD/reg259/big_bird_manuscript/ancestry_estimates/whole_data_new/admixture/finches.fam", sep="\\s+", header=None)
data.index = fam[1]

# add the admixture info to the master dataframe
data.rename(columns={0: "anc_fort", 5: "anc_scan", 4: "anc_fuli"}, inplace=True)
for anc in ["anc_fort", "anc_scan", "anc_fuli"]:
    master[anc] = master.index.map(data[anc])

# remove any that were not sequenced
master.dropna(inplace=True)

# select training and test data
training = master[master["First.year"] != 2022].copy()
test = master[master["First.year"] == 2022].copy()
features = ["PC(Beak Size)", "PC(Beak Pointedness)", "weight", "anc_fort", "anc_scan", "anc_fuli"]

# function to find the outliers
def find_outliers(data):
    Q1 = data.quantile(0.25)
    Q3 = data.quantile(0.75)
    IQR = Q3 - Q1
    lower = Q1 - 1.5 * IQR
    upper = Q3 + 1.5 * IQR
    return(lower, upper)

# remove outliers of each species (except hybrids) and for each feature
inds_to_drop = []
for spp in ["scandens", "fortis", "fuliginosa"]:
    subset = training[(training["Species"] == spp)]
    for feature in features:
        lower, upper = find_outliers(subset[feature])
        outliers = subset[(subset[feature] < lower) | (subset[feature] > upper)].index.to_list()
        for outlier in outliers:
            inds_to_drop.append(outlier)
training.drop(inds_to_drop, inplace=True)

# split into features and labels
X_all = training[features].values
Y_all = training["Species"].values
X_test = test[features].values
Y_test = test["Species"].values

# evaluate the performace of the model
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
accuracies = []
for fold, (train_idx, val_idx) in enumerate(skf.split(X_all, Y_all), 1):
    # split into training and validation sets
    X_train, X_val = X_all[train_idx], X_all[val_idx]
    Y_train, Y_val = Y_all[train_idx], Y_all[val_idx]

    # standardize features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_val_scaled = scaler.transform(X_val)

    # fit Random Forest and predict on validation data
    clf = RandomForestClassifier(class_weight="balanced", random_state=seed)
    clf.fit(X_train_scaled, Y_train)
    Y_pred = clf.predict(X_val_scaled)

    # report the accuracy
    acc = accuracy_score(Y_val, Y_pred)
    accuracies.append(acc)
    print(f"Fold {fold} Accuracy: {round(acc * 100, 2)}%")

print(f"\nAverage Accuracy Across 5 Folds: {round(np.mean(accuracies) * 100, 2)}%")

# standardize features
scaler = StandardScaler()
X_all_scaled = scaler.fit_transform(X_all)
X_test_scaled = scaler.transform(X_test)

# fit Random Forest and predict on the test data
clf = RandomForestClassifier(class_weight="balanced", random_state=seed)
clf.fit(X_all_scaled, Y_all)
test["prediction"] = clf.predict(X_test_scaled)
test["same"] = (test["Species"] == test["prediction"])

####################################################################################
# Post-ML Plots

# plot the pca data again but with predictions for 2022
different = test[~test["same"]]
plt.figure(figsize=(6.4, 4.8))
sns.scatterplot(data=test[test["same"]], x="PC(Beak Size)", y="PC(Beak Pointedness)", hue="Species", palette=color_map_hex, s=30, alpha=0.6, zorder=2)
for spp, color in color_map_hex.items():
    train_sub = training[training["Species"] == spp]
    if not train_sub.empty:
        plt.scatter(train_sub["PC(Beak Size)"], train_sub["PC(Beak Pointedness)"], facecolors="none", edgecolors=color, marker="o", linewidths=0.8, s=20, alpha=0.5, zorder=1)
    test_sub = different[different["prediction"] == spp]
    if not test_sub.empty:
        plt.scatter(test_sub["PC(Beak Size)"], test_sub["PC(Beak Pointedness)"], c=test_sub["Species"].map(color_map_hex), edgecolors=color, marker="o", linewidths=1.5, s=50, alpha=0.75, zorder=3)
train_handle = mlines.Line2D([], [], color="gray", marker="o", linestyle="None", markersize=6, markerfacecolor="None", markeredgewidth=1, label="training")
test_handle = mlines.Line2D([], [], color="gray", marker="o", linestyle="None", markersize=6, markerfacecolor="gray", label="2022")
plt.xlabel("PC1, Beak Size")
plt.ylabel("PC2, Beak Pointedness")
leg = plt.legend(title=None, handles=[train_handle, test_handle], loc="upper left")
sns.despine(top=True, right=True, left=True, bottom=True)
plt.tight_layout()
plt.savefig("pca_rfc_ml.png", dpi=450, bbox_inches="tight", pad_inches=0.05)
plt.close()

# now make a series of ancestry vs. PC1 plots
labels = [r"$\mathit{fortis}$ ancestry", r"$\mathit{scandens}$ ancestry", r"$\mathit{fuliginosa}$ ancestry"]
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(23, 6))
axs = axs.flatten()
for i, anc, label in zip([0, 1, 2], ["anc_fort", "anc_scan", "anc_fuli"], labels):
    sns.scatterplot(data=test[test["same"]], x=anc, y="PC(Beak Size)", hue="Species", palette=color_map_hex, s=30, alpha=0.6, zorder=2, ax=axs[i])
    for spp, color in color_map_hex.items():
        train_sub = training[training["Species"] == spp]
        if not train_sub.empty:
            axs[i].scatter(train_sub[anc], train_sub["PC(Beak Size)"], facecolors="none", edgecolors=color, marker="o",
                           linewidths=0.8, s=20, alpha=0.5, zorder=1)
        test_sub = different[different["prediction"] == spp]
        if not test_sub.empty:
            axs[i].scatter(test_sub[anc], test_sub["PC(Beak Size)"], c=test_sub["Species"].map(color_map_hex), edgecolors=color, marker="o",
                           linewidths=1.5, s=50, alpha=0.75, zorder=3)
    axs[i].set_xlabel(label)
    axs[i].set_ylabel("PC1, Beak Size")
    leg = axs[i].legend(title=None)
    for text in leg.get_texts():
        if text.get_text() in SPECIES_NAMES:
            text.set_fontstyle("italic")
sns.despine(top=True, right=True, left=True, bottom=True)
plt.tight_layout()
plt.savefig("anc_vs_PC1.png", dpi=450, bbox_inches="tight", pad_inches=0.05)
plt.close()

# now make a series of ancestry vs. PC2 plots
fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(23, 6))
axs = axs.flatten()
for i, anc, label in zip([0, 1, 2], ["anc_fort", "anc_scan", "anc_fuli"], labels):
    sns.scatterplot(data=test[test["same"]], x=anc, y="PC(Beak Pointedness)", hue="Species", palette=color_map_hex, s=30, alpha=0.6, zorder=2, ax=axs[i])
    for spp, color in color_map_hex.items():
        train_sub = training[training["Species"] == spp]
        if not train_sub.empty:
            axs[i].scatter(train_sub[anc], train_sub["PC(Beak Pointedness)"], facecolors="none", edgecolors=color, marker="o",
                           linewidths=0.8, s=20, alpha=0.5, zorder=1)
        test_sub = different[different["prediction"] == spp]
        if not test_sub.empty:
            axs[i].scatter(test_sub[anc], test_sub["PC(Beak Pointedness)"], c=test_sub["Species"].map(color_map_hex), edgecolors=color, marker="o",
                           linewidths=1.5, s=50, alpha=0.75, zorder=3)
    axs[i].set_xlabel(label)
    axs[i].set_ylabel("PC2, Beak Pointedness")
    leg = axs[i].legend(title=None)
    for text in leg.get_texts():
        if text.get_text() in SPECIES_NAMES:
            text.set_fontstyle("italic")
sns.despine(top=True, right=True, left=True, bottom=True)
plt.tight_layout()
plt.savefig("anc_vs_PC2.png", dpi=450, bbox_inches="tight", pad_inches=0.05)
plt.close()

# print out the ones that changed
print(different)