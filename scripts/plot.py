# Plots benchmark results from two CSVs (own + reference Wang et al.)
# side by side. Computes a merged dataframe internally — no manual merge
# needed beforehand.
#
# Usage:
#   python3 scripts/plot.py \
#       --mine    results/results_mine.csv \
#       --ref     results/results_ref.csv \
#       --output  results/
#
# Author: Jakov

import argparse
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── CLI ──────────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser(
    description="Plot Bamboo Filter benchmark results, "
                "our implementation vs Wang et al. reference"
)
parser.add_argument("--mine", required=True,
                    help="CSV from ./benchmark_runner (our implementation)")
parser.add_argument("--ref", required=False, default=None,
                    help="CSV from ./reference_benchmark (Wang et al.). "
                         "Optional — if omitted, only our results are plotted.")
parser.add_argument("--output", required=True,
                    help="Directory for output PNG files")
args = parser.parse_args()

os.makedirs(args.output, exist_ok=True)

# ── Load and merge ───────────────────────────────────────────────────────────

df_mine = pd.read_csv(args.mine)
if args.ref:
    df_ref = pd.read_csv(args.ref)
    df = pd.concat([df_mine, df_ref], ignore_index=True)
else:
    df = df_mine.copy()

filters = sorted(df["filter"].unique())
print(f"Loaded {len(df)} rows. Filters: {filters}")

synthetic = df[df["data_source"] == "synthetic"].copy()
ecoli     = df[df["data_source"] == "ecoli"].copy()

# ── Shared style ─────────────────────────────────────────────────────────────

STYLE = {
    "figure.dpi":         150,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.grid":          True,
    "grid.alpha":         0.3,
    "font.size":          11,
    "lines.linewidth":    2,
    "lines.markersize":   7,
}
plt.rcParams.update(STYLE)

FILTER_COLORS = {
    "Bamboo":          "#1f77b4",
    "BambooReference": "#d62728",
}
FILTER_LABELS = {
    "Bamboo":          "Our implementation",
    "BambooReference": "Wang et al. reference",
}
FILTER_MARKERS = {
    "Bamboo":          "o",
    "BambooReference": "s",
}

def color_for(f):  return FILTER_COLORS.get(f, "#999999")
def label_for(f):  return FILTER_LABELS.get(f, f)
def marker_for(f): return FILTER_MARKERS.get(f, "^")

# ── Plot 1: Insert time vs sequence length (one panel per k) ─────────────────

ks_synth = sorted(synthetic["k"].unique())
fig, axes = plt.subplots(1, len(ks_synth), figsize=(4.2 * len(ks_synth), 4.5),
                         sharey=True)
if len(ks_synth) == 1:
    axes = [axes]

for ax, k in zip(axes, ks_synth):
    for f in filters:
        sub = synthetic[(synthetic["k"] == k) & (synthetic["filter"] == f)] \
            .sort_values("seq_length")
        if sub.empty:
            continue
        ax.plot(sub["seq_length"], sub["insert_ms"],
                marker=marker_for(f), color=color_for(f),
                label=label_for(f))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Sequence length [bp]")
    ax.set_title(f"k = {k}")
    ax.legend(fontsize=8)

axes[0].set_ylabel("Insert time [ms]")
fig.suptitle("Insert time vs sequence length — synthetic DNA", fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(args.output, "plot_insert_time.png"))
plt.close(fig)
print("Saved: plot_insert_time.png")

# ── Plot 2: FPR vs k ─────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
for f in filters:
    sub = synthetic[synthetic["filter"] == f].sort_values("k")
    grouped = sub.groupby("k")["false_positive_rate"].mean().reset_index()
    if grouped.empty:
        continue
    ax.plot(grouped["k"], grouped["false_positive_rate"],
            marker=marker_for(f), color=color_for(f), label=label_for(f))
ax.set_xlabel("k-mer size")
ax.set_ylabel("False positive rate")
ax.set_title("Synthetic DNA (averaged over sequence lengths)")
ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=2))
ax.legend(fontsize=9)

ax = axes[1]
ks_ecoli = sorted(ecoli["k"].unique())
if ks_ecoli:
    x = range(len(ks_ecoli))
    width = 0.4 if len(filters) > 1 else 0.6
    for i, f in enumerate(filters):
        sub = ecoli[ecoli["filter"] == f].set_index("k").reindex(ks_ecoli)
        offset = (i - (len(filters) - 1) / 2) * width
        ax.bar([xi + offset for xi in x], sub["false_positive_rate"],
               width, label=label_for(f), color=color_for(f), alpha=0.85)
    ax.set_xticks(list(x))
    ax.set_xticklabels([str(k) for k in ks_ecoli])
    ax.set_xlabel("k-mer size")
    ax.set_ylabel("False positive rate")
    ax.set_title("E. coli K-12 genome")
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=2))
    ax.legend(fontsize=9)
else:
    ax.set_title("E. coli K-12 genome (no data)")

fig.suptitle("False positive rate vs k-mer size", fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(args.output, "plot_fpr.png"))
plt.close(fig)
print("Saved: plot_fpr.png")

# ── Plot 3: Memory ───────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, src, title in zip(axes,
                          ["synthetic", "ecoli"],
                          ["Synthetic DNA", "E. coli K-12 genome"]):
    data = df[df["data_source"] == src]
    if data.empty:
        ax.set_title(f"{title} (no data)")
        continue

    for f in filters:
        sub = data[data["filter"] == f]
        grouped = sub.groupby("k")["bits_per_item"].mean().reset_index() \
                                                   .sort_values("k")
        if grouped.empty:
            continue
        ax.plot(grouped["k"], grouped["bits_per_item"],
                marker=marker_for(f), color=color_for(f), label=label_for(f))

    ax.axhline(y=12.0, linestyle="--", color="gray", alpha=0.6,
               label="Theoretical min (12-bit tag)")
    ax.set_xlabel("k-mer size")
    ax.set_ylabel("Bits per item")
    ax.set_title(title)
    if src == "ecoli":
        ax.set_yscale("log")
    ax.legend(fontsize=9)

fig.suptitle("Memory efficiency (bits per item)", fontsize=13)
fig.tight_layout()
fig.savefig(os.path.join(args.output, "plot_memory.png"))
plt.close(fig)
print("Saved: plot_memory.png")

# ── Plot 4: Lookup time on E. coli ───────────────────────────────────────────

if not ecoli.empty:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

    for ax, col, title in zip(
            axes,
            ["lookup_pos_ms", "lookup_neg_ms"],
            ["Positive lookups (inserted k-mers)",
             "Negative lookups (random k-mers)"]):

        ks_ecoli = sorted(ecoli["k"].unique())
        x = range(len(ks_ecoli))
        width = 0.4 if len(filters) > 1 else 0.6
        for i, f in enumerate(filters):
            sub = ecoli[ecoli["filter"] == f].set_index("k").reindex(ks_ecoli)
            offset = (i - (len(filters) - 1) / 2) * width
            ax.bar([xi + offset for xi in x], sub[col],
                   width, label=label_for(f), color=color_for(f), alpha=0.85)
        ax.set_xticks(list(x))
        ax.set_xticklabels([f"k={k}" for k in ks_ecoli])
        ax.set_ylabel("Lookup time [ms]")
        ax.set_title(title)
        ax.set_yscale("log")
        ax.legend(fontsize=9)

    fig.suptitle("Lookup time on E. coli K-12 genome", fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(args.output, "plot_lookup_ecoli.png"))
    plt.close(fig)
    print("Saved: plot_lookup_ecoli.png")
else:
    print("Skipping plot 4 — no ecoli data")

# ── Plot 5: Ratio plot — ours / reference ────────────────────────────────────

if "Bamboo" in filters and "BambooReference" in filters:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, col, title, ylabel in zip(
            axes,
            ["insert_ms", "memory_bytes"],
            ["Insert time ratio (ours / reference)",
             "Memory ratio (ours / reference)"],
            ["Time ratio", "Memory ratio"]):

        for src, marker in [("synthetic", "o"), ("ecoli", "s")]:
            data = df[df["data_source"] == src]
            mine = data[data["filter"] == "Bamboo"][["k", "seq_length", col]]
            ref  = data[data["filter"] == "BambooReference"][["k", "seq_length", col]]
            merged = mine.merge(ref, on=["k", "seq_length"],
                                suffixes=("_mine", "_ref"))
            if merged.empty:
                continue
            merged["ratio"] = merged[f"{col}_mine"] / merged[f"{col}_ref"]
            grouped = merged.groupby("k")["ratio"].mean().reset_index() \
                                                  .sort_values("k")
            ax.plot(grouped["k"], grouped["ratio"],
                    marker=marker, label=src)

        ax.axhline(y=1.0, linestyle=":",  color="green",  alpha=0.7,
                   label="parity (1.0)")
        ax.axhline(y=2.0, linestyle="--", color="orange", alpha=0.7,
                   label="2x (−10 bod threshold)")
        ax.axhline(y=3.0, linestyle="--", color="red",    alpha=0.7,
                   label="3x (−15 bod threshold)")

        ax.set_xlabel("k-mer size")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_yscale("log")
        ax.legend(fontsize=8)

    fig.suptitle("Performance ratio: our implementation vs Wang et al. reference",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(args.output, "plot_ratio.png"))
    plt.close(fig)
    print("Saved: plot_ratio.png")

print(f"\nAll plots saved to: {args.output}")