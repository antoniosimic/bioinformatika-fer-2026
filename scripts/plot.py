# Plots benchmark results from benchmark_runner CSV output.
# Usage: python3 scripts/plot_results.py --input results/results.csv --output results/
# Author: Jakov

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── CLI ──────────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser(description="Plot Bamboo Filter benchmark results")
parser.add_argument("--input",  required=True, help="Path to benchmark CSV file")
parser.add_argument("--output", required=True, help="Directory for output PNG files")
args = parser.parse_args()

os.makedirs(args.output, exist_ok=True)

# ── Load data ────────────────────────────────────────────────────────────────

df = pd.read_csv(args.input)

# Separate synthetic and ecoli data
synthetic = df[df["data_source"] == "synthetic"].copy()
ecoli     = df[df["data_source"] == "ecoli"].copy()

# ── Shared style ─────────────────────────────────────────────────────────────

STYLE = {
    "figure.dpi":       150,
    "axes.spines.top":  False,
    "axes.spines.right":False,
    "axes.grid":        True,
    "grid.alpha":       0.3,
    "font.size":        11,
}
plt.rcParams.update(STYLE)
COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]

# ── Plot 1: Insert time vs sequence length (synthetic, one line per k) ───────

fig, ax = plt.subplots(figsize=(8, 5))
for i, k in enumerate(sorted(synthetic["k"].unique())):
    sub = synthetic[synthetic["k"] == k].sort_values("seq_length")
    ax.plot(sub["seq_length"], sub["insert_ms"],
            marker="o", label=f"k={k}", color=COLORS[i % len(COLORS)])

ax.set_xscale("log")
ax.set_xlabel("Sequence length (bp)")
ax.set_ylabel("Insert time (ms)")
ax.set_title("Bamboo Filter — Insert time vs sequence length (synthetic DNA)")
ax.legend(title="k-mer size")
fig.tight_layout()
fig.savefig(os.path.join(args.output, "plot_insert_time.png"))
plt.close(fig)
print("Saved: plot_insert_time.png")

# ── Plot 2: False positive rate vs k (synthetic + ecoli side by side) ────────

fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for ax, data, title in zip(
        axes,
        [synthetic, ecoli],
        ["Synthetic DNA", "E. coli K-12 genome"]):

    if data.empty:
        ax.set_title(f"{title}\n(no data)")
        continue

    # For synthetic: one line per seq_length; for ecoli: single bar chart
    if title == "Synthetic DNA":
        for i, length in enumerate(sorted(data["seq_length"].unique())):
            sub = data[data["seq_length"] == length].sort_values("k")
            ax.plot(sub["k"], sub["false_positive_rate"],
                    marker="o", label=f"len={length:,}",
                    color=COLORS[i % len(COLORS)])
        ax.legend(title="Sequence length", fontsize=9)
    else:
        sub = data.sort_values("k")
        ax.bar([str(k) for k in sub["k"]], sub["false_positive_rate"],
               color=COLORS[0], alpha=0.75)

    ax.set_xlabel("k-mer size")
    ax.set_ylabel("False positive rate")
    ax.set_title(f"False positive rate — {title}")
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=2))

fig.tight_layout()
fig.savefig(os.path.join(args.output, "plot_fpr.png"))
plt.close(fig)
print("Saved: plot_fpr.png")

# ── Plot 3: Memory (bits per item) vs k ──────────────────────────────────────

fig, ax = plt.subplots(figsize=(8, 5))

# Theoretical minimum for a 12-bit tag filter
ax.axhline(y=12.0, linestyle="--", color="gray", alpha=0.6,
           label="Theoretical min (12-bit tag)")

for i, source in enumerate(df["data_source"].unique()):
    sub = df[df["data_source"] == source].sort_values("k")
    # Average bits_per_item across all lengths for this source
    grouped = sub.groupby("k")["bits_per_item"].mean().reset_index()
    ax.plot(grouped["k"], grouped["bits_per_item"],
            marker="s", label=source, color=COLORS[i % len(COLORS)])

ax.set_xlabel("k-mer size")
ax.set_ylabel("Bits per item")
ax.set_title("Bamboo Filter — Memory efficiency (bits per item)")
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join(args.output, "plot_memory.png"))
plt.close(fig)
print("Saved: plot_memory.png")

# ── Plot 4: Lookup time — positive vs negative queries (ecoli) ───────────────

if not ecoli.empty:
    fig, ax = plt.subplots(figsize=(8, 5))
    x = range(len(ecoli))
    width = 0.35
    ks_label = [f"k={k}" for k in ecoli.sort_values("k")["k"]]
    ecoli_sorted = ecoli.sort_values("k")

    bars1 = ax.bar([i - width/2 for i in x],
                   ecoli_sorted["lookup_pos_ms"],
                   width, label="Positive (inserted k-mers)",
                   color=COLORS[0], alpha=0.8)
    bars2 = ax.bar([i + width/2 for i in x],
                   ecoli_sorted["lookup_neg_ms"],
                   width, label="Negative (random k-mers)",
                   color=COLORS[1], alpha=0.8)

    ax.set_xticks(list(x))
    ax.set_xticklabels(ks_label)
    ax.set_ylabel("Lookup time (ms)")
    ax.set_title("Bamboo Filter — Lookup time on E. coli genome")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(args.output, "plot_lookup_ecoli.png"))
    plt.close(fig)
    print("Saved: plot_lookup_ecoli.png")
else:
    print("Skipping plot 4 — no ecoli data in CSV (run with --fasta to add it)")

print("\nDone. All plots saved to:", args.output)