# Plots benchmark results from one or two CSVs (our implementation and,
# optionally, the Wang et al. reference). Produces four PNG figures:
#
#   plot_time.png    insert + lookup time, synthetic and E. coli
#   plot_fpr.png     false positive rate vs k
#   plot_memory.png  bits per item vs k
#   plot_ratio.png   ours / reference summary (only when --ref is given)
#
# Usage:
#   python3 scripts/plot.py \
#       --mine    results/results_mine.csv \
#       --ref     results/results_ref.csv \
#       --output  results/
#
# Author: Jakov Malić

import argparse
import os
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

# ── Load ─────────────────────────────────────────────────────────────────────

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

plt.rcParams.update({
    "figure.dpi":         150,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.grid":          True,
    "grid.alpha":         0.3,
    "font.size":          11,
    "lines.linewidth":    2,
    "lines.markersize":   7,
})

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

def style(f):
    """Return (color, label, marker) for a filter name."""
    return (FILTER_COLORS.get(f, "#999999"),
            FILTER_LABELS.get(f, f),
            FILTER_MARKERS.get(f, "^"))


def plot_lines_per_k(ax, data, ycol, ylabel, title, logy=True):
    """Line plot of `ycol` vs k for each filter, mean over seq_length."""
    for f in filters:
        sub = data[data["filter"] == f]
        grouped = sub.groupby("k")[ycol].mean().reset_index().sort_values("k")
        if grouped.empty:
            continue
        color, label, marker = style(f)
        ax.plot(grouped["k"], grouped[ycol],
                marker=marker, color=color, label=label)
    ax.set_xlabel("k-mer size")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if logy:
        ax.set_yscale("log")
    ax.legend(fontsize=9)


# ── Plot 1: Time (insert + lookup, synthetic + ecoli) ────────────────────────
#
# 2x2 grid: rows = phase (insert, lookup_pos), cols = data source.
# Y-axis time in ms (log scale). X-axis k-mer size.

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
for row, (metric, label) in enumerate([
        ("insert_ms",     "Insert time [ms]"),
        ("lookup_pos_ms", "Lookup time [ms]")]):
    for col, (src, title_src, data) in enumerate([
            ("synthetic", "Synthetic DNA", synthetic),
            ("ecoli",     "E. coli K-12 genome", ecoli)]):
        ax = axes[row][col]
        if data.empty:
            ax.set_title(f"{title_src} (no data)")
            continue
        plot_lines_per_k(ax, data, metric, label,
                         f"{title_src} — {label.split(' [')[0].lower()}")

fig.suptitle("Time performance — ours vs reference", fontsize=14)
fig.tight_layout()
fig.savefig(os.path.join(args.output, "plot_time.png"))
plt.close(fig)
print("Saved: plot_time.png")

# ── Plot 2: False positive rate ──────────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, (src, title, data) in zip(axes, [
        ("synthetic", "Synthetic DNA", synthetic),
        ("ecoli",     "E. coli K-12 genome", ecoli)]):
    if data.empty:
        ax.set_title(f"{title} (no data)")
        continue
    plot_lines_per_k(ax, data, "false_positive_rate",
                     "False positive rate", title, logy=False)
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(xmax=1, decimals=2))

fig.suptitle("False positive rate vs k-mer size", fontsize=14)
fig.tight_layout()
fig.savefig(os.path.join(args.output, "plot_fpr.png"))
plt.close(fig)
print("Saved: plot_fpr.png")

# ── Plot 3: Memory ───────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, (src, title, data) in zip(axes, [
        ("synthetic", "Synthetic DNA", synthetic),
        ("ecoli",     "E. coli K-12 genome", ecoli)]):
    if data.empty:
        ax.set_title(f"{title} (no data)")
        continue
    plot_lines_per_k(ax, data, "bits_per_item",
                     "Bits per item", title, logy=(src == "ecoli"))
    # Theoretical floor: 12 bits per stored tag. Real cost is higher
    # because of the load factor (slots reserved but unused).
    ax.axhline(y=12.0, linestyle="--", color="gray", alpha=0.5,
               label="12-bit tag floor")
    ax.legend(fontsize=9)

fig.suptitle("Memory efficiency (bits per item)", fontsize=14)
fig.tight_layout()
fig.savefig(os.path.join(args.output, "plot_memory.png"))
plt.close(fig)
print("Saved: plot_memory.png")

# ── Plot 4: Total memory in bytes ────────────────────────────────────────────
#
# Same panels as plot_memory.png but in absolute KB / MB instead of the
# normalised bits-per-item — useful to see the raw footprint of the filter
# at each k, not just how efficient it is per stored item.

def fmt_kb(x, _pos):
    """Tick label: render a KB value as '<n> KB' or '<n> MB' (no 10^x)."""
    if x <= 0:
        return ""
    if x >= 1024:
        return f"{x / 1024:g} MB"
    return f"{x:g} KB"

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax, (src, title, data) in zip(axes, [
        ("synthetic", "Synthetic DNA", synthetic),
        ("ecoli",     "E. coli K-12 genome", ecoli)]):
    if data.empty:
        ax.set_title(f"{title} (no data)")
        continue
    for f in filters:
        sub = data[data["filter"] == f]
        grouped = sub.groupby("k")["memory_bytes"].mean().reset_index() \
                                                  .sort_values("k")
        if grouped.empty:
            continue
        color, label, marker = style(f)
        # Convert bytes → KB for readability; tick formatter promotes to MB.
        ax.plot(grouped["k"], grouped["memory_bytes"] / 1024.0,
                marker=marker, color=color, label=label)
    ax.set_xlabel("k-mer size")
    ax.set_ylabel("Total memory")
    ax.set_title(title)
    ax.set_yscale("log")
    # Default log axis only labels powers of 10 (sparse — for our range
    # we'd see just "1 MB" with everything else clipped). Promote 2- and
    # 5-subticks to major so we get readable 1, 2, 5, 10, 20, 50, …
    # spacing across the whole range.
    ax.yaxis.set_major_locator(
        ticker.LogLocator(base=10.0, subs=(1.0, 2.0, 5.0), numticks=20))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(fmt_kb))
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.legend(fontsize=9)

fig.suptitle("Total memory used by the filter", fontsize=14)
fig.tight_layout()
fig.savefig(os.path.join(args.output, "plot_memory_total.png"))
plt.close(fig)
print("Saved: plot_memory_total.png")

# ── Plot 5: Ratio summary (only if reference data is present) ────────────────
#
# Shows ours / reference for both insert time and memory, per k, with one
# line per data source. A horizontal reference at 1.0 marks parity — below
# means we're better (faster / smaller), above means worse.

if "Bamboo" in filters and "BambooReference" in filters:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax, (col, title, ylabel) in zip(axes, [
            ("insert_ms",    "Insert time ratio",  "Time ratio (ours / reference)"),
            ("memory_bytes", "Memory ratio",       "Memory ratio (ours / reference)")]):
        for src, marker in [("synthetic", "o"), ("ecoli", "s")]:
            data = df[df["data_source"] == src]
            mine = data[data["filter"] == "Bamboo"][["k", "seq_length", col]]
            ref  = data[data["filter"] == "BambooReference"][["k", "seq_length", col]]
            merged = mine.merge(ref, on=["k", "seq_length"],
                                suffixes=("_mine", "_ref"))
            if merged.empty:
                continue
            merged = merged[merged[f"{col}_ref"] > 0]      # skip bad rows
            merged["ratio"] = merged[f"{col}_mine"] / merged[f"{col}_ref"]
            grouped = merged.groupby("k")["ratio"].mean().reset_index() \
                                                  .sort_values("k")
            ax.plot(grouped["k"], grouped["ratio"],
                    marker=marker, label=src)

        ax.axhline(y=1.0, linestyle=":", color="gray", alpha=0.7,
                   label="parity")
        ax.set_xlabel("k-mer size")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_yscale("log")
        ax.legend(fontsize=9)

    fig.suptitle("Performance ratio: our implementation vs Wang et al. reference",
                 fontsize=14)
    fig.tight_layout()
    fig.savefig(os.path.join(args.output, "plot_ratio.png"))
    plt.close(fig)
    print("Saved: plot_ratio.png")

print(f"\nAll plots saved to: {args.output}")
