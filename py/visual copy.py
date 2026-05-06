import pandas as pd

df = pd.read_csv(
    "/home/sanshou/project/toolkit/output/HYS_sample1_sorted.bed.gz",
    sep="\t",
    header=None,
    names=["chr", "start", "end", "name", "score"]
)
df["length"] = df["end"] - df["start"]

df = df[df["length"] > 0]
import matplotlib.pyplot as plt


def plot_fragment_distribution(df, out):
    plt.figure()
    plt.hist(df["length"], bins=100)
    plt.xlabel("Fragment size")
    plt.ylabel("Count")
    plt.title("Fragment Size Distribution")
    plt.savefig(out)
    plt.close()

def plot_chr_distribution(df, out):
    counts = df["chr"].value_counts()

    plt.figure()
    counts.plot(kind="bar")
    plt.xlabel("Chromosome")
    plt.ylabel("Reads")
    plt.title("Chromosome Distribution")
    plt.savefig(out)
    plt.close()

def plot_coverage(df, out):
    df["mid"] = (df["start"] + df["end"]) // 2

    plt.figure()
    plt.hist(df["mid"], bins=200)
    plt.xlabel("Genome position")
    plt.ylabel("Reads")
    plt.title("Coverage Distribution")
    plt.savefig(out)
    plt.close()

plot_fragment_distribution(df, "fragment_distribution.png")
plot_chr_distribution(df, "chr_distribution.png")
plot_coverage(df, "coverage_distribution.png")

