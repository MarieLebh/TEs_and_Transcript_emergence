from collections import defaultdict
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

def get_genome_len(fasta):
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    genome_len = 0
    for value in seq_dict.values():
        genome_len += len(value.seq)
    return genome_len

def get_category(attr):

    m = re.search(r"\((.*?)\)", attr)

    if not m:
        return "Unknown"

    parts = [x.strip() for x in m.group(1).split(",")]

    if len(parts) == 3:
        family, order, superclass = parts
        if order in ["TIR", "MITE"]:
            order = "DNA"
        return order

    elif len(parts) == 2:
        family, superclass = parts
        if family in ["TIR", "MITE"]:
            family = "DNA"
        return family

    if parts[0] in ["TIR", "MITE"]:
            te = "DNA"
    return te

def main(BED, GENOME_SIZE, pop):

    records = []

    with open(BED) as fh:
        for line in fh:
            if not line.strip():
                continue

            chrom, start, end, attr, score, strand = line.rstrip().split("\t")

            records.append(
                (
                    chrom,
                    int(start),
                    int(end),
                    get_category(attr)
                )
            )

    records.sort(key=lambda x: (x[0], x[1], x[2]))

    claimed = defaultdict(list)

    class_bp = defaultdict(int)
    overlap_bp = 0

    for chrom, start, end, category in records:

        length = end - start  # BED is half-open

        overlap = 0

        for s, e in claimed[chrom]:

            if e <= start:
                continue

            if s >= end:
                break

            overlap += max(0, min(end, e) - max(start, s))

        owned = max(0, length - overlap)

        class_bp[category] += owned
        overlap_bp += overlap

        claimed[chrom].append((start, end))
        claimed[chrom].sort()

    masked_bp = sum(class_bp.values())
    unmasked_bp = GENOME_SIZE - masked_bp

    print(f"\n{pop}\nGenome composition\n")
    print(f"Genome size : {GENOME_SIZE:,}")
    print(f"Masked bp   : {masked_bp:,}")
    print(f"Unmasked bp : {unmasked_bp:,}")
    print(f"Overlap bp  : {overlap_bp:,}\n")

    print("Category\tbp\tpercent")

    for cat, bp in sorted(class_bp.items(),
                        key=lambda x: x[1],
                        reverse=True):

        pct = bp / GENOME_SIZE * 100

        print(f"{cat}\t{bp}\t{pct:.2f}")

    print(f"Unmasked\t{unmasked_bp}\t{100*unmasked_bp/GENOME_SIZE:.2f}")


    # Add unmasked bases
    plot_data = dict(class_bp)
    plot_data["Unmasked"] = unmasked_bp

    # Remove small categories (<0.1% genome)
    plot_data = {
        k: v
        for k, v in plot_data.items()
        if v / GENOME_SIZE * 100 > 0.1
    }

    labels = list(plot_data.keys())
    sizes = list(plot_data.values())

    plt.figure(figsize=(8, 8))


    TE_COLORS = {
        "LTR": "darkblue",
        "ClassI": "blue",
        "LINE": "#588bc7",
        "Non-LTR": "#588bc7",
        "nLTR": "#588bc7",
        "SINE": "#d2d2f7",
        "TIR": "#d99d1c",
        "MITE": "#d99d1c",
        "DNA":"orange",
        "ClassII":"orange",
        "HELITRON": "#e8d258",
        "Helitron": "#e8d258",
        "RC": "#e8d258",
        "Unknown": "grey",
        "unknown": "grey",
        "Other": "#a28ea3",
        "Low_complexity": "purple",
        "Simple_repeat": "pink",
        "Unmasked": "black",
        "Satellite": "green"
    }

    colors = [TE_COLORS.get(cat, "#cccccc")for cat in labels] 

    plt.pie(
        sizes,
        labels=labels,
        #autopct=lambda p: f"{p:.1f}%",
        startangle=90,
        colors = colors
    )


    plt.title(f"TransposonUltimate {pop}")
    plt.tight_layout()
    plt.savefig(f"{pop}TE_genomic_TransposonUltimate.png", dpi=300)
    #plt.show()


BED = "/global/group/research/m_lebh01/TEpaperCorrection/TEannotationBed/TransposonUltimate/DK.bed"
pops = ["DK", "ES", "FI", "SE", "UA", "TR", "ZI"]

for pop in pops:
    GENOME_SIZE = get_genome_len(f"/global/students/research/m_lebh01/DESwoMAN_Oct2025/TEData_circos_MCHelper/MCHelperGenomes/Dmel{pop}.fa")
    BED = f"/global/group/research/m_lebh01/TEpaperCorrection/TEannotationBed/TransposonUltimate/{pop}.bed"
    main(BED,GENOME_SIZE, pop)
