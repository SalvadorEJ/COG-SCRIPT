#!/usr/bin/env python3
# ============================================================
# COG functional profile from any eggNOG-mapper output file
# ============================================================
#
# This script reads an eggNOG-mapper annotation file and:
#   1. Finds the "COG_category" column.
#   2. Counts COG functional categories.
#   3. Correctly handles multiple categories such as FP, KL, CH.
#   4. Calculates percentages.
#   5. Saves a CSV summary table.
#   6. Creates a publication-quality PNG bar plot.
#
# Supported input formats:
#   - Excel: .xlsx, .xls
#   - Tab-separated text: .tsv, .txt
#
# The input file only needs to contain a column named:
#       COG_category
#
# Requirements:
#   pip install pandas openpyxl matplotlib
#
# Usage:
#   1. Put this script in the same folder as your eggNOG file.
#   2. Change INPUT_FILE below.
#   3. Run:
#
#          python COG_profile_eggNOG.py
#
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import sys


# ------------------------------------------------------------
# USER SETTINGS
# ------------------------------------------------------------

# Name of the eggNOG-mapper output file.
# Change this to your own file.
INPUT_FILE = "eggnog_annotations.xlsx"

# Output files will automatically use the input filename.
# Example:
#   genome.emapper.annotations.xlsx
# becomes:
#   genome.emapper.annotations_COG_categories.csv
#   genome.emapper.annotations_COG_categories.png


# ------------------------------------------------------------
# STANDARD COG FUNCTIONAL CATEGORIES
# ------------------------------------------------------------

COG_NAMES = {
    "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
}


# ------------------------------------------------------------
# READ INPUT FILE
# ------------------------------------------------------------

input_path = Path(INPUT_FILE)

if not input_path.exists():
    sys.exit(
        f"\nERROR: Input file was not found:\n{input_path.resolve()}\n"
        "\nPlease change INPUT_FILE at the beginning of the script.\n"
    )

suffix = input_path.suffix.lower()

try:
    if suffix in [".xlsx", ".xls"]:
        # Read the first worksheet.
        df = pd.read_excel(input_path)

    elif suffix in [".tsv", ".txt"]:
        # eggNOG-mapper tabular output is usually tab-separated.
        df = pd.read_csv(
            input_path,
            sep="\t",
            comment="#",
            dtype=str
        )

    else:
        sys.exit(
            "ERROR: Unsupported file format.\n"
            "Use .xlsx, .xls, .tsv or .txt"
        )

except Exception as e:
    sys.exit(f"ERROR while reading the input file:\n{e}")


# ------------------------------------------------------------
# FIND THE COG_CATEGORY COLUMN
# ------------------------------------------------------------

# First look for the exact column name.
if "COG_category" in df.columns:
    cog_column = "COG_category"

else:
    # Try a case-insensitive search.
    possible_columns = [
        col for col in df.columns
        if str(col).lower() == "cog_category"
    ]

    if possible_columns:
        cog_column = possible_columns[0]

    else:
        print("\nERROR: The file does not contain a 'COG_category' column.")
        print("\nAvailable columns:")
        for col in df.columns:
            print(f"  - {col}")

        sys.exit(
            "\nPlease verify that the file is an eggNOG-mapper annotation "
            "file containing COG_category.\n"
        )


# ------------------------------------------------------------
# CLEAN COG CATEGORY DATA
# ------------------------------------------------------------

cog = (
    df[cog_column]
    .fillna("")
    .astype(str)
    .str.strip()
)

# Values such as "-" indicate that no COG category was assigned.
cog = cog.replace(
    ["-", "NA", "nan", "*", ""],
    ""
)

# Remove whitespace.
cog = cog.str.replace(r"\s+", "", regex=True)


# ------------------------------------------------------------
# SPLIT MULTIPLE COG CATEGORIES
# ------------------------------------------------------------

# eggNOG can assign more than one COG category to a protein.
#
# Examples:
#
#   F   -> F
#   FP  -> F + P
#   KL  -> K + L
#   CH  -> C + H
#
# Each letter is therefore counted as an individual COG
# functional category.

all_categories = []

for value in cog:

    if value == "":
        continue

    for letter in value:

        if letter in COG_NAMES:
            all_categories.append(letter)


# ------------------------------------------------------------
# COUNT CATEGORIES
# ------------------------------------------------------------

counts = pd.Series(all_categories).value_counts()

# Force the standard COG order, including categories with zero.
result = pd.DataFrame({
    "COG": list(COG_NAMES.keys()),
    "Category": list(COG_NAMES.values()),
})

result["Proteins"] = result["COG"].map(counts).fillna(0).astype(int)


# ------------------------------------------------------------
# CALCULATE PERCENTAGES
# ------------------------------------------------------------

total_assignments = result["Proteins"].sum()

if total_assignments > 0:
    result["Percent"] = (
        result["Proteins"] / total_assignments * 100
    ).round(2)

else:
    result["Percent"] = 0.0


# ------------------------------------------------------------
# CREATE OUTPUT FILENAMES
# ------------------------------------------------------------

base_name = input_path.stem

output_csv = input_path.parent / (
    f"{base_name}_COG_categories.csv"
)

output_png = input_path.parent / (
    f"{base_name}_COG_categories.png"
)


# ------------------------------------------------------------
# SAVE SUMMARY TABLE
# ------------------------------------------------------------

result.to_csv(
    output_csv,
    index=False
)


# ------------------------------------------------------------
# CREATE BAR PLOT
# ------------------------------------------------------------

# Use a large figure suitable for publication/manuscript work.
fig, ax = plt.subplots(figsize=(12, 7))

bars = ax.bar(
    result["COG"],
    result["Proteins"]
)

ax.set_xlabel(
    "COG functional category",
    fontsize=12
)

ax.set_ylabel(
    "Number of proteins",
    fontsize=12
)

ax.set_title(
    "COG functional profile",
    fontsize=14
)

# Add the number of proteins above each bar.
for bar, value in zip(bars, result["Proteins"]):

    if value > 0:

        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value,
            str(value),
            ha="center",
            va="bottom",
            fontsize=8
        )


# Keep the plot clean.
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.xticks(
    rotation=0,
    fontsize=10
)

plt.yticks(
    fontsize=10
)

plt.tight_layout()


# Save at 300 dpi.
fig.savefig(
    output_png,
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)


# ------------------------------------------------------------
# PRINT SUMMARY
# ------------------------------------------------------------

total_proteins = len(df)

proteins_with_cog = (cog != "").sum()
proteins_without_cog = (cog == "").sum()

print("\n============================================")
print("COG functional profile")
print("============================================\n")

print(f"Input file: {input_path}")
print(f"Total proteins/rows: {total_proteins}")
print(f"Proteins with COG category: {proteins_with_cog}")
print(f"Proteins without COG category: {proteins_without_cog}")
print(f"Total COG category assignments: {total_assignments}")

print("\nCOG categories:")
print(result.to_string(index=False))

print("\nOutput files:")
print(f"  {output_csv}")
print(f"  {output_png}")

print("\nAnalysis completed successfully.")
