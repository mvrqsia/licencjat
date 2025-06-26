import numpy as np
import pandas as pd
from Bio import SeqIO
from pathlib import Path

# Nazwy plików zliczeń RT po filtracji
files = {
    "plus": "/Users/mariachmielorz/Desktop/RTS_counts/plus2_sense.tab",
    "minus": "/Users/mariachmielorz/Desktop/RTS_counts/minus1_sense.tab"
}

# Przetwarzanie danych: konwersja do TSV
for label, path in files.items():
    with open(path) as f:
        lines = f.readlines()
    lines = [l.strip() for l in lines if l.strip()]
    table = ["id\tbase\tcount\n"]
    for l in lines:
        if l.startswith("NC_"):
            id = l
        else:
            base, count = l.split()[:2]
            table.append(f"{id}\t{base}\t{count}\n")
    with open(f"/Users/mariachmielorz/Desktop/RTS_counts/{label}.tsv", "w") as f:
        f.writelines(table)

# Wczytanie TSV do Pandas
plus = pd.read_csv("/Users/mariachmielorz/Desktop/RTS_counts/plus.tsv", sep="\t")
minus = pd.read_csv("/Users/mariachmielorz/Desktop/RTS_counts/minus.tsv", sep="\t")

# Wczytanie sekwencji referencyjnej
ref = list(SeqIO.parse("/Users/mariachmielorz/Desktop/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna", "fasta"))

# Generowanie sygnału DMS
Path("/Users/mariachmielorz/Desktop/DMS_signal/").mkdir(parents=True, exist_ok=True)

# Przypisanie odpowiedniego chromosomu 
for r in ref:
    chr_id = r.id  
    plus_chr = plus[plus["id"] == chr_id]
    minus_chr = minus[minus["id"] == chr_id]

    if plus_chr.empty or minus_chr.empty:
        continue

    
    assert (plus_chr["base"].values == minus_chr["base"].values).all()

    count_plus = plus_chr["count"].astype(int).values
    count_minus = minus_chr["count"].astype(int).values

    # Różnica sygnału (bez wartości ujemnych)
    count = np.maximum(count_plus - count_minus, 0)

    # Normalizacja (2–8% górnych sygnałów)
    top10 = np.sort(count)[-int(0.1 * len(count)):]
    signal = count / np.mean(top10[int(0.2 * len(top10)):])

    df = pd.DataFrame({
    "pos": np.arange(1, len(signal) + 1),  # dodaje numer pozycji w chromosomie
    "base": plus_chr["base"].values,
    "signal": signal
    })


    # Wyzeruj sygnały na G i T
    df.loc[df["base"] == "G", "signal"] = 0
    df.loc[df["base"] == "T", "signal"] = 0

    df.to_csv(f"/Users/mariachmielorz/Desktop/DMS_signal/{chr_id}.tsv", sep="\t", index=False)
