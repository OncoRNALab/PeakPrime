#!/usr/bin/env python3
import pandas as pd
import argparse

def load_profile(path):
    df = pd.read_csv(path, sep="\t", comment="#")
    # Drop the first two columns ("bin labels" and "bins")
    df = df.drop(columns=df.columns[:2])
    return df

def compute_peak(df, method="max", window=None):
    """Compute peak signal from profile dataframe"""
    if method == "max":
        return df.max(axis=1)   # max across bins for each sample
    elif method == "mean":
        if window is None:
            raise ValueError("Must specify window=(start,end) for mean method")
        return df.iloc[:, window[0]:window[1]].mean(axis=1)
    else:
        raise ValueError("method must be 'max' or 'mean'")

def main():
    parser = argparse.ArgumentParser(description="Compute TES/TSS ratio from deepTools profile tables")
    parser.add_argument("--tss", required=True, help="Profile file for TSS (from plotProfile --outFileNameData)")
    parser.add_argument("--tes", required=True, help="Profile file for TES (from plotProfile --outFileNameData)")
    parser.add_argument("-o", "--out", default="tes_tss_ratios.txt", help="Output table")

    args = parser.parse_args()

    tss = load_profile(args.tss)
    tes = load_profile(args.tes)

    tss_peak = compute_peak(tss, method="max")
    tes_peak = compute_peak(tes, method="max")

    ratio = tes_peak / tss_peak

    out = pd.DataFrame({
        "TSS_peak": tss_peak,
        "TES_peak": tes_peak,
        "TES_TSS_ratio": ratio
    })
    out.to_csv(args.out, sep="\t", index=False)
    print(f"✅ Wrote ratios to {args.out}")

if __name__ == "__main__":
    main()
