import os
import pandas as pd

PATIENTS = ["NC1", "NC6", "NC8", "NC83", "NC94", "NC183"]
BASE_DIR = "~/WeavePop_Chronic_Crypto"
RESULTS_DIR = os.path.join(BASE_DIR, "results_chronic_1/01.Samples/cnv")
MAPPING_FILE = os.path.join(BASE_DIR, "sample_name_mapping.csv")

for patient in PATIENTS:
    # 1. Read patient samples file
    samples_file = os.path.join(BASE_DIR, f"config/{patient}_samples.txt")
    if not os.path.exists(samples_file):
        print(f"Missing {samples_file}, skipping {patient}")
        continue
    samples_df = pd.read_csv(samples_file, header=None, names=["sample", "doi"])
    all_cnv = []
    for sample in samples_df["sample"]:
        cnv_file = os.path.join(RESULTS_DIR, sample, "cnv_calls.tsv")
        if not os.path.exists(cnv_file):
            print(f"Missing {cnv_file}, skipping sample {sample}")
            continue
        cnv_df = pd.read_csv(cnv_file, sep="\t")
        cnv_df = cnv_df[cnv_df["repeat_fraction"] < 0.5]
        # Drop specified columns if they exist
        drop_cols = ["depth", "norm_depth", "overlap_bp", "smooth_depth"]
        cnv_df = cnv_df.drop(columns=[col for col in drop_cols if col in cnv_df.columns])
        all_cnv.append(cnv_df)
    if not all_cnv:
        print(f"No CNV data for {patient}")
        continue
    concat_df = pd.concat(all_cnv, ignore_index=True)
    # 2. Merge DOI
    concat_df = concat_df.merge(samples_df, left_on="sample", right_on="sample", how="left")
    concat_df = concat_df.rename(columns={"doi": "DOI"})
    # 3. Save concatenated file
    out_file = os.path.join(BASE_DIR, f"pipeline_test/{patient}_concatenated_sv_genes.tsv")
    concat_df.to_csv(out_file, sep="\t", index=False)

    # 4. Merge with sample renaming
    mapping_df = pd.read_csv(MAPPING_FILE)
    merged_df = concat_df.merge(mapping_df, left_on="sample", right_on="seq_id", how="left")
    merged_df = merged_df.drop(columns=["sample", "seq_id"])
    merged_df = merged_df.fillna("n/a")
    # 5. Collapse identical rows
    collapsed_df = merged_df.groupby(
        merged_df.columns.difference(['strain_id', 'DOI']).tolist()
    ).agg({'strain_id': ', '.join, 'DOI': 'min'}).reset_index()
    # 6. Reorder columns for output
    desired_order = [
        "accession", "start", "end", "region_size", "repeat_fraction",
        "cnv", "strain_id", "DOI", "feature_id"
    ]
    final_cols = [col for col in desired_order if col in collapsed_df.columns]
    collapsed_df = collapsed_df[final_cols]
    # Sort rows by accession, start, end
    sort_cols = [col for col in ["accession", "start", "end"] if col in collapsed_df.columns]
    collapsed_df = collapsed_df.sort_values(by=sort_cols)
    # 7. Save collapsed file
    collapsed_file = os.path.join(BASE_DIR, f"pipeline_test/{patient}_collapsed_sv_genes.tsv")
    collapsed_df.to_csv(collapsed_file, sep="\t", index=False)
    print(f"Processed {patient}: {len(all_cnv)} samples, saved to {collapsed_file}")
