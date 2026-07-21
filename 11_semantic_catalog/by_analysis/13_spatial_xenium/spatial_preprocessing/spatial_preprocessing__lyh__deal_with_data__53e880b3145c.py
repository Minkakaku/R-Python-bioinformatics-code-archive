import sys
import pandas as pd
import numpy as np
import scanpy as sc
import os
import warnings
import argparse

# Suppress warnings
warnings.filterwarnings("ignore")

def parse_arguments():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(description="Process Xenium data with ROI filtering and formatting.")

    parser.add_argument("-i", "--input_dir", type=str, required=True,
                        help="Path to the Xenium output folder (containing cell_feature_matrix and cells.parquet)")
    
    parser.add_argument("-o", "--output_dir", type=str, required=True,
                        help="Directory to save the output files")
    
    parser.add_argument("-r", "--roi_dir", type=str, required=True,
                        help="Directory containing ROI CSV files (e.g., ./67610)")
    
    parser.add_argument("-g", "--group_file", type=str, default=None,
                        help="(Optional) Path to the group mapping CSV file")
    
    parser.add_argument("-c", "--chip_id", type=str, default="Unknown_Chip",
                        help="(Optional) Chip ID string to add to metadata")

    return parser.parse_args()

def normalize_path(path_value):
    """
    Resolve a path to an absolute path for stable access after cwd changes.
    """
    if path_value is None:
        return None
    return os.path.abspath(path_value)

def get_sample_name(filename):
    """
    Strip only the final extension from an ROI filename.
    """
    return os.path.splitext(filename)[0]

def read_xenium(input_dir):
    """
    Reads Xenium data into an AnnData object.
    """
    print(f"Reading Xenium data from: {input_dir}")
    
    # Path construction with validation
    tmp_path = os.path.join(input_dir, "cell_feature_matrix")
    meta_path = os.path.join(input_dir, 'cells.parquet')
    
    if not os.path.exists(tmp_path):
        raise FileNotFoundError(f"10x matrix directory not found: {tmp_path}")
    if not os.path.exists(meta_path):
        raise FileNotFoundError(f"Metadata parquet file not found: {meta_path}")
    
    # Read expression matrix
    adata = sc.read_10x_mtx(
        path=tmp_path,
        var_names="gene_symbols",
        make_unique=True 
    )
    
    # Read and merge metadata efficiently
    meta = pd.read_parquet(meta_path).set_index('cell_id')
    adata.obs = adata.obs.join(meta, how='left')
    
    # Store spatial coordinates (convert to numpy array)
    adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy()
    # Flip Y axis for correct visualization
    adata.obsm['spatial'][:, 1] = -adata.obsm['spatial'][:, 1]
    
    # Preserve raw counts
    adata.layers["counts"] = adata.X.copy()
    
    return adata

def main():
    # 1. Parse arguments
    args = parse_arguments()
    args.input_dir = normalize_path(args.input_dir)
    args.output_dir = normalize_path(args.output_dir)
    args.roi_dir = normalize_path(args.roi_dir)
    args.group_file = normalize_path(args.group_file)

    # 2. Setup Output Directory
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        print(f"Created output directory: {args.output_dir}")
    
    # Change working directory to output directory (optional, but follows original logic)
    os.chdir(args.output_dir)
    print("Current working directory set to: " + os.getcwd())

    # 3. Load Xenium Data
    try:
        adata = read_xenium(args.input_dir)
    except Exception as e:
        print(f"Error loading Xenium data: {e}")
        sys.exit(1)

    # 4. Process ROI CSV Files
    print(f"Processing ROI files from: {args.roi_dir}")
    if not os.path.exists(args.roi_dir):
        print(f"Error: ROI directory {args.roi_dir} does not exist.")
        sys.exit(1)

    csv_files = [f for f in os.listdir(args.roi_dir) if f.endswith('.csv')]
    
    if not csv_files:
        print("Warning: No CSV files found in the ROI directory.")
    
    # Generate group list file based on filenames
    cleaned_list = [get_sample_name(csv_file) for csv_file in csv_files]
    df_groups = pd.DataFrame(cleaned_list, columns=['sample'])
    df_groups.to_csv("group.csv", index=False)
    print("Saved group.csv generated from filenames.")

    # Read and concatenate ROI CSVs
    df_list = []
    for filename in csv_files:
        file_path = os.path.join(args.roi_dir, filename)
        try:
            tmp = pd.read_csv(file_path)
            tmp.columns = [str(col).strip() for col in tmp.columns]
            sample_name = get_sample_name(filename)
            cell_id_column = next(
                (col for col in tmp.columns if col.strip().lower() == "cell id"),
                None
            )

            if cell_id_column is None:
                print(f"Warning: 'Cell ID' column not found in {filename}, skipping.")
                continue

            tmp = tmp[[cell_id_column]].copy()
            tmp = tmp.rename(columns={cell_id_column: "Cell ID"})
            tmp["Cell ID"] = tmp["Cell ID"].astype(str).str.strip()
            tmp = tmp[tmp["Cell ID"] != ""]
            tmp["sample"] = sample_name
            df_list.append(tmp)

        except Exception as e:
            print(f"Error reading {filename}: {e}")

    if not df_list:
        print("Error: No valid ROI data found. Exiting.")
        sys.exit(1)

    df = pd.concat(df_list, ignore_index=True)
    df = df.drop_duplicates(subset=['Cell ID'])

    # 5. Filter AnnData based on ROIs
    print("Filtering AnnData object based on ROIs...")
    adata.obs["Cell ID"] = adata.obs.index.astype(str)
    mask = adata.obs['Cell ID'].isin(df['Cell ID'])
    adata = adata[mask, :].copy()
    
    # Merge ROI info into adata.obs
    adata.obs = adata.obs.merge(df, on='Cell ID', how='left', validate='one_to_one')
    adata.obs = adata.obs.set_index("Cell ID")

    # 6. Process Group File (Optional)
    if args.group_file and os.path.exists(args.group_file):
        print(f"Merging external group file: {args.group_file}")
        group_df = pd.read_csv(args.group_file)
        if "sample" not in group_df.columns:
            print("Warning: Group file does not contain a 'sample' column. Skipping merge.")
        else:
            group_df["sample"] = group_df["sample"].astype(str).str.strip()
        
            # Reset index to merge
            adata.obs["Cell ID"] = adata.obs.index
            adata.obs = adata.obs.merge(
                group_df,
                how="left",
                on="sample",
                validate="many_to_one",
                suffixes=("", "_group")
            )
            adata.obs = adata.obs.set_index("Cell ID")
    elif args.group_file:
        print(f"Warning: Provided group file path does not exist: {args.group_file}")

    # 7. Save Outputs
    print("Saving outputs...")
    
    # Add Chip ID
    adata.obs["chip_id"] = args.chip_id
    # Save raw h5ad
    adata.write_h5ad("./raw.h5ad")

    print("Processing complete.")

if __name__ == "__main__":
    main()
