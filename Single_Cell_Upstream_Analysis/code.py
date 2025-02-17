from scipy.io import mmread
import pandas as pd

# Paths to the files
matrix_file = "/Volumes/Jagannath/Projects/Single_Cell_Upstream/single_cell/simpleaf_quant/af_quant/alevin/quants_mat.mtx"
rows_file = "/Volumes/Jagannath/Projects/Single_Cell_Upstream/single_cell/simpleaf_quant/af_quant/alevin/quants_mat_rows.txt"
cols_file = "/Volumes/Jagannath/Projects/Single_Cell_Upstream/single_cell/simpleaf_quant/af_quant/alevin/quants_mat_cols.txt"

# Load the matrix
count_matrix = mmread(matrix_file).tocsc()  # Load as a sparse matrix

# Load rows (genes) and columns (cell barcodes)
genes = pd.read_csv(rows_file, header=None)[0].tolist()
cell_barcodes = pd.read_csv(cols_file, header=None)[0].tolist()

# Convert sparse matrix to a dense DataFrame (optional, for analysis)
count_df = pd.DataFrame(count_matrix.toarray(), index=genes, columns=cell_barcodes)

# Display the data
print(count_df.head())

# Save as a CSV (optional)
count_df.to_csv("gene_count_matrix.csv")
