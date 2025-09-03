import h5py
import numpy as np
import pandas as pd

## QUESTION 1 AND 2 

SRR = ("SRR28910565", "SRR28910566", "SRR28910567", "SRR28910568")

def print_structure(name, obj):
    print(name)

for srr in SRR:
    file_path = f"/7_SepalAI/3_Project_brunner/Tanycytes/data/{srr}/outs/filtered_feature_bc_matrix.h5"
    print(f"\nProcessing file: {file_path}")

    try:
        with h5py.File(file_path, "r") as f:
            matrix_data = np.array(f["matrix/data"])
            features = np.array(f["matrix/features/id"])
            barcodes = np.array(f["matrix/barcodes"])
            print("\nMatrix data shape:", matrix_data.shape)
            print("Number of features:", len(features))
            print("Number of barcodes:", len(barcodes))
    except KeyError as e:
        print(f"Could not find expected dataset: {e}")
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"An error occurred while processing {srr}: {e}")

## QUESTION 3
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

srr = "SRR28910565"
file_path = f"/7_SepalAI/3_Project_brunner/Tanycytes/data/{srr}/outs/filtered_feature_bc_matrix.h5"

with h5py.File(file_path, "r") as f:
    data = f["matrix/data"][:]
    indices = f["matrix/indices"][:]
    indptr = f["matrix/indptr"][:]
    shape = f["matrix/shape"][:]

    matrix = sp.csc_matrix((data, indices, indptr), shape=shape).T  

    features = f["matrix/features/name"][:].astype(str)
    barcodes = f["matrix/barcodes"][:].astype(str)

avg_expression = np.asarray(matrix.mean(axis=0)).ravel()

max_idx = avg_expression.argmax()
top_feature = features[max_idx]
top_value = avg_expression[max_idx]

print(f"Top feature in {srr}: {top_feature} (average expression = {top_value:.2f})")


# QUESTION 4 
import h5py
import numpy as np
import scipy.sparse as sp

SRR = ("SRR28910565", "SRR28910566", "SRR28910567", "SRR28910568")

average_nfeatures = {}

for srr in SRR:
    file_path = f"/7_SepalAI/3_Project_brunner/Tanycytes/data/{srr}/outs/filtered_feature_bc_matrix.h5"
    try:
        with h5py.File(file_path, "r") as f:
            data = f["matrix/data"][:]
            indices = f["matrix/indices"][:]
            indptr = f["matrix/indptr"][:]
            shape = f["matrix/shape"][:]
            matrix = sp.csc_matrix((data, indices, indptr), shape=shape)
            nfeatures_per_cell = np.array((matrix > 0).sum(axis=0)).ravel()
            average_nfeatures[srr] = nfeatures_per_cell.mean()
    
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"Error processing {srr}: {e}")

sorted_samples = sorted(average_nfeatures.items(), key=lambda x: x[1], reverse=True)

for srr, avg in sorted_samples:
    print(f"{srr}: average nFeatures per cell = {avg:.2f}")


# QUESTION 5
import h5py
import numpy as np
import scipy.sparse as sp
import pandas as pd

srr = "SRR28910567"
file_path = f"/7_SepalAI/3_Project_brunner/Tanycytes/data/{srr}/outs/filtered_feature_bc_matrix.h5"

with h5py.File(file_path, "r") as f:
    data = f["matrix/data"][:]
    indices = f["matrix/indices"][:]
    indptr = f["matrix/indptr"][:]
    shape = f["matrix/shape"][:]

    matrix = sp.csc_matrix((data, indices, indptr), shape=shape)
    features = f["matrix/features/name"][:].astype(str)

matrix = matrix.tocsc()
mean_per_feature = np.array(matrix.mean(axis=1)).ravel()  # mean across cells
squared_sum = np.array(matrix.power(2).mean(axis=1)).ravel()
var_per_feature = squared_sum - mean_per_feature**2

dispersion = var_per_feature / (mean_per_feature + 1e-8)  # avoid division by zero

df = pd.DataFrame({
    "feature": features,
    "mean": mean_per_feature,
    "variance": var_per_feature,
    "dispersion": dispersion
})

df_sorted = df.sort_values("dispersion", ascending=False)
top10 = df_sorted.head(10)
print(top10)
