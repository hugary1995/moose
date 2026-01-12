import pandas as pd

threshold = 5e8  # 0.5 seconds in nanoseconds

data = pd.read_csv("tmp.csv")

# Sort by compile_time descending and filter by threshold
data_sorted = data.sort_values(by="compile_time", ascending=False)
data_sorted = data_sorted[data_sorted["compile_time"] > threshold]

# Find the relative directory relative to MOOSE root
MOOSE_ROOT = "/home/gary/projects/moose"
data_sorted["relative_path"] = data_sorted["header_file"].apply(
    lambda x: x.replace(MOOSE_ROOT + "/", "").rsplit("/", 1)[0]
)

# Print out files exceeding the threshold grouped by directory
grouped = data_sorted.groupby(data_sorted["relative_path"])
for dir_path, group in grouped:
    print(f"Directory: {dir_path}")
    for _, row in group.iterrows():
        print(f"  {row['header_file']}: {row['compile_time']/1e9:.6f} seconds")
    print()
