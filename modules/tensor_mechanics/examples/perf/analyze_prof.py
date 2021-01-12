import pandas as pd
pd.set_option("display.max_rows", None, "display.max_columns", None)

# profiles to analyze
nonAD_profile = 'reg'
AD_profile = 'ad'
AD_optimized_profile = 'ad_opt'
comparison_filename = 'comparison.csv'

# read in the timings
def read_prof(filename):
    with open(filename) as f:
        lines = f.readlines()
        table = pd.DataFrame(columns=['function', filename])
        for row in range(5, len(lines)):
            line_stripped = lines[row].strip('\n')
            entries = line_stripped.split(None, 5)
            table.loc[len(table)] = [entries[5].replace('AD',''), float(entries[0].strip('s'))]
        return table


nonAD_data = read_prof(nonAD_profile)
AD_data = read_prof(AD_profile)
AD_optimized_data = read_prof(AD_optimized_profile)

data = nonAD_data.set_index('function').join(AD_data.set_index('function'), how='outer')
data = data.join(AD_optimized_data.set_index('function'), how='outer')
data = data.fillna(0)
data = data.sort_values(by=AD_profile, ascending=False)
data.to_csv(comparison_filename)
