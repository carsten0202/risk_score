# Script for outputting the rankings based on the frequencies of Klitz 2003

import pandas as pd

# Sample Data: Table A
data_A = {
    'key':   ["03:0X 03:02", "01:02 06:02", "05:01 02:01",
			  "02:01 02:02", "05:05 03:01", "01:0X 05:01", 
			  "03:0X 03:01", "01:03 06:03", "02:01 03:03",
			  "04:01 04:02", "01:0X 05:03", "03:02 03:03",
			  "01:02 06:09", "01:03 06:01"],
    'value': [0.154, 0.1437, 0.1314,
			  0.1108, 0.1106, 0.1088,
			  0.0648, 0.0567, 0.0366,
			  0.0227, 0.0208, 0.0079,
			  0.0071, 0.0066],
}


table_A = pd.DataFrame(data_A)

# Sample Data: Table B
data_B = {
	'key1': ["05:01 02:01", "03:0X 03:01", "03:0X 03:02",
             "03:0X 03:01", "05:01 02:01", "03:0X 03:01",
    		 "04:01 04:02", "05:01 02:01", "05:01 02:01",
             "05:01 02:01", "01:02 06:02", "05:01 02:01",
    		 "05:01 02:01", "03:0X 03:02", "03:0X 03:01",
			 "03:02 03:03", "03:0X 03:01", "05:01 02:01"],
    'key2': ["03:0X 03:02", "01:0X 05:01", "03:0X 03:02",
			 "03:0X 03:02", "03:0X 03:01", "02:01 02:02",
			 "02:01 02:02", "01:0X 05:01", "02:01 02:02",
			 "05:01 02:01", "01:02 06:02", "05:05 03:01",
			 "01:03 06:03", "01:03 06:03", "03:0X 03:01",
			 "03:0X 03:02", "05:05 03:01", "01:02 06:02",],
 	'beta': [3.63, 3.14, 2.31,
			 2.16, 1.08, 1.06, 
			 0.61, 0.17, 0.17,
			 0.05, -0.24, -0.51,
			 -0.65, -0.69, -0.87,
			 -1.15, -1.47, -1.94],
}
table_B = pd.DataFrame(data_B)

# Merge Table A with Table B to get values corresponding to key1
merged_1 = pd.merge(table_B, table_A, left_on='key1', right_on='key')
merged_1 = merged_1.rename(columns={'value': 'value1'}).drop('key', axis=1)

# Merge again to get values corresponding to key2
merged_2 = pd.merge(merged_1, table_A, left_on='key2', right_on='key')
merged_2 = merged_2.rename(columns={'value': 'value2'}).drop('key', axis=1)

# Calculate the product of values for each pair
merged_2['product'] = merged_2['value1'] * merged_2['value2']

# Resulting DataFrame with products
result = merged_2[['key1', 'key2', 'product', 'beta']]

# Sort the result by the 'product' column in ascending order
result_sorted = result.sort_values(by='product', ascending=False)

print(result_sorted)
