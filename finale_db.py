"""
Retrieve data from the Finaledb database and save it to a CSV file.
"""
import requests
import pandas as pd

# List of IDs to query
ids = []
for i in range(85754, 86003, 1):
    ids.append('EE' + str(i))

# Output CSV file name
output_file = 'jiang.csv'

# Initialize a list to store all result dicts
info_table = pd.DataFrame()

for seq_id in ids:
    url = f'http://finaledb.research.cchmc.org/api/v1/seqrun?id={seq_id}'
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()

        if 'results' in data:
            # If data['results'] is a list, we take the first element, else we handle it as a single result
            try:
                if not isinstance(data['results'], list):
                    result = data['results']
                else:
                    result = data['results'][0]
            except:
                print(f"Error processing results for ID {seq_id}: {data['results']}")
            
            # Convert the result to a DataFrame
            result_dict = {
                'ID': 'EE' + str(result.get('id')),
                'sample_name': result['sample']['name'],
                'assay': result['assay'],
                'readlen': result['seqConfig']['readlen'],
                'instrument': result['seqConfig']['instrument'],
                'seq_layout': result['seqConfig']['seq_layout'],
                'fragNum': result['fragNum']['hg38'],
                'publication_name': result['publication']['citeShort'],
                'publication_doi': result['publication']['identifiers']['doi'],
                'tissue': result['sample']['tissue'],
                'disease': result['sample']['disease'],
            }
            # Append the result DataFrame to the info_table DataFrame
            info_table = pd.concat([info_table, pd.DataFrame([result_dict])], ignore_index=True)

        else:
            print(f"No results found for ID {seq_id}")

    except requests.exceptions.RequestException as e:
        print(f"Request failed for ID {seq_id}: {e}")

info_table.to_csv(output_file, index=False)
print(f"Data saved to {output_file}")
