uv run fetch_uniprot_proteins.py --query "(reviewed:true) AND (organism_id:9606)" --output reviewed_human_proteins.csv
pip install esm@git+https://github.com/Biohub/esm.git@main

curl -L "https://data.omim.org/downloads/<API Key>/mimTitles.txt" -o mimTitles.txt