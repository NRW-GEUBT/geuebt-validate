#!/usr/bin/env bash
set -Eeu


DB_URL=https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
DB_PATH="$HOME/.nrw-geuebt/taxdump"
FILENAME="taxdump"

echo "Downloading taxdump database to ${DB_PATH} and extracting"

# Create Subdirectories
mkdir -p ${DB_PATH}
cd ${DB_PATH}
# Download
target="$FILENAME.tar.gz"
if [[ -f $target ]]; then
  echo "The file $target already exists. Skipping download."
else
  wget --output-document $target $DB_URL
  [[ $? -eq 0 ]] && [[ -s $target ]] && download_hash=$(openssl dgst -r -sha256 $target) && download_success=1
fi
[[ -n $download_hash ]] && echo "$download_hash" > $FILENAME.sha256
date --iso-8601='minutes' >> $FILENAME.timestamp
echo $DB_URL > $FILENAME.source
[[ "$download_success" == 1 ]] && tar -xzvf $FILENAME.tar.gz && rm $FILENAME.tar.gz
echo "Downlaod complete"