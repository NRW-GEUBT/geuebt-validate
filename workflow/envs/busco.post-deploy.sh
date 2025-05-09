#!/usr/bin/env bash
set -Eeu

# URL of the .tar archive
# If changing this path need to change the db name in `rules/common.smk` !
url="https://busco-data.ezlab.org/v5/data/lineages/bacteria_odb10.2024-01-08.tar.gz"

# Local directory to save the file
local_dir="${CONDA_PREFIX}/busco/"

# Name of the downloaded file
file_name="bacteria_odb10"

# Remote file (includes .tar.gz extension)
remote_name=$(basename "$url")

# Create local directory if it doesn't exist
mkdir -p "$local_dir"

# Check if the file or directory already exists locally
if [[ -e "$local_dir/$file_name" ]]; then

  echo "The file already exists, skipping download"

else

  # Download the file and save it to the local directory
  echo "Downloading file..."
  curl -L -o "$local_dir/$remote_name" "$url"

  # Calculate the hash of the downloaded file
  download_hash=$(openssl dgst -r -sha256 "$local_dir/$remote_name")
  
  echo "$download_hash" > "$local_dir/$file_name.sha256"
  date --iso-8601='minutes' > "$local_dir/$file_name.timestamp"
  echo "$url" > "$local_dir/$file_name.source"
  tar -xzvf "$local_dir$remote_name" -C "$local_dir/" && rm "$local_dir/$remote_name"

fi
