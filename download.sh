#!/bin/bash
# Download Data
curl -L "https://www.dropbox.com/scl/fi/4qsos5hqfavs09f4b8e2s/visiumhd_nanocourse.zip?rlkey=r3rfprfwrqc6dnvwyuu5zzho9&st=1shp9z0w&dl=1" -o visiumhd_nanocourse.zip

# Move files
unzip visiumhd_nanocourse.zip
mv visiumhd_nanocourse/data_processed lessons/
rm -rf visiumhd_nanocourse
rm -rf __MACOSX
