#!/bin/bash
mkdir -p resources/vep_cache
wget -c http://ftp.ensembl.org/pub/release-104/variation/indexed_vep_cache/homo_sapiens_vep_104_GRCh38.tar.gz -O resources/vep_cache/cache.tar.gz
if [ $? -eq 0 ]; then
  tar -xzf resources/vep_cache/cache.tar.gz -C resources/vep_cache/
  rm resources/vep_cache/cache.tar.gz
fi
