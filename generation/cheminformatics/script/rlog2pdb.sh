#!/usr/bin/env bash 

set -ex 

dir1="$1"
dir2="$2" 
mkdir -p "$dir2"
dir2="$(cd "$dir2" && pwd)"

cd "$dir1" || exit 1
for file in *.log; do
  [ -f "$file" ] || continue 
  fname="${file%.log}" 
  obabel -i log "$file" -O "${fname}.pdb" 
  mv -n "${fname}.pdb" "$dir2"/
done 
