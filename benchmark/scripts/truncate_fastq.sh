#!/bin/bash
# Truncate every read to BM_READLEN nt. in.fq.gz -> out.fq.gz
set -euo pipefail
source "$(dirname "$0")/settings.sh"
IN=$1; OUT=$2
pigz -dc -p 4 "$IN" \
| mawk -v L="$BM_READLEN" 'NR%4==2||NR%4==0 { $0 = substr($0, 1, L) } { print }' \
| pigz -c -p 6 > "$OUT"
