#!/bin/bash
# puts the one julia script of all jobs into $BUDDY (it takes L theta r job eps nev set as arguments)
POL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cp "$POL_DIR/polfed_bands.jl" "$BUDDY/polfed_bands.jl"
echo "cat_jl: polfed_bands.jl -> \$BUDDY"
