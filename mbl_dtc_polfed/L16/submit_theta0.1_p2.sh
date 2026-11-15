#!/bin/bash
# V2.1 L=16, theta=0.1, part 2/3 (10000 jobs); run prepare.sh 0.1 once before
bash "$(dirname "$0")/submit_part.sh" 0.1 2
