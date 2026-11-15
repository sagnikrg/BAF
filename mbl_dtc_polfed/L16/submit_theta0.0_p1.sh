#!/bin/bash
# V2.1 L=16, theta=0.0, part 1/3 (10000 jobs); run prepare.sh 0.0 once before
bash "$(dirname "$0")/submit_part.sh" 0.0 1
