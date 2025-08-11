#!/bin/bash

ID=$1
BASE_DIR=$2

echo "Running $ID at $time" 

python "${BASE_DIR}/code/py/filer_alignment.py" "${BASE_DIR}/data/BUSCO/fasta_BUSCO_concat_lowmiss/mq.call.${ID}.fa" > "${BASE_DIR}/data/BUSCO/phylo_BUSCO_lowmiss/mq.call.${ID}.seq.ex.fa" 

clipkit "${BASE_DIR}/data/BUSCO/phylo_BUSCO_lowmiss/mq.call.${ID}.seq.ex.fa"  -s nt -gc "\-*xXN" -m gappy -g 0.5 -of phylip_relaxed -o "${BASE_DIR}/data/BUSCO/phylo_BUSCO_lowmiss/mq.call.${ID}.seq.ex.phy"

/data/genomics/apps/iqtree-2.2.2.7-Linux/bin/iqtree2 -redo -s "${BASE_DIR}/data/BUSCO/phylo_BUSCO_lowmiss/mq.call.${ID}.seq.ex.phy" -B 1000 -alrt 1000 --boot-trees -pre "${BASE_DIR}/data/BUSCO/phylo_BUSCO_lowmiss/iqtree/mq.call.${ID}.seq.ex.tree"
