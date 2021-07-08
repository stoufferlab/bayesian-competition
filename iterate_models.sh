#!/bin/bash
#$ -j y
#$ -cwd

PATH=/share/apps/R/R-3.6.2/bin/:/share/apps/gcc-8.3/bin:$PATH
LD_LIBRARY_PATH=/share/apps/gcc-8.3/lib64:$LD_LIBRARY_PATH




Rscript code/iterate_models_cluster.R $SGE_TASK_ID




