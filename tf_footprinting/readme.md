## tf_footprinting

This directory contains processing and analysis steps for running transcription factor footprinting with the ChromBPNet framework.
Execution broadly follows the order:

 - **chrombpnet-smk-pipeline**: executes the core training and evaluation pipelines for motif instance identification in cell types and HSCs per developmental week
 - **postprocess**: scripts used to aggregate and format some pipeline outputs to deliverable formats
 - **analysis**: analysis scripts (e.g., generating manuscript figures)
 - **tobias-smk-pipeline**: pipeline for generating images of aggregrated bias-corrected coverage plots over putative TF binding sites

