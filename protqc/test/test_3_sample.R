devtools::load_all()

library(protqc)

setwd("~/mark/FD/QDP本地化/Package_renew/Quartet-Protein-QC/protqc/test")

## example data
example_prot_data_path <- './proteomics_pipeline_data_template.csv'
example_prot_metadata_path <- './proteomics_pipeline_meta_template.csv'

example_prot_data_path <- system.file("extdata","proteomics_pipeline_data_template.csv",package = "protqc")
example_prot_metadata_path <- system.file("extdata","proteomics_pipeline_meta_template.csv",package = "protqc")


output_dir <- './'
prot_result <- protqc::qc_conclusion(example_prot_data_path, example_prot_metadata_path, output_dir, plot=TRUE)

prot_result

doc_file_path_example <- system.file("extdata", "Quartet_temp.docx", package = "protqc")

GenerateProteinReport(Prot_result=prot_result, doc_file_path=doc_file_path_example, output_path=output_dir)
