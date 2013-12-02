/srv/gs1/software/R/R-2.15.1/bin/R --no-save --no-restore --silent --args 2 3 10 11 CDS control over < ~/EMI_RIBOSEQ/RiboSeq_Analysis_Pipeline.R > OUTFILE_c_o_rib
/srv/gs1/software/R/R-2.15.1/bin/R --no-save --no-restore --silent --args 2 3 12 13 CDS control siRNA < ~/EMI_RIBOSEQ/RiboSeq_Analysis_Pipeline.R > OUTFILE_c_s_rib
/srv/gs1/software/R/R-2.15.1/bin/R --no-save --no-restore --silent --args 10 11 12 13 CDS over siRNA < ~/EMI_RIBOSEQ/RiboSeq_Analysis_Pipeline.R > OUTFILE_o_s_rib
/srv/gs1/software/R/R-2.15.1/bin/R --no-save --no-restore --silent --args 4 5 6 7  EX control over < ~/EMI_RIBOSEQ/RiboSeq_Analysis_Pipeline.R > OUTFILE_c_o_rna
/srv/gs1/software/R/R-2.15.1/bin/R --no-save --no-restore --silent --args 4 5 8 9 EX control siRNA < ~/EMI_RIBOSEQ/RiboSeq_Analysis_Pipeline.R > OUTFILE_rna
/srv/gs1/software/R/R-2.15.1/bin/R --no-save --no-restore --silent --args 6 7 8 9 EX over  siRNA < ~/EMI_RIBOSEQ/RiboSeq_Analysis_Pipeline.R > OUTFILE_rna