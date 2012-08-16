exp_ = read.table("/Volumes/agabi/proteomics/MARINA/PTX/raw_datasets/pilot/TOPPAS/TOPPAS_out/023-IDMapper/TMT-mapped.tsv", header=TRUE)

pep_iso = data.frame(accession = exp_$prot.accession, peptide= exp_$prot.name, spectrum=1:length(exp_[,1]), search.engine="?", X126_mass=126.12, X127_mass=127.12, X128_mass=128.13, X129_mass=129.13, X130_mass=130.14, X131_mass=131.13, X126_ions=exp_$i126, X127_ions=exp_$i127, X128_ions=exp_$i128, X129_ions=exp_$i129, X130_ions=exp_$i130, X131_ions=exp_$i131)
ib_pep <- new("TMT6plexSpectra", pep_iso)


library(ggplot2)
maplot(ib_pep,channel1="126",channel2="129",ylim=c(0.5,2),main="Ratio 126:129")
