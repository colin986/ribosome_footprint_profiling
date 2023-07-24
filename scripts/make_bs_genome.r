

try(system("mkdir -p reference_genome/BSgenome/sequence/", intern = TRUE))
try(system("cp reference_genome/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna reference_genome/BSgenome/sequence/genome.fna", intern = TRUE))
try(system("../bin/kentUtils/bin/linux.x86_64/faToTwoBit reference_genome/BSgenome/sequence/genome.fna reference_genome/BSgenome/sequence/genome.2bit", 
intern = TRUE))

seed_info <- c("Package: BSgenome.cgr.ncbi",
             "Title: Limited genome sequence for Chinese hamser (NCBI version PICRH)",
             "Description: Limited genome sequence for Chinese hamster",
             "Version: 1.0.0",
             "genome: BSgenome.PICRH",
             "provider: NCBI",
             "release_date: sept 2019",
             "organism_biocview: cgr",
             "BSgenomeObjname: NCBI_cgr",
             "organism: Cricetulus griseus",
             "common_name: Chinese hamster",
             "seqs_srcdir: reference_genome/BSgenome/sequence",
             "seqfile_name: genome.2bit",
             "circ_seqs: \"MT\"")
cat(seed_info, file = "reference_genome/BSgenome/picrh.seed", sep = "\n")

forgeBSgenomeDataPkg("reference_genome/BSgenome/picrh.seed",destdir = "reference_genome/BSgenome")


try(system("cd reference_genome/BSgenome && R CMD build BSgenome.cgr.ncbi/ && cd ..", intern = FALSE))

try(system("R CMD check reference_genome/BSgenome/BSgenome.cgr.ncbi_1.0.0.tar.gz", intern = FALSE))
try(system("R CMD INSTALL reference_genome/BSgenome/BSgenome.cgr.ncbi_1.0.0.tar.gz", intern = FALSE))  