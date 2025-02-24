#!/usr/bin/env sh

OUTPUT_PATH="mmseqs2_output"
LOG_FILE="${OUTPUT_PATH}/mmseqs2.log"
SENSITIVITY=7.5
THREADS=64

# Function to run a command and write the output to a log file
run_command() {
    local RED="\033[31m"
    local BLUE="\033[34m"
    local DIM="\033[2m"
    local RESET="\033[0m"
    local timestamp="[$(date '+%Y-%m-%d %H:%M:%S')]"
    local formatted_msg="${DIM}${timestamp}${RESET} ${BLUE}Running:${RESET} $*"
    echo "${formatted_msg}"
    "$@" >> "${LOG_FILE}" 2>&1
    local exit_code=$?
    if [ $exit_code -ne 0 ]; then
        local error_timestamp="[$(date '+%Y-%m-%d %H:%M:%S')]"
        echo "${DIM}${error_timestamp}${RESET} ${RED}Error:${RESET} Command failed with exit code ${exit_code}${RESET}"
        echo "${error_timestamp} Error: Command failed with exit code ${exit_code}. Please find the execution log in ${LOG_FILE}." >> "${LOG_FILE}"
        exit $exit_code
    fi
}

# Check if the output directory already exists
if [ -d "${OUTPUT_PATH}" ]; then
    echo "\033[31mError: ${OUTPUT_PATH} directory already exists.\033[0m"
    echo "\033[2mPlease remove it manually before running this script: rm -rf ${OUTPUT_PATH}\033[0m"
    exit 1
fi
mkdir "${OUTPUT_PATH}"

# Create a database of the query sequences and extract the ORFs
run_command mmseqs createdb dataset_challenge.fasta ${OUTPUT_PATH}/query --dbtype 0 --shuffle 1 --createdb-mode 1 --write-lookup 0 --id-offset 0 --compressed 0 -v 3
run_command mmseqs extractorfs ${OUTPUT_PATH}/query ${OUTPUT_PATH}/orfs --min-length 30 --max-length 32734 --max-gaps 2147483647 --contig-start-mode 2 --contig-end-mode 2 --orf-start-mode 1 --forward-frames 1,2,3 --reverse-frames 1,2,3 --translation-table 1 --translate 1 --use-all-table-starts 0 --id-offset 0 --create-lookup 0 --threads ${THREADS} --compressed 0 -v 3

# Perform a quick alignment of the ORFs to the target database and thiscard the ORFs that are either not the longest one in the contig or that had no hit
run_command mmseqs prefilter ${OUTPUT_PATH}/orfs ictv_nr_db/ictv_nr_db ${OUTPUT_PATH}/orfs_pref --sub-mat 'aa:blosum62.out,nucl:nucleotide.out' --seed-sub-mat 'aa:VTML80.out,nucl:nucleotide.out' -s 2 -k 0 --target-search-mode 0 --k-score seq:2147483647,prof:2147483647 --alph-size aa:21,nucl:5 --max-seq-len 65535 --max-seqs 1 --split 0 --split-mode 2 --split-memory-limit 0 -c 0 --cov-mode 0 --comp-bias-corr 1 --comp-bias-corr-scale 1 --diag-score 0 --exact-kmer-matching 0 --mask 1 --mask-prob 0.9 --mask-lower-case 0 --mask-n-repeat 0 --min-ungapped-score 3 --add-self-matches 0 --spaced-kmer-mode 1 --db-load-mode 0 --threads ${THREADS} --compressed 0 -v 3
run_command mmseqs rescorediagonal ${OUTPUT_PATH}/orfs ictv_nr_db/ictv_nr_db ${OUTPUT_PATH}/orfs_pref ${OUTPUT_PATH}/orfs_aln --sub-mat 'aa:blosum62.out,nucl:nucleotide.out' --rescore-mode 2 --wrapped-scoring 0 --filter-hits 0 -e 100 -c 0 -a 0 --cov-mode 0 --min-seq-id 0 --min-aln-len 0 --seq-id-mode 0 --add-self-matches 0 --sort-results 0 --db-load-mode 0 --threads ${THREADS} --compressed 0 -v 3
run_command mmseqs recoverlongestorf ${OUTPUT_PATH}/orfs ${OUTPUT_PATH}/orfs_aln ${OUTPUT_PATH}/orfs_aln_recovered.list --threads ${THREADS} -v 3
awk '$3 > 1 { print $1 }' ${OUTPUT_PATH}/orfs_aln.index > ${OUTPUT_PATH}/orfs_aln_remain.list
cat ${OUTPUT_PATH}/orfs_aln_recovered.list ${OUTPUT_PATH}/orfs_aln_remain.list > ${OUTPUT_PATH}/orfs_aln.list
run_command mmseqs createsubdb ${OUTPUT_PATH}/orfs_aln.list ${OUTPUT_PATH}/orfs ${OUTPUT_PATH}/orfs_filter --subdb-mode 1 --id-mode 0 -v 3
run_command mmseqs rmdb ${OUTPUT_PATH}/orfs_filter_h -v 3
run_command mmseqs createsubdb ${OUTPUT_PATH}/orfs_aln.list ${OUTPUT_PATH}/orfs_h ${OUTPUT_PATH}/orfs_filter_h --subdb-mode 1 --id-mode 0 -v 3

# Perform the actual alignment between the selected ORFs and the target database, then assign a taxonomy to each ORF using the 2bLCA approach
run_command mmseqs prefilter ${OUTPUT_PATH}/orfs_filter ictv_nr_db/ictv_nr_db ${OUTPUT_PATH}/pref_0 --sub-mat 'aa:blosum62.out,nucl:nucleotide.out' --seed-sub-mat 'aa:VTML80.out,nucl:nucleotide.out' -k 0 --target-search-mode 0 --k-score seq:2147483647,prof:2147483647 --alph-size aa:21,nucl:5 --max-seq-len 65535 --max-seqs 300 --split 0 --split-mode 2 --split-memory-limit 0 -c 0 --cov-mode 0 --comp-bias-corr 1 --comp-bias-corr-scale 1 --diag-score 1 --exact-kmer-matching 0 --mask 1 --mask-prob 0.9 --mask-lower-case 0 --mask-n-repeat 0 --min-ungapped-score 15 --add-self-matches 0 --spaced-kmer-mode 1 --db-load-mode 0 --threads ${THREADS} --compressed 0 -v 3 -s ${SENSITIVITY}
run_command mmseqs lcaalign ${OUTPUT_PATH}/orfs_filter ictv_nr_db/ictv_nr_db ${OUTPUT_PATH}/pref_0 ${OUTPUT_PATH}/orfs_tax_aln --sub-mat 'aa:blosum62.out,nucl:nucleotide.out' -a 0 --alignment-mode 1 --alignment-output-mode 0 --wrapped-scoring 0 -e 1 --min-seq-id 0 --min-aln-len 0 --seq-id-mode 0 --alt-ali 0 -c 0 --cov-mode 0 --max-seq-len 65535 --comp-bias-corr 1 --comp-bias-corr-scale 1 --max-rejected 5 --max-accept 30 --add-self-matches 0 --db-load-mode 0 --score-bias 0 --realign 0 --realign-score-bias -0.2 --realign-max-seqs 2147483647 --corr-score-weight 0 --gap-open aa:11,nucl:5 --gap-extend aa:1,nucl:2 --zdrop 40 --threads ${THREADS} --compressed 0 -v 3
run_command mmseqs lca ictv_nr_db/ictv_nr_db ${OUTPUT_PATH}/orfs_tax_aln ${OUTPUT_PATH}/orfs_tax --blacklist '' --tax-lineage 0 --compressed 0 --threads ${THREADS} -v 3
run_command mmseqs filterdb ${OUTPUT_PATH}/orfs_tax_aln ${OUTPUT_PATH}/orfs_tax_aln_first --extract-lines 1 --threads ${THREADS} --compressed 0 -v 3

# Generate TSV files:
# 1. `orfs_taxonomy.tsv` contains the taxonomy assigned to each ORF
# 2. `orfs_taxonomy_aln.tsv` contains the alignment information of the best hit of each ORF
run_command mmseqs createtsv ${OUTPUT_PATH}/orfs_filter ${OUTPUT_PATH}/orfs_tax ${OUTPUT_PATH}/orfs_taxonomy.tsv --full-header
perl -pe 's/\t(?=(.*\t){3})/_/g; s/"//g' ${OUTPUT_PATH}/orfs_taxonomy.tsv | sponge ${OUTPUT_PATH}/orfs_taxonomy.tsv
run_command mmseqs convertalis ${OUTPUT_PATH}/orfs_filter ictv_nr_db/ictv_nr_db ${OUTPUT_PATH}/orfs_tax_aln_first ${OUTPUT_PATH}/orfs_taxonomy_aln.tsv --format-output qheader,target,fident,bits,evalue
perl -pe 's/\t(?=(.*\t){4})/_/g' ${OUTPUT_PATH}/orfs_taxonomy_aln.tsv | sponge ${OUTPUT_PATH}/orfs_taxonomy_aln.tsv
