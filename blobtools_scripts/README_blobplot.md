###-----------Fernando's Instructions--------------
###installation with CONDA

# Create a new environment named 'blobtools_env' and install blobtools into it
conda create -n blobtools_env blobtools

##Activate env

conda activate blobtools_env \ 
blobtools --help

###Git in 
/home/frvsbi/bin/blobtools

####Download NCBI taxdump and create nodesdb
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data \
tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp \
./blobtools nodesdb --nodes data/nodes.dmp --names data/names.dmp

##CREATE 
./blobtools create -i example/assembly.fna -b example/mapping_1.sorted.bam -t example/blast.out -o example/test && \
./blobtools view -i example/test.blobDB.json && \
./blobtools plot -i example/test.blobDB.json 

# -----------Workflow I use for creating blobplots------------------
# parse final.contigs.fa with split_fasta.py then run sbatch blastn.sh 
python split_fasta.py final.contigs.fa output_directory 10000 [10000 is number of sequences per file, feel free to change]

# Run sbatch bwa_samtools.sh on final.contigs.fa (alter script to set paths)

##this will give you three sorted and indexed bam files (one for each replicate); you can specify multiple .bam mapping files in blobtools create, however this then seems to generate three separate graphs - I prefer to merge my .bams so that I only have one graph per site

##!sometimes there is an issue with the arrays conflicting in trying to create/use the final.contigs.fa.pac file at the same time and only the third sorted.bam will be generated; if this is the case I just delete that file, leave the files generated from the indexing of the referece genome, comment out that portion of the code, and run the script again 

# merge the three .sorted.bam files
samtools merge output.file input1.bam input2.bam input3.bam 

##I created a new conda env with samtools in it since I was having issues with the spack version; I usually just run this command in sinteractive and it takes a few minutes

# sort and index merged .bam file with sbatch sort_index_merged_bam.sh (alter script to set paths)

##sorting and indexing the merged bam is too intensive for sinteractive, so use an sbatch script for this

# merge the blastn output files
cat *.out > merged_blastn.out

# create initial blobplot 
blobtools create -i example/assembly.fna -b example/merged.sorted.bam -t example/merged_blast.out -o example/test \
blobtools plot -i example/test.blobDB.json 

##usually no-hits are a significant portion (top ranked taxa)
##I usually do this in sinteractive, but there is an sbatch script that Fernando provided as well

# use view function to get table of contigs and associated taxa
blobtools view -i test.blobDB.json

# extract the no-hit contig names using the table generated from view function
awk -F'\t' '{ if ($6 == "no-hit") print $1 >> "no_hit_names.txt"}' blobDB.table.txt

# find the fasta sequences associated with the no-hit contigs - you may need to install seqtk to do it with this command; I just did a conda install
seqtk subseq /shared/rc/microalgae/metagenomes/assemblies/final.contigs.fa no_hit_names.txt >> no_hit.fasta

# parse no-hit.fasta with split_fasta.py
python split_fasta.py no_hit.fasta output_directory 10000 [10000 is number of sequences per file, feel free to change]

# run sbatch nohit_blastx.sh (alter script to set paths)

##I think you will have to install diamond into your home directory - some versions of diamond are compatible with the ncbi database format, however it seems that the one in the spack env is not and the conda install versions do not have this support yet - I chose to install a version of diamond that was compatible rather than trying to install a new database. Here is the website I referenced: https://github.com/bbuchfink/diamond/issues/439 and and the file downloaded: http://github.com/bbuchfink/diamond/releases/download/v2.0.9/diamond-linux64.tar.gz. You may need to use the prepdb command before running the blastx.

**##you may also want to look into formatting options to see if you can pull taxID number directly in this step - I have not been able to, so there are additional steps to make the blastx output usable for blobplot** 

# merge the blastx output files
cat *.tsv > merged_blastx.tsv

# extract the tax accessions and pull their corresponding tax ID from database (I do in sinteractive)
cut -f2 merged_blastx.tsv > acessions.txt \
spack activate default-genomics-x86_64-25032001 \
blastdbcmd -entry_batch /shared/rc/microalgae/metagenomes/assemblies/accessions.txt -db /shared/rc/datasets/genomics/ncbi-blast-nr_2024-10-01/nr -outfmt "%a %T" > blastx_acc2taxid.tsv

# merge blastx_acc2taxid.tsv with merged_blastx.tsv
awk 'FNR==NR {gsub(/\.[0-9]+$/, "", $1); acc2tax[$1]=$2; next} {
    acc=$2; gsub(/\.[0-9]+$/, "", acc); taxid=acc2tax[acc];
    if (taxid != "") {
        print $1, taxid, $NF
    } else {
        print $1, "no-hit", $NF 
    }
}' blastx_acc2taxid.tsv blastx_merged_site6.tsv > blobtools_blastx_input.tsv

# create blobplot incorporating blastx results (see https://blobtools.readme.io/docs/create for all options)
blobtools create -i example/assembly.fna -b example/merged.sorted.bam -t example/merged_blast.out -t example/blobtools_blastx_input.tsv -x "bestsumorder" -o example/test_with_nohit -n \
blobtools plot -i example/test_with_nohit.blobDB.json -x "bestsumorder"
