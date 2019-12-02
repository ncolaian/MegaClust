#! /bin/bash
#
#SBATCH --job-name=hmmsearch
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --time=1-0

#hmm dir is first, fasta path second, output third
> $3;
for i in $1/*;
	do hmmsearch $i $2 >> $3;
done;
