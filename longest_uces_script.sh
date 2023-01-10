#obtendo UCEs de maior comprimento

cd ~/Ocoronatus/alignments/internaltrim/mafft/
ls *.fasta | sed 's/.fasta//g' | parallel "python ~/Ocoronatus/fox_scripts/fasta_to_tabular.py {}.fasta {}.tab 0 1"
sed -i 's/[?-]//g' *.tab
ls *.fasta | sed 's/.fasta//g' | parallel "(perl -ne 'print (\$l = \$_) if (length > length(\$l));' {}.tab | tail -1)  > {}_longest_contig.tab"
ls *.fasta | sed 's/.fasta//g' | parallel "sed -i 's/_.*\t/\t/g' {}_longest_contig.tab"
cat *longest_contig* > longest_uces.tab
python ~/Ocoronatus/fox_scripts/tabular_to_fasta.py longest_uces.tab 1 2 longest_uces.fas
