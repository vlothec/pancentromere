#!/bin/bash
out='/production/rb996/001WGA/001SNPs'
raw='/home/rb996/001ArabidopsisGenus/000Athaliana/000LongReadGenomes/004_Final66/all/66_genomes_datafiles/002WholeGenome'
ref='Col0_Hifi.fasta'
cd $out

################################STEP ONE pairwise alignment to Col0_Hifi. Replace query string with accession of interest##################################
#split by chromosome for memory 

minimap2 -a -x asm5 --cs -r2k -t 16 ${raw}/${ref} ${raw}/query.Chr1to5.fasta > query.sam
samtools sort -@16 -o query.bam query.sam
samtools index query.bam

bcftools mpileup -f ${raw}/${ref} --threads 15 -r Chr1 -o query.Chr1.vcf -A -O v query.bam
bcftools mpileup -f ${raw}/${ref} --threads 15 -r Chr2 -o query.Chr2.vcf -A -O v query.bam
bcftools mpileup -f ${raw}/${ref} --threads 15 -r Chr3 -o query.Chr3.vcf -A -O v query.bam
bcftools mpileup -f ${raw}/${ref} --threads 15 -r Chr4 -o query.Chr4.vcf -A -O v query.bam
bcftools mpileup -f ${raw}/${ref} --threads 15 -r Chr5 -o query.Chr5.vcf -A -O v query.bam

bcftools call -c query.Chr1.vcf -o query.Chr1.called.vcf
bcftools call -c query.Chr2.vcf -o query.Chr2.called.vcf
bcftools call -c query.Chr3.vcf -o query.Chr3.called.vcf
bcftools call -c query.Chr4.vcf -o query.Chr4.called.vcf
bcftools call -c query.Chr5.vcf -o query.Chr5.called.vcf

bgzip -f -@ 16 query.Chr1.called.vcf
bgzip -f -@ 16 query.Chr2.called.vcf
bgzip -f -@ 16 query.Chr3.called.vcf
bgzip -f -@ 16 query.Chr4.called.vcf
bgzip -f -@ 16 query.Chr5.called.vcf


tabix -f -p vcf query.Chr1.called.vcf.gz
tabix -f -p vcf query.Chr2.called.vcf.gz
tabix -f -p vcf query.Chr3.called.vcf.gz
tabix -f -p vcf query.Chr4.called.vcf.gz
tabix -f -p vcf query.Chr5.called.vcf.gz

################################STEP TWO merge to correct population level VCF##################################


bcftools merge -m all --threads 16 -o T2T_Athaliana_66.Chr1.vcf *Chr1.called.vcf.gz
bcftools merge -m all --threads 16 -o T2T_Athaliana_66.Chr2.vcf *Chr2.called.vcf.gz
bcftools merge -m all --threads 16 -o T2T_Athaliana_66.Chr3.vcf *Chr3.called.vcf.gz
bcftools merge -m all --threads 16 -o T2T_Athaliana_66.Chr4.vcf *Chr4.called.vcf.gz
bcftools merge -m all --threads 16 -o T2T_Athaliana_66.Chr5.vcf *Chr5.called.vcf.gz

bgzip -f -@ 16 T2T_Athaliana_66.Chr1.vcf
bgzip -f -@ 16 T2T_Athaliana_66.Chr2.vcf
bgzip -f -@ 16 T2T_Athaliana_66.Chr3.vcf
bgzip -f -@ 16 T2T_Athaliana_66.Chr4.vcf
bgzip -f -@ 16 T2T_Athaliana_66.Chr5.vcf

tabix -f -p vcf T2T_Athaliana_66.Chr1.vcf.gz
tabix -f -p vcf T2T_Athaliana_66.Chr2.vcf.gz
tabix -f -p vcf T2T_Athaliana_66.Chr3.vcf.gz
tabix -f -p vcf T2T_Athaliana_66.Chr4.vcf.gz
tabix -f -p vcf T2T_Athaliana_66.Chr5.vcf.gz

################################STEP THREE Filter for missing sites and heterozygous sites##################################

bcftools view -e 'GT[*]="mis"' T2T_Athaliana_66.Chr1.vcf.gz > T2T_Athaliana_66.Chr1.nomissing.vcf
bcftools view -e 'GT[*]="mis"' T2T_Athaliana_66.Chr2.vcf.gz > T2T_Athaliana_66.Chr2.nomissing.vcf
bcftools view -e 'GT[*]="mis"' T2T_Athaliana_66.Chr3.vcf.gz > T2T_Athaliana_66.Chr3.nomissing.vcf
bcftools view -e 'GT[*]="mis"' T2T_Athaliana_66.Chr4.vcf.gz > T2T_Athaliana_66.Chr4.nomissing.vcf
bcftools view -e 'GT[*]="mis"' T2T_Athaliana_66.Chr5.vcf.gz > T2T_Athaliana_66.Chr5.nomissing.vcf

bcftools view -e 'GT[*]="het"' T2T_Athaliana_66.Chr1.nomissing.vcf > T2T_Athaliana_66.Chr1.nomissing.nohet.vcf
bcftools view -e 'GT[*]="het"' T2T_Athaliana_66.Chr2.nomissing.vcf > T2T_Athaliana_66.Chr2.nomissing.nohet.vcf
bcftools view -e 'GT[*]="het"' T2T_Athaliana_66.Chr3.nomissing.vcf > T2T_Athaliana_66.Chr3.nomissing.nohet.vcf
bcftools view -e 'GT[*]="het"' T2T_Athaliana_66.Chr4.nomissing.vcf > T2T_Athaliana_66.Chr4.nomissing.nohet.vcf
bcftools view -e 'GT[*]="het"' T2T_Athaliana_66.Chr5.nomissing.vcf > T2T_Athaliana_66.Chr5.nomissing.nohet.vcf

################################STEP FOUR Biallelic SNPs and Non Variant sites for Pi##################################

bcftools view -M2 -m2 -v snps T2T_Athaliana_66.Chr1.nomissing.nohet.vcf > T2T_Athaliana_66.Chr1.nomissing.nohet.justSNPs.vcf
bcftools view -M2 -m2 -v snps T2T_Athaliana_66.Chr2.nomissing.nohet.vcf > T2T_Athaliana_66.Chr2.nomissing.nohet.justSNPs.vcf
bcftools view -M2 -m2 -v snps T2T_Athaliana_66.Chr3.nomissing.nohet.vcf > T2T_Athaliana_66.Chr3.nomissing.nohet.justSNPs.vcf
bcftools view -M2 -m2 -v snps T2T_Athaliana_66.Chr4.nomissing.nohet.vcf > T2T_Athaliana_66.Chr4.nomissing.nohet.justSNPs.vcf
bcftools view -M2 -m2 -v snps T2T_Athaliana_66.Chr5.nomissing.nohet.vcf > T2T_Athaliana_66.Chr5.nomissing.nohet.justSNPs.vcf

bcftools view -C0 T2T_Athaliana_66.Chr1.nomissing.nohet.vcf > T2T_Athaliana_66.Chr1.nomissing.nohet.justNV.vcf
bcftools view -C0 T2T_Athaliana_66.Chr2.nomissing.nohet.vcf > T2T_Athaliana_66.Chr2.nomissing.nohet.justNV.vcf
bcftools view -C0 T2T_Athaliana_66.Chr3.nomissing.nohet.vcf > T2T_Athaliana_66.Chr3.nomissing.nohet.justNV.vcf
bcftools view -C0 T2T_Athaliana_66.Chr4.nomissing.nohet.vcf > T2T_Athaliana_66.Chr4.nomissing.nohet.justNV.vcf
bcftools view -C0 T2T_Athaliana_66.Chr5.nomissing.nohet.vcf > T2T_Athaliana_66.Chr5.nomissing.nohet.justNV.vcf

bgzip -@ 16 T2T_Athaliana_66.Chr1.nomissing.nohet.justSNPs.vcf
tabix -p vcf T2T_Athaliana_66.Chr1.nomissing.nohet.justSNPs.vcf.gz

bgzip -@ 16 T2T_Athaliana_66.Chr2.nomissing.nohet.justSNPs.vcf
tabix -p vcf T2T_Athaliana_66.Chr2.nomissing.nohet.justSNPs.vcf.gz

bgzip -@ 16 T2T_Athaliana_66.Chr3.nomissing.nohet.justSNPs.vcf
tabix -p vcf T2T_Athaliana_66.Chr3.nomissing.nohet.justSNPs.vcf.gz

bgzip -@ 16 T2T_Athaliana_66.Chr4.nomissing.nohet.justSNPs.vcf
tabix -p vcf T2T_Athaliana_66.Chr4.nomissing.nohet.justSNPs.vcf.gz

bgzip -@ 16 T2T_Athaliana_66.Chr5.nomissing.nohet.justSNPs.vcf
tabix -p vcf T2T_Athaliana_66.Chr5.nomissing.nohet.justSNPs.vcf.gz



bcftools concat -a --threads 6 -o T2T_Athaliana_66.Chr1.nomissing.nohet.justSNPsNVs.vcf T2T_Athaliana_66.Chr1.nomissing.nohet.just*.vcf.gz
bcftools concat -a --threads 6 -o T2T_Athaliana_66.Chr2.nomissing.nohet.justSNPsNVs.vcf T2T_Athaliana_66.Chr2.nomissing.nohet.just*.vcf.gz
bcftools concat -a --threads 6 -o T2T_Athaliana_66.Chr3.nomissing.nohet.justSNPsNVs.vcf T2T_Athaliana_66.Chr3.nomissing.nohet.just*.vcf.gz
bcftools concat -a --threads 6 -o T2T_Athaliana_66.Chr4.nomissing.nohet.justSNPsNVs.vcf T2T_Athaliana_66.Chr4.nomissing.nohet.just*.vcf.gz
bcftools concat -a --threads 6 -o T2T_Athaliana_66.Chr5.nomissing.nohet.justSNPsNVs.vcf T2T_Athaliana_66.Chr5.nomissing.nohet.just*.vcf.gz


bcftools concat --threads 16 -o T2T_Athaliana_66.Chr1to5.nomissing.nohet.justSNPs.vcf T2T_Athaliana_66.Chr*.nomissing.nohet.justSNPs.vcf

bgzip -@ 16 T2T_Athaliana_66.Chr1.nomissing.nohet.justSNPsNVs.vcf
tabix -p vcf T2T_Athaliana_66.Chr1.nomissing.nohet.justSNPsNVs.vcf.gz

bgzip -@ 16 T2T_Athaliana_66.Chr2.nomissing.nohet.justSNPsNVs.vcf
tabix -p vcf T2T_Athaliana_66.Chr2.nomissing.nohet.justSNPsNVs.vcf.gz

bgzip -@ 16 T2T_Athaliana_66.Chr3.nomissing.nohet.justSNPsNVs.vcf
tabix -p vcf T2T_Athaliana_66.Chr3.nomissing.nohet.justSNPsNVs.vcf.gz

bgzip -@ 16 T2T_Athaliana_66.Chr4.nomissing.nohet.justSNPsNVs.vcf
tabix -p vcf T2T_Athaliana_66.Chr4.nomissing.nohet.justSNPsNVs.vcf.gz

bgzip -@ 16 T2T_Athaliana_66.Chr5.nomissing.nohet.justSNPsNVs.vcf
tabix -p vcf T2T_Athaliana_66.Chr5.nomissing.nohet.justSNPsNVs.vcf.gz

################################STEP FIVE Remove Repeat sites##################################

bcftools view -R Col_0.Chr1to5.TAIR10TECentrDNAmasked.ATGCs.1based.bed -o T2T_Athaliana_66.Chr2.nomissing.nohet.justSNPsNVs.norepeats.vcf T2T_Athaliana_66.Chr2.nomissing.nohet.justSNPsNVs.vcf.gz 
bcftools view -R Col_0.Chr1to5.TAIR10TECentrDNAmasked.ATGCs.1based.bed -o T2T_Athaliana_66.Chr3.nomissing.nohet.justSNPsNVs.norepeats.vcf T2T_Athaliana_66.Chr3.nomissing.nohet.justSNPsNVs.vcf.gz
bcftools view -R Col_0.Chr1to5.TAIR10TECentrDNAmasked.ATGCs.1based.bed -o T2T_Athaliana_66.Chr4.nomissing.nohet.justSNPsNVs.norepeats.vcf T2T_Athaliana_66.Chr4.nomissing.nohet.justSNPsNVs.vcf.gz
bcftools view -R Col_0.Chr1to5.TAIR10TECentrDNAmasked.ATGCs.1based.bed -o T2T_Athaliana_66.Chr5.nomissing.nohet.justSNPsNVs.norepeats.vcf T2T_Athaliana_66.Chr5.nomissing.nohet.justSNPsNVs.vcf.gz

bgzip -@ 16 T2T_Athaliana_66.Chr1.nomissing.nohet.justSNPs.norepeats.vcf
tabix -p vcf T2T_Athaliana_66.Chr1.nomissing.nohet.justSNPs.norepeats.vcf.gz
bgzip -@ 16 T2T_Athaliana_66.Chr2.nomissing.nohet.justSNPs.norepeats.vcf
tabix -p vcf T2T_Athaliana_66.Chr2.nomissing.nohet.justSNPs.norepeats.vcf.gz
bgzip -@ 16 T2T_Athaliana_66.Chr3.nomissing.nohet.justSNPs.norepeats.vcf
tabix -p vcf T2T_Athaliana_66.Chr3.nomissing.nohet.justSNPs.norepeats.vcf.gz
bgzip -@ 16 T2T_Athaliana_66.Chr4.nomissing.nohet.justSNPs.norepeats.vcf
tabix -p vcf T2T_Athaliana_66.Chr4.nomissing.nohet.justSNPs.norepeats.vcf.gz
bgzip -@ 16 T2T_Athaliana_66.Chr5.nomissing.nohet.justSNPs.norepeats.vcf
tabix -p vcf T2T_Athaliana_66.Chr5.nomissing.nohet.justSNPs.norepeats.vcf.gz
bcftools concat -o T2T_Athaliana_66.Chr1to5.nomissing.nohet.justSNPs.norepeats.vcf T2T_Athaliana_66.Chr*.nomissing.nohet.justSNPs.norepeats.vcf.gz

bgzip -@ 16 T2T_Athaliana_66.Chr1.nomissing.nohet.justNV.vcf
bgzip -@ 16 T2T_Athaliana_66.Chr2.nomissing.nohet.justNV.vcf
bgzip -@ 16 T2T_Athaliana_66.Chr3.nomissing.nohet.justNV.vcf
bgzip -@ 16 T2T_Athaliana_66.Chr4.nomissing.nohet.justNV.vcf
bgzip -@ 16 T2T_Athaliana_66.Chr5.nomissing.nohet.justNV.vcf
tabix -p vcf T2T_Athaliana_66.Chr1.nomissing.nohet.justNV.vcf.gz
tabix -p vcf T2T_Athaliana_66.Chr2.nomissing.nohet.justNV.vcf.gz
tabix -p vcf T2T_Athaliana_66.Chr3.nomissing.nohet.justNV.vcf.gz
tabix -p vcf T2T_Athaliana_66.Chr4.nomissing.nohet.justNV.vcf.gz
tabix -p vcf T2T_Athaliana_66.Chr5.nomissing.nohet.justNV.vcf.gz

bcftools view -R Col_0.Chr1to5.TAIR10TECentrDNAmasked.ATGCs.1based.bed -o T2T_Athaliana_66.Chr1.nomissing.nohet.justNV.norepeats.vcf T2T_Athaliana_66.Chr1.nomissing.nohet.justNV.vcf.gz
bcftools view -R Col_0.Chr1to5.TAIR10TECentrDNAmasked.ATGCs.1based.bed -o T2T_Athaliana_66.Chr2.nomissing.nohet.justNV.norepeats.vcf T2T_Athaliana_66.Chr2.nomissing.nohet.justNV.vcf.gz
bcftools view -R Col_0.Chr1to5.TAIR10TECentrDNAmasked.ATGCs.1based.bed -o T2T_Athaliana_66.Chr3.nomissing.nohet.justNV.norepeats.vcf T2T_Athaliana_66.Chr3.nomissing.nohet.justNV.vcf.gz
bcftools view -R Col_0.Chr1to5.TAIR10TECentrDNAmasked.ATGCs.1based.bed -o T2T_Athaliana_66.Chr4.nomissing.nohet.justNV.norepeats.vcf T2T_Athaliana_66.Chr4.nomissing.nohet.justNV.vcf.gz
bcftools view -R Col_0.Chr1to5.TAIR10TECentrDNAmasked.ATGCs.1based.bed -o T2T_Athaliana_66.Chr5.nomissing.nohet.justNV.norepeats.vcf T2T_Athaliana_66.Chr5.nomissing.nohet.justNV.vcf.gz

bgzip -@ 16 T2T_Athaliana_66.Chr1.nomissing.nohet.justNV.norepeats.vcf
bgzip -@ 16 T2T_Athaliana_66.Chr2.nomissing.nohet.justNV.norepeats.vcf
bgzip -@ 16 T2T_Athaliana_66.Chr3.nomissing.nohet.justNV.norepeats.vcf
bgzip -@ 16 T2T_Athaliana_66.Chr4.nomissing.nohet.justNV.norepeats.vcf
bgzip -@ 16 T2T_Athaliana_66.Chr5.nomissing.nohet.justNV.norepeats.vcf
tabix -p vcf T2T_Athaliana_66.Chr1.nomissing.nohet.justNV.norepeats.vcf.gz
tabix -p vcf T2T_Athaliana_66.Chr2.nomissing.nohet.justNV.norepeats.vcf.gz
tabix -p vcf T2T_Athaliana_66.Chr3.nomissing.nohet.justNV.norepeats.vcf.gz
tabix -p vcf T2T_Athaliana_66.Chr4.nomissing.nohet.justNV.norepeats.vcf.gz
tabix -p vcf T2T_Athaliana_66.Chr5.nomissing.nohet.justNV.norepeats.vcf.gz


bcftools concat -a --threads 6 -o T2T_Athaliana_66.Chr1.nomissing.nohet.norepeats.justSNPsNVs.vcf T2T_Athaliana_66.Chr1.nomissing.nohet.justNV.norepeats.vcf.gz T2T_Athaliana_66.Chr1.nomissing.nohet.justSNPs.norepeats.vcf.gz

bcftools concat -a --threads 6 -o T2T_Athaliana_66.Chr2.nomissing.nohet.norepeats.justSNPsNVs.vcf T2T_Athaliana_66.Chr2.nomissing.nohet.justNV.norepeats.vcf.gz T2T_Athaliana_66.Chr2.nomissing.nohet.justSNPs.norepeats.vcf.gz

bcftools concat -a --threads 6 -o T2T_Athaliana_66.Chr3.nomissing.nohet.norepeats.justSNPsNVs.vcf T2T_Athaliana_66.Chr3.nomissing.nohet.justNV.norepeats.vcf.gz T2T_Athaliana_66.Chr3.nomissing.nohet.justSNPs.norepeats.vcf.gz

bcftools concat -a --threads 6 -o T2T_Athaliana_66.Chr4.nomissing.nohet.norepeats.justSNPsNVs.vcf T2T_Athaliana_66.Chr4.nomissing.nohet.justNV.norepeats.vcf.gz T2T_Athaliana_66.Chr4.nomissing.nohet.justSNPs.norepeats.vcf.gz

bcftools concat -a --threads 6 -o T2T_Athaliana_66.Chr5.nomissing.nohet.norepeats.justSNPsNVs.vcf T2T_Athaliana_66.Chr5.nomissing.nohet.justNV.norepeats.vcf.gz T2T_Athaliana_66.Chr5.nomissing.nohet.justSNPs.norepeats.vcf.gz



##DONE!!!