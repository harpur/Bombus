###
# Admixture Within and Among Bombus
###





#BIMP

plink --file BIMP --noweb --thin 0.1 --maf 0.1 --make-bed 
for K in 1 2 3 4 5 6 7 8 ; \
do /home/brock/admixture/admixture  --cv=10 plink.bed $K -j7 | tee log${K}.out; done
grep -h CV log*.out
#K=1 CV error = 2.07322

#BTERR
vcftools --vcf BtertoBimp.indel.dp.bl.dd.maxdp.recode.vcf --plink  --out BTERR
plink --file BTERR --noweb --thin 0.1 --maf 0.1 --make-bed
for K in 1 2 3 4 5 6 ; \
do /home/brock/admixture/admixture  --cv=10 plink.bed $K -j7 | tee log${K}.out; done
#K=1 CV error =  1.99645

#BMEL - problem inverting matrix
vcftools --vcf BmeltoBimp.indel.dp.bl.dd.maxdp.recode.vcf --plink  --out BMEL
plink --file BMEL --noweb --thin 0.1 --maf 0.1 --make-bed
for K in 1 2 ; \
do /home/brock/admixture/admixture  --cv=10 plink.bed $K -j7 | tee log${K}.out; done



#All three together
plink --file BMEL --noweb --thin 0.1 --maf 0.1 --recode12 --out BMEL12
plink --file BTERR --noweb --thin 0.1 --maf 0.1 --recode12 --out BTERR12
plink --file BIMP --noweb --thin 0.1 --maf 0.1 --recode12 --out BIMP12



plink --file BMEL12 --make-bed --out BMEL12
plink --file BTERR12 --make-bed --out BTERR12
plink --file BIMP12 --make-bed --out BIMP12

plink --file BIMP --noweb --thin 0.1 --maf 0.1 --make-bed
for K in 1 2 3 4 5 6 7 8 9 10; \
do /home/brock/admixture/admixture  --cv=10 plink.bed $K -j7 | tee log${K}.out; done











