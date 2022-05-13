# Summary of fixes used for the NGSBioAfrica2022 VM so far

#### To fix lumpyexpress error:

```
conda create --name thelumpyenv python=2.7

conda activate thelumpyenv

conda install pysam

conda install numpy

cd /home/manager/course_data/structural_variation/data/exercise3/

samtools view -bh -F 1294 ERR1015069.bam | samtools sort -O bam -T ERR1015069.temp -o ERR1015069.discordants.bam

samtools index ERR1015069.discordants.bam 

samtools view -h ERR1015069.bam | extractSplitReads_BwaMem -i stdin | samtools view -b - | samtools sort -O bam -T ERR1015069.temp -o ERR1015069.splitters.bam

samtools index ERR1015069.splitters.bam

lumpyexpress -B ERR1015069.bam -S ERR1015069.splitters.bam -D ERR1015069.discordants.bam -o ERR1015069.vcf

ls
# deactivate the environment when done
conda deactivate

```



#### To install snpdists

```
conda install snp-dists
```

#### To install figtree

```
sudo apt install figtree
# the password is manager
```

#### To install sleuth for module 7:
```
# Begin by running R

R

BiocManager::install()

install.packages(“devtools”)

BiocManager::install(“pachterlab/sleuth”)

# continue with the practical from here
```
#### To install salmon for mod 8:
```
conda create -n salmon salmon
conda activate salmon

# run your salmon commands

# deactivate the environment when done

conda deactivate
```
