# Protokół Analizy Danych Z DMS-seq
Protokół powstał na podstawie https://github.com/f1lem0n/mgr. 
## Środowisko
Skrypty należy uruchomić dopiero po ustawieniu odpowiedniego środowiska poprzez program ```Anaconda```: 
```bash
conda env create -f environment.yml
conda activate mgr
```
## Pobieranie danych i genomu referencyjnego 
Genom referencyjny pobrano ręcznie z baz danych Narodowego Centrum Biotechnologii Informacyjnej (NCBI), a także odpowiadający mu plik gtf.
Indyfikator NCBI: GCF_000001735.4

Dane pobrano z Archiwum Odczytów Sekwencyjnych (SRA) NCBI i skonwersowano na pliki typu FASTA.

Skryt bash użyty do pobierania:
```bash
#!/bin/bash

mkdir -p sra_data
cd sra_data || exit
for RUN in SRR933551 SRR933552 SRR933556 SRR933557; do
    echo "Pobieram: $RUN"

    prefetch "$RUN"

    if [ -f "$RUN/$RUN.sra" ]; then
        echo "Konwertuję: $RUN"
        
        fastq-dump --split-files --gzip "$RUN/$RUN.sra"
    fi
done
```
## Kontrola Jakości Odczytów
Do kontroli jakości i przycinania sekwencji adapterowych użyto narzędzia ```fastp```:
```bash
mkdir -p output/reads_trimmed
mkdir -p output/QC_premap
for f in $(ls /data_rzodkiewnik/sra_data/*.fastq.gz); do
    fastp -i $f \
        -o output/reads_trimmed/$(basename $f) \
        -h output/QC_premap/$(basename $f).html
done 
```
##Ondeksowanie Genomu Referencyjnego 
Przed mapowaniem należy zbudować indeks genomu referencyjnego. W tymc elu wykorzystano program ```STAR```:
```bash
STAR \
    --runMode genomeGenerate \
    --genomeSAindexNbases 10 \
    --genomeFastaFiles /datarzodkiewnik/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna \
    --sjdbGTFfile /datarzodkiewnik/GCF_000001735.4/genomic.gtf \
    --sjdbOverhang 99 \
    --genomeDir /datarzodkiewnik/GCF_000001735.4/output/STAR_index
```
##Mapowanie Odczytów
Znowu wykorzystano program ```STAR```:
```bash
mkdir -p output/STAR_alignment

for f in $(ls /Users/mariachmielorz/Desktop/Semestr6/data_rzodkiewnik/output/reads_trimmed/*.fastq.gz); do
    /Users/mariachmielorz/Desktop/Semestr6/STAR_2.7.11b/MacOSX_x86_64/STAR \
        --genomeDir output/STAR_index \
        --readFilesCommand "gunzip -c" \
        --readFilesIn $f \
        --runThreadN 9 \
        --outFileNamePrefix /Users/mariachmielorz/Desktop/Semestr6/data_rzodkiewnik/output/STAR_alignment/$(basename $f .fastq.gz)_ \
        --alignEndsType EndToEnd \
        --outFilterMultimapNmax 10 \
        --seedSearchStartLmax 25 \
        --outSAMtype BAM SortedByCoordinate
 done 
```
##Wykrywanie Sygnału DMS
Na początku żyto narzędzia ```rf-count``` z pakietu ```RNAFramework``` aby zliczyć liczbe zatrzymań odwortnej transkryptazy (RT):
```bash
/data_rzodkiewnik/scripts/RNAFramework/rf-count \
    /data_rzodkiewnik/output/STAR_alignment/SRR933551*.bam \
    -f /data_rzodkiewnik/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna \
    -o output/RTS_counts/ \
    --primary-only
```
Potem użyto ```rf-rctools``` do konwersji na format ```tab```:
```bash
/data_rzodkiewnik/scripts/RNAFramework/rf-count \
    /data_rzodkiewnik/output/STAR_alignment/SRR933551*.bam \
    -f /data_rzodkiewnik/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna \
    -o output/RTS_counts/ \
    --primary-only
```
Na koniec oby obliczyć sygnał od DMS i by wygodnie podzielić nasz genom na odpowiednie chromosomy użyto załączonego pliku ```sygnal_dms.py```.



