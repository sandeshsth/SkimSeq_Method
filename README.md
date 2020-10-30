## demultiplexNextera

USAGE:

```
demultiplexNextera.pl file.R1.fastq file.R2.fastq i1.fastq i2.fastq barcode.txt
```

Files:
1. Read 1 fastq file: `file.R1.fastq`
2. Read 2 fastq file: `file.R2.fastq` 
3. i7 barcode fastq file: `i1.fastq`
4. i5 barcode fastq file: `i2.fastq` 
5. Barcode file: `barcode.txt`

Barcode file (tab separated without header):

```
Sample1 i1_i2
Sample2 AACCACTC_AAGACTGG
```

Citation:
```
Shrestha S. 2019. Demultiplexing dual-indexed paired-end FASTQ reads. GitHub Repository. https://github.com/sandeshsth/demultiplexNextera.
```
