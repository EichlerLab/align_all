# align_all

A general Snakemake-based read alignment pipeline for long-read and short-read sequencing data, including PacBio HiFi, ONT, ONT ultra-long, ONT methylation-aware reads, ONT BAM input, and Illumina paired-end reads.

This workflow provides a shared alignment framework across multiple sequencing technologies. Depending on input type, it supports direct mapping, BAM-to-FASTQ conversion for selected ONT BAM inputs, per-file or split-batch alignment, merging, duplicate marking, and optional CRAM conversion.

## Overview

At a high level, the pipeline:

- reads sample/input definitions from a manifest file
- reads one or more references from a YAML config
- chooses an aligner preset based on data type
- aligns reads per input file or split batch
- sorts and merges alignments into sample-level BAM files
- optionally converts final BAMs to CRAM after marking duplicates

## Repository contents

```text
align_all/
├── align_all.smk
├── map_contigs.yaml
├── aln.tab
├── runsnake
├── runlocal
└── smk_scripts/
```

## Requirements

This workflow assumes:

- Snakemake is available
- a DRMAA-enabled cluster environment is available for cluster submission
- Singularity/Apptainer-compatible execution is available for containerized runs
- standard alignment tools are available through the configured container image

The workflow uses the container image:

```text
docker://eichlerlab/align-basics:0.3
```

## Configuration

The main workflow uses a YAML config file such as `map_contigs.yaml`.

Example:

```yaml
MANIFEST: aln.tab
REF:
  GRCh38: /net/eichler/vol28/eee_shared/assemblies/hg38/no_alt/hg38.no_alt.fa
NBATCHES: 15
```

### Config fields

#### `MANIFEST`
Path to the manifest file.

#### `REF`
A YAML dictionary where keys are reference names and values are paths to reference FASTA files.

#### `NBATCHES` *(optional)*
Number of batches used for split mapping. If not specified, the workflow defaults to `15`.

#### `ALN_PARAMS` *(optional)*
Additional alignment parameters appended to the default command for each technology.

## Manifest format

The manifest is a three-column tab-separated file:

```text
SAMPLE  FOFN  TYPE
```

Example:

```text
HG002_8X  fofn/HG002_PacBio_HiFi_8X.fofn  PacBio_HiFi
HG002_8X  fofn/HG002_ONT_8X.fofn          ONT
HG002_8X  fofn/HG002_ONT_UL_8X.fofn       ONT_UL
```

### Manifest columns

#### `SAMPLE`
Sample name.

#### `FOFN`
Path to a file-of-file-names.

#### `TYPE`
Input data type. Supported values in the current workflow include:

- `CLR`
- `CCS`
- `PacBio_HiFi`
- `PacBio_HiFi_BAM`
- `REVIO`
- `ONT`
- `ONT_UL`
- `Illumina`
- `HiFi_minimap`
- `ONT_methyl`
- `ONT_BAM`
- `ONT_Q20`

## FOFN format

FOFN stands for “file of file names”.

### Long-read FOFN

For ONT and PacBio-style inputs, the FOFN should contain one input file per line.

Example:

```text
/net/eichler/vol28/projects/long_read_archive/nobackups/nhp/Asia_NLE/raw_data/PacBio_HiFi/m54329U_200812_234915.Q20.fastq
/net/eichler/vol28/projects/long_read_archive/nobackups/nhp/Asia_NLE/raw_data/PacBio_HiFi/m54329U_200814_231018.Q20.fastq
```

For `ONT_BAM` and `PacBio_HiFi_BAM`, the input file listed in the FOFN must end with `.bam`. This mode is intended for **unaligned BAM** input and is used when sequence data and relevant read tags should be carried forward into the downstream mapping process.

### Illumina FOFN

For Illumina paired-end reads, each line should contain read 1 and read 2 separated by a space.

Example:

```text
/net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/fastq/ILLUMINA_NA19704_1.fastq.gz /net/eichler/vol26/projects/1kg_high_coverage/nobackups/data/Illumina/fastq/ILLUMINA_NA19704_2.fastq.gz
```

## Alignment backends

The workflow chooses mapping commands by input type.

### PacBio / HiFi

- `CLR` → `pbmm2 align --preset SUBREAD`
- `CCS` → `pbmm2 align --preset CCS`
- `PacBio_HiFi` / `PacBio_HiFi_BAM` / `REVIO` → `pbmm2 align --preset HiFi`
- `HiFi_minimap` → `minimap2 -ax map-hifi`

### ONT

- `ONT` / `ONT_UL` → `minimap2 -ax map-ont`
- `ONT_methyl` / `ONT_BAM` → `minimap2 -ax map-ont -y`
- `ONT_Q20` → `minimap2 -ax lr:hq`

### Illumina

- `Illumina` → `bwa mem -Y -K 100000000`

The workflow uses `samtools` and `sambamba` for BAM/CRAM processing, sorting, merging, indexing, and duplicate marking.

## Running the pipeline

### Option 1: use a `snakesub` alias

This tutorial assumes the following alias is defined in your shell profile:

```bash
alias snakesub='mkdir -p log; snakemake --ri --jobname "{rulename}.{jobid}" --drmaa " -V -cwd -j y -o ./log -e ./log -l h_rt={resources.hrs}:00:00 -l mfree={resources.mem}G -pe serial {threads} -w n -S /bin/bash" -w 60'
```

Run the default workflow with:

```bash
snakesub -s /path/to/align_all/align_all.smk -j {JOBS} -kp --use-singularity --singularity-args "--bind /net:/net"
```

### Option 2: use the bundled wrapper script

You can symlink the wrapper and run it directly:

```bash
ln -s /path/to/align_all/runsnake .
./runsnake 30
```

You may also pass additional Snakemake options through the wrapper, for example a dry run:

```bash
./runsnake 30 -n
```

### Option 3: local execution

For local runs, use:

```bash
./runlocal 8
```

## Default target behavior

By default, the workflow target is defined by `rule all`, which produces merged BAM outputs:

```text
{ref}/{aln}/{sample}.all.sorted.bam
{ref}/{aln}/{sample}.all.sorted.bam.bai
```

This means that running the pipeline without an explicit target will **not** generate CRAM files.

## CRAM generation

CRAM output is implemented as a separate target rule.

If you want final CRAM files, explicitly request `get_crams`.

### Using `snakesub`

```bash
snakesub -s /path/to/align_all/align_all.smk -j {JOBS} -kp --use-singularity --singularity-args "--bind /net:/net" get_crams
```

### Using `runsnake`

```bash
./runsnake 30 get_crams
```

This generates:

```text
{ref}/{aln}/cram/{sample}.final.cram
{ref}/{aln}/cram/{sample}.final.cram.crai
```

CRAM conversion is performed from the duplicate-marked BAM produced by `mark_duplicates`.

## Workflow details

### Reference preparation

The workflow symlinks the requested reference into a local `ref/` directory and creates required indices.

### ONT BAM input handling

For `ONT_BAM`, the workflow first converts an **unaligned BAM** to gzipped FASTQ using `samtools fastq -T "*"`, then builds an index before alignment. This preserves BAM tags by passing all available tags through FASTQ header output during conversion, so tag information can be retained during the handoff from BAM input to the mapping step.

### Split alignment

For selected non-Illumina and non-`PacBio_HiFi_BAM` inputs, the workflow can split reads into batches and align them separately before merge.

During split mapping, the workflow uses `samtools-fqidx-comment` to preserve FASTQ header comments when extracting read subsets.

### Merge and duplicate marking

Per-read or per-batch BAMs are merged into:

```text
{ref}/{aln}/{sample}.all.sorted.bam
```

Duplicate marking produces:

```text
{ref}/{aln}/{sample}.all.sorted.md.bam
```

## Outputs

### Default outputs

The default target produces merged BAM outputs:

```text
{ref}/{aln}/{sample}.all.sorted.bam
{ref}/{aln}/{sample}.all.sorted.bam.bai
```

### Intermediate / downstream outputs

Depending on execution path, the workflow may also produce:

```text
{ref}/{aln}/{sample}.{read}.sorted.bam
{ref}/{aln}/{sample}.all.sorted.md.bam
{ref}/{aln}/{sample}.all.sorted.md.bam.bai
```

### Optional CRAM outputs

When `get_crams` is explicitly requested:

```text
{ref}/{aln}/cram/{sample}.final.cram
{ref}/{aln}/cram/{sample}.final.cram.crai
```

## Example minimal setup

### `map_contigs.yaml`

```yaml
MANIFEST: aln.tab
REF:
  GRCh38: /net/eichler/vol28/eee_shared/assemblies/hg38/no_alt/hg38.no_alt.fa
NBATCHES: 15
```

### `aln.tab`

```text
SAMPLE  FOFN  TYPE
HG002   fofn/HG002.hifi.fofn  PacBio_HiFi
HG002   fofn/HG002.ont.fofn   ONT
HG002   fofn/HG002.ont.ubam.fofn   ONT_BAM
```

### run default BAM target

```bash
./runsnake 30
```

### run CRAM target

```bash
./runsnake 30 get_crams
```

## Notes

- `PacBio_HiFi_BAM` and `Illumina` are handled as whole-file alignments rather than split-batch extraction.
- `ONT_methyl`, `ONT_BAM`, and related ONT modes include explicit read-group handling in alignment parameters.
- The workflow defines both `all` and `get_crams` as local rules.

