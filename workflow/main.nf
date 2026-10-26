#!/usr/bin/env nextflow
// STV_GSE137303 single-end RNA-seq workflow: samplesheet -> FastQC -> [Salmon] -> MultiQC.
// Architecture follows nf-core/rnaseq (FASTQC, SALMON_INDEX/QUANT, MULTIQC modules) at a much smaller scale.
// Salmon runs only with --run_quant true and a reference (approved: SGN SL4.0 genome decoys + ITAG4.0 cDNA, D-001).

include { FASTQC       } from './modules/fastqc'
include { SALMON_INDEX } from './modules/salmon_index'
include { SALMON_QUANT } from './modules/salmon_quant'
include { MULTIQC      } from './modules/multiqc'

params {
    // Samplesheet; relative paths (sheet and fastq_1) resolve against project_dir.
    input: String = 'config/samplesheet.csv'
    // Project root; empty = parent directory of this workflow/ folder.
    project_dir: String = ''

    // Quantification (off by default).
    run_quant: Boolean = false
    salmon_index: String = ''
    transcript_fasta: String = ''
    genome_fasta: String = ''
    // Overrides the per-sample strandedness mapping (e.g. 'A', 'SR'); empty = use samplesheet.
    salmon_libtype: String = ''
    // salmon index -k (Salmon default; docs: k = 31 works well for reads >= 75 bp).
    salmon_kmer: Integer = 31
    // Single-end reads carry no fragment-length information: Salmon uses these as the fragment-length
    // distribution for effective lengths. Library fragment length is not reported (SRA InsertSize = 0),
    // so the Salmon defaults are kept and passed explicitly (docs/salmon_quantification.md).
    salmon_fld_mean: Integer = 250
    salmon_fld_sd: Integer = 25
    // salmon quant --gcBias (D-012; recommended by the DESeq2 vignette). --seqBias is not used.
    salmon_gc_bias: Boolean = true

    multiqc_title: String = 'STV_GSE137303 raw FASTQ QC'

    // Publish locations relative to outputDir (default: results/).
    fastqc_dir: String = 'qc/fastqc_raw'
    multiqc_dir: String = 'qc/multiqc_raw'
    salmon_dir: String = 'quantification/salmon'
    salmon_index_dir: String = 'genome/index'
    pipeline_info_dir: String = 'pipeline_info'
}

workflow {
    main:
    def root_dir = params.project_dir ? file(params.project_dir) : projectDir.parent
    def sheet = resolvePath(params.input, root_dir)
    def samples = readSamplesheet(sheet, root_dir)
    log.info "Samplesheet: ${sheet} | samples: ${samples.size()} | layout: SINGLE | run_quant: ${params.run_quant}"

    ch_reads = channel.fromList(samples)

    FASTQC(ch_reads)
    ch_multiqc = FASTQC.out.zip.map { _meta, zip -> zip }

    ch_salmon = channel.empty()
    ch_index_files = channel.empty()
    if (params.run_quant) {
        if (!params.salmon_index && !params.transcript_fasta) {
            error("--run_quant true requires --salmon_index or --transcript_fasta")
        }
        def ch_index = params.salmon_index ? channel.value(resolvePath(params.salmon_index, root_dir)) : null
        if (!params.salmon_index) {
            def decoys = params.genome_fasta ? resolvePath(params.genome_fasta, root_dir) : []
            SALMON_INDEX(resolvePath(params.transcript_fasta, root_dir), decoys)
            ch_index = SALMON_INDEX.out.index
            ch_index_files = SALMON_INDEX.out.index.mix(SALMON_INDEX.out.decoys, SALMON_INDEX.out.log)
        }
        def ch_quant_in = ch_reads.map { meta, fastq ->
            tuple(meta + [libtype: params.salmon_libtype ?: salmonLibtype(meta.strandedness)], fastq)
        }
        SALMON_QUANT(ch_quant_in, ch_index)
        ch_salmon = SALMON_QUANT.out.results
        ch_multiqc = ch_multiqc.mix(SALMON_QUANT.out.results.map { _meta, dir -> dir })
    }

    MULTIQC(ch_multiqc.collect(), params.multiqc_title)

    ch_versions = channel.topic('versions')
        .map { _process_name, tool, version -> "${tool}\t${version}" }
        .unique()
        .collectFile(name: 'software_versions.tsv', newLine: true, sort: true)

    publish:
    fastqc_raw = FASTQC.out.html.mix(FASTQC.out.zip)
    multiqc_raw = MULTIQC.out.report.mix(MULTIQC.out.data)
    salmon = ch_salmon.map { meta, dir ->
        [sample: meta.id, condition: meta.condition, replicate: meta.replicate, libtype: meta.libtype, quant_dir: dir]
    }
    salmon_index = ch_index_files
    versions = ch_versions
}

output {
    fastqc_raw {
        path params.fastqc_dir
    }
    multiqc_raw {
        path params.multiqc_dir
    }
    salmon {
        path params.salmon_dir
        index {
            path "${params.salmon_dir}/samples.csv"
            header true
        }
    }
    salmon_index {
        path params.salmon_index_dir
    }
    versions {
        path params.pipeline_info_dir
    }
}

// ---------------------------------------------------------------------------
// Samplesheet parsing and validation (fails before any task is submitted)
// ---------------------------------------------------------------------------

def resolvePath(path_str, root_dir) {
    def p = path_str.startsWith('/') ? file(path_str) : root_dir.resolve(path_str)
    if (!p.exists()) {
        error("File not found: ${p}")
    }
    return p
}

def readSamplesheet(sheet, root_dir) {
    def required = ['sample', 'condition', 'replicate', 'run_accession', 'library_layout', 'strandedness', 'fastq_1']
    def rows = sheet.splitCsv(header: true)
    if (rows.isEmpty()) {
        error("Samplesheet is empty: ${sheet}")
    }
    def missing = required.findAll { col -> !rows[0].containsKey(col) }
    if (missing) {
        error("Samplesheet ${sheet} is missing column(s): ${missing.join(', ')}")
    }
    def duplicated = rows.countBy { row -> row.sample }.findAll { _id, n -> n > 1 }.keySet()
    if (duplicated) {
        error("Duplicated sample IDs in ${sheet}: ${duplicated.join(', ')}")
    }
    return rows.collect { row -> sampleFromRow(row, root_dir) }
}

def sampleFromRow(row, root_dir) {
    def id = row.sample
    if (!(id ==~ /[A-Za-z0-9_.-]+/)) {
        error("Invalid sample ID '${id}' (allowed: letters, digits, '_', '.', '-')")
    }
    if (row.library_layout != 'SINGLE') {
        error("Sample ${id}: library_layout '${row.library_layout}' is not supported; this workflow is single-end only")
    }
    if (row.fastq_2) {
        error("Sample ${id}: fastq_2 given but the dataset is single-end")
    }
    if (!(row.strandedness in ['auto', 'unstranded', 'forward', 'reverse'])) {
        error("Sample ${id}: strandedness must be auto, unstranded, forward or reverse (got '${row.strandedness}')")
    }
    if (!row.condition) {
        error("Sample ${id}: empty condition")
    }
    if (!(row.fastq_1 ==~ /.+\.(fastq|fq)\.gz/)) {
        error("Sample ${id}: fastq_1 must be a .fastq.gz/.fq.gz file (got '${row.fastq_1}')")
    }
    def fastq = resolvePath(row.fastq_1, root_dir)
    def meta = [
        id: id,
        condition: row.condition,
        replicate: row.replicate,
        run: row.run_accession,
        single_end: true,
        strandedness: row.strandedness,
    ]
    return tuple(meta, fastq)
}

def salmonLibtype(strandedness) {
    // Single-end Salmon library types (Salmon docs): U, SF, SR; 'A' = automatic detection.
    def types = [auto: 'A', unstranded: 'U', forward: 'SF', reverse: 'SR']
    return types[strandedness]
}
