#!/usr/bin/env nextflow
// Bootstrap test: proves launch, process execution, workflow-output publishing, task cache and -resume.
// Uses generated toy values only (no RNA-seq data). Run from the project root:
//   nextflow run workflow/bootstrap/main.nf -w work/bootstrap
//   nextflow run workflow/bootstrap/main.nf -w work/bootstrap -resume

params {
    n_items: Integer = 3
}

process MAKE_TEXT {
    tag "item_${id}"

    input:
    val id

    output:
    path "item_${id}.txt", emit: txt

    script:
    """
    printf 'bootstrap item %s\\n' ${id} > item_${id}.txt
    """
}

process CHECKSUM {
    input:
    path txts

    output:
    path 'checksums.md5', emit: md5

    script:
    """
    md5sum ${txts} | sort -k2 > checksums.md5
    """
}

workflow {
    main:
    ids = channel.fromList((1..params.n_items).toList())
    MAKE_TEXT(ids)
    CHECKSUM(MAKE_TEXT.out.txt.collect())

    publish:
    items = MAKE_TEXT.out.txt
    checksums = CHECKSUM.out.md5
}

output {
    items {
        path 'items'
    }
    checksums {
        path '.'
    }
}
