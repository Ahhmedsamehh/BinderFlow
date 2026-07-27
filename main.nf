#!/usr/bin/env nextflow

include { rf_diffusion } from './modules/rf_diffusion.nf'

workflow {
    main:
    inputs_channel = channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [ id: row.id, contigs: row.contigs.trim(), hotspots: row.hotspots.trim().replaceAll(';', ',') ]
            tuple(meta, file(row.target_pdb, checkIfExists: true))
        }
    rf_diffusion(inputs_channel)

    publish:
    rf_diffusion_output_backbones = rf_diffusion.out.backbones
    rf_diffusion_output_metadata = rf_diffusion.out.metadata
}

output {
    rf_diffusion_output_backbones {
        path { meta, files -> "${meta.id}/backbones" }
    }
    rf_diffusion_output_metadata {
        path { meta, files -> "${meta.id}/metadata" }
    }
}
