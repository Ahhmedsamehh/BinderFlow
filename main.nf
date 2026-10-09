#!/usr/bin/env nextflow

include { rf_diffusion } from './modules/rf_diffusion.nf'
include { protein_mpnn } from './modules/protein_mpnn.nf'
include { predict_complex } from './modules/af.nf'
include { developability_features; rank_candidates } from './modules/ranker.nf'

workflow {
    main:
    if (!params.input) error 'Provide --input samples.csv or -params-file params.yaml.'
    if (!(params.stop_after in ['rfdiffusion', 'full'])) error '--stop_after must be rfdiffusion or full.'
    if (!workflow.stubRun && !params.models) error 'Provide --models /path/to/weights (or MODELS_PATH).'
    if (workflow.profile.tokenize(',').contains('stub') && !workflow.stubRun) {
        error 'The stub profile requires -stub-run; it does not perform inference.'
    }
    ['num_designs', 'cpus', 'max_forks'].each { name ->
        if (!(params[name].toString() ==~ /[1-9][0-9]*/)) error "--${name} must be a positive integer."
    }
    if (!(params.design_startnum.toString() ==~ /[0-9]+/)) error '--design_startnum must be a nonnegative integer.'
    def sheet = file(params.input, checkIfExists: true)
    // Include imported Python helpers in task cache keys, not only entry scripts.
    def source = file("${projectDir}/bin/*.py").sort { it.name }.collect { it.name + '\n' + it.text }.join('\n')
    def code_fingerprint = java.security.MessageDigest.getInstance('SHA-256').digest(source.getBytes('UTF-8')).encodeHex().toString()
    def model_dir = file(workflow.stubRun ? "${projectDir}/tests/fixtures" : params.models, checkIfExists: true)
    if (!java.nio.file.Files.isDirectory(model_dir)) error '--models must name a directory.'
    def af_models = null
    def ranking_config = null
    if (params.stop_after == 'full') {
        ['sequences_per_backbone', 'sequence_seed'].each { name ->
            if (!(params[name].toString() ==~ /[1-9][0-9]*/)) error "--${name} must be a positive integer."
        }
        if (!workflow.stubRun && !params.colabfold_models) error 'Provide --colabfold_models /path/to/alphafold/weights.'
        af_models = file(workflow.stubRun ? "${projectDir}/tests/fixtures" : params.colabfold_models, checkIfExists: true)
        if (!java.nio.file.Files.isDirectory(af_models)) error '--colabfold_models must name a directory.'
        ranking_config = file(params.ranking_config, checkIfExists: true)
    }

    // Collect first: duplicate IDs or a bad row must fail before any inference.
    inputs_channel = channel.fromPath(sheet)
        .splitCsv(header: true)
        .toList()
        .flatMap { rows ->
            if (!rows) error 'The sample sheet contains no samples.'
            def seen = [] as Set
            rows.collect { row ->
                if (!row.keySet().containsAll(['id', 'target_pdb', 'contigs', 'hotspots'])) {
                    error 'Required CSV columns: id,target_pdb,contigs,hotspots.'
                }
                def id = row.id?.trim()
                if (!id || !(id ==~ /[A-Za-z0-9][A-Za-z0-9_-]*/)) error 'Sample IDs must contain only letters, numbers, underscores or hyphens.'
                if (!seen.add(id)) error "Duplicate sample ID: ${id}"
                if (!row.target_pdb?.trim() || !row.contigs?.trim()) error "Missing target_pdb or contigs for ${id}."
                def target = java.nio.file.Paths.get(row.target_pdb.trim())
                if (!target.isAbsolute()) target = sheet.parent.resolve(target)
                target = file(target.normalize(), checkIfExists: true)
                if (!java.nio.file.Files.isRegularFile(target)) error "Target is not a file: ${target}"
                tuple([id: id, contigs: row.contigs.trim(),
                       hotspots: (row.hotspots ?: '').trim().replace(';', ',')], target)
            }
        }
    rf_diffusion(inputs_channel, model_dir, code_fingerprint)
    if (params.stop_after == 'full') {
        designs = rf_diffusion.out.backbones.flatMap { meta, files ->
            (files instanceof List ? files : [files]).collect { backbone ->
                tuple(meta + [backbone_id: backbone.baseName], backbone)
            }
        }
        protein_mpnn(designs, code_fingerprint)
        predict_complex(protein_mpnn.out.sequences, af_models, code_fingerprint)
        developability_features(predict_complex.out.predictions, code_fingerprint)
        rank_candidates(developability_features.out.features.collect(), ranking_config, code_fingerprint)
    }

    publish:
    backbones = rf_diffusion.out.backbones
    metadata = rf_diffusion.out.metadata
    reports = rf_diffusion.out.reports
    sequences = params.stop_after == 'full' ? protein_mpnn.out.sequences : channel.empty()
    predictions = params.stop_after == 'full' ? predict_complex.out.predictions : channel.empty()
    features = params.stop_after == 'full' ? developability_features.out.features : channel.empty()
    ranked_csv = params.stop_after == 'full' ? rank_candidates.out.csv : channel.empty()
    ranked_json = params.stop_after == 'full' ? rank_candidates.out.json : channel.empty()
    ranked_fasta = params.stop_after == 'full' ? rank_candidates.out.fasta : channel.empty()
    ranking_method = params.stop_after == 'full' ? rank_candidates.out.method : channel.empty()
}

output {
    backbones { path { meta, files ->
        (files instanceof List ? files : [files]).each { f -> f >> "${meta.id}/backbones/${f.name}" }
    } }
    metadata { path { meta, files ->
        (files instanceof List ? files : [files]).each { f -> f >> "${meta.id}/metadata/${f.name}" }
    } }
    reports { path { meta, manifest, run, log ->
        [manifest, run, log].each { f -> f >> "${meta.id}/reports/${f.name}" }
    } }
    sequences { path { meta, folder -> folder >> "${meta.id}/sequences/${meta.backbone_id}" } }
    predictions { path { meta, folder -> folder >> "${meta.id}/predictions/${meta.backbone_id}" } }
    features { path { f -> f >> "features/${f.name}" } }
    ranked_csv { path '.' }
    ranked_json { path '.' }
    ranked_fasta { path '.' }
    ranking_method { path '.' }
}
