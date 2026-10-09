def rankQuote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

process developability_features {
    tag "${meta.backbone_id}"
    container params.use_containers ? params.container_features : null
    conda params.use_containers ? null : "${projectDir}/envs/features.yml"
    containerOptions ''
    cpus params.cpus
    memory params.memory

    input:
    tuple val(meta), path(predictions)
    val code_fingerprint

    output:
    path "${meta.backbone_id}.features.json", emit: features

    // Also executed in stub mode: descriptors are computed on synthetic PDBs.
    script:
    "${rankQuote(params.features_python)} ${params.use_containers ? '\"\$(command -v rank_binders.py)\"' : rankQuote("${projectDir}/bin/rank_binders.py")} features --predictions ${rankQuote(predictions)} --output ${rankQuote(meta.backbone_id + '.features.json')}"
}

process rank_candidates {
    container params.use_containers ? params.container_features : null
    conda params.use_containers ? null : "${projectDir}/envs/features.yml"
    containerOptions ''
    cpus 1
    memory params.memory

    input:
    path features
    path config
    val code_fingerprint

    output:
    path 'ranked_candidates.csv', emit: csv
    path 'ranked_candidates.json', emit: json
    path 'ranked_binders.fasta', emit: fasta
    path 'ranking_method.json', emit: method

    script:
    def files = (features instanceof List ? features : [features]).collect { rankQuote(it) }.join(' ')
    "${rankQuote(params.features_python)} ${params.use_containers ? '\"\$(command -v rank_binders.py)\"' : rankQuote("${projectDir}/bin/rank_binders.py")} rank --config ${rankQuote(config)} --inputs ${files}"
}
