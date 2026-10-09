def afCommand(sequences, models, settings, stub) {
    def q = { v -> "'" + v.toString().replace("'", "'\"'\"'") + "'" }
    def runner = settings.use_containers ? '"\$(command -v run_colabfold.py)"' : q("${projectDir}/bin/run_colabfold.py")
    def args = ['--sequences', sequences, '--models', models,
                '--executable', settings.colabfold_executable, '--num-models', settings.af_num_models,
                '--recycles', settings.af_recycles, '--seed', settings.structure_seed]
    if (stub) args.add('--stub')
    "${q(settings.colabfold_python)} ${runner} " + args.collect(q).join(' ')
}

process predict_complex {
    tag "${meta.backbone_id}"
    container params.use_containers ? params.container_colabfold : null
    conda params.use_containers ? null : "${projectDir}/envs/colabfold.yml"
    cpus params.cpus
    memory params.af_memory
    maxForks params.max_forks

    input:
    tuple val(meta), path(sequences)
    path models
    val code_fingerprint

    output:
    tuple val(meta), path('predictions'), emit: predictions

    script:
    afCommand(sequences, models, params, false)

    stub:
    afCommand(sequences, models, params, true)
}
