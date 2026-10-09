def mpnnCommand(meta, backbone, settings, stub) {
    def q = { v -> "'" + v.toString().replace("'", "'\"'\"'") + "'" }
    def runner = settings.use_containers ? '"\$(command -v run_mpnn.py)"' : q("${projectDir}/bin/run_mpnn.py")
    def args = ['--backbone', backbone, '--sample-id', meta.id, '--backbone-id', meta.backbone_id,
                '--contigs', meta.contigs, '--script', settings.proteinmpnn_script,
                '--count', settings.sequences_per_backbone, '--temperature', settings.mpnn_temperature,
                '--seed', settings.sequence_seed]
    if (stub) args.add('--stub')
    "${q(settings.proteinmpnn_python)} ${runner} " + args.collect(q).join(' ')
}

process protein_mpnn {
    tag "${meta.backbone_id}"
    container params.use_containers ? params.container_proteinmpnn : null
    conda params.use_containers ? null : "${projectDir}/envs/proteinmpnn.yml"
    cpus params.cpus
    memory params.memory
    maxForks params.max_forks

    input:
    tuple val(meta), path(backbone)
    val code_fingerprint

    output:
    tuple val(meta), path('sequences'), emit: sequences

    script:
    mpnnCommand(meta, backbone, params, false)

    stub:
    mpnnCommand(meta, backbone, params, true)
}
