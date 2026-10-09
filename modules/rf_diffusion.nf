def inferenceCommand(meta, target, models, settings, isStub) {
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def runner = settings.use_containers ? '"\$(command -v run_rfdiffusion.py)"' : quote("${projectDir}/bin/run_rfdiffusion.py")
    def arguments = [
        '--sample-id', meta.id, '--target', target, '--contigs', meta.contigs,
        '--hotspots', meta.hotspots, '--models', models,
        '--script', settings.rfdiffusion_script,
        '--num-designs', settings.num_designs, '--design-startnum', settings.design_startnum,
        '--deterministic', settings.deterministic,
        '--noise-scale-ca', settings.noise_scale_ca,
        '--noise-scale-frame', settings.noise_scale_frame
    ]
    if (isStub) arguments.add('--stub')
    // The container image pip-installs RFdiffusion and sets DGLBACKEND itself. The conda
    // envs supply only the dependencies, so point Python at the checkout --script names.
    def prelude = ''
    if (!settings.use_containers && !isStub) {
        def repository = new File(settings.rfdiffusion_script.toString()).parentFile.parent
        prelude = 'export DGLBACKEND=pytorch\n' +
                  "export PYTHONPATH=\"${repository}:${repository}/env/SE3Transformer\${PYTHONPATH:+:\$PYTHONPATH}\"\n"
    }
    prelude + "${quote(settings.python)} ${runner} " + arguments.collect(quote).join(' ')
}

process rf_diffusion {
    tag "${meta.id}"
    container params.use_containers ? params.container_rfdiffusion : null
    conda params.use_containers ? null : "${projectDir}/envs/rfdiffusion.yml"
    cpus params.cpus
    memory params.memory
    maxForks params.max_forks

    input:
    tuple val(meta), path(target_pdb, stageAs: 'target.pdb')
    path models
    val code_fingerprint

    output:
    tuple val(meta), path("outputs/${meta.id}_*.pdb"), emit: backbones
    tuple val(meta), path("outputs/${meta.id}_*.trb"), emit: metadata
    tuple val(meta), path('outputs/manifest.csv'), path('outputs/run.json'), path('outputs/inference.log'), emit: reports

    script:
    inferenceCommand(meta, target_pdb, models, params, false)

    stub:
    inferenceCommand(meta, target_pdb, models, params, true)
}
