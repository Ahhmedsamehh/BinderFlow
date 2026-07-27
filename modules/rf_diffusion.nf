process rf_diffusion {
    container 'ahhmedsamehh/rfdiffusion-standalone:latest'

    input:
    tuple val(meta), path(target_pdb)

    output:
    tuple val(meta), path("outputs/${meta.id}_*.pdb"), emit: backbones
    tuple val(meta), path("outputs/${meta.id}_*.trb"), emit: metadata

    script:
    """
    python3.9 /app/RFdiffusion/scripts/run_inference.py \\
        inference.input_pdb=${target_pdb} \\
        inference.output_prefix=outputs/${meta.id} \\
        inference.model_directory_path=\$MODELS_PATH \\
        inference.num_designs=${params.num_designs} \\
        inference.write_trajectory=False \\
        'contigmap.contigs=[${meta.contigs}]' \\
        'ppi.hotspot_res=[${meta.hotspots}]' \\
        denoiser.noise_scale_ca=${params.noise_scale_ca} \\
        denoiser.noise_scale_frame=${params.noise_scale_frame}
    """

}
