process AGGREGATE_RESULTS_KFOLD_CV {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir_prefix}/ref_${params.subsample_ref}/aggregated_results/kfold_cv", mode: 'copy'

    input:
    path cv_metrics

    output:
    path "**all_metrics.tsv", emit: all_cv_metrics

    script:
    """
    python $projectDir/bin/aggregate_results.py --metrics ${cv_metrics}
    """
}

process PLOT_KFOLD_CV_RESULTS {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir_prefix}/ref_${params.subsample_ref}/plots/kfold_cv", pattern: "**png", mode: 'copy'

    input:
    path all_cv_metrics

    output:
    path "**.png"

    script:
    ref_keys = params.ref_keys.join(" ")
    """
    python $projectDir/bin/plot_results.py --aggregated_results ${all_cv_metrics} \\
        --color_mapping_file ${params.color_mapping_file} \\
        --mapping_file ${params.mapping_file} \\
        --ref_keys ${ref_keys}
    """
}


process AGGREGATE_RESULTS_LOOCV {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir_prefix}/ref_${params.subsample_ref}/aggregated_results/loocv", mode: 'copy'

    input:
    path cv_metrics

    output:
    path "**all_metrics.tsv", emit: all_cv_metrics

    script:
    """
    python $projectDir/bin/aggregate_results.py --metrics ${cv_metrics}
    """
}

process PLOT_LOOCV_RESULTS {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir_prefix}/ref_${params.subsample_ref}/plots/loocv", pattern: "**png", mode: 'copy'

    input:
    path all_cv_metrics

    output:
    path "**.png"

    script:
    ref_keys = params.ref_keys.join(" ")
    """
    python $projectDir/bin/plot_results.py --aggregated_results ${all_cv_metrics} \\
        --color_mapping_file ${params.color_mapping_file} \\
        --mapping_file ${params.mapping_file} \\
        --ref_keys ${ref_keys}
    """
}

workflow {
    // Collect cross-validation results
    Channel
        .fromPath("${params.outdir_prefix}/ref_${params.subsample_ref}/*/kfold_cv/*/*/*/**summary.scores.tsv", type: 'file')
        .flatten()
        .toList()
        .set { kfold_cv_results_files }


    Channel
        .fromPath("${params.outdir_prefix}/ref_${params.subsample_ref}/*/loocv/*/*/*/**summary.scores.tsv", type: 'file')
        .flatten()
        .toList()
        .set { loocv_results_files }

    // Aggregate KFOLD cross-validation results
    AGGREGATE_RESULTS_KFOLD_CV(kfold_cv_results_files)
    .set { all_kfold_cv_metrics }
    // Plot KFOLD cross-validation results
    PLOT_KFOLD_CV_RESULTS(all_kfold_cv_metrics)



    // Aggregate LOOCV results
    AGGREGATE_RESULTS_LOOCV(loocv_results_files)
    .set { all_loocv_metrics }
    // Plot LOOCV results
    PLOT_LOOCV_RESULTS(all_loocv_metrics)
}

workflow.onComplete {
    println "Workflow completed!"
    println "Plots in ${params.outdir_prefix}/ref_${params.subsample_ref}/plots"
}