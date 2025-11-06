#!/usr/bin/env nextflow

process SAVE_PARAMS_TO_FILE {
    publishDir (
        "${params.outdir}",
        mode: "copy"
    )

    output:
    file "params.yaml"

    script:
    """
    cat <<EOF > params.yaml
    organism: ${params.organism}
    census_version: ${params.census_version}
    ref_keys: ${params.ref_keys}
    subsample_ref: ${params.subsample_ref}
    cutoff: ${params.cutoff}
    remove_unknown: ${params.remove_unknown}
    ref_split: ${params.ref_split}
    ref_collections: ${params.ref_collections}
    integration_method: ${params.integration_method}
    dims: ${params.dims}
    max_features: ${params.max_features}
    k_anchor: ${params.k_anchor}
    k_score: ${params.k_score}
    k_weight: ${params.k_weight}
    outdir: ${params.outdir}
    normalization_method: ${params.normalization_method}
    subset_type: ${params.subset_type}
    batch_correct: ${params.batch_correct}
    nmads: ${params.nmads}
    use_gap: ${params.use_gap}
    EOF
    """
}

process RUN_SETUP {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    input:
    val organism
    val census_version

    output:
    path "scvi-${params.organism}-${census_version}/"

    script:
    """
    python $projectDir/bin/setup.py --organism ${organism} --census_version ${census_version}
    """
}

process GET_CENSUS_ADATA {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    publishDir path: "${params.outdir}/refs/scvi_integrated/", pattern: "**.h5ad", mode: 'copy'
    publishDir path: "${params.outdir}/refs/obs/", pattern: "**obs.tsv", mode: 'copy'

    input:
    val ref_collections

    output:
    path "**.h5ad", emit: ref_paths_adata
    path "**obs.tsv", emit: ref_obs
    path "**yaml", emit: ref_region_mapping

    script:
    ref_keys = params.ref_keys.join(' ')
    """
    python $projectDir/bin/get_census_adata.py \\
        --organism ${params.organism} \\
        --census_version ${params.census_version} \\
        --subsample_ref ${params.subsample_ref} \\
        --relabel_path ${params.mapping_file} \\
        --split_column ${params.ref_split} \\
        --ref_collections ${ref_collections} \\
        --ref_keys ${ref_keys} \\
        --seed ${params.seed}
    """
}

process REF_PROCESS_SEURAT {
    conda '/home/rschwartz/anaconda3/envs/r4.3'
    publishDir path: "${params.outdir}/refs/SCT_integrated/", pattern: "**.rds", mode: 'copy'

    input:
    path h5ad_file

    output:
        path "${h5ad_file.getName().toString().replace('.h5ad','.rds')}", emit: ref_paths_seurat
        tuple val(ref_name), path(h5ad_file), path("${h5ad_file.getName().toString().replace('.h5ad','_sct_integrated.h5ad')}"), emit: sct_scvi_integrated_h5ads

    script:
    ref_name = h5ad_file.getName().toString().replace('.h5ad','')
    """
    Rscript $projectDir/bin/ref_preprocessing.R --h5ad_file ${h5ad_file} \\
            --normalization_method ${params.normalization_method} \\
            --dims ${params.dims} \\
            --nfeatures ${params.nfeatures} \\
            --k.score ${params.k_score} \\
            --k.anchor ${params.k_anchor} \\
            --k.weight ${params.k_weight} \\
            ${params.batch_correct ? '--batch_correct' : ''} 
    """
}

process EVALUATE_INTEGRATION {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    // Only use subsample_ref in the output directory, not cutoff
    publishDir path: "${params.outdir}/refs/integration/${ref_name}/", mode: 'copy'

    input:
    tuple val(ref_name), path(scvi_h5ad), path(sct_h5ad)
    output:
    path "**png"
    path "**tsv"
    path "**svg"

    script:
    batch_fields = params.batch_fields.join(" ")
    """
    python $projectDir/bin/evaluate_integration.py --ref_name ${ref_name} \\
        --scvi_h5ad ${scvi_h5ad} \\
        --sct_h5ad ${sct_h5ad} \\
        --batch_fields ${batch_fields}
    """
}
process CV_SEURAT {
    conda '/home/rschwartz/anaconda3/envs/r4.3'
    publishDir path: "${params.outdir}/cutoff_${params.cutoff}/kfold_cv/seurat/${ref_name}/probs", pattern: "**tsv", mode: "copy"

    input:
    val ref_path
    val ref_keys

    output:
    tuple path("*obs.tsv"), val(ref_path), path("*prediction_scores.tsv"), emit: pred_scores_channel 

    script:
    ref_name = ref_path.getName().toString().replace('.rds','')
    """
    Rscript $projectDir/bin/KFOLD_CV_seurat.R --ref_path ${ref_path} \\
        --ref_keys ${ref_keys} \\
        --integration_method ${params.integration_method} \\
        --dims ${params.dims} \\
        --max.features ${params.max_features} \\
        --k.score ${params.k_score} \\
        --k.anchor ${params.k_anchor} \\
        --k.weight ${params.k_weight} \\
        --normalization_method ${params.normalization_method} \\
        --n_folds ${params.n_folds}
    """
}

process CV_SCVI {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir path: "${params.outdir}/cutoff_${params.cutoff}/kfold_cv/scvi/${ref_name}/probs", pattern: "**tsv", mode: "copy"

    input:
    val ref_path
    val ref_keys

    output:
    tuple path("*obs.tsv"), val(ref_path), path("probs/*tsv"), emit: probs_channel

    script:
    ref_name = ref_path.getName().toString().replace('.h5ad','')
    """
    python $projectDir/bin/KFOLD_CV_scvi.py --ref_path ${ref_path} \\
                 --ref_keys ${ref_keys}  \\
                 --n_folds ${params.n_folds}     
    """
}

process CLASSIFY_KFOLD_CV {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir path: "${params.outdir}/cutoff_${params.cutoff}/kfold_cv/${method}/${ref_name}/label_transfer_metrics", pattern: "label_transfer_metrics**", mode: 'copy'
    publishDir path: "${params.outdir}/cutoff_${params.cutoff}/kfold_cv/${method}/${ref_name}/confusion", pattern: "confusion**", mode: 'copy'
    publishDir path: "${params.outdir}/cutoff_${params.cutoff}/kfold_cv/${method}/${ref_name}/predicted_meta", pattern: "predicted_meta**", mode: 'copy'

    input:
    val ref_keys
    tuple val(method), path(cv_obs), path(ref_path), path(probs_path)
    val ref_region_mapping

    output:
    path("**summary.scores.tsv"), emit: summary_score_channel
    path "confusion/**"
    tuple val(method), path("${ref_path}"), path("predicted_meta/**tsv"), emit: predicted_meta_channel

    script:
    ref_name = ref_path.getName().split('\\.')[0]
    """
    python $projectDir/bin/classify.py \\
        --obs ${cv_obs} \\
        --ref_name ${ref_name} \\
        --ref_keys ${ref_keys} \\
        --cutoff ${params.cutoff} \\
        --probs ${probs_path} \\
        --mapping_file ${params.mapping_file} \\
        --ref_region_mapping ${ref_region_mapping} \\
        --method ${method} \\
        ${params.use_gap ? '--use_gap' : ''}
    """ 
}



process SEURAT_PAIR_PREDICT {
    conda '/home/rschwartz/anaconda3/envs/r4.3'
    publishDir path: "${params.outdir}/cutoff_${params.cutoff}/pairwise/seurat/${query_name}/${ref_name}", pattern: "**tsv", mode: "copy"

    input:
    tuple path(ref_path), path(query_path)
    val ref_keys

    output:
    tuple val("seurat"), val(ref_name), val(query_name), path("**obs.tsv"), path("**prediction_scores.tsv"), emit: predicted_probs_seurat

    script:
    ref_name = ref_path.getName().split('\\.')[0]
    query_name = query_path.getName().split('\\.')[0]
    """
    Rscript $projectDir/bin/predict_seurat.R --query_path ${query_path} \\
        --ref_path ${ref_path} \\
        --ref_keys ${ref_keys} \\
        --integration_method ${params.integration_method} \\
        --dims ${params.dims} \\
        --max.features ${params.max_features} \\
        --k.score ${params.k_score} \\
        --k.anchor ${params.k_anchor} \\
        --k.weight ${params.k_weight} \\
        --normalization_method ${params.normalization_method}
    """
}

process SCVI_PAIR_PREDICT {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir path: "${params.outdir}/cutoff_${params.cutoff}/pairwise/scvi/${query_name}/${ref_name}", pattern: "**tsv", mode: "copy"

    input:
    tuple path(ref_adata), path(query_adata)
    val ref_keys

    output:
    tuple val("scvi"), val(ref_name), val(query_name), path("**obs.tsv"), path("**prob.df.tsv"), emit: predicted_probs_scvi


    script:
    ref_name = ref_adata.getName().split('\\.')[0]
    query_name = query_adata.getName().split('\\.')[0]
    """
    python $projectDir/bin/predict_scvi.py --ref_path ${ref_adata} \\
        --query_path ${query_adata} \\
        --ref_keys ${ref_keys}
    """
}

// classify pairwise predictions with classify.py
process CLASSIFY_PAIRWISE {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir path: "${params.outdir}/cutoff_${params.cutoff}/pairwise/${method}/${query_name}/${ref_name}", pattern: "label_transfer_metrics**", mode: 'copy'
    publishDir path: "${params.outdir}/cutoff_${params.cutoff}/pairwise/${method}/${query_name}/${ref_name}", pattern: "confusion**", mode: 'copy'
    publishDir path: "${params.outdir}/cutoff_${params.cutoff}/pairwise/${method}/${query_name}/${ref_name}", pattern: "predicted_meta**", mode: 'copy'

    input:
    tuple val(method), val(ref_name), val(query_name), path(pairwise_obs), path(probs_path)
    val ref_keys
    val ref_region_mapping


    output:
    path("**summary.scores.tsv")
    path "confusion/**"
    path "predicted_meta/**tsv"

    script:
    """
    python $projectDir/bin/classify.py \\
        --obs ${pairwise_obs} \\
        --ref_name ${ref_name} \\
        --ref_keys ${ref_keys} \\
        --cutoff ${params.cutoff} \\
        --probs ${probs_path} \\
        --mapping_file ${params.mapping_file} \\
        --ref_region_mapping ${ref_region_mapping} \\
        --method ${method} \\
        ${params.use_gap ? '--use_gap' : ''}
    """ 
}

process LOOCV_SEURAT {
    conda '/home/rschwartz/anaconda3/envs/r4.3'
    publishDir path: "${params.outdir}/cutoff_${params.cutoff}/loocv/seurat/${ref_name}/", pattern: "**tsv", mode: "copy"

    input:
    path ref_path
    val ref_keys

    output:
    path "**loocv.obs.tsv", emit: loocv_obs_seurat
    path "**loocv.prediction.scores.tsv",  emit: loocv_probs_seurat

    script:
    ref_name = ref_path.getName().toString().replace('.rds','')
    """
    Rscript $projectDir/bin/LOOCV_seurat.R --ref_path ${ref_path} \\
        --ref_keys ${ref_keys} \\
        --integration_method ${params.integration_method} \\
        --dims ${params.dims} \\
        --max.features ${params.max_features} \\
        --k.score ${params.k_score} \\
        --k.anchor ${params.k_anchor} \\
        --k.weight ${params.k_weight} \\
        --normalization_method ${params.normalization_method}
    """
}

process LOOCV_SCVI {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir path: "${params.outdir}/cutoff_${params.cutoff}/loocv/scvi/${ref_name}/", pattern: "**tsv", mode: "copy"

    input:
    path ref_path
    val ref_keys

    output:
    path "**loocv.obs.tsv", emit: loocv_obs_scvi
    path "**loocv.prob.df.tsv", emit: loocv_probs_scvi

    script:
    ref_name = ref_path.getName().toString().replace('.h5ad','')
    """
    python $projectDir/bin/LOOCV_scvi.py --ref_path ${ref_path} \\
        --ref_keys ${ref_keys}
    """
}

process CLASSIFY_LOOCV {
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir path: "${params.outdir}/cutoff_${params.cutoff}/loocv/${method}/${held_out_name}/label_transfer_metrics", pattern: "label_transfer_metrics**", mode: 'copy'
    publishDir path: "${params.outdir}/cutoff_${params.cutoff}/loocv/${method}/${held_out_name}/confusion", pattern: "confusion**", mode: 'copy'
    publishDir path: "${params.outdir}/cutoff_${params.cutoff}/loocv/${method}/${held_out_name}/predicted_meta", pattern: "predicted_meta**", mode: 'copy'

    input:
    tuple val(held_out_name), val(method), path(loocv_probs), path(loocv_obs)
    val ref_keys
    val ref_region_mapping

    output:
    path("**summary.scores.tsv")
    path "confusion/**"
    path "predicted_meta/**tsv"

    script:
    """
    python $projectDir/bin/classify.py \\
        --query_name ${held_out_name} \\
        --obs ${loocv_obs} \\
        --ref_name whole_cortex \\
        --ref_keys ${ref_keys} \\
        --cutoff ${params.cutoff} \\
        --probs ${loocv_probs} \\
        --mapping_file ${params.mapping_file} \\
        --ref_region_mapping ${ref_region_mapping} \\
        --method ${method} \\
        ${params.use_gap ? '--use_gap' : ''}
    """ 
}

// Workflow definition
workflow {

    // Save parameters to file
    SAVE_PARAMS_TO_FILE()

    // ----------------------- prepare references -----------------------
    ref_keys = params.ref_keys.join(' ')
    model_path = RUN_SETUP(params.organism, params.census_version)
    ref_collections = params.ref_collections.collect { "\"${it}\"" }.join(' ')

    // get census adata for each reference collection
    GET_CENSUS_ADATA(ref_collections)
    GET_CENSUS_ADATA.out.ref_paths_adata.flatten().set { ref_paths_adata }
    GET_CENSUS_ADATA.out.ref_region_mapping.set { ref_region_mapping }

    // integrate with Seurat SCTransform and IntegrateData
    // saves a PCA embedding of the SCT integrated data in anndata obsm
    // see IntegrateData 
    // also saves the SCT integrated h5ad for evaluation

    REF_PROCESS_SEURAT(ref_paths_adata)
    REF_PROCESS_SEURAT.out.ref_paths_seurat.set { ref_paths_seurat }
    REF_PROCESS_SEURAT.out.sct_scvi_integrated_h5ads.set { sct_scvi_integrated_h5ads }


    // ---------------evaluate integration -----------------------

    EVALUATE_INTEGRATION(sct_scvi_integrated_h5ads) 


    // ----------------------- run k-fold CV on all references

    CV_SEURAT(ref_paths_seurat, ref_keys)
    CV_SEURAT.out.pred_scores_channel.set { pred_scores_channel }

    CV_SCVI(ref_paths_adata, ref_keys)
    CV_SCVI.out.probs_channel.set { probs_channel }

    probs_channel.map { ['scvi', it[0], it[1], it[2]] }
    .set { scvi_probs_channel }

    pred_scores_channel.map { ['seurat', it[0], it[1], it[2]] }
    .set { seurat_pred_scores_channel }

    scvi_probs_channel.concat( seurat_pred_scores_channel )
    .set { combined_probs_channel }

    CLASSIFY_KFOLD_CV(ref_keys, combined_probs_channel, ref_region_mapping)

    CLASSIFY_KFOLD_CV.out.summary_score_channel
    .flatten()
    .toList()
    .set { ct_metrics }


    ref_paths_adata.filter { it.getName().contains("whole_cortex") }
    .set { whole_cortex_adata }

    ref_paths_seurat.filter { it.getName().contains("whole_cortex") }
    .set { whole_cortex_seurat }

    // LOOCV SEURAT ---------------------------------------------

    LOOCV_SEURAT(whole_cortex_seurat, ref_keys)
    LOOCV_SEURAT.out.loocv_obs_seurat
    .set { loocv_obs_seurat }
    LOOCV_SEURAT.out.loocv_probs_seurat
    .set { loocv_probs_seurat }

loocv_probs_seurat
  .map { list -> list.sort { it.getName() } }  // ensures stable order
  .flatMap { list ->
      list.collect { file ->
          held_out_name = file.getName().split('_loocv.prediction.scores.tsv')[0]
          [held_out_name, file]
      }
  }
  .set { seurat_loocv_probs }

loocv_obs_seurat
  .map { list -> list.sort { it.getName() } }
  .flatMap { list ->
      list.collect { file ->
          held_out_name = file.getName().split('_loocv.obs.tsv')[0]
          [held_out_name, file]
      }
  }
  .set { seurat_loocv_obs }


    seurat_loocv_results = seurat_loocv_probs
    .join(seurat_loocv_obs, by: 0)
    .map { held_out_name, probs, obs -> 
        def method = 'seurat'
        [held_out_name, method, probs, obs] 
    }
    // SCVI LOOCV ---------------------------------------------

    LOOCV_SCVI(whole_cortex_adata, ref_keys)
    LOOCV_SCVI.out.loocv_obs_scvi
    .set { loocv_obs_scvi }
    LOOCV_SCVI.out.loocv_probs_scvi
    .set { loocv_probs_scvi }

    loocv_probs_scvi.flatten()
    .map { file ->
        held_out_name = file.getName().split('_loocv.prob.df.tsv')[0]
        [held_out_name, file]
    }.set { scvi_loocv_probs }

    loocv_obs_scvi.flatten()
    .map { file ->
        held_out_name = file.getName().split('_loocv.obs.tsv')[0]
        [held_out_name, file]
    }.set { scvi_loocv_obs }

    scvi_loocv_results = scvi_loocv_probs
    .join(scvi_loocv_obs, by: 0)
    .map { held_out_name, probs, obs -> 
        def method = 'scvi'
        [held_out_name, method, probs, obs] 
    }

    scvi_loocv_results.view()
    // combine results
    scvi_loocv_results.concat(seurat_loocv_results)
    .set { combined_loocv_results }
   // CLASSIFY_LOOCV(combined_loocv_results, ref_keys, ref_region_mapping)


      //// Pairwise prediction processes (requires map_valid_labels) -----------------------

    //ref_paths_seurat.combine(ref_paths_seurat)
    //.filter{ pair -> pair[0] != pair[1] }
        //.set { seurat_pairs }

    //ref_paths_adata.combine(ref_paths_adata)
        //.filter{ pair -> pair[0] != pair[1]}
        //.set { scvi_pairs }

    //SEURAT_PAIR_PREDICT(seurat_pairs, ref_keys)
    //SCVI_PAIR_PREDICT(scvi_pairs, ref_keys)

    //predicted_probs_seurat = SEURAT_PAIR_PREDICT.out.predicted_probs_seurat
    //predicted_probs_scvi = SCVI_PAIR_PREDICT.out.predicted_probs_scvi

    //// concat
    //predicted_probs_scvi.concat(predicted_probs_seurat)
    //.set { combined_pairwise_probs_channel }


    //CLASSIFY_PAIRWISE(combined_pairwise_probs_channel, ref_keys, ref_region_mapping)
}

workflow.onComplete {
    println "Successfully completed"
    println ( workflow.success ? 
    """
    ===============================================================================
    Pipeline execution summary
    -------------------------------------------------------------------------------

    Run as      : ${workflow.commandLine}
    Started at  : ${workflow.start}
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    Config files: ${workflow.configFiles}
    Exit status : ${workflow.exitStatus}
    Output directory : ${params.outdir}

    --------------------------------------------------------------------------------
    ================================================================================
    """.stripIndent() : """
    Failed: ${workflow.errorReport}
    exit status : ${workflow.exitStatus}
    """.stripIndent()
    )
}

workflow.onError = {
    println "Error: something went wrong, check the pipeline log at '.nextflow.log"
}