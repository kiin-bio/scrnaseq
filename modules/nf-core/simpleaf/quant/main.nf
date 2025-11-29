process SIMPLEAF_QUANT {
    tag "${meta.id ?: meta4.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/simpleaf:0.19.4--ha6fb395_0':
        'biocontainers/simpleaf:0.19.4--ha6fb395_0' }"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    tuple val(meta), val(chemistry), val(reads)                         // chemistry and reads (val for lazy download)
    tuple val(meta2), path(index), path(txp2gene)                       // index and t2g mapping
    tuple val(meta3), val(cell_filter), val(number_cb), path(cb_list)   // cell filtering strategy
    val resolution                                                      // UMI resolution
    tuple val(meta4), path(map_dir)                                     // mapping results

    output:
    tuple val(meta), path("${prefix}/af_map")       , emit: map, optional: true // missing if map_dir is provided
    tuple val(meta), path("${prefix}/af_quant")     , emit: quant
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"

    // The second required input is a cell filtering strategy.
    cf_option = cellFilteringArgs(cell_filter, number_cb, cb_list)

    meta = map_dir ? meta4 : meta2 + meta3 + meta
    meta = meta + [ "filtered": cell_filter != "unfiltered-pl" ]

    // Prepare reads for download if provided
    def reads_list = reads instanceof List ? reads : (reads ? [reads] : [])
    def reads_str  = reads_list.collect { "\"${it}\"" }.join(' ')
    def has_reads  = reads_list.size() > 0
    def t2g_arg    = txp2gene ? "--t2g-map ${txp2gene}" : ""

    """
    # Helper function to download files based on URL type
    download_file() {
        local src="\$1"
        local dest="\$2"
        
        if [[ "\$src" == s3://* ]]; then
            aws s3 cp "\$src" "\$dest" --only-show-errors
        elif [[ "\$src" == az://* ]] || [[ "\$src" == https://*.blob.core.windows.net/* ]]; then
            azcopy copy "\$src" "\$dest"
        elif [[ "\$src" == gs://* ]]; then
            gsutil cp "\$src" "\$dest"
        elif [[ "\$src" == http://* ]] || [[ "\$src" == https://* ]]; then
            curl -sL "\$src" -o "\$dest"
        elif [[ -f "\$src" ]]; then
            ln -s "\$src" "\$dest"
        else
            echo "ERROR: Cannot access file: \$src" >&2
            exit 1
        fi
    }

    # Build mapping arguments
    mapping_args=""
    if [[ -n "${map_dir}" && "${map_dir}" != "[]" ]]; then
        mapping_args="--map-dir ${map_dir}"
    elif [[ ${has_reads} == true ]]; then
        # Download reads to local directory preserving order
        mkdir -p fastq_dir
        local_reads=()
        for read_path in ${reads_str}; do
            filename=\$(basename "\$read_path")
            download_file "\$read_path" "fastq_dir/\$filename"
            local_reads+=("fastq_dir/\$filename")
        done

        # Separate forward from reverse pairs (reads come as R1,R2,R1,R2,...)
        forward_reads=""
        reverse_reads=""
        for ((i=0; i<\${#local_reads[@]}; i+=2)); do
            if [ -n "\$forward_reads" ]; then
                forward_reads="\$forward_reads,\${local_reads[\$i]}"
                reverse_reads="\$reverse_reads,\${local_reads[\$((i+1))]}"
            else
                forward_reads="\${local_reads[\$i]}"
                reverse_reads="\${local_reads[\$((i+1))]}"
            fi
        done

        mapping_args="${t2g_arg} --chemistry ${chemistry} --index ${index} --reads1 \$forward_reads --reads2 \$reverse_reads"
    else
        echo "ERROR: Neither reads nor map_dir provided" >&2
        exit 1
    fi

    # export required var
    export ALEVIN_FRY_HOME=.

    # prep simpleaf
    simpleaf set-paths

    # run simpleaf quant
    simpleaf quant \\
        \$mapping_args \\
        --resolution ${resolution} \\
        --output ${prefix} \\
        --threads ${task.cpus} \\
        --anndata-out \\
        ${cf_option} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alevin-fry: \$(alevin-fry --version | sed -e "s/alevin-fry //g")
        piscem: \$(piscem --version | sed -e "s/piscem //g")
        salmon: \$(salmon --version | sed -e "s/salmon //g")
        simpleaf: \$(simpleaf --version | sed -e "s/simpleaf //g")
    END_VERSIONS
    """

    stub:
    prefix    = task.ext.prefix ?: "${meta.id}"
    """
    export ALEVIN_FRY_HOME=.

    mkdir -p ${prefix}/af_map
    mkdir -p ${prefix}/af_quant/alevin

    touch ${prefix}/af_map/map.rad
    touch ${prefix}/af_map/unmapped_bc_count.bin
    touch ${prefix}/af_quant/alevin/quants_mat_rows.txt
    touch ${prefix}/af_quant/map.collated.rad
    touch ${prefix}/af_quant/permit_freq.bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alevin-fry: \$(alevin-fry --version | sed -e "s/alevin-fry //g")
        piscem: \$(piscem --version | sed -e "s/piscem //g")
        salmon: \$(salmon --version | sed -e "s/salmon //g")
        simpleaf: \$(simpleaf --version | sed -e "s/simpleaf //g")
    END_VERSIONS
    """
}

// We have mutual exclusive options for permit list generation.
// 1. 'k' (knee), which is a flag for the knee method and any value provided will be ignored;
// 2. 'f' (forced-cells), which takes an integer indicating the exact number of cells to recover;
// 3. 'e' (expect-cells), which takes an integer indicating the expected number of cells to recover;
// 4. 'x' (explicit-pl), which takes a string indicating the path to a valid permit list;
// 5. 'u' (unfiltered-pl), which takes an empty string (if `chemistry` is defined as "10xv2" or "10xv3"), or a string indicating the path to a valid white list file.
// The difference between (4) and (5) is that (4) contains the exact permit list to filter the observed barcodes, while (5) will use the white list to generate a permit list via barcode correction.

// We have two ways to take these options. `-u` is implied by the presence of the input `whitelist` channel. The options can also be passed as arguments to ext.args. Therefore, we must check two things:
// 1. if there is at least one of the options in the args list, and
// 2. if none of the four options are in the args list, there must be a non-empty whitelist channel.

def cellFilteringArgs(cell_filter_method, number_cb, cb_list) {
    def pl_options = ["knee", "forced-cells", "explicit-pl", "expect-cells", "unfiltered-pl"]

    // try catch unintentional underscore in method name
    def method = cell_filter_method.replaceAll('_','-')

    def number = number_cb
    if (!method) {
        error "No cell filtering method was provided; cannot proceed."
    } else if (! method in pl_options) {
        error "Invalid cell filtering method, '${method}', was provided; cannot proceed. possible options are ${pl_options.join(',')}."
    }

    if (method == "unfiltered-pl") {
        return "--${method} ${cb_list}"
    } else if (method == "explicit-pl") {
        return "--${method} ${cb_list}"
    } else if (method == "knee") {
        return "--${method}"
    } else {
        if (!number) {
            error "Could not find the corresponding 'number' field for the cell filtering method '${method}'; please use the following format: [method:'${method}',number:3000]."
        }
        return "--${method} ${number}"
    }
}

def mappingArgs(chemistry, reads, index, txp2gene, map_dir) {
    if ( map_dir ) {
        if (reads) {
            error "Found both reads and map_dir. Please provide only one of the two."
        }
        return "--map-dir ${map_dir}"
    } else {
        if (!reads) {
            error "Missing read files; could not proceed."
        }
        if (!index) {
            error "Missing index files; could not proceed."
        }
        if (!chemistry) {
            error "Missing chemistry; could not proceed."
        }

        def (forward, reverse) = reads.collate(2).transpose()

        def t2g = txp2gene ? "--t2g-map ${txp2gene}" : ""
        def mapping_args = """${t2g} \\
        --chemistry ${chemistry} \\
        --index ${index} \\
        --reads1 ${forward.join( "," )} \\
        --reads2 ${reverse.join( "," )}"""
        return mapping_args
    }
}
