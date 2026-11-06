execute_process(
    COMMAND "${NTRIAGE}"
        --bed "${SOURCE_DIR}/tests/tiny.bed"
        --fasta "${SOURCE_DIR}/tests/tiny.fa"
        --refs "${SOURCE_DIR}/tests/tiny_refs.fa"
        --colors "${SOURCE_DIR}/tests/tiny_colors.tsv"
        --out "${BINARY_DIR}/tiny_result"
        -k 5 --window 6 --bins 8 -c 2 --min-hits 1
        --min-period 8 --max-period 16 --period-target 12 --boundary-span 8
        --write-fasta
    RESULT_VARIABLE result
    OUTPUT_VARIABLE stdout
    ERROR_VARIABLE stderr
)
if(NOT result EQUAL 0)
    message(FATAL_ERROR "N-TRIAGE integration run failed: ${stderr}")
endif()

file(STRINGS "${BINARY_DIR}/tiny_result.clusters.tsv" cluster_lines)
file(STRINGS "${BINARY_DIR}/tiny_result.features.tsv" feature_lines)
file(STRINGS "${BINARY_DIR}/tiny_result.annotations.bed" annotation_lines)
file(STRINGS "${BINARY_DIR}/tiny_result.acres.fasta" fasta_headers REGEX "^>")
file(STRINGS "${BINARY_DIR}/tiny_result.acres.bed" acre_lines)
file(STRINGS "${BINARY_DIR}/tiny_result.arrays.tsv" array_lines)
file(STRINGS "${BINARY_DIR}/tiny_result.colors.tsv" color_lines)
list(LENGTH cluster_lines cluster_count)
list(LENGTH feature_lines feature_count)
list(LENGTH annotation_lines annotation_count)
list(LENGTH fasta_headers fasta_count)
list(LENGTH acre_lines acre_count)
list(LENGTH array_lines array_count)
list(LENGTH color_lines color_count)
if(NOT cluster_count EQUAL 3 OR NOT feature_count EQUAL 3 OR
   annotation_count LESS 2 OR NOT fasta_count EQUAL 2 OR
   NOT acre_count EQUAL 2 OR NOT array_count EQUAL 3 OR NOT color_count EQUAL 4)
    message(FATAL_ERROR "Unexpected integration output sizes")
endif()

file(READ "${BINARY_DIR}/tiny_result.annotations.bed" annotations)
if(NOT annotations MATCHES "left_component[^\n]*1,2,3" OR
   NOT annotations MATCHES "right_component[^\n]*4,5,6")
    message(FATAL_ERROR "RGB overrides were not written to BED column 9")
endif()

execute_process(
    COMMAND "${NTRIAGE}"
        --bed "${SOURCE_DIR}/tests/deep.bed"
        --fasta "${SOURCE_DIR}/tests/deep.fa"
        --out "${BINARY_DIR}/deep_result"
        -k 5 --window 8 --bins 16 -c 2 --min-hits 1
        --min-period 16 --max-period 16 --period-target 16 --boundary-span 8
        --deep-subfamilies --sketch-size 64 --switch-flank 3 --switch-radius 3
        --pangenome --threads 2
    RESULT_VARIABLE deep_result
    OUTPUT_VARIABLE deep_stdout
    ERROR_VARIABLE deep_stderr
)
if(NOT deep_result EQUAL 0)
    message(FATAL_ERROR "N-TRIAGE deep-subfamily integration run failed: ${deep_stderr}")
endif()

file(STRINGS "${BINARY_DIR}/deep_result.switches.tsv" switch_lines)
file(STRINGS "${BINARY_DIR}/deep_result.switch_boundaries.tsv" boundary_lines)
file(STRINGS "${BINARY_DIR}/deep_result.switch_subfamilies.bed" subfamily_lines)
file(STRINGS "${BINARY_DIR}/deep_result.switch_subfamily_summary.tsv" summary_lines)
file(STRINGS "${BINARY_DIR}/deep_result.switches.bed" switch_bed_lines)
file(STRINGS "${BINARY_DIR}/deep_result.switch_arrays.bed" switch_array_lines)
file(STRINGS "${BINARY_DIR}/deep_result.acres.bed" deep_acre_lines)
file(STRINGS "${BINARY_DIR}/deep_result.arrays.tsv" deep_array_lines)
list(LENGTH switch_lines switch_count)
list(LENGTH boundary_lines boundary_count)
list(LENGTH subfamily_lines subfamily_count)
list(LENGTH summary_lines summary_count)
list(LENGTH switch_bed_lines switch_bed_count)
list(LENGTH switch_array_lines switch_array_count)
list(LENGTH deep_acre_lines deep_acre_count)
list(LENGTH deep_array_lines deep_array_count)
if(NOT switch_count EQUAL 2 OR NOT boundary_count EQUAL 2 OR
   NOT subfamily_count EQUAL 7 OR NOT summary_count EQUAL 3 OR
   NOT switch_bed_count EQUAL 2 OR NOT switch_array_count EQUAL 2 OR
   NOT deep_acre_count EQUAL 6 OR NOT deep_array_count EQUAL 2)
    message(FATAL_ERROR "Unexpected deep-subfamily output sizes: ${deep_stderr}")
endif()
file(READ "${BINARY_DIR}/deep_result.switches.tsv" switches)
file(READ "${BINARY_DIR}/deep_result.switch_subfamilies.bed" subfamilies)
file(READ "${BINARY_DIR}/deep_result.switches.bed" switches_bed)
file(READ "${BINARY_DIR}/deep_result.switch_arrays.bed" switch_arrays_bed)
if(NOT switches MATCHES "ARRAY_000777[^\n]*48" OR
   NOT switches_bed MATCHES "switchchr[	]48[	]49[	]SWITCH_000001" OR
   NOT switch_arrays_bed MATCHES "switchchr[	]0[	]96[	]ARRAY_000777" OR
   NOT subfamilies MATCHES "switch_subfamily_1" OR
   NOT subfamilies MATCHES "switch_subfamily_2" OR
   NOT subfamilies MATCHES "\tleft\t" OR
   NOT subfamilies MATCHES "\tright\t")
    message(FATAL_ERROR "Deep-subfamily calls did not recover the synthetic switch")
endif()
