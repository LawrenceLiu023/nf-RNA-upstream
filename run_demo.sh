#!/usr/bin/env bash
nextflow -log './test/logs/test.log' \
    run main.nf \
    -params-file './test/params/params_test.yaml' \
    -with-report './test/reports/report_test.html'

nextflow -log './logs/20250421-12_cfRNA_single_patient.log' \
    run main.nf \
    -params-file './params/params_20250421-12_cfRNA_single_patient.yaml' \
    -with-report './reports/report_20250421-12_cfRNA_single_patient.html'
