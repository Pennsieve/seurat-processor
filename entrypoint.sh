#!/bin/bash
set -e

if [ -n "$AWS_LAMBDA_RUNTIME_API" ]; then
    echo "Running in AWS Lambda mode"
    exec Rscript /app/processor/main.R
else
    echo "Running in local/container mode"
    exec Rscript /app/processor/main.R
fi
