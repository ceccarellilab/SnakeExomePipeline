#!/bin/bash
rm html_pipeline.zip
rm -r ./report
mkdir report

source pipeline/bin/activate
jupyter nbconvert --no-input --execute --ExecutePreprocessor.timeout=-1 --output-dir='./report' --to html pipeline.ipynb
zip -r html_pipeline.zip report/