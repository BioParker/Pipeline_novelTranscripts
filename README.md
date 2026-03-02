# Pipeline_novelTranscripts

Nextflow pipeline for identifying transcripts that contain novel cassette exons from IsoQuant outputs

## Parameters

All essential parameters can be set in config file.

### Key files

  *	models: Path to transcript models gtf output from IsoQuant
  * model_tpms: Path to model transcript tpms output from IsoQuant
  * refgtf: Path to GENCODE reference gtf

### Paths

  * novelTranscripts_conda: Path to conda environment containing necessary packages and dependecies. See novelTranscripts.yaml and required R packages below.


## **Required R packages**

	•	tidyverse (2.0.0)
	•	rtracklayer (1.70.1)
	•	optparse (1.7.5)
