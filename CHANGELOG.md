### 1.5.0

Update Kraken2 and BUSCO databases to the latest versions

### 1.4.0

Add Fastq MD5 to data model

### 1.3.4

Fix authenticaiton

### 1.3.3

Add API Authentication

### 1.3.2

Fix parsing of metadata code fields, now ensures everything is strings
 
### 1.3.1

Fix value parsing for q30, ref_coverage and seq_depth to allow None in case the samples are from third party

### 1.3.0

Report version to API

### 1.2.2

Fix isolate_id field in sequence db

### 1.2.0

Deployement of dependencies in the conda env path
Fix Json response parsing

### 1.0.0

Rework for integration with Geuebt-API and Metadata V3

### 0.2.7

Allows to bypass Modelchoice and provide another validation model when calling validate_metata.main() directly (for internal use only)

### 0.2.6

Add support for 25.11.2024 date format form excel
Relaxed the fraction_majority_genus cutoff to 0.9 for all species to accomodate isolate with many plasmids

### 0.2.5

Fix missing entries in isolate sheet

### 0.2.4

Fix isoalte sheet creation with new metadata definition

### 0.2.3

Quickfix asssembly qc merging

### 0.2.2

Quickfix invalid field type seq_depth from int to float

### 0.2.1

Remove deprecated param

### 0.2.0

- Implements GEÜBT Metadata v2 (03-2024)
- Changed input validation from a JSON Schema to a more flexible Pydantic Model.

### 0.1.7

Removedhard QC threshold oncontigs > 1kb and GC percent. This are optionnal parameters in the ASU and should be checked in the assembly.

### 0.1.6

Remove unnescessary `isolate_id` field check in assembly QC

### 0.1.5

Relaxed regex for sample naming, subsample number now accepts tow to many digits (eg. 2023-0001254-01 or 2023-0001254-5236)

### 0.1.4

Added 'ringtrial' sample type

### 0.1.3

Hotfix for fasta name in the isolate_sheets. Fastas are intrnally renamed to match the `isolate_id` field in order to ensure unique naming in the filesystem. This is now correctlyrefelcted in the isolate sheets.

### 0.1.2

Checking config file with schema

### 0.1.1

Adding unit tests

### 0.1.0

Workflow is working.
So far only JSON outputs and no documentation.
