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
