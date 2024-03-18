# Guide for users

## Installation

### Conda

Install conda from any distribution, i.e. miniconda.
You can follow the setup guide from the [Bioconda team](https://bioconda.github.io/).

We advise installing the mamaba solver in the base environement to speed up
environments creation.

```bash
conda install mamba -n base -c conda-forge
```

### Run environment

The running environement simply required a recent python version (>= 3.9) and snakemake.
If you followed the steps above just run:

```bash
mamba create -n snakemake snakemake
```

### Install module and databases

Download the [latest realease](https://github.com/NRW-GEUBT/geuebt-validate/releases/latest)
and unpack it.

If you're feeling brave, clone the repository form Github:

```bash
git clone https://github.com/NRW-GEUBT/geuebt-validate
```

All software and databases dependencies will be installed during the first run.
It may take a little time so be patient!

### Manually install databases

If for any reason you want to manually install database it is possible to provide
a custom path in the configuration.

Following databases are required:

#### Kraken2

Preformatted kraken2 databases are available from:
[https://benlangmead.github.io/aws-indexes/k2](https://benlangmead.github.io/aws-indexes/k2)

The workflow uses the Standard-8 (2023-03-14) as default.

#### Busco

The workflow uses the bacteria_odb10 for BUSCO.
If you wish to install an alternatve database, see the [user's guide](https://busco.ezlab.org/busco_userguide.html#offline).

#### Taxonomy

TaxidTools relies on the local availability of the Taxdump files which can be
obtained from teh [NCBI servers](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/).

Using other taxonomy definitions is not supported.

## Configuration

The configuaration can be defined in two ways:

- either edit and locally save the `config/config.yaml` files and provide its path
  to the snakemake command with the `--configfile` argument

- or provide the parameters directly to the snakemake command with
  `--config <ARGNAME>=<VALUE>`

### User defined parameters

Following arguments must be provided for each run:

| Parameter | Type | Description |
| --- | --- | --- |
| `workdir` | path-like string | Path to the ouptut directory |
| `metadata` | path-like string | Path to the metadata file in TSV format |
| `fasta_dir` | path-like string | Path to the folder containing assembly files in fasta format |

### Optional parameters

Following parameters are optional and will revert to defaults if not set:

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `max_threads_per_job` | integer | 1 | Max number of threads assigned to a single job |
| `kraken2_db` | path-like string | Default installation in `~/.nrw-geuebt/` | Path to the kraken2 database folder |
| `busco_db` | path-like string | Default installation in `~/.nrw-geuebt/` | Path to the busco database folder |
| `taxdump` | path-like string | Default installation in `~/.nrw-geuebt/` | Path to the taxdump database folder |
| `min_contig_length` | integer | 500 | Minimal contig length (in bp)<br>to be included in the calculation of assembly<br>statistics |

## Usage

The workflow can be started with:

```bash
snakemake --use-conda --conda-prefix <PATH TO CONDA ENVS> --configfile <PATH TO CONFIG> --cores <N>
```

See the [snakemake documentaiton](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for more information.

## Input formats

### Metadata format

The metadata must be provided as a tab-delimited text data, with strict requirements:

| Field name             | Definition                                                                                            | Required? | Allowed Values                                                                                 |
| ---------------------- | ----------------------------------------------------------------------------------------------------- | --------- | ---------------------------------------------------------------------------------------------- |
| isolate_id             | Eindeutiger Bezeichner der Isolat                                                                     | required  | free text                                                                                      |
| sample_id              | Eindeutiger Bezeichner der Probe                                                                      | required  | free text                                                                                      |
| alt_isolate_id         | Alternativer Bezeichnung der Siolat (z.B. Laborrinternenummer)                                        | optional  | free text                                                                                      |
| organism               | Erwarteter Organismus                                                                                 | required  | oneof "Listeria monocytogenes","Salmonella enterica", "Escherichia coli", "Campylobacter spp." |
| isolation_org          | Identifizierungsname des für die Isolation verantwortlichen Labor                                     | required  | oneof "OWL", "RRW", "WFL", "RLD","MEL", "other"                                                |
| sequencing_org         | Identifizierungsname des für die Sequenzierung verantwortlichen Labor                                 | required  | oneof "OWL", "RRW", "WFL", "RLD","MEL", "other"                                                |
| bioinformatics_org     | Identifizierungsname des für die Primär Analyse verantwortlichen Labor                                | required  | oneof "OWL", "RRW", "WFL", "RLD","MEL", "other"                                                |
| third_party_owner      | Ggfs. Name der Dritteigentümer (z.B. bei Humanem Isolaten oder RV)                                    | optional  | free text                                                                                      |
| extraction_method      | Name des Kits oder Methode für die DNA Extraktion                                                     | optional  | free text                                                                                      |
| library_kit            | Name des Kits für die DNA-Library Erstellung                                                          | optional  | free text                                                                                      |
| sequencing_kit         | Name des Kits für die DNA Sequenzierung                                                               | optional  | free text                                                                                      |
| sequencing_instrument  | Model des Gerät für die DNA Sequenzierung                                                             | optional  | free text                                                                                      |
| assembly_method        | Name und Version des verwendeten Software oder Pipeline für die Assembly                              | optional  | free text                                                                                      |
| sample_type            | Art der Probe (z.B. Lebensmittel, Umfeldprobe, veterinärmedizinisches Material, usw.)                 | required  | oneof "Lebensmittel", "Futtermittel", "Umfeld", "Tiergesundheit", "Human", "unknown"           |
| fasta_name             | Eindeutiger Bezeichner der Sequenzdaten (Dateiname)                                                   | required  | free text                                                                                      |
| fasta_md5              | MD5 Checksum der Sequenzdaten                                                                         | required  | MD5                                                                                            |
| collection_date        | Datum der Probenahme                                                                                  | required  | ISOdate                                                                                        |
| customer               | Einsendende KOB                                                                                       | required  | Kodierung (2-5 buchstaben)                                                                     |
| manufacturer           | Hersteller der Probe                                                                                  | required  | Freitext OR  "unknown"                                                                         |
| collection_place       | Entnahmeort der Probe                                                                                 | required  | AVV / ADV Text                                                                                 |
| collection_place_code  | Entnahmeort der Probe                                                                                 | required  | AVV / ADV Code                                                                                 |
| description            | Bezeichnung laut Probeentnahmeschein bzw. bei Tupfern Bezeichnung der Entnahmeortes und ggfs. Tierart | required  | free text                                                                                      |
| manufacturer_type      | Einordnung des Entnahmeortes nach AVV-Katalog                                                         | optional  | AVV Text                                                                                       |
| manufacturer_type_code | Einordnung des Entnahmeortes nach AVV-Katalog                                                         | optional  | AVV Code                                                                                       |
| matrix                 | Einordnung der Probenmatrix nach AVV-Katalog                                                          | optional  | AVV/TSN Text                                                                                   |
| matrix_code            | Einordnung der Probenmatrix nach AVV-Katalog                                                          | optional  | AVV/TSN Code                                                                                   |
| collection_cause       | Grund für die Probenahme nach AVV-Katalog                                                             | optional  | AVV Text oder Freitext                                                                         |
| collection_cause_code  | Grund für die Probenahme nach AVV-Katalog                                                             | optional  | AVV Code                                                                                       |
| lot_number             | Los-Nr./Chargen-Nr. des Herstellers                                                                   | optional  | free text                                                                                      |
| seq_depth              | Durchschnittliche Sequenzierungstiefe nach Assembly                                                   | required  | integer                                                                                        |
| ref_coverage           | Antei Abdeckung des Referenzgenoms                                                                    | required  | float [0-1]                                                                                    |
| q30                    | Anteil der Basen mit PHRED-Score mind. 30                                                             | required  | float [0-1]                                                                                    |

See an example of metadata table in the `config` folder.

### Sequence files

The assemblies must be provided as fasta file.
Wrapped and unwrapped fastas as well as multifastas are allowed.
There are no special requirements for sequence headers.

## Assembly QC testing

Assemblies are tested to ensure the quality of input data.
Acceptance criteria are defined by expected species as follows:

| Parameter | Listeria monocytogenes | Salmonella enterica | Capylobacter spp. | Escherichia coli |
| --- | --- | --- | --- | --- |
| Sequencing depth | 20 - 200 | 30 - 200 | 20 - 200 | 40 - 200 |
| Q30 | >= 0.80 | >= 0.80 | >= 0.80 | >= 0.80 |
| Number of contigs > 1 kb | 300 | 300 | 300 | 500 |
| Assembly size (bp) | 2,700,000 - 3,200,000 | 4,300,000 - 5,200,000 | 1,500,000 - 1,900,000 | 4,500,000 - 5,900,000 |
| %GC | 33.9 - 41.9 | 48.1 - 56.1 | 26.4 - 35.3 | 46.6 - 54.6 |
| Orthologs found (%) | >=95 % | >=95 % | >=80% | >=95 % |
| Duplicated orthologs (%) | <= 5 % | <= 5 % | <= 5 % | <= 5 % |
| Fraction majority Genus | >= 0.95 | >= 0.95 | >= 0.90 | >= 0.90 |
| Majority genus | Listeria | Salmonella | Campylobacter | Escherichia or Shigella |

## Results

Results to be used for the next steps are located in the `staging` folder in the workdir.

### Status report

A JSON report of the QC checks in the form

```json
{
    "2022-0232977-01": {
        "STATUS": "PASS",
        "MESSAGES": []
    },
    "2022-0232977-02": {
        "STATUS": "FAIL",
        "MESSAGES": [
            "Message checksum",
            "Error message assembly"
        ]
    },
    "2022-0232977-03": {
        "STATUS": "FAIL",
        "MESSAGES": ["Error message checksum"]
    },
    "2022-0232977-04": {
        "STATUS": "FAIL",
        "MESSAGES": [
            "Error message metadata 1",
            "Error message metadata 2"
        ]
    }
}
```

### Isolate datasheets

Isolate information are summarized in single JSON files as well as in a merged JSON.
Note that these are generated only for samples satisfying all filters.
The follow the same structure, here for a single entry:

```json
{
    "isolate_id": "2022-0232977-01",
    "sample_id": "2022-0232977",
    "organism": "Salmonella enterica",
    "isolate_name_alt": null,
    "third_party_flag": "true",
    "third_party_owner": "BfR",
    "sample_type": "unknown",
    "fasta_name": "0232977-01.fasta",
    "fasta_md5": "b9621a5adb8e59cbea8034bfc04fca44",
    "isolation": {
        "org_name": "unknown"
    },
    "sequencing": {
        "org_name": "RRW",
        "extraction_method": null,
        "library_method": null,
        "sequencing_instrument": null
    },
    "bioinformatics": {
        "org_name": "RRW",
        "assembly_method": null
    },
    "epidata": {
        "collection_date": "2023-03-28",
        "municipality": "unknown",
        "country": "DE",
        "cause": "Proficiency Test",
        "collected_by": "unknown",
        "manufacturer": null,
        "designation": "unknown",
        "manufacturer_type": "unknown",
        "sample_description": "unkown",
        "lot_number": null
    },
    "mlst": {
        "sequence_type": "34",
        "scheme": "senterica_achtman_2",
        "alleles": {
            "purE": "5",
            "hemD": "12",
            "thrA": "2",
            "sucA": "9",
            "hisD": "9",
            "dnaN": "19",
            "aroC": "10"
        }
    },
    "qc_metrics": {
        "sequencing_depth": 45.6,
        "ref_coverage": 0.99,
        "q30": 0.9,
        "N50": 284333,
        "L50": 5,
        "n_contigs_1kbp": 53,
        "assembly_size": 4960041,
        "GC_perc": 52.13,
        "orthologs_found": 98.4,
        "duplicated_orthologs": 0.0,
        "majority_genus": "Salmonella",
        "fraction_majority_genus": 0.9971668689656122,
        "majority_species": "Salmonella enterica",
        "fraction_majority_species": 0.9971575751364787
    }
}
```

### Assemblies

Assemblies are renamed so that the fastas are in the form `<ISOLATE ID>.fa` and
copied to the `staging/fastas/`folder of the workdir.
