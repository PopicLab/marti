# ```marti```: a lightweight framework for long-read cDNA artifact analysis

##### Table of Contents
[Overview](#overview)  
[Installation](#install)  
[User Guide](#guide)

<a name="overview"></a>
### Overview 

```marti``` is a lightweight framework for the classification and quantification of artefactual long-read cDNA 
constructs. At a high level, ```marti``` classifies each read based on the presence, absence, and 
location of predefined target oligos (e.g., PCR and sequencing adapters) within the read sequence. 
It takes as input (1) a uBAM file with long cDNA reads and (2) a YAML configuration file with search 
parameters and target sequences of interest. It outputs (1) a uBAM file with annotated cDNA reads 
(per-read annotations include the artifact category, the coordinates of all target oligos and the native 
cDNA sequence, and the structural representation of the read), (2) TSV reports with counts for each 
artifact category and read structure, respectively, and (3) log files with examples of each artifact category, 
which include detailed parsing information and annotation of target oligos in these reads. Optionally, 
```marti``` can be configured to output an additional uBAM file containing only proper (non-artefactual) trimmed reads.

<a name="install"></a>
### Installation

1. Clone the repository: ```$> git clone git@github.com:PopicLab/marti.git```  
2. Create a new conda environment based on your platform:  
  ```$> conda create --prefix envs/marti --file envs/conda-linux-64.lock ``` (Linux)  
  ```$> conda create --prefix envs/marti --file envs/conda-osx-64.lock``` (MacOS)  
3. Activate the environment: ```$> conda activate envs/marti```  
4. Build ```marti```: ```$> build.sh```

<a name="guide"></a>
### User guide

#### Quick start

To run ```marti```: ```$> ./build/bin/marti <path/to/config.yaml>```

#### Recommended workflow
1. Create a new directory for your experiment
2. Create a new YAML config file in this directory (see the provided templates)
3. Populate the YAML config file with the parameters specific to this experiment
4. Run ```marti``` with this YAML file as input. ```marti``` will automatically 
create auxiliary folders and files with results in the directory where the config YAML file is located

#### Inputs

A template config file is provided in ```docs/template.yaml```.

The key ```marti``` parameters include:
* The two PCR adapter sequences (referred to as ```adapterA``` and ```adapterB```)
* The ```TSO_adapters``` and ```RT_adapters``` lists, which assign PCR adapters to the TSO (template-switch oligo) and 
the RT (reverse transcription) primer sites, respectively
* The SLS (switch leader sequence) of the TSO adapter(s) (e.g. ```adapterA_SLS```)
* The length of the cDNA terminal regions (boundaries) within which the TSO, RT primer, and the polyAs are 
expected to be located in a proper read
* The minimum length of the polyA tail  
* The maximum expected sequencing error rate

Full list of parameters:

* ```input_bam``` path to the input uBAM/BAM file
* ```adapterA``` PCR adapter sequence 
* ```adapterB``` PCR adapter sequence 
* ```adapterA_SLS``` SLS sequence for the adapterA TSO
* ```adapterB_SLS``` SLS sequence for the adapterB TSO
* ```adapterA_SPLIT_SLS``` SLS sequence for the adapterA TSO when split from the adapter
* ```adapterB_SPLIT_SLS``` SLS sequence for the adapterB TSO when split from the adapter
* ```TSO_adapters``` list of one or two PCR adapters in the TSO (e.g. ```[adapterA, adapterB]```)
* ```RT_adapters``` list of one or two PCR adapters in the RT (e.g. ```[adapterB]```)
* ```terminal_adapter_search_buffer [default=100]``` distance from each end of the read, within which the adapter sequence is expected to be found (fully contained)
* ```terminal_polyA_search_buffer [default=150]``` distance from each end of the read, within which the polyA/polyT sequence is expected to start
* ```min_polyA_match [default=20]``` minimum length of a polyA match
* ```max_err_rate [default=0.1]``` determines the maximum distance allowed for a target match
* ```max_err_rate_polya [default=0.1]``` determines the maximum distance allowed for a polyA match
* ```min_rq [default=1.0]``` minimum required read quality score for classification (PacBio only)
* ```max_reads_to_process [default=all]``` maximum number of reads to process from the BAM file
* ```output_proper_trimmed [default=false]``` output a separate BAM file with proper reads only (trimmed, adjusted to forward strand) (uBAMs only!)
* ```n_threads [default=1]``` number of threads to use for classification
* ```n_max_reads_in_mem [default=100000]``` maximum number of reads of load into memory at a time

#### Artifact classes 

* ```Proper``` properly constructed non-artefactual reads
* ```TsoTso``` RT artifact caused by the internal priming of the SLS; results in the presence of a TSO at both ends of the read
* ```RtRt``` RT artifact caused by the internal priming of the oligo dT; results in the presence of the RT PCR adapter and a polyT at both ends of the read
* ```InternalPrimingRT``` PCR artifact caused by the internal priming of the RT PCR adapter; results in a missing polyA
* ```InternalPrimingTSO``` PCR artifact caused by the internal priming of the TSO PCR adapter; results in a missing SLS
* ```OnlyPolyA``` artifact class assigned when only As or Ts are found inside the read sequence, likely a result of RNA degradation
* ```MissingAdapter```  sequencing artifact; assigned when the read is missing at least one PCR adapter at either end
* ```InternalAdapter``` sequencing artifact; assigned when an extra adapter is found inside the read (in addition to the two adapters already found at each end)
* ```RetainedSMRTBell``` (PacBio reads only) sequencing artifact; assigned when the PacBio SMRTbell adapter is found inside the read
* ```Unk``` artifact of unknown type; the read structure that does not correspond to a proper read nor any other known artifact class
* ```LowRQ``` reads with low read quality (rq) score
* ```TooShort``` reads that are too short to classify

#### Outputs

##### Classified annotated BAM file (```.classified.bam```)

BAM file containing all the input reads annotated with additional ```marti``` tags.

###### Marti tags

* ```pr [int]``` flag indicating whether the read is proper (1) or not (0)
* ```sf [int]``` (proper reads only) flag indicating whether the original read was on the forward (0) or reverse strand (1)
* ```cd [array]``` (proper reads only) stores the 0-based start and end (inclusive) coordinates of the cDNA segment
* ```lb [string]``` class profile: comma-separated list of one or multiple assigned artifact classes (listed in alphabetical order)
* ```st[string]``` structure profile: ellipsis-separated ordered list of key sequence types found within the read
* ```ch [string]``` comma-separated list of key sequence types found within the read in the format ```type:start:end:lev```, where
  * ```type``` denotes the sequence type (e.g. ```adapterA```/```polyA```/```cDNA```)
  * ```start```, ```end``` (inclusive) are 0-based positions of each subsequence in the read
  * ```lev```: is the Levenshtein distance between the target sequence and its match in the read
* ```th[string]``` comma-separated list of all targets found within the read

##### Proper cDNA uBAM file (if configured) (```.proper.cdna.bam```)
uBAM file containing only proper reads with trimmed adapters (the polyA is kept).
Note: reads that were found in the reverse complemented configuration are flipped to the forward strand.
Note: this functionality is available only for unmapped BAM inputs.

##### Reports (```reports/``` folder)
* ```class_counts.tsv``` TSV file providing the number of reads with each artifact class profile
* ```structure_counts.tsv``` TSV file providing the number of reads with each artifact class and structure profile
* ```<artifact_class_profile>.txt``` (e.g. ```TsoTso.txt```) examples of reads from each artifact category with detailed 
parsing information and target annotation
* ```reports/main.log```: log file with configuration and runtime information
