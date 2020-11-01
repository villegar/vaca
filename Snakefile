####### Libraries #######
from utils import extractFilenames, findLibraries, loadGenome, verifyGenome, which
from utils import expand_list as el

####### Global variables #######
EXTENSION = config["reads"]["extension"]
PREFIX = config["reads"]["prefix"]
READS = config["reads"]["path"]
PAIRED_END = [True if config["reads"]["end_type"] == "pe" else False][0]
try:
    BBDUK_OPTIONS = config["bbduk"]["options"]
except:
    raise ValueError("bbduk > options not found in the configuration file")
if PAIRED_END:
    FORWARD_READ_ID = [config["reads"]["forward_read_id"]]
    REVERSE_READ_ID = [config["reads"]["reverse_read_id"]]
    ENDS = [FORWARD_READ_ID,REVERSE_READ_ID]
    SUFFIX = "_" + FORWARD_READ_ID[0] + "." + EXTENSION
else:
    ENDS = []
    FORWARD_READ_ID = []
    REVERSE_READ_ID = []
    SUFFIX = "." + EXTENSION

LIBS = findLibraries(READS,PREFIX,SUFFIX)

###### Multithread configuration #####
CPUS_FASTQC = 4
CPUS_TRIMMING = 5
CPUS_HISAT2_INDEX = 40
CPUS_ALIGNMENT = 10
CPUS_READCOUNTS = 20

ADAPTER = which("bbduk")

####### Output directories #######
REF_GENOME = "GENOME/"
RAW_FASTQC = "1.QC.RAW/"
TRIMMED_READS = "2.TRIMMED/"
TRIMMED_READS_FASTQC = "3.QC.TRIMMED/"
ALIGNMENT = "4.ALIGNMENT/"
ALIGNMENT_QC = "5.QC.ALIGNMENT/"
COUNTS = "6.COUNTS/"
RMD = "7.RMD/"
REPORTS = "999.REPORTS/"

####### Reference datasets #######
FA,GTF = loadGenome(config["genome"])
GENOME_FILENAMES = {"FA":FA,"GTF":GTF}
verifyGenome(config["genome"],REF_GENOME + FA, REF_GENOME + GTF)

RAW_ENDS = [""]
if PAIRED_END:
    RAW_ENDS = el(["_"],ENDS)

####### Rules #######
rule all:
    input:
        expand(RAW_FASTQC + "{raw_reads}{raw_ends}_fastqc.{format}",
            raw_reads = LIBS, raw_ends = RAW_ENDS, format = ["html","zip"])
        # expand(RAW_FASTQC + "{raw_reads}{raw_ends}_fastqc.{format}",
        #     raw_reads = LIBS, raw_ends = [1, 2], format = ["html","zip"])
    output:
        expand(REPORTS + "Report_{step}.html", step = ["FastQC_Raw"])
    params:
        logs 	= directory(LOGS),
        reports	= directory(REPORTS)
    run:
        shell("multiqc -f -o {params.reports} -n Report_FastQC_Raw.html -d " + RAW_FASTQC)

rule fastqc_raw:
    input:
        reads = READS_PATH + "{raw_reads}{raw_ends}." + EXTENSION
    output:
        html = RAW_FASTQC + "{raw_reads}{raw_ends}_fastqc.html",
        zip  = RAW_FASTQC + "{raw_reads}{raw_ends}_fastqc.zip"
    message:
        "FastQC on raw data"
    log:
        RAW_FASTQC + "{raw_reads}{raw_ends}.log"
    threads:
        CPUS_FASTQC
    shell:
        "fastqc -o " + RAW_FASTQC + " -t {threads} {input.reads} 2> {log}"
