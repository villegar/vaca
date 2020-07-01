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

####### Reference datasets #######
FA,GTF = loadGenome(config["genome"])
GENOME_FILENAMES = {"FA":FA,"GTF":GTF}
verifyGenome(config["genome"],REF_GENOME + FA, REF_GENOME + GTF)

RAW_ENDS = [""]
if PAIRED_END:
    RAW_ENDS = el(["_"],ENDS)
