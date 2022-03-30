from .in_gerp import reader_Gerp, SCHEMA_GERP
from .in_spliceai import reader_SpliceAI, SCHEMA_SPLICE_AI
from .in_dbnsfp4 import reader_dbNSFP4, SCHEMA_DBNSFP_4
from .in_dbsnp import reader_SNP, SCHEMA_dbSNP
from .in_gnomad211 import reader_GNOMAD211, SCHEMA_GNOMAD_2_1
from .in_gtf import reader_GTF, SCHEMA_GTF
#========================================
sIngestModes = {
    "Gerp":     [SCHEMA_GERP,       reader_Gerp],
    "gnomAD":   [SCHEMA_GNOMAD_2_1, reader_GNOMAD211],
    "dbNSFP":   [SCHEMA_DBNSFP_4,   reader_dbNSFP4],
    "dbSNP":    [SCHEMA_dbSNP,      reader_SNP],
    "SpliceAI": [SCHEMA_SPLICE_AI,  reader_SpliceAI],
    "gtf":      [SCHEMA_GTF,        reader_GTF]
}

#========================================
def getIngestModeSetup(mode):
    global sIngestModes
    return sIngestModes[mode]
