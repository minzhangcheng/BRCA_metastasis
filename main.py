from GEO import probeExpr
from GEO import geneExpr
from GEO.expression import probeExpr2geneExpr
from GEO.meta import metaQC
from GEO.meta import metaAnalysis
from GEO.meta import meta2Csv
from GEO.validate import meta2Heatmap
from GEO.validate import meta2Forest

if __name__ == '__main__':

    probeExpr('/home/minzhang/workspace/pyMeta/BRCA_metastasis.json')
    geneExpr('/home/minzhang/workspace/pyMeta/BRCA_metastasis.json', True)
    metaQC('/home/minzhang/workspace/pyMeta/BRCA_metastasis.json')
    metaAnalysis('/home/minzhang/workspace/pyMeta/BRCA_metastasis.json')

    meta2Csv('/home/minzhang/workspace/pyMeta/BRCA_metastasis.json')

    meta2Heatmap('/home/minzhang/workspace/pyMeta/BRCA_metastasis.json')
    meta2Forest('/home/minzhang/workspace/pyMeta/BRCA_metastasis.json')
