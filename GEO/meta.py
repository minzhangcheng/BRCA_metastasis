##############################################################################
#
# Copyright (C) 2018 Minzhang Cheng
# Contact: minzhangcheng@gmail.com
#
# GNU Lesser General Public License Usage
# This file may be used under the terms of the GNU Lesser General Public
# License version 3 as published by the Free Software Foundation and
# appearing in the file LICENSE included in the packaging of this file.
# Please review the following information to ensure the GNU Lesser
# General Public License version 3 requirements will be met:
# http://www.gnu.org/licenses/gpl-3.0.html
#
##############################################################################

##############################################################################
# Generate R script for metaQC and meta-analysis.
# Depending on R package: MetaQC, MAMA, DEDS, CONOR, if running in R.
##############################################################################


import os
import logging
import multiprocessing
import random
import numpy as np
import operator
import statistics
import pandas as pd
import json

from GEO.utility import config
from GEO.utility import mkdir


def geneExpr2metaqcTable(gse, gpl, metaqcDir,
                         geneExprDir, clinicalDir):
    expr = dict()
    clinical = dict()
    with open('{}/{}_{}_expr.tsv'.format(geneExprDir, gse, gpl), 'r') as rf:
        lines = rf.readlines()
        gsmList = [i.split('.')[0] for i in lines[0].strip().split('\t')]
        for line in lines[1:]:
            if len(line) > 3:
                gene, *e = [i.strip('\n "') for i in line.split('\t')]
                expr[gene] = e
    with open('{}/{}_{}_clinical.csv'.format(clinicalDir, gse, gpl), 'r') as rf:
        for line in rf.readlines()[1:]:
            gsm, c = [i.strip('\n "') for i in line.split(',')]
            clinical[gsm] = c
    with open('{}/{}_{}.txt'.format(metaqcDir, gse, gpl), 'w') as wf:
        s = '\t{}\n'.format('\t'.join(gsmList))
        s += 'clinic\t{}'.format('\t'.join([clinical[gsm] for gsm in gsmList]))
        for gene in expr:
            s += '\n{}\t{}'.format(gene, '\t'.join(expr[gene]))
        print(s, file=wf, end='')


def metaQC(settingFile):
    settings = config(settingFile)
    metaqcDir = settings['metaQC directory']
    geneExprDir = settings['gene expression table directory']
    clinicalDir = settings['simplified clinical data directory']
    gse_gpl_list = ['{}_{}'.format(info['gse'], info['gpl'])
                    for info in settings['datasets information']
                    if info['metaQC']]
    gmtFile = settings['gmt File']

    p = [(info['gse'], info['gpl'], metaqcDir, geneExprDir, clinicalDir)
         for info in settings['datasets information']
         if info['metaQC']]
    with multiprocessing.Pool(processes=settings['thread']) as pool:
        results = pool.starmap(geneExpr2metaqcTable, p)

    s = 'rm(list = ls())\nlibrary("MetaDE")\nlibrary("MetaQC")\n\n'
    s += 'setwd("{}")\nmemory.limit(16000)\n\n'.format(metaqcDir)
    s += 'study.names <- c({})\n'.format(', '.join(['"%s"' % i
                                                    for i in gse_gpl_list]))
    s += 'raw <- MetaDE.Read(study.names, skip=rep(1, {}), via="txt", matched=T, log=T)\n'.format(len(gse_gpl_list))
    s += 'gc()\n'
    s += 'merged <- MetaDE.merge(raw)\ngc()\n'
    s += 'Data.QC<-list()\nfor(i in 1:%d){\n' % len(gse_gpl_list)
    s += '  colnames(merged[[i]][[1]])<-merged[[i]][[2]]\n'
    s += '  Data.QC[[i]]<-impute.knn(merged[[i]][[1]])$data\n}\n'
    s += 'names(Data.QC)<-names(merged)\n'
    s += 'QC <- MetaQC(Data.QC, "{}", filterGenes=F,verbose=TRUE, isParallel=T, resp.type="Twoclass")\n'.format(gmtFile)
    s += 'gc()\nrunQC(QC, B=1e4, fileForCQCp="{}")\n'.format(gmtFile)
    s += 'png(file="{}/metaQC.png", width=2048, height=2048, bg="transparent")\n'.format(metaqcDir)
    s += 'plot(QC)\ndev.off()\n'
    s += 'sink("{}/metaQC.txt")\nprint(QC)\nsink()\n'.format(metaqcDir)
    s += 'save.image("{}/metaQC.RData")'.format(metaqcDir)

    with open('{}/metaQC.R'.format(metaqcDir), 'w') as wf:
        print(s, file=wf, end='')

    if settings['run R in python']:
        if settings['log']:
            logging.info('Begin QC by R.')
            cmd = 'Rscript {}/metaQC.R >{}/metaQC.R.log'. \
                format(metaqcDir, metaqcDir)
            if os.system(cmd) == 0:
                logging.info('QC Successfully.')
                return True
            else:
                logging.error('QC Failed.')
                return False


def metaAnalysis(settingFile):
    settings = config(settingFile)
    metaDir = settings['meta directory']
    geneExprDir = settings['gene expression table directory']
    exprDir = settings['expression table directory']
    clinicalDir = settings['simplified clinical data directory']
    clinicalTag = settings['clinical phenotype']
    gse_gpl_list = ['{}_{}'.format(info['gse'], info['gpl'])
                    for info in settings['datasets information']
                    if info['meta']]
    metaMethod = settings['meta method']

    if metaMethod == 'MAMA':
        s = 'rm(list = ls())\nlibrary("MAMA")\nlibrary("metaMA")\n'
        s += 'library("affyPLM")\nlibrary("affy")\n\n'
        s += 'memory.limit(16000)\n\n'
        for gse in gse_gpl_list:
            if settings['mapping to gene in meta']:
                s += '{} <- as.matrix(read.table("{}/{}_expr.tsv"))\n' \
                    .format(gse, geneExprDir, gse)
            else:
                s += '{} <- as.matrix(read.table("{}/{}_expr.tsv"))\n' \
                    .format(gse, exprDir, gse)
        s += 'GEDM<-list({})\n'.format(', '.join(gse_gpl_list))
        s += 'rm({})\ngc()\n\n'.format(', '.join(gse_gpl_list))
        s += 'setwd("{}")\n'.format(metaDir)
        for gse in gse_gpl_list:
            s += 'annot_%s <- read.csv("%s/%s_clinical.csv")\n' \
                 % (gse, clinicalDir, gse)
            s += 'row.names(annot_%s) <- annot_%s$CEL_Number\n' \
                 % (gse, gse)
            s += 'annot_%s$%s <- as.factor(annot_%s$%s)\n' \
                 % (gse, clinicalTag, gse, clinicalTag)
        s += 'Clinical <- list(%s)\n' \
             % ', '.join(['annot_%s' % gse for gse in gse_gpl_list])
        s += 'rm(%s)\ngc()\n' \
             % ', '.join(['annot_%s' % gse for gse in gse_gpl_list])
        s += 'datanames <- c(%s)\n\n' \
             % ', '.join(['"%s"' % gse for gse in gse_gpl_list])
        s += 'setClass("esets",slots=list(GEDM="list",clinical="list",datanames="character"),package = "MAMA")\n'
        s += 'HCC_meta_data <- new("esets",GEDM=GEDM,clinical=Clinical,datanames=datanames)\n'
        s += 'rm(Clinical,datanames)\ngc()\n\n'

        s += 'pval <- metaMA(HCC_meta_data,"%s",which="pval")\n' % clinicalTag
        s += 'gc()\n'
        s += 'es2 <- ES.GeneMeta(HCC_meta_data,"%s",nperm=1000)\n\n' % clinicalTag

        s += 'results1 <- join.results(pval,type=1,genenames=rownames(GEDM(HCC_meta_data)[[1]]))\n'
        s += 'p_value <- as.data.frame(results1)\n'
        s += 'rawpval = 2 * (1 - pnorm(abs(pval$TestStatistic)))\n'
        s += 'FDR_pval <- p.adjust(rawpval, method="BY", n=length(rawpval))\n'
        s += 'p_value$c_pval <- rawpval\n'
        s += 'p_value$FDR <- FDR_pval\n'
        s += 'rm(rawpval, FDR_pval)\ngc()\n'
        s += 'es2_theScores <- es2$theScores\n'
        s += 'es2_ScoresFDR <- es2$ScoresFDR\n'
        s += 'es2_ScoresFDR <- es2_ScoresFDR$two.sided\n'
        s += 'write.table(p_value, file="p_value", sep="\\t", col.names = T)\n'
        s += 'write.table(es2_ScoresFDR, file="es", sep="\\t", col.names = T)\n'

        with open('{}/{}.R'.format(metaDir, metaMethod), 'w') as wf:
            print(s, file=wf)

        if settings['run R in python']:
            if settings['log']:
                logging.info('Begin meta-analysis by R.')
                cmd = 'Rscript {}/{}.R >{}/{}.R.log'. \
                    format(metaDir, metaMethod, metaDir, metaMethod)
                if os.system(cmd) == 0:
                    logging.info('meta-analysis Successfully.')
                    return True
                else:
                    logging.error('meta-analysis Failed.')
                    return False

    elif metaMethod == 'mDEDS':
        with open('{}/clinic.csv'.format(metaDir), 'w') as wf:
            print('gsm,{}'.format(clinicalTag), file=wf, end='')
            for dataset in gse_gpl_list:
                with open('{}/{}_clinical.csv'.format(clinicalDir, dataset), 'r') as rf:
                    s = '\n'
                    for line in rf.readlines()[1:]:
                        s += line
                    print(s, file=wf, end='')
        s = 'rm(list = ls())\nlibrary("MAMA")\nlibrary("metaMA")\n'
        s += 'library("affyPLM")\nlibrary("affy")\n'
        s += 'library("CONOR")\nlibrary("DEDS")\n\n'
        s += 'memory.limit(16000)\n\n'

        gse = gse_gpl_list[0]
        if settings['mapping to gene in meta']:
            s += 'a <- read.table("%s/%s_expr.tsv")\n' % (geneExprDir, gse)
        else:
            s += 'a <- read.table("%s/%s_expr.tsv")\n' % (exprDir, gse)
        for gse in gse_gpl_list[1:]:
            if settings['mapping to gene in meta']:
                s += 'b <- read.table("%s/%s_expr.tsv")\n' % (geneExprDir, gse)
            else:
                s += 'b <- read.table("%s/%s_expr.tsv")\n' % (exprDir, gse)
            s += 'm <- xpn(a, b, iterations=10)\na <- m$x\nb <- m$y\n'
            s += 'm <- merge(a, b, by="row.names")\n'
            s += 'row.names(m) <- m$Row.names\na <- m[, -1]\n'
            s += 'a <- a[, order(colnames(a))]\na <- as.matrix(a)\n'
            s += 'rm(b, m)\ngc()\n\n'
        s += 'merged.expr <- a\n'

        s += 'save(merged.expr, file="%s/merged.expr.RData")\n' % metaDir
        s += 'rm(a)\ngc()\n\n'

        s += 'annot <- read.csv("{}/clinic.csv", row.names=1, header=T)\n'.format(metaDir)
        s += 'merged <- as.data.frame(merged.expr)\n'
        s += 'annot <- annot[order(row.names(annot)), ]\n'
        s += 'merged <- merged[, order(colnames(merged))]\n'
        s += 'merged <- as.matrix(merged)\n'
        # s += 'merged <- 2 ^ merged\n'
        s += 'rm(merged.expr)\ngc()\n'
        s += 'deds <- deds.stat.linkC(merged, annot, B=1300)\n'
        s += 'save(deds, file="{}/deds.RData")\n'.format(metaDir)
        s += 'gc()\n\n'

        for i in settings['max different expression gene included']:
            s += 'r <- topgenes(deds, number=%d, sort.by="fc")\n' % i
            s += 'write.csv(r, file="%s/mDEDS_%d.csv")\n' % (metaDir, i)

        s += 'write.table(rownames(merged), file="{}/geneOrder.tsv", sep="\t")'.format(metaDir)

        with open('{}/{}.R'.format(metaDir, metaMethod), 'w') as wf:
            print(s, file=wf)

        if settings['run R in python']:
            if settings['log']:
                logging.info('Begin meta-analysis by R.')
            cmd = 'Rscript {}/{}.R >{}/{}.R.log'. \
                format(metaDir, metaMethod, metaDir, metaMethod)
            if os.system(cmd) == 0:
                if settings['log']:
                    logging.info('meta-analysis Successfully.')
                return True
            else:
                if settings['log']:
                    logging.error('meta-analysis Failed.')
                return False


def meta2Csv(settingFile):
    settings = config(settingFile)
    metaDir = settings['meta directory']
    metaMethod = settings['meta method']

    if not settings['mapping to gene in meta'] \
            and settings['probe mapping to gene']:
        annotDir = settings['probe annotation directory']
        gpls = [ds['gpl']
                for ds in settings['datasets information'] if ds['meta']]
        gpls = set(gpls)
        rfs = [open('{}/{}_probe2symbol.json'.format(annotDir, gpl), 'r') for gpl in gpls]
        probe2symbols = [json.load(rf) for rf in rfs]
        probe2symbol = dict()
        for i in probe2symbols:
            probe2symbol.update(i)
        for rf in rfs:
            rf.close()

    if metaMethod == 'MAMA':
        pass

    elif metaMethod == 'mDEDS':
        geneCount = max(settings['max different expression gene included'])
        geneOrderFile = '{}/geneOrder.tsv'.format(metaDir)
        with open(geneOrderFile, 'r') as rf:
            lines = rf.readlines()[1:]
            table = [[i.replace('"', '').strip() for i in line.split('\t')]
                     for line in lines if len(line.split('\t')) > 1]
            geneOrder = {
                i[0]: i[1] for i in table
            }

        with open('{}/mDEDS_{}.csv'.format(metaDir, geneCount), 'r') as rf:
            lines = rf.readlines()
            title = [i.replace('"', '').strip() for i in lines[0].split(',')]
            table = [[i.replace('"', '').strip() for i in line.split(',')]
                     for line in lines[1:] if len(line.split(',')) > 1]
            data = dict()
            for i in range(len(title)):
                t = title[i]
                if t.upper() == 'geneOrder'.upper():
                    if settings['mapping to gene in meta']:
                        data['input order'] = [row[i] for row in table]
                        data['symbol'] = [geneOrder[j] for j in data['input order']]
                    else:
                        data['input order'] = [row[i] for row in table]
                        data['probe'] = [geneOrder[j] for j in data['input order']]
                        if settings['probe mapping to gene']:
                            data['symbol'] = list()
                            for p in data['probe']:
                                if p in probe2symbol:
                                    data['symbol'].append(probe2symbol[p])
                                else:
                                    data['symbol'].append('')
                elif t == '':
                    data['order'] = [row[i] for row in table]
                elif t.lower() == 'fc'.lower():
                    data['fc'] = [row[i] for row in table]
                    data['abs(fc)'] = [-abs(float(row[i])) for row in table]
                else:
                    data[t.lower()] = [row[i] for row in table]

            df = pd.DataFrame(data).sort_values(by=['deds', 'abs(fc)'])
            df.to_csv(
                '{}/__meta_results.csv'.format(metaDir),
                index=False
            )

