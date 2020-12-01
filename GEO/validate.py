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
# Generate R script for validating meta-analysis results.
# Depending on R package: meta, pheatmap, gplots, RColorBrewer, if running in R.
##############################################################################


import os
import logging
import multiprocessing
import random
import numpy as np
import operator
import statistics
import json
import pandas as pd

from GEO.utility import config
from GEO.utility import mkdir


"""
def sampling(count, population=None, groupList=None, times=1):
    results = list()
    for i in range(0, times):
        if groupList:
            spl = list()
            for j in groupList:
                c = int(len(j) * count / sum([len(j) for j in groupList]))
                spl.extend(random.sample(j, c))
            results.append(spl)
        else:
            results.append(random.sample(population, count))
    if times == 1:
        return results[0]
    else:
        return results
"""


def sampling(count, gsmGroupDict, stratification='balance', times=1):
    results = list()
    for i in range(0, times):
        if stratification == 'random':
            gsmList = random.sample(list(gsmGroupDict.keys()), count)
            results.append({gsm: gsmGroupDict[gsm] for gsm in gsmList})
        else:
            values = set(gsmGroupDict.values())
            grouped = {value: [gsm for gsm in gsmGroupDict
                               if gsmGroupDict[gsm] == value]
                       for value in values}
            counts = {value: len(grouped[value]) for value in values}
            sum_counts = sum([counts[v] for v in values])
            if sum([counts[v] for v in counts]) <= count:
                results.append(gsmGroupDict)
            elif stratification == 'balance':
                c = sorted(counts.items(), key=lambda item:item[1])
                sortedValues = [i[0] for i in c]
                if c[0][1] >= count / len(values):
                    selected = {v: random.sample(grouped[v], int(count/len(values)))
                                for v in values}
                    gsmList = list()
                    for v in values:
                        gsmList.extend(selected[v])
                    results.append({gsm: gsmGroupDict[gsm] for gsm in gsmList})
                else:
                    c = {v: count / len(sortedValues) for v in sortedValues}
                    for i in range(0, len(sortedValues)):
                        if c[sortedValues[i]] <= counts[sortedValues[i]]:
                            break
                        else:
                            c[sortedValues[i]] = counts[sortedValues[i]]
                            for j in range(i+1, len(sortedValues)):
                                c[sortedValues[j]] += (count / len(values) - c[sortedValues[i]]) / (len(values) - 1 - i)
                    selected = {v: random.sample(grouped[v], int(c[v])) for v in
                                values}
                    gsmList = list()
                    for v in values:
                        gsmList.extend(selected[v])
                    results.append({gsm: gsmGroupDict[gsm] for gsm in gsmList})
            else:
                c = {v: counts[v] / sum_counts * count for v in values}
                selected = {v: random.sample(grouped[v], int(c[v])) for v in values}
                gsmList = list()
                for v in values:
                    gsmList.extend(selected[v])
                results.append({gsm: gsmGroupDict[gsm] for gsm in gsmList})
    return results


"""
def samplingGse(gse, gpl, count, clinicalDir, stratification=False, times=1):
    if stratification:
        gsmGroupedList = dict()
        with open('{}/{}_{}_clinical.csv'.format(clinicalDir, gse, gpl), 'r') as rf:
            lines = rf.readlines()
            for line in lines[1:]:
                line = line.split(',')
                if len(line) > 1:
                    gsm = line[0].strip('\t \n\'"')
                    clinic = line[1].strip('\t \n')
                    if clinic not in gsmGroupedList:
                        gsmGroupedList[clinic] = list()
                    gsmGroupedList[clinic].append(gsm)
        gsmGroupedList = [gsmGroupedList[clinic] for clinic in gsmGroupedList]
        if sum([len(i) for i in gsmGroupedList]) > count:
            return sampling(count, groupList=gsmGroupedList, times=times)
        else:
            return samplingGse(gse, gpl, count, clinicalDir, False, times)
    else:
        gsmList = list()
        with open('{}/{}_{}_clinical.csv'.format(clinicalDir, gse, gpl), 'r') as rf:
            lines = rf.readlines()
            for line in lines[1:]:
                gsmList.append(line.split(',')[0].strip('\t \n\'"'))
        if len(gsmList) > count:
            return sampling(count, population=gsmList, times=times)
        else:
            return [gsmList for i in range(0, times)]
"""


def samplingGse(gse, gpl, count, clinicalDir, stratification="", times=1):
    gsmGroupDict = dict()
    with open('{}/{}_{}_clinical.csv'.format(clinicalDir, gse, gpl), 'r') as rf:
        lines = rf.readlines()
        for line in lines[1:]:
            line = line.split(',')
            if len(line) > 1:
                gsm = line[0].strip('\t \n\'"')
                clinic = line[1].strip('\t \n')
                gsmGroupDict[gsm] = clinic
    return sampling(count, gsmGroupDict, stratification, times)


def heatMapExpr(expr, gse, gpl, gsmList, probeList,
                outputFile, outlierReplace=True):
    exp = {probe: [expr[probe][gsm] for gsm in gsmList]
           for probe in probeList}
    if outlierReplace:
        q = {probe: list(np.percentile(exp[probe], [25, 50, 75]))
             for probe in exp}
        n = 1.5
        cutoff = {probe: (q[probe][1] - n * (q[probe][2] - q[probe][0]),
                          q[probe][1] + n * (q[probe][2] - q[probe][0]),
                          q[probe][1])
                  for probe in exp}
        exp = {probe: {gsm: expr[probe][gsm] for gsm in gsmList}
               for probe in probeList}
        e = dict()
        for probe in exp:
            e[probe] = dict()
            for gsm in gsmList:
                if exp[probe][gsm] < cutoff[probe][0]:
                    e[probe][gsm] = -1
                elif exp[probe][gsm] > cutoff[probe][1]:
                    e[probe][gsm] = 1
                else:
                    e[probe][gsm] = (exp[probe][gsm] - cutoff[probe][2]) / (cutoff[probe][1] - cutoff[probe][2]) * n
    else:
        e = exp
    with open(outputFile, 'w') as wf:
        s = '\t{}'.format('\t'.join(gsmList))
        for probe in probeList:
            s += '\n{}\t{}'.format(probe,
                                   '\t'.join([str(e[probe][gsm])
                                              for gsm in gsmList]))
        print(s, file=wf, end='')
    return outputFile


def heatMapGroup(gse, gpl, gsmGroupDict, outputFile):
    with open(outputFile, 'w') as wf:
        s = '\tgroup\n'
        s += '\n'.join(['{}\t{}'.format(gsm, gsmGroupDict[gsm])
                        for gsm in gsmGroupDict])
        print(s, file=wf, end='')
    return outputFile


def exprGSE(gse, gpl, exprDir):
    expr = dict()
    with open('{}/{}_{}_expr.tsv'.format(exprDir, gse, gpl), 'r') as rf:
        lines = rf.readlines()
        gsmList = lines[0].strip('\t\n ').split('\t')
        gsmList = [i.strip('\'"').split('.')[0] for i in gsmList]
        for line in lines[1:]:
            probe, *exp = [i.strip('\'"')
                           for i in line.strip('\n ').split('\t')]
            expr[probe] = dict(zip(gsmList, [float(i) for i in exp]))
    return expr


def MAMA2Heatmap(gse_gpl_list, metaDir, exprDir, clinicalDir, clinicalTag,
                 pCut, esCut, stratification, ramdomTimes, heatmapSampleCount,
                 outlierReplace=True, runningR=True, log=False):
    pVal = dict()
    es = dict()
    with open('{}/p_value'.format(metaDir), 'r') as rf:
        lines = rf.readlines()
        title = [i.strip('"\'\n \t') for i in lines[0].split('\t')]
        col_id = title.index('c_pval') + 1
        for line in lines[1:]:
            s = [i.strip('"\n\t \'') for i in line.split('\t')]
            if len(s) > col_id:
                pVal[s[0]] = float(s[col_id])
    with open('{}/es'.format(metaDir), 'r') as rf:
        lines = rf.readlines()
        title = [i.strip('"\'\n \t') for i in lines[0].split('\t')]
        col_id = title.index('zSco') + 1
        for line in lines[1:]:
            s = [i.strip('"\n\t \'') for i in line.split('\t')]
            if len(s) > col_id:
                es[s[0]] = float(s[col_id])
    with open('{}/pVal_es.tsv'.format(metaDir), 'w') as wf:
        s = '\tp\tes'
        for probe in pVal:
            if probe in es:
                s += '{}\t{}\t{}'.format(probe, pVal[probe], es[probe])
        print(s, file=wf, end='')
    rScript = 'library(gplots)\nlibrary(pheatmap)\nlibrary(RColorBrewer)\n\n'
    for p in pCut:
        for e in esCut:
            dir = '{}/{}_{}'.format(metaDir, p, e)
            mkdir(dir)
            selectedProbe = dict()
            for probe in pVal:
                if probe in es:
                    if pVal[probe] < p and abs(es[probe]) > e:
                        selectedProbe[probe] = es[probe]
            selectedProbe = sorted(selectedProbe.items(),
                                   key=operator.itemgetter(1), reverse=True)
            selectedProbe = [i[0] for i in selectedProbe]
            with open('{}/pVal_es.tsv'.format(dir), 'w') as wf:
                s = '\tp\tes\n'
                s += '\n'.join(['{}\t{}\t{}'.format(probe, pVal[probe], es[probe])
                                for probe in selectedProbe])
                print(s, file=wf, end='')
            for dataset in gse_gpl_list:
                gse = dataset[0]
                gpl = dataset[1]
                mkdir('{}/{}_{}'.format(dir, gse, gpl))
                expr = exprGSE(gse, gpl, exprDir)
                for count in heatmapSampleCount:
                    gsmGroupList = samplingGse(gse, gpl, count,
                                                  clinicalDir,
                                                  stratification=stratification,
                                                  times=ramdomTimes)
                    selectedSamples = [d.keys() for d in gsmGroupList]
                    for i in range(0, ramdomTimes):
                        if count < len(selectedProbe):
                            filteredProbe = selectedProbe[0:int(count/2)]
                            filteredProbe.extend(selectedProbe[int(len(selectedProbe)-count/2):])
                            heatMapExpr(expr, gse, gpl, selectedSamples[i],
                                        filteredProbe,
                                        '{}/{}_{}/{}_{}_expr.tsv'.format(dir, gse, gpl, count, i),
                                        outlierReplace)
                            heatMapGroup(gse, gpl, gsmGroupList[i],
                                         '{}/{}_{}/{}_{}_group.tsv'.format(dir, gse, gpl, count, i))
                        else:
                            heatMapExpr(expr, gse, gpl, selectedSamples[i],
                                        selectedProbe,
                                        '{}/{}_{}/{}_{}_expr.tsv'.format(dir, gse, gpl, count, i),
                                        outlierReplace)
                            heatMapGroup(gse, gpl, gsmGroupList[i],
                                         '{}/{}_{}/{}_{}_group.tsv'.format(dir, gse, gpl, count, i))
                        rScript += '\nexpr <- read.table("{}/{}_{}/{}_{}_expr.tsv")\n'.format(dir, gse, gpl, count, i)
                        rScript += 'x <- as.matrix(expr)\n'
                        rScript += 'png(file="{}/{}_{}/1_{}_{}.png", width=1024, height=1024, bg="transparent")\n'\
                            .format(dir, gse, gpl, count, i)
                        rScript += 'pheatmap(x, cutree_rows=2, cutree_cols=2, color=greenred(75), border_color=NA)\n'
                        rScript += 'dev.off()\n'
                        rScript += 'png(file="{}/{}_{}/2_{}_{}.png", width=1024, height=1024, bg="transparent")\n'\
                            .format(dir, gse, gpl, count, i)
                        rScript += 'heatmap.2(x, col=greenred, scale="row", trace="none")\n'
                        rScript += 'dev.off()\n'
    with open('{}/heatmap.R'.format(metaDir), 'w') as wf:
        print(rScript, file=wf)
    if runningR:
        if log:
            logging.info('Begin to calculate expression by R.')
        cmd = 'Rscript {}/heatmap.R > {}/heatmap.R.log'. \
            format(metaDir, metaDir)
        r = os.system(cmd)
        if log:
            if r == 0:
                logging.info('Calculate expression Successfully.')
            else:
                logging.error('Calculate expression Failed.')


def mDEDS2Heatmap(gse_gpl_list, metaDir, exprDir, clinicalDir, clinicalTag,
                  differentGeneCount, stratification, ramdomTimes,
                  heatmapSampleCount, outlierReplace=True,
                  runningR=True, log=False):
    fc = dict()
    gene = dict()
    with open('{}/geneOrder.tsv'.format(metaDir), 'r') as rf:
        for line in rf.readlines():
            s = line.split('\t')
            if len(s) > 1:
                s = [i.strip('"\' \n') for i in s]
                gene[s[0]] = s[1]
    with open('{}/mDEDS_{}.csv'.format(metaDir, max(differentGeneCount))) as rf:
        lines = rf.readlines()
        title = [i.strip('"\'\n \t').lower() for i in lines[0].split(',')]
        gene_col_id = title.index('geneorder')
        fc_col_id = title.index('fc')
        for line in lines[1:]:
            s = [i.strip('"\n\t \'') for i in line.split(',')]
            if len(s) > fc_col_id:
                fc[gene[s[gene_col_id]]] = float(s[fc_col_id])
    selectedProbe = sorted(fc.items(),
                           key=operator.itemgetter(1), reverse=True)
    selectedProbe = [i[0] for i in selectedProbe]
    selectedProbe = ["KRT80", "SAPCD2", "FHOD1", "CCNB2", "ANLN", "COL4A1", "LOC100506119", "CDKN3", "IER5L", "MMP11", "FAM83D", "LOC286052", "AEBP1", "PCDH17", "SLITRK5", "BUB1B", "CCNB1", "SLC7A5", "DONSON", "SUGCT", "ST8SIA6-AS1", "AVEN", "PCAT6", "KIF14", "MEX3A", "ACKR1", "P2RY12", "WDR78", "MDH1B", "RP11-53O19.3", "C1orf21", "AGBL2", "RLN2", "CCDC176", "HAUS1", "DYNLRB2", "MED13L", "FCGBP", "KIAA1551", "UBXN10", "SUSD3", "PNRC2", "CASC1", "FAM120AOS", "CCDC170", "STC2", "CCR6", "FAM161B", "RP11-28F1.2", "LINC00472"]
    rScript = 'library(gplots)\nlibrary(pheatmap)\nlibrary(RColorBrewer)\n\n'
    for dataset in gse_gpl_list:
        gse = dataset[0]
        gpl = dataset[1]
        mkdir('{}/{}_{}'.format(metaDir, gse, gpl))
        expr = exprGSE(gse, gpl, exprDir)
        for count in heatmapSampleCount:
            gsmGroupList = samplingGse(gse, gpl, count,
                                       clinicalDir,
                                       stratification=stratification,
                                       times=ramdomTimes)
            selectedSamples = [d.keys() for d in gsmGroupList]
            for i in range(0, ramdomTimes):
                if count < len(selectedProbe):
                    filteredProbe = selectedProbe[0:int(count / 2)]
                    filteredProbe.extend(selectedProbe[int(len(selectedProbe) - count / 2):])
                    heatMapExpr(expr, gse, gpl, selectedSamples[i],
                                filteredProbe,
                                '{}/{}_{}/{}_{}_expr.tsv'.format(metaDir, gse, gpl, count, i),
                                outlierReplace)
                    heatMapGroup(gse, gpl, gsmGroupList[i],
                                 '{}/{}_{}/{}_{}_group.tsv'.format(metaDir, gse, gpl, count, i))
                else:
                    heatMapExpr(expr, gse, gpl, selectedSamples[i],
                                selectedProbe,
                                '{}/{}_{}/{}_{}_expr.tsv'.format(metaDir, gse, gpl, count, i),
                                outlierReplace)
                    heatMapGroup(gse, gpl, gsmGroupList[i],
                                 '{}/{}_{}/{}_{}_group.tsv'.format(metaDir, gse, gpl, count, i))
                rScript += '\nexpr <- read.table("{}/{}_{}/{}_{}_expr.tsv")\n'.format(metaDir, gse, gpl, count, i)
                rScript += 'x <- as.matrix(expr)\n'
                rScript += 'png(file="{}/{}_{}/1_{}_{}.png", width=1024, height=1024, bg="transparent")\n'\
                    .format(metaDir, gse, gpl, count, i)
                rScript += 'pheatmap(x, cutree_rows=2, cutree_cols=2, color=greenred(75), border_color=NA)\n'
                rScript += 'dev.off()\n'
                rScript += 'png(file="{}/{}_{}/2_{}_{}.png", width=1024, height=1024, bg="transparent")\n'\
                    .format(metaDir, gse, gpl, count, i)
                rScript += 'heatmap.2(x, col=greenred, scale="row", trace="none")\n'
                rScript += 'dev.off()\n'
                rScript += '\nclinic <- read.table("{}/{}_{}/{}_{}_group.tsv", head=T, row.names=1)\n'.format(metaDir, gse, gpl, count, i)
                rScript += 'c <- clinic\n'
                rScript += 'annotation_c <- data.frame(c)\n'
                rScript += 'rownames(annotation_c) <- colnames(x)\n'
                for j in ['correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski']:
                    rScript += 'png(file="{}/{}_{}/3_{}_{}_{}.png", width=1024, height=1024, bg="transparent")\n'.format(metaDir, gse, gpl, count, i, j)
                    rScript += 'pheatmap(as.matrix(x), annotation_col=annotation_c, color=bluered(200), border_color=NA, cutree_rows=2, cutree_cols=2, clustering_distance_cols="{}", scale="column")\n'.format(j)
                    rScript += 'dev.off()\n'

    with open('{}/heatmap.R'.format(metaDir), 'w') as wf:
        print(rScript, file=wf)
    if runningR:
        if log:
            logging.info('Begin to calculate expression by R.')
        cmd = 'Rscript {}/heatmap.R > {}/heatmap.R.log'. \
            format(metaDir, metaDir)
        r = os.system(cmd)
        if log:
            if r == 0:
                logging.info('Calculate expression Successfully.')
            else:
                logging.error('Calculate expression Failed.')


def meta2Heatmap(settingFile):
    settings = config(settingFile)
    metaDir = settings['meta directory']
    if settings['mapping to gene in meta']:
        exprDir = settings['gene expression table directory']
    else:
        exprDir = settings['expression table directory']
    clinicalDir = settings['simplified clinical data directory']
    clinicalTag = settings['clinical phenotype']
    pCut = settings['p-value cutoff']
    esCut = settings['effect size cutoff']
    gse_gpl_list = [(info['gse'], info['gpl'])
                    for info in settings['datasets information']
                    if info['meta']]
    ramdomTimes = settings['randomized times']
    heatmapSampleCount = settings['max sample in heatmap']
    outlierReplace = settings['outliers replace']
    stratification = settings['randomized sampling']
    runningR = settings['run R in python']
    differentGeneCount = settings['max different expression gene included']
    log = settings['log']

    if settings['meta method'] == 'MAMA':
        return MAMA2Heatmap(gse_gpl_list, metaDir, exprDir, clinicalDir,
                               clinicalTag, pCut, esCut, stratification,
                               ramdomTimes, heatmapSampleCount,
                               outlierReplace, runningR, log)
    elif settings['meta method'] == 'mDEDS':
        return mDEDS2Heatmap(gse_gpl_list, metaDir, exprDir, clinicalDir,
                               clinicalTag, differentGeneCount, stratification,
                               ramdomTimes, heatmapSampleCount, outlierReplace,
                               runningR, log)


def forest(geneList, gse_gplList, exprDir, clinicalDir, outDir, clinicalTag,
           outlierReplace, runningR=True, log=False):
    geneList = [gene for gene in geneList if '/' not in gene]
    data = {gene: dict() for gene in geneList}
    for g in gse_gplList:
        gse = '{}_{}'.format(g[0], g[1])
        with open('{}/{}_clinical.csv'.format(clinicalDir, gse), 'r') as rf:
            lines = rf.readlines()[1:]
            lines = [line.strip().split(',') for line in lines if len(line) > 3]
            lines = [[item.strip(' \n\t"\'') for item in line] for line in lines]
            clinic = {line[0]: line[1] for line in lines}
            control = [gsm.upper().split('.')[0]
                       for gsm in clinic
                       if clinic[gsm] == '0']
            treat = [gsm.upper().split('.')[0]
                     for gsm in clinic
                     if clinic[gsm] == '1']
        with open('{}/{}_expr.tsv'.format(exprDir, gse), 'r') as rf:
            lines = rf.readlines()
            gsmList = lines[0].strip().split('\t')
            gsmList = [gsm.strip(' \t\n\'"').upper().split('.')[0] for gsm in gsmList]
            for line in lines[1:]:
                if len(line) > 3:
                    line = [item.strip('\t\n \'"') for item in line.split('\t')]
                    if line[0] in geneList:
                        gene = line[0]
                        expr = [float(i) for i in line[1:]]
                        if outlierReplace:
                            q = list(np.percentile(expr, [25, 50, 75]))
                            n = 1.5
                            cutoff = (q[1] - n * (q[2] - q[0]),
                                      q[1] + n * (q[2] - q[0]),
                                      q[1])
                            exp = list()
                            for e in expr:
                                if e < cutoff[0]:
                                    e = -1
                                elif e > cutoff[1]:
                                    e = 1
                                else:
                                    e = (e - cutoff[2]) / (cutoff[1] - cutoff[2]) * n
                                exp.append(e)
                            expr = exp
                        expr = dict(zip(gsmList, expr))
                        controlExpr = [float(expr[gsm]) for gsm in control if gsm in expr]
                        treatExpr = [float(expr[gsm]) for gsm in treat if gsm in expr]
                        if gse not in data[gene]:
                            data[gene][gse] = dict()
                        data[gene][gse]['NContrl'] = len(control)
                        data[gene][gse]['MeanControl'] = statistics.mean(controlExpr)
                        data[gene][gse]['StdControl'] = statistics.stdev(controlExpr)
                        data[gene][gse]['NTreat'] = len(treat)
                        data[gene][gse]['MeanTreat'] = statistics.mean(treatExpr)
                        data[gene][gse]['StdTreat'] = statistics.stdev(treatExpr)

    mkdir('{}/tmp'.format(outDir))
    rGeneList = list()
    for i in range(0, len(geneList)):
        gene = geneList[i]
        if len(data[gene]) > 1:
            with open('{}/tmp/{}_{}.csv'.format(outDir, i, gene), 'w') as wf:
                s = 'study,n1,mean1,sd1,n2,mean2,sd2\n'
                s += '\n'.join(
                    ['{},{},{},{},{},{},{}'.format(gse,
                                                   data[gene][gse]['NTreat'], data[gene][gse]['MeanTreat'], data[gene][gse]['StdTreat'],
                                                   data[gene][gse]['NContrl'], data[gene][gse]['MeanControl'], data[gene][gse]['StdControl'])
                    for gse in data[gene]]
                )
                rGeneList.append('{}_{}'.format(i, gene))
                print(s, file=wf, end='')
    rScript = 'library(parallel)\n'
    rScript += 'f <- function(gene) {\n    library(meta)\n'
    rScript += '    a <- read.table(paste("{}/tmp/", gene, ".csv", sep=""), sep=",", header=T, stringsAsFactors=FALSE)\n'.format(outDir)
    rScript += '    metarsmd=metacont(n1,mean1,sd1,n2,mean2,sd2,data=a,sm="SMD",label.c="Non-metastasis",label.e="metastasis",comb.fixed=FALSE,comb.random=TRUE,studlab=study)\n'
    rScript += '    png(file=paste("{}/", gene, ".png", sep=""), width=12, height=6, units="in", res=600, bg="transparent")\n'.format(outDir)
    rScript += '    forest(metarsmd)\n    dev.off()\n'
    rScript += '    c(gene, metarsmd$TE.random, metarsmd$lower.random, metarsmd$upper.random, metarsmd$k)\n}\n\n'
    rScript += 'genes <- c("{}")\n'.format('", "'.join(rGeneList))
    rScript += 'detectCores(logical = F)\ncl <- makeCluster(getOption("cl.cores", 4))\n'
    rScript += 'ci <- parLapply(cl, genes, f)\nstopCluster(cl)\n'
    rScript += 'write.table(ci, file="{}/ci.csv", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)'.format(outDir)

    with open('{}/forest.R'.format(outDir), 'w') as wf:
        print(rScript, file=wf, end='')

    if runningR:
        if log:
            logging.info('Begin to calculate expression by R.')
        cmd = 'Rscript {}/forest.R > {}/forest.R.log'.format(outDir, outDir)
        r = os.system(cmd)
        if log:
            if r == 0:
                logging.info('Forest plot Successfully.')
            else:
                logging.error('Forest plot Failed.')


def MAMA2Forest(gse_gpl_list, metaDir, exprDir, clinicalDir, clinicalTag,
                pCut, esCut, thread=4, gplMapFile='',
                outlierReplace=False, runningR=True, log=False):
    paraList = list()
    for p in pCut:
        for es in esCut:
            outputDir = '{}/{}_{}/forest'.format(metaDir, p, es)
            mkdir(outputDir)
            with open('{}/{}_{}/pVal_es.tsv'.format(metaDir, p, es), 'r') as rf:
                lines = rf.readlines()[1:]
                geneList = [line.split('\t')[0] for line in lines]
            if gplMapFile:
                with open(gplMapFile, 'r') as rf:
                    annot = json.load(rf)
                geneListMapped = list()
                for gene in geneList:
                    if gene in annot:
                        mapped = annot[gene]
                        if '///' in mapped:
                            """
                            mappedList = [i.strip() for i in mapped.split('///')]
                            for g in mappedList:
                                if g not in geneListMapped:
                                    geneListMapped.append(g)
                            """
                            pass
                        else:
                            if mapped not in geneListMapped:
                                geneListMapped.append(mapped)
                geneList = geneListMapped
                with open('{}/geneListMapped.json'.format(outputDir), 'w') as wf:
                    json.dump(geneList, wf, indent=2)
            paraList.append((geneList, gse_gpl_list, exprDir, clinicalDir,
                             outputDir, clinicalTag,
                             outlierReplace, runningR, log))
    with multiprocessing.Pool(processes=thread) as pool:
        pool.starmap(forest, paraList)


def mDEDS2Forest(gse_gpl_list, metaDir, exprDir, clinicalDir, clinicalTag,
                 outlierReplace=False, runningR=True, log=False):
    geneList = list()
    geneSet = set()
    df = pd.read_csv('{}/__meta_results.csv'.format(metaDir))
    symbol = list(df['symbol'])
    for s in symbol:
        s = str(s)
        if s != 'nan' and '/' not in s:
            if s not in geneSet:
                geneList.append(s)
                geneSet.add(s)
    outDir = '{}/forest'.format(metaDir)
    mkdir(outDir)
    forest(geneList, gse_gpl_list, exprDir, clinicalDir, outDir, clinicalTag,
           outlierReplace, runningR, log)


def meta2Forest(settingFile):
    settings = config(settingFile)
    metaDir = settings['meta directory']
    if settings['mapping to gene in valid']:
        exprDir = settings['gene expression table directory']
    else:
        exprDir = settings['expression table directory']
    clinicalDir = settings['simplified clinical data directory']
    clinicalTag = settings['clinical phenotype']
    thread = settings['thread']
    outlierReplace = settings['outliers replace']
    gse_gpl_list = [(info['gse'], info['gpl'])
                    for info in settings['datasets information']]
    gpl = [info['gpl']
           for info in settings['datasets information']
           if info['meta']][0]
    log = settings['log']
    runningR = settings['run R in python']
    gplMapFile = ''
    if settings['probe mapping to gene'] \
            and not settings['mapping to gene in meta']\
            and settings['mapping to gene in valid']:
        gplMapFile = '{}/{}_probe2{}.json'.format(
            settings['probe annotation directory'],
            gpl,
            settings['probe mapping to gene keyword']
        )
    if settings['meta method'] == 'MAMA':
        pCut = settings['p-value cutoff']
        esCut = settings['effect size cutoff']
        return MAMA2Forest(gse_gpl_list, metaDir, exprDir, clinicalDir,
                           clinicalTag, pCut, esCut, thread,
                           gplMapFile, outlierReplace, runningR, log)
    elif settings['meta method'] == 'mDEDS':
        maxDEDSGeneCount = max(settings['max different expression gene included'])
        return mDEDS2Forest(gse_gpl_list, metaDir, exprDir, clinicalDir,
                            clinicalTag, outlierReplace, runningR, log)

