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
# Download original CEL files, prepare clinical data, and generate the R
#  scriptcalculating gene expressions.
#
# Only works on Python 3.3+
##############################################################################


import logging
import multiprocessing
import multiprocessing.dummy
import json
import requests
import bs4
import shutil
import tempfile
import tarfile
import os
import statistics
import sqlite3

from GEO.utility import extractGzipFile
from GEO.utility import mkdir
from GEO.utility import listFtpFiles
from GEO.utility import download
from GEO.utility import transpose
from GEO.utility import config


def extractClinical(datasetsInfo, gseClinicalCsvDir, outputClinicalDir,
                    clinicalTag='label', log=True):
    clinicalDataAll = dict()
    gsmSeparated = dict()
    gpl2gsm = dict()
    for info in datasetsInfo:
        gse = info['gse']
        if isinstance(info['gpl'], str):
            gplList = [info['gpl']]
        elif isinstance(info['gpl'], list):
            gplList = info['gpl']
        else:
            raise(Exception)
        clinicalData = dict()
        gpl = info['gpl']
        with open('{}/{}_{}_clinical.tsv'.format(gseClinicalCsvDir, gse, gpl), 'r') as rf:
            lines = rf.readlines()
            if 'col_name' in info:
                name = lines[0].split('\t')[info['col_id']].strip()
                if name != info['col_name']:
                    logging.error('Column name not matched in %s!\n' % gse)
            for line in lines[1:]:
                gsm = line.split('\t')[0].strip().upper()
                val = line.split('\t')[info['col_id']].strip()
                if val == '':
                    if 'dropNA' in info and not info['dropNA']:
                        val = 'NA'
                    else:
                        continue
                if info['type'] == 'continuous':
                    try:
                        val = float(val)
                    except:
                        if 'dropNA' not in info or info['dropNA']:
                            continue
                        else:
                            val = 'NA'
                    if 'coefficient' in info:
                        if type(val) == type(1.0):
                            val *= info['coefficient']
                elif info['type'] == 'category':
                    if 'map' in info:
                        if val in info['map']:
                            val = info['map'][val]
                        else:
                            if '_default_' in info['map']:
                                val = info['map']['_default_']
                            else:
                                if 'dropNA' not in info or info['dropNA']:
                                    logging.warning('Value of {gsm} in {gse} is not valid.'.format(gsm=gsm, gse=gse))
                                    continue
                elif info['type'] == 'classed':
                    # todo: classify the continuous data into different types
                    pass
                clinicalData[gsm] = val
                clinicalDataAll[gsm] = val
                if gpl not in gpl2gsm:
                    gpl2gsm[gpl] = list()
                gpl2gsm[gpl].append(gsm)
                if (gse, gpl) not in gsmSeparated:
                    gsmSeparated[(gse, gpl)] = list()
                gsmSeparated[(gse, gpl)].append(gsm.upper())
            with open('{}/{}_{}_clinical.csv'.format(outputClinicalDir, gse, gpl), 'w') as wf:
                s = 'gsm,{}\n'.format(clinicalTag)
                s += '\n'.join(['{},{}'.format(gsm, clinicalData[gsm]) for gsm in sorted(clinicalData.keys())])
                print(s, file=wf, end='')
    with open('{}/all_clinical.csv'.format(outputClinicalDir), 'w') as wf:
        s = 'gsm,{}\n'.format(clinicalTag)
        s += '\n'.join(['{},{}'.format(gsm, clinicalDataAll[gsm]) for gsm in sorted(clinicalDataAll.keys())])
        print(s, file=wf, end='')
    for gpl in gpl2gsm:
        with open('{}/{}_clinical.csv'.format(outputClinicalDir, gpl), 'w') as wf:
            s = 'gsm,{}\n'.format(clinicalTag)
            s += '\n'.join(['{},{}'.format(gsm, clinicalDataAll[gsm]) for gsm in sorted(gpl2gsm[gpl])])
            print(s, file=wf, end='')
    return gsmSeparated


def getGseRaw(gse, directory, log=True,
              maxTrial=4, downloadMethod='wget'):
    # dowloadMethod in ['wget', 'ftplib']
    filename = ''
    ftpUrl = 'ftp.ncbi.nlm.nih.gov'
    if len(gse) > 6:
        wd = 'geo/series/GSE%snnn/%s/suppl/' % (gse[3:-3], gse)
    else:
        wd = 'geo/series/GSEnnn/%s/suppl/' % gse
    filename = '%s_RAW.tar' % gse
    url = 'ftp://{}/{}{}'.format(ftpUrl, wd, filename)
    if download(url, filename, directory, maxTrial, downloadMethod, log, False):
        if log:
            logging.info('Download %s successfully' % gse)
        return True, gse, filename
    else:
        msg = 'Download {} failed and terminated: directory: {}, filename {}'. \
            format(gse, wd, filename)
        logging.error(msg)
        return False, gse, filename


def getGseRawList(gseList, directory, log=True,
                  saveFilenameList=True, maxTrial=4, thread=2,
                  downloadMethod='wget'):
    if thread > 1:
        p = set(gseList)
        p = [(gse, directory, log, maxTrial, downloadMethod)
             for gse in p]
        with multiprocessing.dummy.Pool(processes=thread) as pool:
            results = pool.starmap(getGseRaw, p)
    else:
        results = [getGseRaw(gse, directory, log,
                                   maxTrial, downloadMethod)
                   for gse in gseList]
    r = dict()
    for i in results:
        if i[0]:
            if i[1] in r:
                r[i[1]].append(i[2])
            else:
                r[i[1]] = [i[2]]
    if saveFilenameList:
        with open('{}/download_gse_file_list.json'.format(directory), 'w') as wf:
            json.dump(r, wf, indent=2)
    return r


def getGplSoft(gpl, directory, log=True,
               maxTrial=4, downloadMethod='wget'):
    if len(gpl) > 6:
        wd = 'geo/platforms/GPL%snnn/%s/' % (gpl[3:-3], gpl)
    else:
        wd = 'geo/platforms/GPLnnn/%s/' % gpl
    filename = '%s.annot.gz' % gpl
    ftpUrl = 'ftp.ncbi.nlm.nih.gov'
    url = 'ftp://{}/{}'.format(ftpUrl, wd)
    if 'annot' in listFtpFiles(url):
        url = '{}/annot/filename'.format(url)
        if download(url, filename, directory, maxTrial,
                    downloadMethod, log, False):
            if log:
                logging.info('Download %s successfully' % gpl)
            return True, '%s' % gpl, filename
        else:
            if log:
                msg = 'Download {} failed: directory: {}, filename {}'. \
                    format(gpl, wd, filename)
                logging.error(msg)
            return False, '%s' % gpl, filename
    else:
        url = '{}/suppl/'.format(url)
        fileList = [i for i in listFtpFiles(url) if i.find('annot') >= 0]
        success = True
        for filename in fileList:
            url_file = '{}/{}'.format(url, filename)
            if not download(url, filename, directory,
                            maxTrial, downloadMethod, log, False):
                success = False
        if success:
            if log:
                logging.info('Download %s successfully' % gpl)
            return True, '%s' % gpl, ';'.join(fileList)
        else:
            if log:
                msg = 'Download {} failed: directory: {}, filename {}'. \
                    format(gpl, wd, filename)
                logging.error(msg)
            return False, '%s' % gpl, filename


def getGplSoftList(gplList, directory, log=True,
                   saveFilenameList=True, maxTrial=4,
                   thread=4, downloadMethod='wget'):
    if thread > 1:
        p = [(id, directory, log, maxTrial, downloadMethod) for id in gplList]
        with multiprocessing.dummy.Pool(processes=thread) as pool:
            results = pool.starmap(getGplSoft, p)
    else:
        results = []
        for id in gplList:
            results.append(getGplSoft(id, directory, log,
                                            maxTrial, downloadMethod))
    r = {i[1]: i[2] for i in results if i[0]}
    if saveFilenameList:
        with open('{}/download_gpl_file_list.json'.format(directory), 'w') as wf:
            json.dump(r, wf, indent=2)
    return r


def getGseMatrix(gse, gpl, directory, log=True,
                 maxTrial=4, downloadMethod='wget'):
    filename = '{}-{}_series_matrix.txt.gz'.format(gse, gpl)
    ftpUrl = 'ftp.ncbi.nlm.nih.gov'
    if len(gse) > 6:
        wd = 'geo/series/GSE%snnn/%s/matrix' % (gse[3:-3], gse)
    else:
        wd = 'geo/series/GSEnnn/%s/matrix' % gse
    url = '{}/{}'.format(ftpUrl, wd)
    if len(listFtpFiles(url)) > 1:
        file = '{}-{}_series_matrix.txt.gz'.format(gse, gpl)
    else:
        file = '{}_series_matrix.txt.gz'.format(gse)
    url = '{}/{}'.format(url, file)
    if download(url, filename, directory, maxTrial, downloadMethod, log, False):
        if log:
            logging.info('Download {}_{} successfully'.format(gse, gpl))
        return True, gse, gpl, filename
    else:
        if log:
            msg = 'Download %s_%s failed and terminated, directory %s, filename %s, url %s' \
                  % (gse, gpl, wd, filename, url)
            logging.error(msg)
        return False, gse, gpl, filename


def getGseMatrixList(gse_gpl_List, directory, log=True, saveFilenameList=True,
                     maxTrial=4, thread=4, downloadMethod='wget'):
    """
    gse_gpl_List = [('GSE1234', 'GPL1234'),
                    ('GSE1234', 'GPL1234'),
                    ('GSE1234', 'GPL1234')]
    """
    gse_gpl_List = set(gse_gpl_List)
    if thread > 1:
        p = list()
        for id in gse_gpl_List:
            p.append((id[0], id[1], directory, log, maxTrial, downloadMethod))
        with multiprocessing.dummy.Pool(processes=thread) as pool:
            results = pool.starmap(getGseMatrix, p)
    else:
        results = [getGseMatrix(id[0], id[1], directory, log,
                                      maxTrial, downloadMethod)
                   for id in gse_gpl_List]
    r = dict()
    for i in results:
        if i[0]:
            if i[1] not in r:
                r[i[1]] = dict()
            r[i[1]][i[2]] = i[3]
    if saveFilenameList:
        with open('{}/download_matrix_file_list.json'.format(directory), 'w') as wf:
            json.dump(r, wf, indent=2)
    return r


def matrix2expr(matrixFile):
    with open(matrixFile, 'r') as rf:
        exprStr = ''
        exprMatBegin = False
        samples = list()
        # probes = list()
        expr = dict()
        gsmLine = False
        for line in rf.readlines():
            if exprMatBegin:
                if not line.startswith('!'):
                    line = line.split('\t')
                    line = [i.strip('"\' \n') for i in line]
                    # probes.append(line[0])
                    expr[line[0]] = dict(zip(samples, line[1:]))
            elif gsmLine:
                s = line.split('\t')[1:]
                samples = [i.strip('"\' \n').upper() for i in s]
                exprMatBegin = True
            elif line.startswith('!series_matrix_table_begin'):
                gsmLine = True
    return expr


def getGplAnnotTable(gpl, directory, log=True, maxTrial=4):
    url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi'
    try:
        r = requests.get(url, {'view': 'data', 'acc': gpl})
        soup = bs4.BeautifulSoup(r.text, 'lxml')
        r = soup.get_text()
    except:
        if log:
            logging.error('Failed in downloading GPL annotation table for {} !'.format(gpl))
        if maxTrial > 1:
            return getGplAnnotTable(gpl, directory,
                                       log, maxTrial=maxTrial - 1)
        else:
            if log:
                logging.error('Failed in downloading GPL annotation table for {} !\n'
                              'Max trial reached, and give up downloading this file!'.format(gpl))
            return False, gpl
    with open('{}/{}.annot'.format(directory, gpl), 'w') as wf:
        print(r, file=wf, end='')
    if log:
        logging.info('Download GPL annotation table for {} Successfully!'.format(gpl))
    return True, gpl


def getGplAnnotTableList(gplList, directory, log=True, maxTrial=4, thread=4):
    gplList = set(gplList)
    if thread > 1:
        p = list()
        p = [(gpl, directory, log, maxTrial) for gpl in gplList]
        with multiprocessing.dummy.Pool(processes=thread) as pool:
            results = pool.starmap(getGplAnnotTable, p)
    else:
        results = [getGplAnnotTable(gpl, directory, log, maxTrial)
                   for gpl in gplList]
    return [i[1] for i in results if i[0]]


def annot2Map(gpl, gplAnnotDir, key=None, log=True, tag=None):

    if not key:
        key = 'symbol'

    _tag = {
        'symbol': ['symbol'],
        'unigene': ['unigene', 'uni_gene', 'uni gene'],
        'genbank': ['gb', 'genbank']
    }
    if not tag:
        tag = _tag

    def _title(s):
        for k in tag:
            for i in tag[k]:
                if i.lower() in s:
                    return k
        return None

    def _title_match(s, k):
        for i in tag[k]:
            if i.lower() in s.lower():
                return True
        return False

    with open('{}/{}.annot'.format(gplAnnotDir, gpl), 'r') as rf:
        map = dict()
        begin = False
        col_id = False
        for line in rf.readlines():
            if line.startswith('#') or line.startswith(' '):
                continue
            elif not col_id:
                titles = [i.strip() for i in line.split('\t')]
                for i in range(0, len(titles)):
                    if _title_match(titles[i], key):
                        col_id = i
            else:
                item = [i.strip() for i in line.split('\t')]
                if col_id and len(item) > col_id and item[col_id]:
                    map[item[0]] = item[col_id]

    m = dict()
    for probe in map:
        s = map[probe]
        if '///' in s:
            s = [i.strip() for i in s.split('///')]
        else:
            s = [s]
        for i in s:
            if i not in m:
                m[i] = dict()
                m[i][probe] = len(s)

    with open('{}/{}_{}2probe.json'.format(gplAnnotDir, gpl, key), 'w') as wf:
        json.dump(m, wf, indent=2)
    with open('{}/{}_probe2{}.json'.format(gplAnnotDir, gpl, key), 'w') as wf:
        json.dump(map, wf, indent=2)
    if log:
        logging.warning('Cannot map probe to {} in {}!'.format(key, gpl))

    return gpl


def getExprTable(gse, gpl, exprDir, gsmList=None, matrixOnly=False,
                log=True, maxTrial=4, downloadMethod='wget',
                runningR=True, rScriptDir=None,
                gseRawDir=None, celDir=None, matrixDir=None):

    with tempfile.TemporaryDirectory() as tmpDir:

        if matrixOnly:
            if not matrixDir:
                matrixDir = '{}/matrix'.format(tmpDir)
            mkdir(matrixDir)
            r = getGseMatrix(gse, gpl, matrixDir, log, maxTrial, downloadMethod)
            if r[0]:
                filename = '{}/{}'.format(matrixDir, r[3])
                matrixFile = extractGzipFile(filename)
                expr = matrix2expr(matrixFile)
                if gsmList:
                    gsmList = {gsm.upper() for gsm in gsmList}
                    expr = {oligo: {sample: expr[oligo][sample]
                                    for sample in expr[oligo]
                                    if sample in gsmList}
                            for oligo in expr}
                gsmList = sorted(expr[list(expr.keys())[0]].keys())
                s = '\t'.join(gsmList)
                s += '\n'
                s += '\n'.join(['{}\t{}'.format(oligo, '\t'.
                                                join([expr[oligo][sample]
                                                      for sample in gsmList]))
                                for oligo in sorted(expr.keys())])
                with open('{}/{}_{}_expr.tsv'.format(exprDir, gse, gpl), 'w') as wf:
                    print(s, end='', file=wf)
                return True, gse, gpl, '{}/{}_{}_expr.tsv'.format(exprDir, gse, gpl)
            else:
                return False, gse, gpl, ''

        else:
            if not gseRawDir:
                gseRawDir = '{}/gseRaw'.format(tmpDir)
            mkdir(gseRawDir)
            if not celDir:
                celDir = '{}/cel'.format(tmpDir)
            mkdir(gseRawDir)
            celDir = '{}/{}_{}'.format(celDir, gse, gpl)
            mkdir(celDir)

            r = getGseRaw(gse, gseRawDir, log, maxTrial, downloadMethod)
            if r[0]:
                filename = '{}/{}'.format(gseRawDir, r[2])
                with tarfile.open(filename, 'r') as tf:
                    def is_within_directory(directory, target):
                        
                        abs_directory = os.path.abspath(directory)
                        abs_target = os.path.abspath(target)
                    
                        prefix = os.path.commonprefix([abs_directory, abs_target])
                        
                        return prefix == abs_directory
                    
                    def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
                    
                        for member in tar.getmembers():
                            member_path = os.path.join(path, member.name)
                            if not is_within_directory(path, member_path):
                                raise Exception("Attempted Path Traversal in Tar File")
                    
                        tar.extractall(path, members, numeric_owner=numeric_owner) 
                        
                    
                    safe_extract(tf, celDir)
                for file in os.listdir(celDir):
                    gsm = file.split('.')[0].split('_')[0].split('-')[0].upper()
                    if gsmList:
                        if gsm not in gsmList:
                            os.remove('{}/{}'.format(celDir, file))
                            continue
                    f = extractGzipFile('{}/{}'.format(celDir, file), rm=True)
                    os.rename(f, '{}/{}.CEL'.format(celDir, gsm))

                if not rScriptDir:
                    rScriptDir = '{}/rScript'.format(exprDir)
                mkdir(rScriptDir)
                R = 'rm(list = ls())\nlibrary("affyPLM")\nlibrary("affy")\n\n'
                R += 'memory.limit(16000)\n\n'
                R += 'setwd("{}")\n'.format(celDir)
                R += 'raw.set <- ReadAffy()\nrma.data <- rma(raw.set)\n'
                R += '%s <- exprs(rma.data)\n' % gse
                R += 'write.table({}, "{}/{}_{}_expr.tsv", sep="\\t")\n'.\
                    format(gse, exprDir, gse, gpl)
                R += 'rm(raw.set, rma.data)\ngc()\n\n'
                with open('{}/{}_{}_expr.R'.format(rScriptDir, gse, gpl), 'w') \
                        as wf:
                    print(R, file=wf)
                if runningR:
                    if log:
                        logging.info('Begin to calculate expression by R.')
                        cmd = 'Rscript {}/{}_{}_expr.R > {}/{}_{}_expr.R.log'.\
                            format(rScriptDir, gse, gpl, rScriptDir, gse, gpl)
                        if os.system(cmd) == 0:
                            logging.info('Calculate expression Successfully.')
                            return True, gse, gpl, '{}/{}_{}_expr.tsv'.\
                                format(exprDir, gse, gpl)
                        else:
                            logging.error('Calculate expression Failed.')
                            return False, gse, gpl, ''
            else:
                return False, gse, gpl, ''


def getExprTables(gse_gpl_gsmList_matrixOnly_List, exprDir,
                  log=True, maxTrial=4, thread=4, downloadMethod='wget',
                  runningR=True, rScriptDir=None,
                  gseRawDir=None, celDir=None, matrixDir=None):
    p = [(item['gse'], item['gpl'], exprDir, item['gsmList'], item['matrixOnly'],
          log, maxTrial, downloadMethod, runningR, rScriptDir,
          gseRawDir, celDir, matrixDir)
         for item in gse_gpl_gsmList_matrixOnly_List]
    with multiprocessing.Pool(processes=thread) as pool:
        results = pool.starmap(getExprTable, p)
    return results


def probeExpr(settingFile):
    settings = config(settingFile)
    gsmList = extractClinical(settings['datasets information'],
                              settings['directory containing clinical csv files'],
                              settings['simplified clinical data directory'],
                              settings['clinical phenotype'],
                              settings['log'])
    for item in settings['datasets information']:
        item['gsmList'] = gsmList[(item['gse'], item['gpl'])]
        if 'matrixOnly' not in item:
            item['matrixOnly'] = False
    with open('{}/probeExpr.config.json'.format(settings['expression table directory']), 'w') as wf:
        json.dump(settings, wf, indent=2)
    return getExprTables(gse_gpl_gsmList_matrixOnly_List=settings['datasets information'],
                         exprDir=settings['expression table directory'],
                         log=settings['log'],
                         maxTrial=settings['max trail in downloading files'],
                         thread=settings['thread'],
                         downloadMethod=settings['download method'],
                         runningR=settings['run R in python'],
                         rScriptDir=settings['R script directory'],
                         gseRawDir=settings['GSE raw directory'],
                         celDir=settings['CEL file directory'],
                         matrixDir=settings['GSE matrix directory'])


def probeExpr2geneExpr(gse, gpl, exprDir, annotDir, geneExprDir,
                   key='symbol', method='average', multipleProbe='drop',
                   key_mapping=None, dropUnexpercted='sample',
                   log=True, annotMapped=True, maxTrial=4):
    if not annotMapped:
        getGplAnnotTable(gpl, annotDir, log, maxTrial)
        annot2Map(gpl, annotDir, key, log, key_mapping)
    with open('{}/{}_{}_expr.tsv'.format(exprDir, gse, gpl), 'r') as rf:
        lines = rf.readlines()
        samples = [i.strip().strip('\'"')
                   for i in lines[0].split('\t')]
        expr = dict()
        unexperctedValues = list()
        for line in lines[1:]:
            items = line.strip('\n ').split('\t')
            try:
                expr[items[0].strip('"\'')] = [float(i) for i in items[1:]]
            except:
                probe = items[0].strip('"\'')
                expr[probe] = list()
                for i in range(1, len(items[1:])):
                    try:
                        expr[probe].append(float(items[i]))
                    except:
                        unexperctedValues.append((probe, i-1))
                        expr[probe].append('')
        unexperctedSamples = list({i[1] for i in unexperctedValues})
        unexperctedProbes = list({i[0] for i in unexperctedValues})
        if dropUnexpercted == 'sample':
            samples = [samples[i]
                       for i in range(0, len(samples))
                       if i not in unexperctedSamples]
            expr = {probe: [expr[probe][i]
                     for i in range(0, len(expr[probe]))
                     if i not in unexperctedSamples]
                    for probe in expr}
        elif dropUnexpercted == 'probe':
            expr = {probe: expr[probe]
                    for probe in expr
                    if probe not in unexperctedProbes}

    with open('{}/{}_{}2probe.json'.format(annotDir, gpl, key), 'r') as rf:
        gene2probe = json.load(rf)
    geneExpr = dict()
    for gene in gene2probe:
        probes = gene2probe[gene]
        if multipleProbe == 'drop':
            probes = {probe: probes[probe]
                      for probe in probes if probes[probe] == 1}
        elif multipleProbe == 'ignore':
            probes = {probe: 1 for probe in probes}
        if len(probes) == 1:
            probe = list(probes.keys())[0]
            geneExpr[gene] = [i / probes[probe] for i in expr[probe]]
        elif len(probes) > 1:
            probesExpr = [[i / probes[probe] for i in expr[probe]]
                          for probe in probes]
            probesExpr = transpose(probeExpr)
            if method == 'average':
                geneExpr[gene] = [statistics.mean(i) for i in probesExpr]
            if method == 'median':
                geneExpr[gene] = [statistics.median(i) for i in probesExpr]
            if method == 'maximum':
                geneExpr[gene] = [max(i) for i in probesExpr]
    with open('{}/{}_{}_expr.tsv'.format(geneExprDir, gse, gpl), 'w') as wf:
        s = '\t{}'.format('\t'.join(samples))
        s += '\n'
        s += '\n'.join(['{}\t{}'.format(gene, '\t'.join([str(i) for i in geneExpr[gene]]))
                        for gene in geneExpr])
        print(s, file=wf, end='')
    return gse, gpl, '{}/{}_{}_expr.tsv'.format(geneExprDir, gse, gpl)


def geneExpr(settingFile, probeExprDone=False):
    settings = config(settingFile)
    if not probeExprDone:
        probeExpr(settingFile)
    if settings['probe mapping to gene']:
        gplList = {info['gpl'] for info in settings['datasets information']}
        getGplAnnotTableList(gplList,
                             directory=settings['probe annotation directory'],
                             log=settings['log'],
                             maxTrial=settings['max trail in downloading files'],
                             thread=settings['thread'])
        for gpl in gplList:
            annot2Map(gpl, settings['probe annotation directory'],
                      settings['probe mapping to gene keyword'], settings['log'])
        p = [(info['gse'], info['gpl'], settings['expression table directory'],
              settings['probe annotation directory'],
              settings['gene expression table directory'],
              settings['probe mapping to gene keyword'],
              settings['probe mapping to gene method'],
              settings['multiple mapped probe'])
             for info in settings['datasets information']]
        with multiprocessing.Pool(processes=settings['thread']) as pool:
            results = pool.starmap(probeExpr2geneExpr, p)
