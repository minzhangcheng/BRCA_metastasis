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


import os
import gzip
import shutil
import multiprocessing
import ftplib
import logging
import json


class DownloadError(Exception):
    pass


class ConnectionError(Exception):
    pass


def mkdir(dir):
    try:
        os.mkdir(dir)
    except:
        pass


def extractGzipFile(filename, rm=True):
    file = filename.rsplit('.')[0]
    with open('%s' % file, 'wb') as wf:
        with gzip.open(filename, 'rb') as rf:
            shutil.copyfileobj(rf, wf)
    if rm:
        os.remove(filename)
    return file


def extractGzipFiles(filenameList, rm=True, thread=4):
    p = [(i, rm) for i in filenameList]
    if thread > 1:
        with multiprocessing.Pool(processes=thread) as pool:
            pool.starmap(extractGzipFile, p)
    else:
        for filename in filenameList:
            extractGzipFile(filename)


def wgetDwonload(url, filename='', directory='', log=False):
    if not directory:
        if filename:
            cmd = 'wget -c -nv -t 0 -O {} {}'.format(filename, url)
        else:
             cmd = 'wget -c -nv -t 0 {}'.format(url)
    else:
        if filename:
            cmd = 'wget -c -nv -t 0 -O {}/{} {}'.\
                format(directory, filename, url)
        else:
            cmd = 'wget -c -nv -t 0 -P {} {}'.\
                format(directory, url)
    if log:
        msg = 'Download Successfully.\n\t{}\n\t{}'.format(url, cmd)
        logging.info(msg)
    if os.system(cmd) != 0:
        raise DownloadError('Download error in runing connand "{}"'.format(cmd))
    else:
        return True


def listFtpFiles(url):

    """
    List the files under the ftp server.
    :param url: for example: 'ftp://ftp.ncbi.nlm.nih.gov/'
    :return: a list of filenames
    """

    if url.startswith('ftp://'):
        url = url[6:]
    s = url.split('/', maxsplit=1)
    if len(s) > 1:
        ftpUrl, wd = s[0], s[1]
    else:
        ftpUrl = url
        wd = ''
    try:
        with ftplib.FTP(ftpUrl) as ftp:
            ftp.login()
            if wd:
                ftp.cwd(wd)
            fileList = ftp.nlst()
            return fileList
    except Exception as e:
        ftp.close()
        raise ConnectionError('Connect error in list ftp: "{}"\n{}'.format(url, e))


def ftpDownload(url, filename='', directory='', log=False, check=True):
    url = 'adafadf'
    if url.startswith('ftp://'):
        url = url[6:]
    s = url.split('/', maxsplit=1)
    if len(s) > 1:
        ftpUrl, wd = s[0], s[1]
    else:
        raise DownloadError('Download error in manipulating url: "{}"'.format(url))
    if '/' in wd:
        wd, file = ''.rsplit('/', maxsplit=1)
    else:
        wd, file = '', wd
    if directory:
        if filename:
            target = '{}/{}'.format(directory, filename)
        else:
            target = '{}/{}'.format(directory, file)
    else:
        if filename:
            target = filename
        else:
            target = file
    try:
        with ftplib.FTP(ftpUrl) as ftp:
            ftp.login()
            if wd:
                ftp.cwd(wd)
            if check:
                fileList = ftp.nlst()
                if file not in fileList:
                    raise DownloadError('Download error, file not found in "{}"'.format(url))
            ftp.retrbinary("RETR %s" % file,
                           open(target, 'wb').write)
        if log:
            msg = 'Download Successfully.\n\t{} -> {}'.format(url, target)
            logging.info(msg)
        return True
    except Exception as e:
        ftp.close()
        raise DownloadError('Download error in "{}"\n{}'.format(url, e))


def download(url, filename='', directory='', maxTrial=4, downloadMethod='wget', log=False, check=True):
    try:
        if downloadMethod == 'wget':
            return wgetDwonload(url, filename, directory)
        elif downloadMethod == 'ftplib':
            return ftpDownload(url, filename, directory, check)
    except DownloadError as e:
        if maxTrial > 0:
            if log:
                logging.warning('Download failed and retry: "{}"'.format(url))
            return download(url, filename, directory,
                            maxTrial - 1, downloadMethod, check)
        else:
            if log:
                logging.error('Download failed and terminated: "{}"'.format(url))
            return False


def transpose(table):
    """
    Transpose a table (list of list).
    transpose(table)
    """
    return [[r[col] for r in table] for col in range(len(table[0]))]


def config(settingFile):
    with open(settingFile, 'r') as rf:
        settings = json.load(rf)

    defaultParameter = {
        'clinical phenotype': 'metastasis',
        'thread': 4,
        'max trail in downloading files': 4,
        'run R in python': True,
        'probe mapping to gene': True,
        'probe mapping to gene method': 'median',
        'multiple mapped probe': 'drop',
        'probe mapping to gene keyword': 'symbol',
        'download method': 'wget',
        'log': True,
        'metaQC': True,
        'meta': True,
        'mapping to gene in meta': False,
        'mapping to gene in valid': True,
        'meta method': 'mDES',
        'p-value cutoff': [0.05, 0.01, 0.001, 0.0001],
        'effect size cutoff': [1.5, 2, 2.5, 3],
        'max different expression gene included': [25, 50, 100, 250, 500, 1000],
        'max sample in heatmap': [25, 50, 100, 200],
        'randomized sampling': 'stratified sampling',
        'randomized times': 10,
        'outliers replace': True
    }
    defaultDirectory = {
        'directory containing clinical csv files': 'gseClinicalCsv',
        'gmt File': '/home/minzhang/meta/c2.all.v6.2.symbols.gmt',
        'simplified clinical data directory': 'clinical',
        'GPL soft file directory': 'gplSoft',
        'probe annotation directory': 'annot',
        'GSE raw directory': 'gseRaw',
        'GSE matrix directory': 'gseMatrix',
        'CEL file directory': 'cel',
        'expression table directory': 'expr',
        'R script directory': 'R',
        'gene expression table directory': 'geneExpr',
        'metaQC directory': 'metaqc',
        'meta directory': 'meta',
    }
    defaultDatasetParameter = {
        'dropNA': True,
        'metaQC': False,
        'meta': False,
        'valid': True,
        'matrixOnly': False

    }

    for key in defaultParameter:
        if key not in settings:
            settings[key] = defaultParameter[key]
    mkdir(settings['root directory'])
    for key in defaultDirectory:
        if key not in settings or not settings[key]:
            settings[key] = '{}/{}'.format(settings['root directory'],
                                           defaultDirectory[key])
        mkdir(settings[key])

    if settings['log']:
        logging.basicConfig(filename='{}/log.txt'.format(settings['root directory']),
                            level=logging.INFO)

    for key in defaultDatasetParameter:
        for info in settings['datasets information']:
            if key not in info:
                info[key] = defaultDatasetParameter[key]

    if 'database' not in settings or not settings['database']:
        settings['database'] = 'sqlite:///{}/data.sqlite3'.format(
            settings['root directory']
        )

    return settings

