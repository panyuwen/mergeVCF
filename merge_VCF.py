# -*- coding: utf-8 -*-

import numpy as np
import pandas
import gzip
import argparse
import itertools
import socket
import os
import time
import sys

def check(filename):
    with gzip.open(filename) as f:
        headerline = 0
        line = f.readline().decode("utf-8")
        while line.startswith == "##":
            headerline += 1
            line = f.readline().decode("utf-8")

    samplelist = line.strip().split('\t')[9:]
    sitelist = pandas.read_csv(filename,sep='\t',skiprows=range(headerline),usecols=['#CHROM','POS','REF','ALT'],dtype={'#CHROM':'category','POS':'int32','REF':'category','ALT':'category'})
    sitelist.drop_duplicates(['POS'],inplace=True,keep=False)
    sitelist = list(sitelist.apply(lambda x: '{}:{}:{}:{}'.format(str(x['#CHROM']),str(x['POS']),str(x['REF']),str(x['ALT'])),axis=1))

    return headerline, samplelist, sitelist

def merge_sample(filelist,output):
    headerline = {}; samplelist = {}; sitelist = {}
    
    for filename in filelist:
        headerline[filename], samplelist[filename], sitelist[filename] = check(filename)

    '''
    filepairlist = list(itertools.combinations(filelist,2))
    for f1, f2 in filepairlist:
        overlap = list(set(samplelist[f1]) & set(samplelist[f2]))
        if len(overlap) >0:
            samplelist[f2] = list(set(samplelist[f2]) - set(overlap))
        else:
            pass
    '''
    
    site = sitelist[filelist[0]]
    for filename in filelist[1:]:
        site = list(set(site) & set(sitelist[filename]))

    for filename in filelist:
        datatype = dict(zip(['#CHROM','POS','REF','ALT','QUAL','FILTER','INFO','FORMAT']+samplelist[filename], ['category','int32']+['category']*(len(samplelist[filename])+6)))
        if filename == filelist[0]:
            data = pandas.read_csv(filename,sep='\t',dtype=datatype,skiprows=range(headerline[filename]))
            data['tmp'] = data.apply(lambda x: '{}:{}:{}:{}'.format(str(x['#CHROM']),str(x['POS']),str(x['REF']),str(x['ALT'])),axis=1)
            data = data[data['tmp'].isin(site)]
            data.index = data['tmp']
            site = list(data.index)
            data.drop('tmp',axis=1,inplace=True)
            datalist = [data]
        else:
            tmp = pandas.read_csv(filename,sep='\t',dtype=datatype,skiprows=range(headerline[filename]))
            tmp.index = tmp.apply(lambda x: '{}:{}:{}:{}'.format(str(x['#CHROM']),str(x['POS']),str(x['REF']),str(x['ALT'])),axis=1)
            tmp = tmp.loc[site]
            tmp = tmp.iloc[:,9:]
            datalist += [tmp]
    
    data = pandas.concat(datalist,axis=1)

    data.to_csv(output,sep='\t',index=None,compression='gzip',mode='a')

def merge_site(filelist,output):
    headerline = {}; samplelist = {}; sitelist = {}
    
    for filename in filelist:
        headerline[filename], samplelist[filename], sitelist[filename] = check(filename)

    sample = samplelist[filelist[0]]
    for filename in filelist[1:]:
        sample = list(set(sample) & set(samplelist[filename]))

    datatype = dict(zip(['#CHROM','POS','REF','ALT','QUAL','FILTER','INFO','FORMAT']+sample, ['category','int32']+['category']*(len(sample)+6)))
    for filename in filelist:
        if filename == filelist[0]:
            data = pandas.read_csv(filename,sep='\t',usecols=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']+sample,dtype=datatype,skiprows=range(headerline[filename]))
            orderlist = list(data.columns)
            data.to_csv(output,sep='\t',index=None,compression='gzip',mode='a')
        else:
            data = pandas.read_csv(filename,sep='\t',usecols=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']+sample,dtype=datatype,skiprows=range(headerline[filename]))
            data = data[orderlist]
            data.to_csv(output,sep='\t',index=None,header=None,compression='gzip',mode='a')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--filelist", type=str, required = True, \
                        help="list of /path/to/phased.vcf.gz, headers in the 1st file will be kept in the output")
    parser.add_argument("--merge", type=str, required = True, choices=['sample','site'], \
                        help="merge samples (same sites, different samples) or sites (same samples, different sites)")
    parser.add_argument("--out", type=str, required = False, default='out.vcf.gz', \
                        help="/out/file/name")
    args = parser.parse_args()

    with open(args.out+'.logfile','w') as log:
        log.write('python {}\n'.format(sys.argv[0]))
        log.write('{}--filelist  {}\n'.format(' '*8, args.filelist))
        log.write('{}--merge     {}\n'.format(' '*8, args.merge))
        log.write('{}--out       {}\n\n'.format(' '*8, args.out))
        
        log.write('Hostname: '+socket.gethostname()+'\n')
        log.write('Working directory: '+os.getcwd()+'\n')
        log.write('Start time: '+time.strftime("%Y-%m-%d %X",time.localtime())+'\n\n')

    files = list(pandas.read_csv(args.filelist,header=None)[0])

    with gzip.open(args.out,'wb') as fout:
        with gzip.open(files[0]) as fin:
            line = fin.readline().decode("utf-8")
            while line.startswith == "##":
                fout.write(line.encode())
                line = fin.readline().decode("utf-8")

    if args.merge == 'sample':
        merge_sample(files, args.out)
    else:
        merge_site(files, args.out)

    with open(args.out+'.logfile','a') as log:
        log.write("Done.\n")
        log.write('End time: '+time.strftime("%Y-%m-%d %X",time.localtime())+'\n\n')

if __name__ == '__main__':
    main()
