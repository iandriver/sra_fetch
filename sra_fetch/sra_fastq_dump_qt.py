import fnmatch
import os
import subprocess
import GEOparse
import sys
import json
import re
import pandas as pd
import utils
import platform
from collections import OrderedDict
import boto3
import botocore
from Bio import Entrez
from urllib.error import HTTPError
import time


def read_qt5_input(qt5_data):
    gs_text = qt5_data.sra.text()

    if gs_text[0:2] =='GS' and gs_text[3:].isdigit() and len(gs_text[3:]) in [5,6]:
        if gs_text[2] == 'E':
            series = True
        else:
            series =False
    else:
        print(gs_text)
        sys.exit('Please enter a valid GEO series or sample ID (GSE or GSM followed by 6 digits)')
    s3_text = qt5_data.s3.text()
    if s3_text != '':
        if s3_text[0:5] != "s3://":
            sys.exit('Please enter a valid s3 bucket begining with s3://')
    directory_path = qt5_data.directoryLabel.text()
    if directory_path == '':
        cwd = os.getcwd()
        output_path = os.path.join(cwd, gs_text)
    else:
        output_path = os.path.join(directory_path, gs_text)
    email = qt5_data.email.text()
    if email == '':
        sys.exit('Please enter a valid email.')
    local_files_only = qt5_data.local_manifest
    s3_files_only = qt5_data.s3_manifest
    if local_files_only and s3_files_only:
        sys.exit('Please build your manifest from local or s3 not both.')
    return gs_text, series, s3_text, output_path, email, local_files_only, s3_files_only

def read_args_input(args):
    gs_text = args.sra

    if gs_text[0:2] =='GS' and gs_text[3:].isdigit() and len(gs_text[3:]) in [5,6]:
        if gs_text[2] == 'E':
            series = True
        else:
            series =False
    else:
        print(gs_text)
        sys.exit('Please enter a valid GEO series or sample ID (GSE or GSM followed by 6 digits)')
    s3_text = args.s3
    if s3_text != '':
        if s3_text[0:5] != "s3://":
            sys.exit('Please enter a valid s3 bucket begining with s3://')
    directory_path = args.directoryLabel
    if directory_path == '':
        cwd = os.getcwd()
        output_path = os.path.join(cwd, gs_text)
    else:
        output_path = os.path.join(directory_path, gs_text)
    email = args.email
    if email == '':
        sys.exit('Please enter a valid email.')
    local_files_only = args.local_manifest
    s3_files_only = args.s3_manifest
    if local_files_only and s3_files_only:
        sys.exit('Please build your manifest from local or s3 not both.')
    process_num = args.processes
    split_num = str(args.threads)
    return split_num, process_num, gs_text, series, s3_text, output_path, email, local_files_only, s3_files_only


def make_sra_sub_dir(gsm, directory):
    # make the directory
    if platform.system() == "Windows":
        name_regex = r'[\s\*\?\(\),\.\:\%\|\"\<\>]'
    else:
        name_regex = r'[\s\*\?\(\),\.;]'
    directory_path = os.path.abspath(os.path.join(directory, "%s_%s" % (gsm.get_accession(),re.sub(name_regex, '_', gsm.metadata['title'][0]) # the directory name cannot contain many of the signs
    )))
    return directory_path

def gsm_query_to_df(query, email):
    # retrieve IDs for given SRX
    # check if the e-mail is more or less not a total crap
    Entrez.email = email
    if not (Entrez.email is not None and '@' in email and email != '' and '.' in email):
        raise Exception('You have to provide valid e-mail')
    searchdata = Entrez.esearch(db='sra', term=query, usehistory='y', retmode='json')
    answer = json.loads(searchdata.read())
    ids = answer["esearchresult"]["idlist"]
    assert len(ids) == 1, "There should be one and only one ID per SRX"

    # using ID fetch the info
    number_of_trials = 10
    wait_time = 30
    for trial in range(number_of_trials):
        try:
            results = Entrez.efetch(db="sra", id=ids[0], rettype="runinfo", retmode="text").read()
            break
        except HTTPError as httperr:
            if "502" in str(httperr) or "429" in str(httperr):
                sys.stderr.write("Error: %s, trial %i out of %i, waiting for %i seconds." % (str(httperr),
                                                                                             trial,
                                                                                             number_of_trials,
                                                                                             wait_time))
                time.sleep(wait_time)
            else:
                raise httperr
    df = pd.DataFrame([i.split(',') for i in results.split('\n') if i != ''][1:], columns = [i.split(',') for i in results.split('\n') if i != ''][0])
    return(df)

def download_SRA(gsm, queries, split_num, email, metadata_key='auto', directory='./', filetype='sra', aspera=False, keep_sra=False):
    """Download RAW data as SRA file to the sample directory created ad hoc
    or the directory specified by the parameter. The sample has to come from
    sequencing eg. mRNA-seq, CLIP etc.

    An important parameter is a download_type. By default an SRA is accessed by FTP and
    such file is downloaded. This does not require additional libraries. However in order
    to produce FASTA of FASTQ files one would need to use SRA-Toolkit. Thus, it is assumed
    that this library is already installed or it will be installed in the near future. One
    can immediately specify the download type to fasta or fastq.

    :param email: an email (any) - required by NCBI for access
    :param directory: The directory to which download the data. By default current directory is used
    :param filetype: can be sra, fasta, or fastq - for fasta or fastq SRA-Toolkit need to be installed
    :param aspera: bool - use Aspera to download samples, defaults to False
    :param keep_sra: bool - keep SRA files after download, defaults to False

    """

    # Check download filetype
    filetype = filetype.lower()
    if filetype not in ["sra", "fastq", "fasta"]:
        raise Exception("Unknown type to downlod: %s. Use sra, fastq or fasta." % filetype)

    # Setup the query
    ftpaddres = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/{range_subdir}/{record_dir}/{file_dir}/{file_dir}.sra"


    # check if the e-mail is more or less not a total crap
    Entrez.email = email
    if not (Entrez.email is not None and '@' in email and email != '' and '.' in email):
        raise Exception('You have to provide valid e-mail')


    for query in queries:
        df = gsm_query_to_df(query, email)
        # check it first
        try:
            df['download_path']
        except KeyError as e:
            sys.stderr.write('KeyError: ' + str(e) + '\n')
            sys.stderr.write(str(results) + '\n')

        directory_path = make_sra_sub_dir(gsm, directory)
        utils.mkdir_p(os.path.abspath(directory_path))
        count = 0
        for path in df['download_path']:
            count+=1
            print(count)
            if df['LibraryStrategy'][0] == 'RNA-Seq':
                sra_run = path.split("/")[-1]
                print("Analysing %s" % sra_run)
                url = ftpaddres.format(range_subdir=query[:6],
                                           record_dir=query,
                                           file_dir=sra_run)
                filepath = os.path.abspath(os.path.join(directory_path, "%s.sra" % sra_run))
                df.to_csv(os.path.abspath(os.path.join(directory_path, "%s_metadata.txt" % sra_run)), sep='\t')

                if filetype in ["fasta", "fastq","sra"]:
                    ftype = ""
                    if filetype == "fasta":
                        ftype = " --fasta "
                    cmd = "parallel-fastq-dump --sra-id %s --threads %s --outdir %s --split-files --gzip"
                    cmd = cmd % (sra_run, split_num, directory_path)
                    return (cmd, sra_run, directory_path)

def make_manifest(gs_text, series, s3_text, output_path, email, response, local_files_only, s3_files_only):
    # check if the e-mail is more or less not a total crap
    Entrez.email = email
    if not (Entrez.email is not None and '@' in email and email != '' and '.' in email):
        raise Exception('You have to provide valid e-mail')
    gs_object = GEOparse.get_GEO(geo=gs_text)
    gsms_to_use = gs_object.gsms.values()
    gsm_names = [gsm.name for gsm in gsms_to_use]
    sra_dict_by_gsm = dict(zip(gsm_names,[gsm.relations['SRA'][0].split("=")[-1] for gsm in gsms_to_use]))

    filetype='sra'
    sra_already_done = []
    data_manifest_dict = OrderedDict()
    if local_files_only:
        for root, dirnames, filenames in os.walk(output_path):
            metadata_exits = False
            for f in filenames:
                if '_metadata' in f:
                    metadata_exits = True
                    metadata_path = os.path.join(root,f)
            if metadata_exits:
                metadata_df = pd.read_csv(metadata_path, sep='\t', header=0)
                srr_num = metadata_df['Run'][0]
                if metadata_df['LibraryLayout'][0] == 'PAIRED':
                    fastq_1 = srr_num+'_1.fastq.gz'
                    fastq_2 = srr_num+'_2.fastq.gz'
                    fastqs_to_check = [fastq_1,fastq_2]
                else:
                    fastqs_to_check = [srr_num+'_1.fastq.gz']
                data_manifest_dict[sra_dict_by_gsm[metadata_df['SampleName'][0]]] = fastqs_to_check

    elif s3_files_only:
        data_manifest_dict = OrderedDict()
        bucket_name = s3_text.split('/')[2]
        folder_path = '/'.join(s3_text.split('/')[3:])
        s3 = boto3.resource('s3')
        queries = []
        for gsm in gsms_to_use:
            try:
                for sra in gsm.relations['SRA']:
                    query = sra.split("=")[-1]
                    assert 'SRX' in query, "Sample looks like it is not SRA: %s" % query
                    print("Query: %s" % query)
                    sys.stderr.write("Downloading %s files for %s series\n" % (filetype, gsm.name))
                    queries.append(query)
            except KeyError:
                raise NoSRARelationException('No relation called SRA for %s' % gsm.get_accession())
        for query in queries:
            metadata_df = gsm_query_to_df(query, email)
            if metadata_df['LibraryStrategy'][0] == 'RNA-Seq':
                srr_num = metadata_df['Run'][0]
                if metadata_df['LibraryLayout'][0] == 'PAIRED':
                    fastq_1 = srr_num+'_1.fastq.gz'
                    fastq_2 = srr_num+'_2.fastq.gz'
                    fastqs_to_check = [fastq_1,fastq_2]
                else:
                    fastqs_to_check = [srr_num+'_1.fastq.gz']
                for f in fastqs_to_check:
                    try:
                        s3.Object(bucket_name, folder_path+'/'+f).load()
                        file_exists = True
                        print(f, 'found in s3')
                    except botocore.exceptions.ClientError as e:
                        if e.response['Error']['Code'] == "404":
                            file_exists = False
                            pass
                if file_exists:
                    data_manifest_dict[sra_dict_by_gsm[metadata_df['SampleName'][0]]] = fastqs_to_check
    data_manifest_df = pd.DataFrame.from_dict(data_manifest_dict, orient='index')
    if len(data_manifest_df.columns) == 2:
        data_manifest_df.rename(columns={0:'read1', 1:'read2'}, inplace=True)
    else:
        data_manifest_df.rename(columns={0:'read1'}, inplace=True)
    data_manifest_name = gs_text+'_data_manifest.txt'
    data_manifest_path = output_path+'/'+data_manifest_name
    data_manifest_df.to_csv(data_manifest_path, sep='\t')

    if s3_text != '':
        bucket_name = s3_text.split('/')[2]
        folder_path = '/'.join(s3_text.split('/')[3:])
        s3 = boto3.resource('s3')
        s3.Bucket(bucket_name).upload_file(data_manifest_path, folder_path+'/'+data_manifest_name)

def sra_series_fetch(split_num, process_num, gs_text, series, s3_text, output_path, email, response, local_files_only, s3_files_only):
    gs_object = GEOparse.get_GEO(geo=gs_text)
    gsms_to_use = gs_object.gsms.values()
    gsm_names = [gsm.name for gsm in gsms_to_use]
    sra_dict_by_gsm = dict(zip(gsm_names,[gsm.relations['SRA'][0].split("=")[-1] for gsm in gsms_to_use]))

    filetype='sra'
    sra_already_done = []
    if response == 'Y':
        for root, dirnames, filenames in os.walk(output_path):
            fastq_exits = False
            metadata_exits = False
            for f in filenames:
                if 'fastq.gz' in f:
                    fastq_exits = True
                if '_metadata' in f:
                    metadata_exits = True
                    metadata_path = os.path.join(root,f)
            if fastq_exits:
                dir_done = os.path.basename(root).split('_')[1]
                if sra_dict_by_gsm[dir_done] not in sra_already_done:
                    sra_already_done.append(sra_dict_by_gsm[dir_done])
                    print('Already downloaded:', sra_already_done)

    cmd_list = []

    for gsm in gsms_to_use:
        queries = []
        try:
            for sra in gsm.relations['SRA']:
                query = sra.split("=")[-1]
                assert 'SRX' in query, "Sample looks like it is not SRA: %s" % query
                print("Query: %s" % query)


                queries.append(query)
        except KeyError:
            raise NoSRARelationException('No relation called SRA for %s' % gsm.get_accession())
        download_sra_cmd = download_SRA(gsm, queries, split_num, email=email, filetype=filetype, directory=output_path)
        cmd_list.append(download_sra_cmd)
    processes = set()
    max_processes = min(process_num, len(cmd_list))
    for c in cmd_list:
        cmd = c[0]
        print(cmd)
        name = cmd.split(' ')[2]
        ouput_path = cmd.split(' ')[6]
        folder_name = os.path.basename(cmd.split(' ')[6])

        s3_path = s3_text+'/'+folder_name
        if s3_text != '':
            bucket_name = s3_text.split('/')[2]
            folder_path = os.path.basename(c[2].strip('/'))
            s3 = boto3.resource('s3')


            fastq_1 = c[1]+'_1.fastq.gz'
            fastq_2 = c[1]+'_2.fastq.gz'
            fastqs_to_check = [fastq_1,fastq_2]
            file_found = 0
            for f in fastqs_to_check:
                try:
                    s3.Object(bucket_name, folder_path+'/'+f).load()
                    file_found +=1
                    print(os.path.join(bucket_name,folder_path,f), 'found')
                except botocore.exceptions.ClientError as e:
                    if e.response['Error']['Code'] == "404":
                        pass
        if file_found > 0:
            print(c[2], 'Found in s3')
            run=False
        else:
            if s3_text:
                c_run = cmd +' && aws s3 sync %s %s && rm -rf %s' %(cmd.split(' ')[6],s3_path, cmd.split(' ')[6])
            else:
                c_run = cmd
            print("Downloading %s \n" %(name))
            print(c_run)
            run=True
            processes.add(subprocess.Popen([c_run, name], shell=True))
            if len(processes) >= max_processes:
                os.wait()
                processes.difference_update([
                p for p in processes if p.poll() is not None])


    if run:
        for p in processes:
            if p.poll() is None:
                p.wait()

    if not local_files_only and not s3_files_only:
        if s3_text == '':
            local_files_only=True
            make_manifest(gs_text, series, s3_text, output_path, email, response, local_files_only, s3_files_only)
        else:
            s3_files_only=True
            make_manifest(gs_text, series, s3_text, output_path, email, response, local_files_only, s3_files_only)
