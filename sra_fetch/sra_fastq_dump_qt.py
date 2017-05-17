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
    return gs_text, series, s3_text, output_path, email

def download_SRA(gsm, email, metadata_key='auto', directory='./', filetype='sra', aspera=False, keep_sra=False):
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
    from Bio import Entrez
    # Check download filetype
    filetype = filetype.lower()
    if filetype not in ["sra", "fastq", "fasta"]:
        raise Exception("Unknown type to downlod: %s. Use sra, fastq or fasta." % filetype)

    # Setup the query
    ftpaddres = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/{range_subdir}/{record_dir}/{file_dir}/{file_dir}.sra"
    queries = []
    try:
        for sra in gsm.relations['SRA']:
            query = sra.split("=")[-1]
            assert 'SRX' in query, "Sample looks like it is not SRA: %s" % query
            print("Query: %s" % query)
            queries.append(query)
    except KeyError:
        raise NoSRARelationException('No relation called SRA for %s' % gsm.get_accession())

    # check if the e-mail is more or less not a total crap
    Entrez.email = email
    if not (Entrez.email is not None and '@' in email and email != '' and '.' in email):
        raise Exception('You have to provide valid e-mail')

    for query in queries:
        # retrieve IDs for given SRX
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
                if "502" in str(httperr):
                    sys.stderr.write("Error: %s, trial %i out of %i, waiting for %i seconds." % (str(httperr),
                                                                                                 trial,
                                                                                                 number_of_trials,
                                                                                                 wait_time))
                    time.sleep(wait_time)
                else:
                    raise httperr
        df = pd.DataFrame([i.split(',') for i in results.split('\n') if i != ''][1:], columns = [i.split(',') for i in results.split('\n') if i != ''][0])

        # check it first
        try:
            df['download_path']
        except KeyError as e:
            sys.stderr.write('KeyError: ' + str(e) + '\n')
            sys.stderr.write(str(results) + '\n')

        # make the directory
        if platform.system() == "Windows":
            name_regex = r'[\s\*\?\(\),\.\:\%\|\"\<\>]'
        else:
            name_regex = r'[\s\*\?\(\),\.;]'
        directory_path = os.path.abspath(os.path.join(directory, "%s_%s_%s" % ('Supp',gsm.get_accession(),re.sub(name_regex, '_', gsm.metadata['title'][0]) # the directory name cannot contain many of the signs
        )))
        utils.mkdir_p(os.path.abspath(directory_path))

        for path in df['download_path']:
            sra_run = path.split("/")[-1]
            print("Analysing %s" % sra_run)
            url = ftpaddres.format(range_subdir=query[:6],
                                       record_dir=query,
                                       file_dir=sra_run)
            filepath = os.path.abspath(os.path.join(directory_path, "%s.sra" % sra_run))
            df.to_csv(os.path.abspath(os.path.join(directory_path, "%s_metadata.txt" % sra_run)), sep='\t')
            utils.download_from_url(url, filepath, aspera=aspera)

            if filetype in ["fasta", "fastq","sra"]:
                ftype = ""
                if filetype == "fasta":
                    ftype = " --fasta "
                cmd = "fastq-dump --split-files --gzip %s --outdir %s %s"
                cmd = cmd % (ftype, directory_path, filepath)

                process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                sys.stderr.write("Converting to %s/%s_*.%s.gz\n" % (
                    directory_path, sra_run, filetype))
                pout, perr = process.communicate()
                if "command not found" in perr:
                    raise NoSRAToolkitException("fastq-dump command not found")
                else:
                    print(pout)
                    if not keep_sra:
                        # Delete sra file
                        os.unlink(filepath)

def sra_series_fetch(gs_text, series, s3_text, output_path, email):
    gs_object = GEOparse.get_GEO(geo=gs_text)
    gsms_to_use = gs_object.gsms.values()
    filetype='sra'
    for gsm in gsms_to_use:
        sys.stderr.write("Downloading %s files for %s series\n" % (filetype, gsm.name))
        download_SRA(gsm, email=email, filetype=filetype, directory=output_path)
