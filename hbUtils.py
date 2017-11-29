def plot_gene_exp(df, gene, groups = None, by_diet = True):
    '''Function to plot the expression values of a particular gene using the output
    from analyze_repeats.pl in the HOMER package. df is the output of analyze_repeats. 
    Currently function will try and guess the best groups

    In the future you can set groups yourself. groups should be a list with labels that group the columns into subgroups for
    plotting.
    '''
    import pandas as pd
    import re
    import seaborn as sns

    anno_col = 8 # currently 8 annotation columns in output
    n_samples = df.shape[1] - anno_col # number of samples
    samples = pd.Series(df.columns[anno_col:])\
            .str.split(' ').str[0]\
            .str.split('/').str[-1]

    pattern = r'(^[A-Z]{2,6})_([A-Z]+).*_([A-Z]+)Diet.*_([A-Z0-9]+)Week'
    discrim = pd.DataFrame(samples.str.findall(pattern, flags=re.IGNORECASE).str[0].tolist())

    # drop uninformative columns
    nunique = discrim.apply(pd.Series.nunique)
    cols_to_drop = nunique[nunique == 1].index
    discrim = discrim.drop(cols_to_drop, axis=1)

    # convert to labels
    subj_groups = discrim.apply('_'.join, axis=1)

    # create df for plotting
    df.loc[:, 'gene'] = df.loc[:, 'Annotation/Divergence'].str.split('|').str[0]
    row_index = df.loc[:, 'Annotation/Divergence'].str.contains(gene+'\|', flags=re.IGNORECASE)
    df2 = df.loc[row_index, :]
    n_results = df2.shape[0]
    df2 = df2.iloc[:, anno_col:]
    df2.columns = subj_groups.append(pd.Series('Gene'))
    df2 = df2.melt(id_vars = 'Gene', var_name='Group', value_name='Expression')


    if by_diet:
        df2['diet'] = df2['Group'].str.extract('(AMLN|Control)', flags=re.IGNORECASE, expand=True)
        df2['strain'] = df2['Group'].str.split('_').str[0]
        ax = sns.factorplot(x="strain", y="Expression", hue="diet",
                   data=df2, kind="bar", col='Gene', capsize=0.1)
    else:
        # plot the data
        ax = sns.factorplot(x="Group", y="Expression",
                   data=df2, kind="bar", col='Gene', capsize=0.1)
        ax.set_xticklabels(rotation=30)


    return ax

def pull_json_db(url):
    ''' pulls json db and returns it as a json object'''
    import requests, sys, json
    
    r = requests.get(url, headers={ "Accept" : "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    return json.loads(r.text)

def pull_go_terms(term, taxon_id = None, product_type = None):
    ''' this function pulls go gene IDs from the QuickGo API
    interface. Currently simply pulls all genes directly associatd
    with the term and returns a list of the gene id, gene name, 
    associated organism, and gene product type

    plan to add filtering options to filter by animal (simple using the
    request url) AND add pulling from descendents (complicated)

    term is a string including the "GO:" OR a string without the GO: header OR an int
    without the GO header.
    '''

    # import req packages
    import requests, sys, json

    # clean input
    if type(term) == str:
        term = term.replace('GO:', '')

    if type(term) == int:
        term = str(term)

    requestURL = "https://www.ebi.ac.uk/QuickGO/services/annotation/search?goId=GO%3A" + term 
    if taxon_id != None: requestURL += "&taxonId=" + str(taxon_id)
    if product_type != None: requestURl += "&geneProductType=" + str(product_type)
    requestURL += "&limit=100"

    db = pull_json_db(requestURL)

    pagenum = db['pageInfo']['total']
    gene_list = [[x['geneProductId'], x['symbol'], x['taxonId']] for x in db['results']]
    
    if pagenum == 1:
        return gene_list
    
    else:
        for page in range(2, pagenum + 1): # want to start at the second page since we initialized
            db = pull_json_db(requestURL + '&page=' + str(pagenum))
            gene_list.extend([[x['geneProductId'], x['symbol'], x['taxonId']] for x in db['results']])
    
    return gene_list

def find_file_path(fileDir, fileID):
    '''Finds file on server given path to search and file ID
    
    I originaly made this to search the epiglass archives for
    files based on more general IDs'''
    
    import subprocess
    
    # find the file paths in given fileDir
    p = subprocess.Popen(["find", fileDir, "-iname", fileID + '*'],
    stdout=subprocess.PIPE, stderr = subprocess.PIPE)
    out, err = p.communicate()
    
    # alert if no files found
    if out == '':
        print('No files found for ID:' + fileID)
    
    # format bash output into list
    out = out.split('\n')
    
    # remove blanks that come in from the stdout
    try:
        out.remove('')
    except ValueError:
        pass # do nothing
    
    return out

def summarize_fastqc(fastqcDir, fileName):
    '''summarize results inline for a fastqc report given the directory'''
    f = open(fastqcDir + '/summary.txt')
    fastqcSummary = f.read()
    if fastqcSummary.count('FAIL') > 0:
        print(fileName.split('/')[-1] + " failed " + str(fastqcSummary.count('FAIL')) + " fastqc test(s)")
    if fastqcSummary.count('WARN') > 0:
        print(fileName.split('/')[-1] + " warned " + str(fastqcSummary.count('WARN')) + " fastqc test(s)")
    if fastqcSummary.count('WARN') + fastqcSummary.count('FAIL') == 0:
        print(fileName.split('/')[-1] + " passed all fastqc tests, cheers!")

def run_fastqc(archiveDir, outdir, *identifiers):
    '''Python wrapper for fastqc
    
    Takes a list of identifying tag directory names, finds them within
    the archive directory passed, and funs fastqc on them
    '''
    import subprocess  
    
    paths = []
    for i in identifiers:
        # Use find_file_path to get the file paths of the identifiers passed
        paths.extend(find_file_path(archiveDir, i))
    
    for i in paths:
        # run fastqc and unzip the results
        print('\n' + '-' * 40)
        print("trying fastqc on file:" + i)
        subprocess.call(["fastqc", i, "-o", outdir])
        zipPath = outdir + i.split('/')[-1].replace('.fastq.gz', '_fastqc.zip')
        subprocess.call(["unzip", "-o", zipPath, "-d", outdir])
        dirPath = zipPath.replace('.zip', '')
        #summarize fastqc
        summarize_fastqc(dirPath, i)

    # return paths of fastqc reports
    return [outdir + i.split('/')[-1].replace('.fastq.gz', '_fastqc') for i in paths]
    
    # Next the program needs to summarize the data
    # Maybe have it return the paths to the fastqc reports?
