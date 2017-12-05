def pathway_gene_heatmap(df, gene_list):
    '''This function maps the log2fc for each gene within a
    given pathway and displays them as a heatmap for each
    strain within the data frame

    df = analyzeRepeats data frame
    gene_list = python list of gene names'''

    import pandas as pd
    import seaborn as sns
    import numpy as np
    import re

    df = df.copy(deep=True)
    df = df.set_index(df['Annotation/Divergence'].str.split('|').str[0])

    # process analyzeRepeats final
    anno_col = 8 # currently 8 annotation columns in output
    n_samples = df.shape[1] - anno_col # number of samples
    samples = pd.Series(df.columns[anno_col:]).str.split(' ').str[0].str.split('/').str[-1]

    pattern = r'(^[A-Z]{2,6})_([A-Z]+).*_([A-Z]+)Diet.*_([A-Z0-9]+)Week'
    discrim = pd.DataFrame(samples.str.findall(pattern, flags=re.IGNORECASE).str[0].tolist())

    # drop uninformative columns and label the
    # resulting columns as 
    nunique = discrim.apply(pd.Series.nunique)
    cols_to_drop = nunique[nunique == 1].index
    discrim = discrim.drop(cols_to_drop, axis=1)
    discrim.columns = ['strain', 'diet']

    # convert to labels
    subj_groups = discrim.apply('_'.join, axis=1)

    # subset only genes within the pathway, then remove the annotation
    row_index = df.loc[:, 'Annotation/Divergence'].str.match('\||'.join(gene_list.iloc[:,0]), flags=re.IGNORECASE)
    df2 = df.loc[row_index, :]
    n_results = df2.shape[0]
    df2 = df2.iloc[:, anno_col:]

    # calculate log fold changes
    fc_mat = {}
    for strain in discrim.strain.drop_duplicates():
        s_amln = df2.columns.str.contains(strain + '.*AMLN', flags=re.IGNORECASE)
        s_cntrl = df2.columns.str.contains(strain + '.*Control', flags=re.IGNORECASE)
        log2fc = np.log2((df2.loc[:, s_amln].mean(axis=1)+1)/(df2.loc[:, s_cntrl].mean(axis=1)+1))
        fc_mat[strain + '_log2fc'] = log2fc
    cg = sns.clustermap(pd.DataFrame(fc_mat).T, figsize=(15,3), cmap='BrBG', vmin=-3, vmax=3)
    cg.ax_row_dendrogram.set_visible(False)
    cg.ax_col_dendrogram.set_visible(False)
    cg.cax.set_visible(False)

    return cg

def plot_gene_exp(df, gene, plot_style='bar', groups = None, by_refseq = False, by_diet = True):
    '''Function to plot the expression values of a particular gene using the output
    from analyze_repeats.pl in the HOMER package. df is the output of analyze_repeats. 
    Currently function will try and guess the best groups

    In the future you can set groups yourself. groups should be a list with labels that group the columns into subgroups for
    plotting.
    '''
    import pandas as pd
    import re
    import seaborn as sns

    df = df.copy(deep=True)
    df.columns = pd.Series(df.columns)\
        .str.split(r" \(").str[0]\
        .str.split('/.*/').str[-1]
    anno_col = 8 # currently 8 annotation columns in output
    n_samples = df.shape[1] - anno_col # number of samples
    samples = pd.Series(df.columns[anno_col:])\
            .str.split(' ').str[0]\
            .str.split('/').str[-1]
    pattern = r'(^[A-Z]{2,6})[_-]([A-Z]+).*[_-]([A-Z]+)Diet.*[_-]([A-Z0-9]+)Week'
    discrim = pd.DataFrame(samples.str.findall(pattern, flags=re.IGNORECASE).str[0].tolist())

    # drop uninformative columns
    nunique = discrim.apply(pd.Series.nunique)
    cols_to_drop = nunique[nunique == 1].index
    discrim = discrim.drop(cols_to_drop, axis=1)

    # convert to labels
    subj_groups = discrim.apply('_'.join, axis=1)

    # create df for plotting
    df.loc[:, 'Gene'] = df.loc[:, 'Annotation/Divergence'].str.split('|').str[0]
    if by_refseq:
        row_index = df.loc[:, 'Transcript/RepeatID'].str.contains(gene, flags=re.IGNORECASE)
        plot_col = 'Transcript/RepeatID'
    else:
        row_index = df.loc[:, 'Annotation/Divergence'].str.contains(gene+'\|', flags=re.IGNORECASE)
        plot_col = 'Gene'
    df2 = df.loc[row_index, :]
    n_results = df2.shape[0]
    df3 = df2.iloc[:, anno_col:(anno_col + n_samples)]
    # df3.columns = subj_groups.append(pd.Series(plot_col))
    df3.loc[:, plot_col] = pd.Series(df2.loc[:, plot_col], index = df2.index)
    df3 = df3.melt(id_vars = plot_col, var_name='Group', value_name='Expression')


    if by_diet:
        df3['diet'] = df3['Group'].str.extract('(AMLN|Control)', flags=re.IGNORECASE, expand=True)
        df3['strain'] = df3['Group'].str.split('_').str[0]
        ax = sns.factorplot(x="strain", y="Expression", hue="diet",
                   data=df3, kind=plot_style, col=plot_col, # capsize=0.1,
                   sharey=False, col_wrap=3)
    else:
        # plot the data
        ax = sns.factorplot(x="Group", y="Expression",
                   data=df3, kind=plot_style, col=plot_col, # capsize=0.1,
                   sharey=False)
        ax.set_xticklabels(rotation=30)

    return ax

def run_annotateRepeats(dirpaths, genome = 'mm10', countType = 'genes',
                        raw = False, print_out = False, norm = None):
    '''Python wrapper for analyzeRepeates.pl
    Right now a basic implementation, evnetually add more
    functionality for variables to be passed'''
    import subprocess
    aRCall = ['analyzeRepeats.pl', 'rna', genome, '-count', countType]
    if raw:
        aRCall.extend(['-raw'])
    elif norm !=None:
        aRCall.extend(['-' + norm])
    else:
        aRCall.extend(['-norm', '1e7'])
    aRCall.extend(['-d'])
    aRCall.extend(dirpaths)
    if print_out: print(aRCall)
    p = subprocess.Popen(aRCall, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out, err = p.communicate()
    
    return out, err

def run_getDiffExpression(dirpaths, groups, batch = None, genome = 'mm10', countType = 'genes', print_out = False):
    '''Python wrapper for getDiffExpression.pl
    groups should be the experimental groups for analysis
    batch can be used to control for experimental date
    see documentation for getDiffExprssion.pl for more info
    
    eventually implement option to control the algorithm used'''
    import subprocess
    
    # Run annotate repeats with no adjustments
    out, err = run_annotateRepeats(dirpaths, countType = countType, raw = True)
    f = open(outdir + '/tmpARForgDE.txt', 'w')
    f.write(out.decode('utf-8'))
    f.close()
    
    # Run diff exp
    gDECall = ['getDiffExpression.pl', outdir + '/tmpARForgDE.txt']
    gDECall.extend(groups)
    # gDECall.extend(['-edgeR'])
    if batch is not None:
        gDECall.extend(['-batch'])
        gDECall.extend(batch)
    if print_out: print(gDECall)
    p = subprocess.Popen(gDECall, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    gDEOut, gDEErr = p.communicate()
    
    subprocess.call(['rm', outdir + '/tmpARForgDE.txt'])
    return gDEOut, gDEErr

def df_col_subset(x, id_str):
    '''Returns pandas dataframe with only columns containing id_str'''
    import pandas as pd
    import re
    return x[x.columns[pd.Series(x.columns).str.contains(id_str, flags=re.IGNORECASE)]]

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
    if taxon_id != None: requestURL += "&taxonId=" + str(taxon)
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
