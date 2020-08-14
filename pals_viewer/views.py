import json
from django.conf import settings
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from pals.pimp_tools import get_pimp_API_token_from_env, PIMP_HOST, download_from_pimp
from pals.feature_extraction import DataSource
from pals.loader import GNPSLoader
from pals.PLAGE import PLAGE
from pals.ORA import ORA
from pals.GSEA import GSEA
from pals.common import *
from bioservices.kegg import KEGG

import requests

import pandas as pd
import zipfile
import seaborn as sns
import datetime

# Create your views here.
def index(request):
    databases = [
        DATABASE_PIMP_KEGG,
        DATABASE_REACTOME_KEGG,
        DATABASE_REACTOME_CHEBI,
        DATABASE_REACTOME_UNIPROT,
        DATABASE_REACTOME_ENSEMBL,
        DATABASE_GNPS_MOLECULAR_FAMILY,
        DATABASE_GNPS_MS2LDA,
    ]

    reactome_species = [
        REACTOME_SPECIES_ARABIDOPSIS_THALIANA,
        REACTOME_SPECIES_BOS_TAURUS,
        REACTOME_SPECIES_CAENORHABDITIS_ELEGANS,
        REACTOME_SPECIES_CANIS_LUPUS_FAMILIARIS,
        REACTOME_SPECIES_DANIO_RERIO,
        REACTOME_SPECIES_DICTYOSTELIUM_DISCOIDEUM,
        REACTOME_SPECIES_DROSOPHILA_MELANOGASTER,
        REACTOME_SPECIES_GALLUS_GALLUS,
        REACTOME_SPECIES_HOMO_SAPIENS,
        REACTOME_SPECIES_MUS_MUSCULUS,
        REACTOME_SPECIES_ORYZA_SATIVA,
        REACTOME_SPECIES_RATTUS_NORVEGICUS,
        REACTOME_SPECIES_SACCHAROMYCES_CEREVISIAE,
        REACTOME_SPECIES_SUS_SCROFA,
    ]

    content_dict = {}
    content_dict['databases'] = databases
    content_dict['KEGG'] = DATABASE_PIMP_KEGG
    content_dict['reactome_species'] = reactome_species

    try:
        content_dict['local_token'] = get_pimp_API_token_from_env().strip()
    except:
        print("not found pimp token from environment variable")
        content_dict['local_token'] = ''

    return render(request,'pals_viewer/pathway_index.html', content_dict)

def gnps_index(request):
    databases = [
        DATABASE_PIMP_KEGG,
        DATABASE_REACTOME_KEGG,
        DATABASE_REACTOME_CHEBI,
        DATABASE_REACTOME_UNIPROT,
        DATABASE_REACTOME_ENSEMBL,
        DATABASE_GNPS_MOLECULAR_FAMILY,
        DATABASE_GNPS_MS2LDA,
    ]

    reactome_species = [
        REACTOME_SPECIES_ARABIDOPSIS_THALIANA,
        REACTOME_SPECIES_BOS_TAURUS,
        REACTOME_SPECIES_CAENORHABDITIS_ELEGANS,
        REACTOME_SPECIES_CANIS_LUPUS_FAMILIARIS,
        REACTOME_SPECIES_DANIO_RERIO,
        REACTOME_SPECIES_DICTYOSTELIUM_DISCOIDEUM,
        REACTOME_SPECIES_DROSOPHILA_MELANOGASTER,
        REACTOME_SPECIES_GALLUS_GALLUS,
        REACTOME_SPECIES_HOMO_SAPIENS,
        REACTOME_SPECIES_MUS_MUSCULUS,
        REACTOME_SPECIES_ORYZA_SATIVA,
        REACTOME_SPECIES_RATTUS_NORVEGICUS,
        REACTOME_SPECIES_SACCHAROMYCES_CEREVISIAE,
        REACTOME_SPECIES_SUS_SCROFA,
    ]

    content_dict = {}
    content_dict['databases'] = databases
    content_dict['KEGG'] = DATABASE_PIMP_KEGG
    content_dict['reactome_species'] = reactome_species

    return render(request,'pals_viewer/gnps_index.html', content_dict)

def ms2lda_index(request):
    databases = [
        DATABASE_PIMP_KEGG,
        DATABASE_REACTOME_KEGG,
        DATABASE_REACTOME_CHEBI,
        DATABASE_REACTOME_UNIPROT,
        DATABASE_REACTOME_ENSEMBL,
        DATABASE_GNPS_MOLECULAR_FAMILY,
        DATABASE_GNPS_MS2LDA,
    ]

    reactome_species = [
        REACTOME_SPECIES_ARABIDOPSIS_THALIANA,
        REACTOME_SPECIES_BOS_TAURUS,
        REACTOME_SPECIES_CAENORHABDITIS_ELEGANS,
        REACTOME_SPECIES_CANIS_LUPUS_FAMILIARIS,
        REACTOME_SPECIES_DANIO_RERIO,
        REACTOME_SPECIES_DICTYOSTELIUM_DISCOIDEUM,
        REACTOME_SPECIES_DROSOPHILA_MELANOGASTER,
        REACTOME_SPECIES_GALLUS_GALLUS,
        REACTOME_SPECIES_HOMO_SAPIENS,
        REACTOME_SPECIES_MUS_MUSCULUS,
        REACTOME_SPECIES_ORYZA_SATIVA,
        REACTOME_SPECIES_RATTUS_NORVEGICUS,
        REACTOME_SPECIES_SACCHAROMYCES_CEREVISIAE,
        REACTOME_SPECIES_SUS_SCROFA,
    ]

    content_dict = {}
    content_dict['databases'] = databases
    content_dict['KEGG'] = DATABASE_PIMP_KEGG
    content_dict['reactome_species'] = reactome_species

    return render(request,'pals_viewer/ms2lda_index.html', content_dict)

def keypath_get_data(request):
    """
    If data type is token or account, connect to remote host and download data to local.
    If data type is file, just save them to local.

    Return the local file name, column names.
    """
    if request.method != 'POST':
        return None
    
    # get data from POST
    token=request.POST.get('token')
    analysis_id=int(request.POST.get('analysis_id'))
    data_type=request.POST.get('data_type') # data type can be token, account and file

    print('token=', token)
    print('analysis_id=', analysis_id)
    print('data_type=', data_type)

    if data_type == 'token':
        # get data by token
        # download data from PiMP
        int_df, annotation_df, experimental_design = download_from_pimp(token, PIMP_HOST, analysis_id, 'kegg')

        # save data to local
        int_df_filename = save_dataframe_to_csv(int_df, 'int_df.csv')
        annotation_df_filename = save_dataframe_to_csv(annotation_df, 'annotation_df.csv')

        # get column names
        int_df_columns = int_df.columns.to_list()
        annotation_df_columns = annotation_df.columns.to_list()

        # TODO:save token to local

        result = {'status':'success', 'message':'Load data done!', 'data':{
            'int_df': {'filename': int_df_filename, 'columns': int_df_columns},
            'annotation_df': {'filename': annotation_df_filename, 'columns': annotation_df_columns},
            'experimental_design': experimental_design
        }}
        return HttpResponse(json.dumps(result))
    elif data_type == 'account':
        pass
    elif data_type == 'upload_file':
        int_df=request.FILES.get('int_csv')
        annotation_df=request.FILES.get('annotation_csv')
        experimental_design=""

        # save data to local
        int_df_filename = save_upload_file_to_csv(int_df, 'int_df.csv')
        annotation_df_filename = save_upload_file_to_csv(annotation_df, 'annotation_df.csv')

        # read data and transfer them to dataframe
        int_df = pd.read_csv(settings.MEDIA_ROOT+'/'+int_df_filename)
        int_df.set_index('row_id', inplace=True)
        annotation_df = pd.read_csv(settings.MEDIA_ROOT+'/'+annotation_df_filename)
        annotation_df.set_index('row_id', inplace=True)

        # get column names
        int_df_columns = int_df.columns.to_list()
        annotation_df_columns = annotation_df.columns.to_list()

        result = {'status':'success', 'message':'Load data done!', 'data':{
            'int_df': {'filename': int_df_filename, 'columns': int_df_columns},
            'annotation_df': {'filename': annotation_df_filename, 'columns': annotation_df_columns},
            'experimental_design': experimental_design
        }}
        return HttpResponse(json.dumps(result))
    else:
        result = {'status':'error', 'message':'Illegal data type!', 'data':{}}
        return HttpResponse(json.dumps(result))

def save_dataframe_to_csv(df, oldfilename):
    # save data to local
    now_time = str(datetime.datetime.now().strftime('%Y%m%d%H%M%S'))
    newfilename = oldfilename.replace(".", "_{}.".format(now_time))
    newfilepath = settings.MEDIA_ROOT+'/'+newfilename
    df.to_csv(newfilepath)
    print('saved {} to local:{}'.format(oldfilename, newfilepath))
    return newfilename

def save_upload_file_to_csv(upload_file, oldfilename):
    # save data to local
    now_time = str(datetime.datetime.now().strftime('%Y%m%d%H%M%S'))
    newfilename = oldfilename.replace(".", "_{}.".format(now_time))
    newfilepath = settings.MEDIA_ROOT+'/'+newfilename
    with open(newfilepath, 'wb+') as destination:
        for chunk in upload_file.chunks():
            destination.write(chunk)
    print('saved {} to local:{}'.format(oldfilename, newfilepath))

    return newfilename

def analysis(request):

    if request.method != 'POST':
        return None

    # get data from POST
    int_df_filename=request.POST.get('keypath_int_df_filename')
    annotation_df_filename=request.POST.get('keypath_annotation_df_filename')
    experimental_design=json.loads(request.POST.get('experimental_design'))
    pathway_analysis_method=request.POST.get('pathway_analysis_method')
    database=request.POST.get('database')
    reactome_species=request.POST.get('reactome_species')
    reactome_metabolic_pathway_only=request.POST.get('reactome_metabolic_pathway_only')
    reactome_query=request.POST.get('reactome_query')

    print('int_df_filename=', int_df_filename)
    print('annotation_df_filename=', annotation_df_filename)
    print('experimental_design=', experimental_design)
    print('pathway_analysis_method=', pathway_analysis_method)
    print('database=', database)
    print('reactome_species=', reactome_species)
    print('reactome_metabolic_pathway_only=', reactome_metabolic_pathway_only)
    print('reactome_query=', reactome_query)

    # load data from csv and reset first row as index
    # we must set first row as index, because the original table is it
    int_df = pd.read_csv(settings.MEDIA_ROOT+'/'+int_df_filename)
    int_df.set_index('row_id', inplace=True)
    annotation_df = pd.read_csv(settings.MEDIA_ROOT+'/'+annotation_df_filename)
    annotation_df.set_index('row_id', inplace=True)

    # set data source
    ds = DataSource(int_df, annotation_df, experimental_design, database, reactome_species=reactome_species, reactome_metabolic_pathway_only=reactome_metabolic_pathway_only, reactome_query=reactome_query)

    # analysis by different methods
    sorted_by = database+' beer1/beer2 comb_p'
    if pathway_analysis_method == 'PLAGE':
        plage = PLAGE(ds)
        pathway_df = plage.get_pathway_df()
        pathway_df.sort_values(sorted_by, ascending=True, inplace=True)
    elif pathway_analysis_method == 'ORA':
        ora = ORA(ds)
        pathway_df = ora.get_pathway_df()
        pathway_df.sort_values(sorted_by, ascending=True, inplace=True)
    elif pathway_analysis_method == 'GSEA':
        gsea = GSEA(ds)
        pathway_df = gsea.get_pathway_df()
        pathway_df.sort_values(sorted_by, ascending=True, inplace=True)

    pathway_json = json.loads(pathway_df.to_json(orient='split'))
    headers = pathway_json['columns']
    headers.insert(0, '')
    
    rows = []
    for index, row in zip(pathway_json['index'], pathway_json['data']):
        row.insert(0, index)
        rows.append(row)

    table = {
        'headers': headers,
        'rows': rows
    }

    result = {'message':'Analysis done!', 'data':{'table':table}}
    
    return HttpResponse(json.dumps(result))

def show_kegg_diagram(request):
    if request.method != 'POST':
        return None
    
    content_dict = {}
    details = []
    
    # get data from POST
    id=request.POST.get('id') # pathway index of row of table
    pathway_name=request.POST.get('pathway_name') # pathway name of selected row of table
    row=json.loads(request.POST.get('row')) # the select row of table
    stId=id
    print('row=',row)

    # calculate information
    label = '%s: %s' % (stId, pathway_name)
    info_url = 'https://www.genome.jp/dbget-bin/www_bget?%s' % stId
    header = '<h3>{} [<a href="{}" target="_blank" rel="noopener noreferrer">Info</a>]</h3>'.format(label, info_url)
    details.append(header)

    p_value = row['p-value']
    num_hits = row['Formula Hits']
    details.append('<h4>p-value: %.6f</h4>' % p_value)
    details.append('<h4>Formula Hits: %d</h4>' % (num_hits))
    details.append('<h4>Summary:</h4>')

    dict_data = get_kegg_info(stId)
    for k, v in dict_data.items():
        # if k in ['CLASS', 'MODULE', 'DISEASE', 'REL_PATHWAY']:
        if k in ['CLASS']:
            details.append('<p>{}: {}</p>'.format(k, v))

    image_url = 'https://www.genome.jp/kegg/pathway/map/%s.png' % stId
    logger.debug('image_url = %s' % image_url)
    details.append('<img src="{}" class="img-fluid" alt="" srcset="">'.format(image_url))

    result = {'status':'success', 'message':'Analysis done!', 'data':{'details':details}}

    return HttpResponse(json.dumps(result))

def show_reactome_diagram(request):
    if request.method != 'POST':
        return None
    
    content_dict = {}
    details = []
    
    # get data from POST
    token=request.POST.get('token', None)
    id=request.POST.get('id') # pathway index of row of table
    pathway_name=request.POST.get('pathway_name') # pathway name of selected row of table
    row=json.loads(request.POST.get('row')) # the select row of table
    stId=id
    print('row=',row)

    # get reactome info
    status_code, json_response = get_reactome_info(stId)
    if status_code != 200:
        result = {'status':'error', 'message':'Get reactome info fail!', 'data':{'details':[]}}
        return HttpResponse(json.dumps(result))

    # calculate information
    label = '%s: %s' % (stId, pathway_name)
    info_url = 'https://reactome.org/content/detail/%s' % stId

    if token is None:
        header = '<h3>{} [<a href="{}" target="_blank" rel="noopener noreferrer">Info</a>]</h3>'.format(label, info_url)
    else:
        # TODO: the token always show invalid, so we cancel it.
        header = '<h3>{} [<a href="{}" target="_blank" rel="noopener noreferrer">Info</a>]</h3>'.format(label, info_url)
        # viewer_url = 'https://reactome.org/PathwayBrowser/#/%s&DTAB=AN&ANALYSIS=%s' % (stId, token)
        # header = '<h3>{} [<a href="{}" target="_blank" rel="noopener noreferrer">Info</a>] [<a href="{}" target="_blank" rel="noopener noreferrer">Viewer</a>]</h3>'.format(label, info_url, viewer_url)
    details.append(header)

    for summation in json_response['summation']:
        details.append('<p>{}</p>'.format(summation['text']))

    image_url = 'https://reactome.org/ContentService/exporter/diagram/%s.png?quality=8&diagramProfile=standard&analysisProfile=strosobar' % stId
    if token is not None:
        # TODO: the token always show invalid, so we cancel it.
        # image_url += '&token=%s&resource=TOTAL&expColumn=0' % token
        pass

    logger.debug('image_url = %s' % image_url)
    details.append('<img src="{}" class="img-fluid" alt="" srcset="">'.format(image_url))

    result = {'status':'success', 'message':'Analysis done!', 'data':{'details':details}}

    return HttpResponse(json.dumps(result))

def get_kegg_info(stId):
    k = KEGG()
    data = k.get(stId)
    dict_data = k.parse(data)
    return dict_data

def get_reactome_info(stId):
    # refer to https://reactome.org/dev/content-service
    url = 'https://reactome.org/ContentService/data/query/%s' % stId
    logger.debug('Reactome URL: ' + url)

    # make a GET request to Reactome Content service
    response = requests.get(url)
    logger.debug('Response status code = %d' % response.status_code)

    status_code = response.status_code
    if status_code == 200:
        json_response = json.loads(response.text)
    else:
        json_response = None
    return status_code, json_response


def gnps_get_data(request):
    """
    If data type is token or account, connect to remote host and download data to local.
    If data type is file, just save them to local.

    Return the local file name, column names.
    """
    if request.method != 'POST':
        return None
    
    # get data from POST
    metadata_df=request.FILES.get('metadata_csv')

    # save data to local
    metadata_df_filename = save_upload_file_to_csv(metadata_df, 'metadata_df.csv')

    # read data and transfer them to dataframe
    metadata_df = pd.read_csv(settings.MEDIA_ROOT+'/'+metadata_df_filename)
    print('metadata columns:', metadata_df.columns.to_list())

    # get column names
    metadata_df_columns = metadata_df.columns.to_list()

    # get group values
    group_values = metadata_df['group'].drop_duplicates().values.tolist()


    result = {'status':'success', 'message':'Load data done!', 'data':{
        'metadata_df': {'filename': metadata_df_filename, 'columns': metadata_df_columns, 'groups': group_values}
    }}
    return HttpResponse(json.dumps(result))


def gnps_analysis(request):
    if request.method != 'POST':
        return None

    # get data from POST
    metadata_df_filename=request.POST.get('gnps_metadata_df_filename')
    gnps_url=request.POST.get('gnps_url')
    experimental_design=json.loads(request.POST.get('experimental_design'))
    comparisons = experimental_design['comparisons']
    database = DATABASE_GNPS_MOLECULAR_FAMILY

    print('metadata_df_filename=', metadata_df_filename)
    print('gnps_url=', gnps_url)
    print('experimental_design=', experimental_design)

    # load data from csv
    metadata_df = pd.read_csv(settings.MEDIA_ROOT+'/'+metadata_df_filename)

    # load data from remote
    loader = GNPSLoader(database, gnps_url, metadata_df, comparisons)
    database = loader.load_data()

    measurement_df = database.extra_data['measurement_df']
    annotation_df = database.extra_data['annotation_df']
    experimental_design = database.extra_data['experimental_design']

    # set data source
    gnps_ds = DataSource(measurement_df, annotation_df, experimental_design, None, database=database, min_replace=SMALL)

    plage = PLAGE(gnps_ds, num_resamples=1000)
    pathway_df = plage.get_pathway_df(standardize=True)

    pathway_json = json.loads(pathway_df.to_json(orient='split'))
    headers = pathway_json['columns']
    headers.insert(0, '')
    
    rows = []
    for index, row in zip(pathway_json['index'], pathway_json['data']):
        row.insert(0, index)
        rows.append(row)

    table = {
        'headers': headers,
        'rows': rows
    }

    result = {'message':'Analysis done!', 'data':{'table':table}}
    
    return HttpResponse(json.dumps(result))

def ms2lda_analysis(request):
    pass