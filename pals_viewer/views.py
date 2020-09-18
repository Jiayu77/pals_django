import json
from django.conf import settings
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect, JsonResponse
from pals.pimp_tools import get_pimp_API_token_from_env, PIMP_HOST, download_from_pimp
from pals.feature_extraction import DataSource
from pals.loader import GNPSLoader
from pals.PLAGE import PLAGE
from pals.ORA import ORA
from pals.GSEA import GSEA
from pals.common import *
from bioservices.kegg import KEGG

import requests
import pickle
import pandas as pd
import zipfile
import seaborn as sns
import datetime
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

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
    
    # add active nav id
    content_dict['active_nav_id'] = '#nav-pathway'

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
    # add active nav id
    content_dict['active_nav_id'] = '#nav-gnps'

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
    # add active nav id
    content_dict['active_nav_id'] = '#nav-ms2lda'

    return render(request,'pals_viewer/ms_two_lda_index.html', content_dict)

def keypath_get_data(request):
    """
    If data type is token or account, connect to remote host and download data to local.
    If data type is file, just save them to local.

    Return the local file name, column names.
    """
    if request.method != 'POST':
        return None
    
    # get data from POST
    data_type=request.POST.get('data_type') # data type can be token, account and file

    print('data_type=', data_type)

    if data_type == 'token':
        token=request.POST.get('token')
        analysis_id=int(request.POST.get('analysis_id'))
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
        return JsonResponse(result)
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
        return JsonResponse(result)
    else:
        result = {'status':'error', 'message':'Illegal data type!', 'data':{}}
        return JsonResponse(result)

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

def save_object_to_pkl(object_to_save, oldfilename):
    # save object to local
    now_time = str(datetime.datetime.now().strftime('%Y%m%d%H%M%S'))
    newfilename = oldfilename.replace(".", "_{}.".format(now_time))
    newfilepath = settings.MEDIA_ROOT+'/'+newfilename
    pickle_file = open(newfilepath,'wb')
    pickle.dump(object_to_save, pickle_file)
    pickle_file.close()
    print('saved object {} to local:{}'.format(oldfilename, newfilepath))
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
    int_df = pd.read_csv(settings.MEDIA_ROOT+'/'+int_df_filename, skiprows=[1], index_col=0)

    annotation_df = pd.read_csv(settings.MEDIA_ROOT+'/'+annotation_df_filename)
    annotation_df.set_index('row_id', inplace=True)

    # set data source
    ds = DataSource(int_df, annotation_df, experimental_design, database, reactome_species=reactome_species, reactome_metabolic_pathway_only=reactome_metabolic_pathway_only, reactome_query=reactome_query)

    # analysis by different methods
    if pathway_analysis_method == 'PLAGE':
        plage = PLAGE(ds)
        pathway_df = plage.get_pathway_df()
    elif pathway_analysis_method == 'ORA':
        ora = ORA(ds)
        pathway_df = ora.get_pathway_df()
    elif pathway_analysis_method == 'GSEA':
        gsea = GSEA(ds)
        pathway_df = gsea.get_pathway_df()

    pathway_json = json.loads(pathway_df.to_json(orient='split'))
    headers = pathway_json['columns']
    headers.insert(0, '')
    
    rows = []
    for index, row in zip(pathway_json['index'], pathway_json['data']):
        row.insert(0, index)
        rows.append(row)

    # default display pathway name, p-values and F values
    display_columns = ["pw_name", "unq_pw_F", "tot_ds_F", "F_coverage"]
    for comparison in experimental_design['comparisons']:
        comparison_name = comparison['name']
        display_columns.append("{} p-value".format(comparison_name))

    table = {
        'headers': headers,
        'rows': rows,
        'display_columns': display_columns
    }

    result = {'message':'Analysis done!', 'data':{'table':table}}
    
    return JsonResponse(result)

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

    return JsonResponse(result)

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
        return JsonResponse(result)

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

    result = {'status':'success', 'message':'Analysis done!', 'data':{'details':details}}

    return JsonResponse(result)

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
    return JsonResponse(result)


def gnps_analysis(request):
    if request.method != 'POST':
        return None

    # get data from POST
    metadata_df_filename=request.POST.get('gnps_metadata_df_filename')
    gnps_url=request.POST.get('gnps_url')
    comparisons=json.loads(request.POST.get('comparisons'))
    database = DATABASE_GNPS_MOLECULAR_FAMILY

    print('metadata_df_filename=', metadata_df_filename)
    print('gnps_url=', gnps_url)
    print('comparisons=', comparisons)

    # load data from csv
    metadata_df = pd.read_csv(settings.MEDIA_ROOT+'/'+metadata_df_filename)
    # metadata_df.reset_index()

    loader = GNPSLoader(database, gnps_url, metadata_df, comparisons)
    database = loader.load_data()

    # cache datasource
    gnps_load_data_filename = save_object_to_pkl(database, 'gnps_load_data.pkl')

    measurement_df = database.extra_data['measurement_df']
    annotation_df = database.extra_data['annotation_df']
    experimental_design = database.extra_data['experimental_design']

    # set data source
    ds = DataSource(measurement_df, annotation_df, experimental_design, None, database=database, min_replace=SMALL)
    plage = PLAGE(ds)
    df = plage.get_pathway_df()

    df_json = json.loads(df.to_json(orient='split'))
    headers = df_json['columns']
    headers.insert(0, '')

    rows = []
    for index, row in zip(df_json['index'], df_json['data']):
        row.insert(0, index)
        rows.append(row)
    
    # default display pathway name, p-values and F values
    display_columns = ["pw_name", "unq_pw_F", "tot_ds_F", "F_coverage"]
    for comparison in experimental_design['comparisons']:
        comparison_name = comparison['name']
        display_columns.append("{} p-value".format(comparison_name))

    table = {
        'headers': headers,
        'rows': rows,
        'display_columns': display_columns,
    }

    result = {
        'message':'Analysis done!', 
        'data':{
            'table':table,
            'gnps_load_data_filename': gnps_load_data_filename,
        }
    }
    
    return JsonResponse(result)


def gnps_show_details(request):
    if request.method != 'POST':
        return None
    
    content_dict = {}
    details = []

    # get data from POST
    row=json.loads(request.POST.get('row')) # the select row of table
    gnps_load_data_filename=request.POST.get('gnps_load_data_filename')
    print('row=', row)

    # load database from local
    with open(settings.MEDIA_ROOT+'/'+gnps_load_data_filename, "rb") as gnps_load_data_file:
        database = pickle.load(gnps_load_data_file)

    measurement_df = database.extra_data['measurement_df']
    annotation_df = database.extra_data['annotation_df']
    experimental_design = database.extra_data['experimental_design']

    # convert to a DataSource object that can be used by PLAGE
    ds = DataSource(measurement_df, annotation_df, experimental_design, None, database=database, min_replace=SMALL)

    # run PLAGE decomposition on the ds
    df = PLAGE_decomposition(ds)

    # we assume comparisons has only one item
    comparisons = experimental_design['comparisons']
    case = comparisons[0]['case']
    control = comparisons[0]['control']

    significant_column = '%s p-value' % (comparisons[0]['name'])
    df_filtered = process_gnps_results(df, significant_column)

    all_groups, all_samples, entity_dict, intensities_df, dataset_pathways_to_row_ids = get_plot_data(ds, case, control)
    results = {
        'df': df_filtered,
        'all_groups': all_groups,
        'all_samples': all_samples,
        'entity_dict': entity_dict,
        'intensities_df': intensities_df,
        'dataset_pathways_to_row_ids': dataset_pathways_to_row_ids,
        'database_name': database.database_name,
    }

    # quick hack to pass the motif db urls out for plotting/displaying
    if 'motifdb_urls' in database.extra_data:
        results['motifdb_urls'] = database.extra_data['motifdb_urls']

    for key in results:
        print("{}: type is {}".format(key, type(results[key])))

    # if this is molecular family analysis, the selected name will be e.g. 'Molecular Family #123 (...)'
    # since we split by space, we take the token at position 2 and remove the first character to get the index
    # otherwise if this is e.g. MS2LDA analysis, the selected name will be e.g. 'motif_123 (...), so
    # we just take the first token as the index

    # update row name
    # 'pw_name': 'Components',
    # significant_column: 'p-value',
    # 'tot_ds_F': 'No. of members',
    row['Components'] = row['pw_name']
    row['p-value'] = row[significant_column]
    row['No. of members'] = row['tot_ds_F']

    pw_name = row['Components']
    p_value = row['p-value']
    no_members = row['No. of members']
    selected = '%s (p-value=%.6e, members=%d)' % (pw_name, p_value, no_members)

    tokens = selected.split(' ')
    database_name = results['database_name']
    idx = tokens[2][1:] if database_name == DATABASE_GNPS_MOLECULAR_FAMILY else tokens[0]
    # display the selected row
    # st.write(row)
    print('row=',row)

    # write the link to MotifDB if available
    motifdb_link = get_motifdb_link(results, row)
    if motifdb_link is not None:
        # display link
        print('Link to MotifDB: %s' % motifdb_link)

    members = dataset_pathways_to_row_ids[idx]
    member_df = get_member_df(entity_dict, members)
    plot_heatmap_image_filename = plot_heatmap(all_groups, all_samples, intensities_df, member_df, members, row)
    # detail_table = display_member_df(member_df)

    if 'link' in member_df.columns:
        member_df['link'] = member_df['link'].apply(make_clickable)

    print('plot_heatmap_image_filename=',plot_heatmap_image_filename)

    details = []

    # current row
    details.append('<h2>Component Browser</h2>')
    details.append('<h4>Components: {}</h4>'.format(no_members))
    details.append('<h4>p-value:  {}</h4>'.format(p_value))
    details.append('<h4>o. of members: {}</h4>'.format(pw_name))

    # heatmap
    image_url = settings.STATIC_URL+plot_heatmap_image_filename
    logger.debug('image_url = %s' % image_url)
    details.append('<img src="{}" class="img-fluid" alt="" srcset="">'.format(image_url))

    # members table
    details.append('<h2>Members</h2>')
    details.append(member_df.to_html(escape=False))

    result = {'status':'success', 'message':'Analysis done!', 'data':{'details':details}}
    return JsonResponse(result)

def ms2lda_get_data(request):
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
        'metadata_df': {
            'filename': metadata_df_filename,
            'columns': metadata_df_columns,
            'groups': group_values
        }
    }}
    return JsonResponse(result)

def ms2lda_analysis(request):
    # get data from POST
    metadata_df_filename=request.POST.get('ms2lda_metadata_df_filename')
    ms2lda_url=request.POST.get('motif_data_url')
    gnps_url=request.POST.get('peak_data_url')
    comparisons=json.loads(request.POST.get('comparisons'))
    data_type=request.POST.get('data_type')
    database_name = DATABASE_GNPS_MS2LDA

    print('metadata_df_filename=', metadata_df_filename)
    print('gnps_url=', gnps_url)
    print('comparisons=', comparisons)

    # load data from csv
    metadata_df = pd.read_csv(settings.MEDIA_ROOT+'/'+metadata_df_filename)

    if data_type == 'upload_file':
        peak_table_df=request.FILES.get('peak_data_csv')
        peak_table_df_filename = save_upload_file_to_csv(peak_table_df, 'peak_table_df.csv')
        peak_table_df = pd.read_csv(settings.MEDIA_ROOT+'/'+peak_table_df_filename)

        loader = GNPSLoader(database_name, gnps_url, metadata_df, 
            comparisons, gnps_ms2lda_url=ms2lda_url, peak_table_df=peak_table_df)
    else:
        loader = GNPSLoader(database_name, gnps_url, metadata_df, 
            comparisons, gnps_ms2lda_url=ms2lda_url)

    database = loader.load_data()
    # cache datasource
    ms2lda_load_data_filename = save_object_to_pkl(database, 'ms2lda_load_data.pkl')

    measurement_df = database.extra_data['measurement_df']
    annotation_df = database.extra_data['annotation_df']
    experimental_design = database.extra_data['experimental_design']

    # set data source
    gnps_ds = DataSource(measurement_df, annotation_df, experimental_design, None, database=database, min_replace=SMALL)

    plage = PLAGE(gnps_ds)
    df = plage.get_pathway_df()

    df_json = json.loads(df.to_json(orient='split'))
    headers = df_json['columns']
    headers.insert(0, '')
    
    rows = []
    for index, row in zip(df_json['index'], df_json['data']):
        row.insert(0, index)
        rows.append(row)

    # default display pathway name, p-values and F values
    display_columns = ["pw_name", "unq_pw_F", "tot_ds_F", "F_coverage"]
    for comparison in experimental_design['comparisons']:
        comparison_name = comparison['name']
        display_columns.append("{} p-value".format(comparison_name))

    table = {
        'headers': headers,
        'rows': rows,
        'display_columns': display_columns
    }

    result = {
        'message':'Analysis done!', 
        'data':{
            'table':table,
            'ms2lda_load_data_filename': ms2lda_load_data_filename
        }
    }
    
    return JsonResponse(result)

def ms2lda_show_details(request):
    if request.method != 'POST':
        return None
    
    content_dict = {}
    details = []

    # get data from POST
    ms2lda_load_data_filename=request.POST.get('ms2lda_load_data_filename')
    row=json.loads(request.POST.get('row')) # the select row of table
    print('row=', row)

    # load database from local
    with open(settings.MEDIA_ROOT+'/'+ms2lda_load_data_filename, "rb") as ms2lda_load_data_file:
        database = pickle.load(ms2lda_load_data_file)

    measurement_df = database.extra_data['measurement_df']
    annotation_df = database.extra_data['annotation_df']
    experimental_design = database.extra_data['experimental_design']

    # convert to a DataSource object that can be used by PLAGE
    ds = DataSource(measurement_df, annotation_df, experimental_design, None, database=database, min_replace=SMALL)

    # run PLAGE decomposition on the ds
    df = PLAGE_decomposition(ds)

    # we assume comparisons has only one item
    comparisons = experimental_design['comparisons']
    case = comparisons[0]['case']
    control = comparisons[0]['control']

    significant_column = '%s p-value' % (comparisons[0]['name'])
    df = process_gnps_results(df, significant_column)

    all_groups, all_samples, entity_dict, intensities_df, dataset_pathways_to_row_ids = get_plot_data(ds, case, control)
    results = {
        'df': df,
        'all_groups': all_groups,
        'all_samples': all_samples,
        'entity_dict': entity_dict,
        'intensities_df': intensities_df,
        'dataset_pathways_to_row_ids': dataset_pathways_to_row_ids,
        'database_name': database.database_name,
    }

    # quick hack to pass the motif db urls out for plotting/displaying
    if 'motifdb_urls' in database.extra_data:
        results['motifdb_urls'] = database.extra_data['motifdb_urls']
    
    # print('results=', results)

    # if this is molecular family analysis, the selected name will be e.g. 'Molecular Family #123 (...)'
    # since we split by space, we take the token at position 2 and remove the first character to get the index
    # otherwise if this is e.g. MS2LDA analysis, the selected name will be e.g. 'motif_123 (...), so
    # we just take the first token as the index

    # update row name
    # 'pw_name': 'Components',
    # significant_column: 'p-value',
    # 'tot_ds_F': 'No. of members',
    row['Components'] = row['pw_name']
    row['p-value'] = row[significant_column]
    row['No. of members'] = row['tot_ds_F']

    pw_name = row['Components']
    p_value = row['p-value']
    no_members = row['No. of members']
    selected = '%s (p-value=%.6e, members=%d)' % (pw_name, p_value, no_members)

    tokens = selected.split(' ')
    database_name = results['database_name']
    idx = tokens[2][1:] if database_name == DATABASE_GNPS_MOLECULAR_FAMILY else tokens[0]
    # display the selected row
    # st.write(row)
    print('row=',row)

    # write the link to MotifDB if available
    motifdb_link = get_motifdb_link(results, row)
    if motifdb_link is not None:
        # display link
        print('Link to MotifDB: %s' % motifdb_link)

    members = dataset_pathways_to_row_ids[idx]
    member_df = get_member_df(entity_dict, members)
    plot_heatmap_image_filename = plot_heatmap(all_groups, all_samples, intensities_df, member_df, members, row)
    # detail_table = display_member_df(member_df)

    if 'link' in member_df.columns:
        member_df['link'] = member_df['link'].apply(make_clickable)

    print('plot_heatmap_image_filename=',plot_heatmap_image_filename)

    details = []

    # current row
    details.append('<h2>Component Browser</h2>')
    details.append('<h4>Components: {}</h4>'.format(no_members))
    details.append('<h4>p-value:  {}</h4>'.format(p_value))
    details.append('<h4>o. of members: {}</h4>'.format(pw_name))

    # heatmap
    image_url = settings.STATIC_URL+plot_heatmap_image_filename
    logger.debug('image_url = %s' % image_url)
    details.append('<img src="{}" class="img-fluid" alt="" srcset="">'.format(image_url))

    # members table
    details.append('<h2>Members</h2>')
    details.append(member_df.to_html(escape=False))

    result = {'status':'success', 'message':'Analysis done!', 'data':{'details':details}}
    return JsonResponse(result)


def PLAGE_decomposition(ds):
    method = PLAGE(ds)
    df = method.get_pathway_df()
    return df

def get_plot_data(gnps_ds, case, control):
    experimental_design = gnps_ds.get_experimental_design()
    all_samples = []
    all_groups = []
    for group in experimental_design['groups']:
        if group == case or group == control:
            samples = experimental_design['groups'][group]
            all_samples.extend(samples)
            all_groups.extend([group] * len(samples))
    entity_dict = gnps_ds.entity_dict
    intensities_df = gnps_ds.standardize_intensity_df()
    dataset_pathways_to_row_ids = gnps_ds.dataset_pathways_to_row_ids
    return all_groups, all_samples, entity_dict, intensities_df, dataset_pathways_to_row_ids


def make_clickable(link):
    # target _blank to open new window
    # extract clickable text to display for your link
    text = 'GNPSLinkout_Network'
    return f'<a target="_blank" href="{link}">{text}</a>'

def process_gnps_results(df, significant_column):
    # filter results to show only the columns we want
    try:
        df = df.drop(columns=['sf', 'exp_F', 'Ex_Cov', 'unq_pw_F', 'F_coverage'])
    except KeyError:
        pass
    df = df[df.columns.drop(list(df.filter(regex='comb_p')))]

    # sort column
    count_col = 'tot_ds_F'
    df = df.sort_values([significant_column, count_col], ascending=[True, False])

    # reorder and rename columns
    df = df[['pw_name', significant_column, 'tot_ds_F']]

    df = df.rename(columns={
        'pw_name': 'Components',
        significant_column: 'p-value',
        'tot_ds_F': 'No. of members',
    })
    return df

def get_motifdb_link(results, row):
    component = row['Components']
    motifdb_link = None
    if 'motifdb_urls' in results:
        motifdb_url = results['motifdb_urls'][component]
        try:
            float(motifdb_url)  # will throw ValueError if it actually contains a string, otherwise it's a nan
        except ValueError:
            motifdb_link = motifdb_url
    return motifdb_link

def get_member_df(entity_dict, members):
    # get group info
    # print('%s p-value=%.4f' % (pw_name, p_value))
    data = []
    from_gnps = True
    for member in members:
        member_info = entity_dict[member]
        unique_id = member_info['unique_id']
        mz = member_info['mass']
        rt = member_info['RT']
        try:
            intensity = member_info['SumPeakIntensity']
            library_id = member_info['LibraryID']
            gnps_linkout_network = member_info['GNPSLinkout_Network']
            no_spectra = member_info['number of spectra']
            temp = [unique_id, library_id, mz, rt, intensity, no_spectra, gnps_linkout_network]
        except KeyError:
            from_gnps = False
            temp = [unique_id, mz, rt]
        data.append(temp)

    if from_gnps:
        member_df = pd.DataFrame(data, columns=['id', 'LibraryID', 'Precursor m/z', 'RTConsensus', 'PrecursorInt',
                                                'no_spectra', 'link']).set_index('id')
        member_df = member_df.rename(columns={
            'Precursor m/z': 'mass',
            'RTConsensus': 'RT'
        })
    else:
        member_df = pd.DataFrame(data, columns=['id', 'mass', 'RT']).set_index('id')
    return member_df

def plot_heatmap(all_groups, all_samples, intensities_df, member_df, members, row):
    # Create a categorical palette to identify the networks
    used_groups = list(set(all_groups))
    group_pal = sns.husl_palette(len(used_groups), s=.90)
    group_lut = dict(zip(map(str, used_groups), group_pal))

    # Convert the palette to vectors that will be drawn on the side of the matrix
    group_intensities = intensities_df.loc[members][all_samples]
    group_colours = pd.Series(all_groups, index=group_intensities.columns).map(group_lut)
    group_colours.name = 'groups'

    # plot heatmap
    g = sns.clustermap(group_intensities, center=0, cmap='vlag', col_colors=group_colours,
                       col_cluster=False, linewidths=0.75, cbar_pos=(0.1, 0.05, 0.05, 0.18))
    pw_name = row['Components']
    plt.suptitle('%s' % (pw_name), fontsize=18, y=0.89)

    # draw group legend
    for group in used_groups:
        g.ax_col_dendrogram.bar(0, 0, color=group_lut[group], label=group, linewidth=0)
    g.ax_col_dendrogram.legend(loc="right")

    # make the annotated peaks to have labels in bold
    try:
        annotated_df = member_df[member_df['LibraryID'].notnull()]
        annotated_peaks = annotated_df.index.values
        for label in g.ax_heatmap.get_yticklabels():
            if label.get_text() in annotated_peaks:
                label.set_weight("bold")
                label.set_color("green")
    except KeyError:
        pass

    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0)

    # render plot
    # generate image
    now_time = str(datetime.datetime.now().strftime('%Y%m%d%H%M%S'))
    plog_image_name = 'plot_heatmap_{}.png'.format(now_time)
    plt.savefig(settings.STATICFILES_DIRS[0] +'/'+plog_image_name)
    return plog_image_name

