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


def analysis(request):

    if request.method != 'POST':
        return None

    # get data from POST
    token=request.POST.get('token')
    analysis_id=int(request.POST.get('analysis_id'))
    experimental_design=request.POST.get('experimental_design')
    pathway_analysis_method=request.POST.get('pathway_analysis_method')
    database=request.POST.get('database')
    reactome_species=request.POST.get('reactome_species')
    reactome_metabolic_pathway_only=request.POST.get('reactome_metabolic_pathway_only')
    reactome_query=request.POST.get('reactome_query')

    print('token=', token)
    print('analysis_id=', analysis_id)
    print('experimental_design=', experimental_design)
    print('pathway_analysis_method=', pathway_analysis_method)
    print('database=', database)
    print('reactome_species=', reactome_species)
    print('reactome_metabolic_pathway_only=', reactome_metabolic_pathway_only)
    print('reactome_query=', reactome_query)

    # download data from PiMP
    int_df, annotation_df, experimental_design = download_from_pimp(token, PIMP_HOST, analysis_id, 'kegg')
    print('experimental design:', experimental_design)

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

def gnps_analysis(request):
    if request.method != 'POST':
        return None

    # get data from POST
    gnps_url=request.POST.get('gnps_url')
    metadata=request.FILES.get('metadata')
    # comparisons=request.POST.get('comparisons')
    pathway_analysis_method=request.POST.get('pathway_analysis_method')
    database=request.POST.get('database')
    reactome_species=request.POST.get('reactome_species')
    reactome_metabolic_pathway_only=request.POST.get('reactome_metabolic_pathway_only')
    reactome_query=request.POST.get('reactome_query')

    print('gnps_url=', gnps_url)
    print('metadata=', metadata, ',type=', type(metadata))
    # print('comparisons=', comparisons)
    print('pathway_analysis_method=', pathway_analysis_method)
    print('database=', database)
    print('reactome_species=', reactome_species)
    print('reactome_metabolic_pathway_only=', reactome_metabolic_pathway_only)
    print('reactome_query=', reactome_query)

    case = 'More than 30'
    control = 'Less than 10'
    comp_name = 'more_plants/no_plants'
    comparisons = [{'case': case, 'control': control, 'name': comp_name },]

    # save and read metadata csv file
    # reset new file name
    now_time = str(datetime.datetime.now().strftime('%Y%m%d%H%M%S'))
    oldfilename = metadata.name
    newfilename = oldfilename.replace(".", "_{}.".format(now_time))
    newfilepath = settings.MEDIA_ROOT+'/'+newfilename
    with open(newfilepath, 'wb+') as destination:
        for chunk in metadata.chunks():
            destination.write(chunk)
    
    print('saved {} to local:{}'.format(oldfilename, newfilepath))
    metadata_file = newfilepath
    metadata_df = pd.read_csv(metadata_file)

    loader = GNPSLoader(database, gnps_url, metadata_df, comparisons)
    database = loader.load_data()

    measurement_df = database.extra_data['measurement_df']
    annotation_df = database.extra_data['annotation_df']
    experimental_design = database.extra_data['experimental_design']

    # set data source
    gnps_ds = DataSource(measurement_df, annotation_df, experimental_design, None, database=database, min_replace=SMALL)

    plage = PLAGE(gnps_ds, num_resamples=1000)
    pathway_df = plage.get_pathway_df(standardize=True)

    # analysis by different methods
    p_value_col = '%s p-value' % comp_name
    count_col = 'unq_pw_F'
    pathway_df.sort_values([p_value_col, count_col], ascending=[True, False], inplace=True)

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