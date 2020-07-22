import json
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from pals.pimp_tools import get_pimp_API_token_from_env, PIMP_HOST, download_from_pimp
from pals.feature_extraction import DataSource
from pals.PLAGE import PLAGE
from pals.ORA import ORA
from pals.GSEA import GSEA
from pals.common import *

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

    return render(request,'pals_viewer/index.html', content_dict)

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