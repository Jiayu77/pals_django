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
    return render(request,'pals_viewer/index.html')

def analysis(request):

    if request.method != 'POST':
        return None

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

    int_df, annotation_df, experimental_design = download_from_pimp(token, PIMP_HOST, analysis_id, 'kegg')
    print('experimental design:', experimental_design)

    table=[]

    result = {'message':'Analysis done!', 'data':{'table':table}}

    import time
    time.sleep(1)
    
    return HttpResponse(json.dumps(result))