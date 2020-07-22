import json
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect

# Create your views here.
def index(request):
    return render(request,'pals/index.html')

def analysis(request):

    if request.method != 'POST':
        return None

    token=request.POST.get('token')
    analysis_id=request.POST.get('analysis_id')
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

    table=[]

    result = {'message':'Analysis done!', 'data':{'table':table}}
    
    return HttpResponse(json.dumps(result))