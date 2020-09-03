from django.urls import path
from . import views

app_name = 'pals_viewer'

urlpatterns = [
    # pathway analysis
    path('', views.index, name='index'),
    path('pathway/analysis', views.analysis, name='analysis'),
    path('pathway/get_data', views.keypath_get_data, name='keypath_get_data'),
    path('pathway/show_kegg_diagram', views.show_kegg_diagram, name='show_kegg_diagram'),
    path('pathway/show_reactome_diagram', views.show_reactome_diagram, name='show_reactome_diagram'),

    # GNPS analysis
    path('gnps/get_data', views.gnps_get_data, name='gnps_get_data'),
    path('gnps/index', views.gnps_index, name='gnps_index'),
    path('gnps/analysis', views.gnps_analysis, name='gnps_analysis'),
    path('gnps/details', views.gnps_show_details, name='gnps_show_details'),

    # MS2LDA analysis
    path('ms2lda/index', views.ms2lda_index, name='ms2lda_index'),
    path('ms2lda/get_data', views.ms2lda_get_data, name='ms2lda_get_data'),
    path('ms2lda/analysis', views.ms2lda_analysis, name='ms2lda_analysis'),
    path('ms2lda/details', views.ms2lda_show_details, name='ms2lda_show_details'),
]