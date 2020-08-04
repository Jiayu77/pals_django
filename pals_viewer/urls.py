from django.urls import path
from . import views

app_name = 'pals_viewer'

urlpatterns = [
    # pathway analysis
    path('', views.index, name='index'),
    path('analysis', views.analysis, name='analysis'),
    path('pathway/show_kegg_diagram', views.show_kegg_diagram, name='show_kegg_diagram'),
    path('pathway/show_reactome_diagram', views.show_reactome_diagram, name='show_reactome_diagram'),

    # GNPS analysis
    path('gnps/index', views.gnps_index, name='gnps_index'),
    path('gnps/analysis', views.gnps_analysis, name='gnps_analysis'),

    # MS2LDA analysis
    path('ms2lda/index', views.ms2lda_index, name='ms2lda_index'),
    path('ms2lda/analysis', views.ms2lda_analysis, name='ms2lda_analysis'),
]