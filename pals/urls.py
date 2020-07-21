from django.urls import path
from . import views

app_name = 'pals'

urlpatterns = [
    path('', views.index, name='index')
]