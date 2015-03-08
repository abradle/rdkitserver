from django.conf.urls import patterns, url

from testproject import views
from django.conf.urls import patterns, include, url

urlpatterns = patterns('',
    url(r'^$', views.index, name='index'),
    url(r'^rdkit_cluster/', include('rdkit_cluster.urls',namespace="rdkit_cluster")),
    url(r'^rdkit_screen/', include('rdkit_screen.urls',namespace="rdkit_screen")),


)
