from django.conf.urls import patterns, url

from rdkit_cluster import views

urlpatterns = patterns('',
    url(r'^$', views.index, name='index'),
    url(r'^cluster/$', views.cluster, name='cluster'),
    url(r'^cluster_mol_body/$', views.cluster_mol_body, name='cluster_mol_body'),

)

