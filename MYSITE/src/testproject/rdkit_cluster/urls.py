from django.conf.urls import patterns, url

from rdkit_cluster import views

urlpatterns = patterns('',
    url(r'^$', views.index, name='index'),
    url(r'^cluster_simple/$', views.cluster_simple, name='cluster_simple'),
    url(r'^return_sd/$', views.return_sd, name='return_sd'),
)

