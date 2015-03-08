from django.conf.urls import patterns, url

from rdkit_cluster import views

urlpatterns = patterns('',
    url(r'^$', views.index, name='index'),
    url(r'^cluster/$', views.cluster, name='cluster'),
)

