from django.conf.urls import patterns, url

from rdkit_screen import views

urlpatterns = patterns('',
    url(r'^$', views.index, name='index'),
    url(r'^screen/$', views.screen, name='screen'),
    url(r'^screen_mol_body/$', views.screen_mol_body, name='screen_mol_body'),
    url(r'^screen_simple/$', views.screen_simple, name='screen_simple'),

)
