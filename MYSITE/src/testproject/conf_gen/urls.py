from django.conf.urls import patterns, url
from conf_gen import views

urlpatterns = patterns('',
    url(r'^gen_confs/$', views.gen_confs, name='gen_confs'),
    url(r'^get_moments/$', views.get_moments, name='get_moments'),
    url(r'^gen_moments/$', views.gen_moments, name='gen_moments'),
    url(r'^get_progress/$', views.get_progress, name='get_progress'),
    url(r'^get_confs/$', views.get_confs, name='get_confs'),
)
