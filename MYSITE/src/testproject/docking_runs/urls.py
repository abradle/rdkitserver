from django.conf.urls import patterns, url

from docking_runs import views

urlpatterns = patterns('',
    url(r'^start_dock/$', views.start_dock, name='start_dock'),
    url(r'^check_dock/$', views.check_dock, name='check_dock'),
    url(r'^libs_avail/$', views.libs_avail, name='libs_avail'),
    url(r'^docks_avail/$', views.docks_avail, name='docks_avail'),
    url(r'^get_output/$', views.get_output, name='get_output'),
)

