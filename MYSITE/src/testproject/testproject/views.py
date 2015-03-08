
from django.shortcuts import render

def index(request):
    context = {"projects": None}
    return render(request, 'index.html', context)

