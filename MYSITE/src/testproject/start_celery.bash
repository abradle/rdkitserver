export C_FORCE_ROOT=true && celery -A testproject worker -l info > celery.log &
