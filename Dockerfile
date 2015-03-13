FROM ubuntu:14.04
MAINTAINER Anthony Bradley
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y wget build-essential libpq-dev python-dev python-numpy git python-setuptools vim flex bison cmake sqlite3 libsqlite3-dev libboost-dev openbabel libboost-python-dev libboost-regex-dev python-matplotlib python-openbabel python-pip nginx libjpeg8 libjpeg62-dev libfreetype6 libfreetype6-dev unzip python-scipy
RUN apt-get install -y python-psycopg2
RUN apt-get install postgresql postgresql-contrib -y
RUN easy_install ipython Django==1.5.0 tornado
RUN pip install gunicorn django-jfu djutils south
RUN mkdir /RDKit && cd /RDKit && git clone https://github.com/rdkit/rdkit.git
ADD bashrc /root/.bashrc
ADD make_rdkit.bash /make_rdkit.bash
RUN /bin/bash make_rdkit.bash
ADD run_gunicorn.bash /run_gunicorn.bash
ADD nginx.conf /etc/nginx/sites-available/Mysite
RUN ln -s /etc/nginx/sites-available/Mysite /etc/nginx/sites-enabled/Mysite
RUN rm /etc/nginx/sites-enabled/default
ADD set_nginx.bash /set_nginx.bash
ENV SERVER_NAME localhost
ENV SERVER_PORT 9000
ADD MYSITE /MYSITE
RUN mkdir /MYSITE/logs/
RUN python /MYSITE/src/testproject/manage.py collectstatic --noinput
CMD cd /MYSITE && bash /set_nginx.bash && service nginx restart && bash /run_gunicorn.bash
