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
ADD USRCAT /USRCAT
RUN cd /USRCAT && python setup.py install
RUN mkdir /MYSITE/logs/
RUN apt-get install -y rabbitmq-server
RUN pip install celery boto django-celery Closeablequeue
RUN python /MYSITE/src/testproject/manage.py collectstatic --noinput
# Now make the databse
# Run the rest of the commands as the ``postgres`` user created by the ``postgres-9.3`` package when it was ``apt-get installed``
USER postgres
RUN /etc/init.d/postgresql start &&\
    psql --command "ALTER USER postgres WITH PASSWORD 'postgres';" &&\
    createdb -O postgres mydb
USER root
# Hack to get around some AUFS weriedness
RUN mkdir /etc/ssl/private-copy; mv /etc/ssl/private/* /etc/ssl/private-copy/; rm -r /etc/ssl/private; mv /etc/ssl/private-copy /etc/ssl/private; chmod -R 0700 /etc/ssl/private; chown -R postgres /etc/ssl/private
RUN rabbitmq-server & sleep 10 && rabbitmqctl add_user myuser mypassword && rabbitmqctl add_vhost myvhost && rabbitmqctl set_permissions -p myvhost myuser ".*" ".*" ".*" && rabbitmq-plugins enable rabbitmq_management
RUN service postgresql start && python /MYSITE/src/testproject/manage.py syncdb --noinput
CMD rabbitmq-server -detached && cd /MYSITE && bash /set_nginx.bash && service postgresql start && service nginx restart && bash /run_gunicorn.bash
