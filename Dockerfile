FROM python:3.6

RUN mkdir /opt/vtool3
RUN pip install Pybuilder
ADD ./ /opt/vtool3/
RUN cd /opt/vtool3 && pyb install_dependencies && pyb -x run_unit_tests package && pip install target/dist/vtool3-1.0.dev0
RUN useradd -m vtool

USER vtool
WORKDIR /home/vtool
VOLUME ['/home/vtool']

