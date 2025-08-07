FROM tomcat:9.0
# FROM tomcat:9.0.54-jdk11-openjdk
LABEL maintainer="rohdebe1@gmail.com"

RUN apt-get -y update && apt-get install -y

RUN apt-get -y install build-essential

RUN apt-get -y install ant

RUN apt -y install curl

# Add repository
# RUN  echo deb http://deb.debian.org/debian sid main >>/etc/apt/sources.list

# Fix certificate issues
RUN apt-get update && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;

# Install OpenJDK-11
RUN apt-get update && \
    apt-get install -y openjdk-11-jdk && \
    apt-get install -y ant && \
    apt-get clean;

# Setup JAVA_HOME -- useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-11-openjdk-amd64/
RUN export JAVA_HOME

# Copy the needed files to the build directory
COPY ./src /usr/src/ava-formake/src
COPY ./lib /usr/src/ava-formake/lib
COPY makefile /usr/src/ava-formake
COPY makefile.defs /usr/src/ava-formake
COPY build.xml /usr/src/ava-formake

WORKDIR /usr/src/ava-formake

# Build the artifacts
RUN make all

# Set up the web app in Tomcat
# RUN cp target/tomcat_lib/*.* /usr/local/tomcat/shared/lib
RUN cp lib/jstl-1.2.jar /usr/local/tomcat/lib
RUN cp lib/servlet-api-2.4.jar /usr/local/tomcat/lib
RUN cp target/tomcat_lib/*.* /usr/local/tomcat/lib
RUN cp target/archives/*.war /usr/local/tomcat/webapps

# ADD sample.war /usr/local/tomcat/webapps/

EXPOSE 8080
CMD ["catalina.sh", "run"]

# Enable this debugging entry if you need to inspect the running container
# CMD ["sh"]
