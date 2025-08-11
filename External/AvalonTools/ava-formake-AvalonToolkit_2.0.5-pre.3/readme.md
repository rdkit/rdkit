# Avalon Toolkit Repo

This repository contains the public parts of the Avalon toolkit. Building it with <code>make all</code> will create a number of artifacts
for the running platform.

The file <code>makefile.defs</code> contains definitions that may be specific for the respective platform:

## Pointers to JAR archives to be used on the classpath

| Variable | Description |
| -------- | ------- |
| <code>LOG4J_URL</code> | Points to the JAR file that implements LOG4J |
| <code>COMMONS_URL</code> |  |
| <code>CORS_URL</code> |  |
| <code>PROP_UTILS_URL</code> |  |
| <code>XOM_URL</code> |  |
| <code>JSTL_URL</code> |  |
| <code>SERVLET_URL</code> |  |

## Home directories of standard tools

These variable should normally be set and exported already by build environment. The <code>makefile.defs</code> file contains example settings
if the environment doesn't provide it.

| Directory | Tool description |
| -------- | ------- |
| <code>JAVA_HOME</code> | Main root of the Java SDK tools |
| <code>ANT_HOME</code> | Main root of Jakarta <code>Ant</code>. <code>Ant</code> is used to create the JAR artifacts ot the repo. |
