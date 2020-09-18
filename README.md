# pals_django
This project is a ui system for PALS project developed by Django.

Installation instruction

Operating environment installation

The environment to be installed in this project is divided into two parts, one is the basic environment, and the other is the database environment.
Basic environment
The project is developed based on python3.8, and all python libraries are configured using pipenv. So we only need to follow the steps below:
1. Install python3.8 version.
2. Install pipenv.
3. The configuration file of pipenv already exists in the project directory. Use the configuration file to start pipenv directly in the project directory to enter the configured pipenv environment.


Database environment

This project needs to use the reactome database, and needs:
1. Install the neo4j database.
2. Import the reaction gene library into neo4j.
Startup project
1. Start the neo4j database to ensure normal access to the reactome database
2. Run python manage.py runserver in the project directory to start the project
