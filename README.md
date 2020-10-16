# Mordred-Descriptors

REST Microservice that communicates with  jaqpot-api and calculates the respective mordred descriptors.

**How to run mordred-descriptors (Docker)**

After downloading the source code from github, navigate till the app path.

`docker build -t mordred .`  
`docker run --name mordredcontainer -p 8000:8000 mordred`

You can navigate to  http://127.0.0.1:8000/docs to interact with the Mordred Descriptors' API

Additional documentation can be found at  http://127.0.0.1:8000/redoc
 
**Learn more about Mordred Descriptors** 

http://mordred-descriptor.github.io/documentation/master/index.html
