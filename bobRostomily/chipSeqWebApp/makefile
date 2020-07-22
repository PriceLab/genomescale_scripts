ID = 21jul2020
APP_IMAGE_NAME = rosto_webapp_$(ID)

build:
	docker build -t $(APP_IMAGE_NAME) -f Dockerfile .

bash:
	docker run -it \
           -p 3838:3838 \
           -v ~/github/LC_Systemsbiology/webApps/lcApp/datavol:/datavol \
           --entrypoint /bin/bash \
           $(APP_IMAGE_NAME)

