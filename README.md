# To start container

    sudo docker rm -f emdemo; sudo docker run --restart always --name emdemo -d -p 3846:3838 -t cannin/emdemo shiny-server