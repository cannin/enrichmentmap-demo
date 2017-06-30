# Get Docker Container
    docker pull cannin/emdemo

# Start Container
    docker run --restart always --name emdemo -d -p 3838:3838 -t cannin/emdemo shiny-server

# Stop/Remove Running Container
    docker rm -f emdemo

# Run Shiny with an Interactive Shell
## Command
    docker run -i --name emdemo -p 3838:3838 -t cannin/emdemo bash

## Then run shiny-server
    shiny-server

## Get another interactive shell
    docker exec -i -t emdemo bash
