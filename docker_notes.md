# Rudimentary notes on Docker

```sh
# Set up docker env (only needed if you use command line docker-machine)
docker-machine start
eval $(docker-machine env default)

# Create container
docker build -t varikn_test_container .
# Log in container
docker container run --interactive --tty  varikn_test_container
# Inside container
mkdir build; (cd build; cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=1
  && make && ctest --verbose)

# clean up
docker image rm varikn_test_container
docker machine stop
```

Notes for mac docker setup

```sh
# Install needed Homebrew packages
brew install docker docker-machine
brew cask install virtualbox
# Create VM
docker-machine create default
```
