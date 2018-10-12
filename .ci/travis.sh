#!/bin/bash -x

echo "executing travis.sh"
if [[ $TEST_LANG == paper ]]
then
  if git diff --name-only $TRAVIS_COMMIT_RANGE | grep 'paper/'
  then
    echo "Building the paper..."
    export CELERITE_BUILDING_PAPER=true
    source "$( dirname "${BASH_SOURCE[0]}" )"/setup-tectonic.sh
    return
  fi
  export CELERITE_BUILDING_PAPER=false
  return
fi
