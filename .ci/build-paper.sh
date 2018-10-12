#!/bin/bash -x

echo "executing build-paper.sh"
#if git diff --name-only $TRAVIS_COMMIT_RANGE | grep 'paper/'
#then
echo "Building the paper..."
export CELERITE_BUILDING_PAPER=true
source "$( dirname "${BASH_SOURCE[0]}" )"/setup-tectonic.sh
#return
#fi

# Build the paper
cd paper
make

# Push to GitHub
if [ -n "$GH_TOKEN" ]; then
  cd $TRAVIS_BUILD_DIR
  git checkout --orphan $TRAVIS_BRANCH-pdf
  # next line very important 
  git rm -rf .
  ls
  git add -f paper/paper.pdf
  git -c user.name='travis' -c user.email='travis' commit -m "building the paper"
  git push -q -f https://$GH_TOKEN@github.com/$TRAVIS_REPO_SLUG $TRAVIS_BRANCH-pdf
fi
