#!/bin/bash -x

# Build the paper
cd paper
make

# Push to GitHub
if [ -n "$GH_TOKEN" ]; then
  cd $TRAVIS_BUILD_DIR
  git checkout --orphan $TRAVIS_BRANCH-pdf
  git rm -rf .
  git add -f paper/paper.pdf
  git -c user.name='travis' -c user.email='travis' commit -m "building the paper"
  git push -q -f https://$GITHUB_USER:$GH_TOKEN@github.com/$TRAVIS_REPO_SLUG $TRAVIS_BRANCH
fi
