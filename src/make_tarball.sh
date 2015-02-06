#! /bin/sh

repo=lare3d
cur=`pwd`
dir=$(mktemp -d -t$repo)

git init -q $dir/$repo
git push -q --tags $dir/$repo HEAD:tmp
cd $dir/$repo
git checkout -q tmp
git submodule update --init --recursive
cstring=$(git describe --abbrev=0 --match v[0-9]* HEAD | cut -c2-)
fullstring=$(git describe --match v[0-9]* HEAD | cut -c2-)

# Append date if this is not a tagged release
if [ "$cstring"x != "$fullstring"x ]; then
  # Only append date if there are actually changes
  git diff --quiet "$(git describe --abbrev=0 --match v[0-9]* HEAD)"
  if [ $? -ne 0 ]; then
    dt=$(git log --pretty=format:%ci -1 HEAD | cut -f1 -d' ')
    cstring="${cstring}-$dt"
  fi
fi

(cd SDF/VisIt
/bin/sh gen_commit_string)
/bin/sh src/gen_commit_string
rm -rf .git
cd $dir
mv $repo $repo-$cstring
tar -cf - $repo-$cstring | gzip -c > $repo-$cstring.tar.gz
cd $cur
mv $dir/$repo-$cstring.tar.gz .
rm -rf $dir
