
for f in reference/summary/*
do
  diff $f ${f/"reference"/"out"} > diffs/${f/"reference/summary"/}
done
