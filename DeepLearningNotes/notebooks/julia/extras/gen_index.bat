tree -H '.' -FL 1 --noreport --charset utf-8 --matchdirs -P "*.ipynb" | grep -v extras > index.html
