tree -H '.' -T "Jupyter notebooks <br> (use right click and save to download)" -FL 1 --noreport --charset utf-8 --matchdirs -P "*.ipynb" | grep -v extras > index.html
