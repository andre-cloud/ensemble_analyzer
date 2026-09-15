import re

text = "! optts %maxcore 4000 %eprnmr nuclei H = all {shift} end end \n! Freq \n%pal nprocs_group 4 end %scf MaxIter 200 end"
parts = re.split(r'(?=[%!])', text)
for i, p in enumerate(parts):
    print(f"[{i}]: {repr(p)}")
