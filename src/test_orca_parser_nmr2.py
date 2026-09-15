from enan_calculators._orca import _split_post_blocks

text = "! NMR\n %eprnmr nuclei H = all {shift} end end\n%geom end"
pre, post = _split_post_blocks(text)
print("PRE:")
print(pre)
print("POST:")
print(post)
