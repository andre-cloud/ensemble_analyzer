from enan_calculators._orca import _split_post_blocks

text = """%frag
  fragment definitions...
end
%scf
  MaxIter 200
end
"""
pre, post = _split_post_blocks(text)
print("PRE:")
print(pre)
print("POST:")
print(post)
