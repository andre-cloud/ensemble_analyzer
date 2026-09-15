from enan_calculators._orca import _split_post_blocks

text1 = """! NMR
%eprnmr nuclei H = all {shift} end end
%geom end
"""

text2 = """! NMR
%eprnmr
  nuclei H = all {shift} end
end
%geom end
"""

text3 = """! NMR
%eprnmr
  nuclei H = all {shift}
end end
%geom end
"""

for i, t in enumerate([text1, text2, text3]):
    print(f"--- TEXT {i+1} ---")
    pre, post = _split_post_blocks(t)
    print("PRE:\n" + pre)
    print("POST:\n" + post)
    print()
