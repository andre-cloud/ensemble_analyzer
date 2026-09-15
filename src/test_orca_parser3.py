_POST_COORDS_KEYWORDS = frozenset({'%frag', '%eprnmr', '%nmr', '%rel', '%epr'})

def _split_post_blocks_new(text: str) -> tuple[str, str]:
    if not text.strip():
        return text, ""
    
    lines = text.split('\n')
    pre = []
    post = []
    
    i = 0
    in_post_block = False
    current_block = []
    
    while i < len(lines):
        line = lines[i]
        no_comment = line.split('#')[0].strip().lower()
        
        if not in_post_block:
            kw_match = next((kw for kw in _POST_COORDS_KEYWORDS if no_comment.startswith(kw)), None)
            if kw_match:
                # Check if it ends with ' end' or is literally exactly just the keyword + end
                if no_comment.endswith(" end") or no_comment == kw_match + " end":
                    post.append(line)
                else:
                    in_post_block = True
                    current_block.append(line)
            else:
                pre.append(line)
        else:
            current_block.append(line)
            if no_comment == 'end':
                post.extend(current_block)
                current_block = []
                in_post_block = False
                
        i += 1
        
    if current_block:
        post.extend(current_block)
        
    return '\n'.join(pre), '\n'.join(post)

text1 = """%frag
  fragment definitions...
end
%scf
  MaxIter 200
end
"""

text2 = """%scf
  MaxIter 200
end
%eprnmr gtensor true end # comment
%geom
end
"""

pre, post = _split_post_blocks_new(text1)
print("T1 PRE:\n" + pre + "\nT1 POST:\n" + post)

pre, post = _split_post_blocks_new(text2)
print("T2 PRE:\n" + pre + "\nT2 POST:\n" + post)
