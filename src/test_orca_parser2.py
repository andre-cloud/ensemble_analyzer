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
        stripped = line.strip().lower()
        
        if not in_post_block:
            # Check if this line starts a post-coords block
            kw_match = next((kw for kw in _POST_COORDS_KEYWORDS if stripped.startswith(kw)), None)
            if kw_match:
                in_post_block = True
                current_block.append(line)
                # If the block is somehow single-line "end" (rare but possible like %frag ... end)
                if stripped.endswith(" end") or stripped == "end":
                    # But usually %block starts a multiline block
                    pass
            else:
                pre.append(line)
        else:
            current_block.append(line)
            # Check if this line ends the block
            # In ORCA, 'end' on a line closes the block
            if stripped == 'end' or stripped.split('#')[0].strip() == 'end':
                post.extend(current_block)
                current_block = []
                in_post_block = False
                
        i += 1
        
    # If we reached EOF and were still in a post block (missing 'end'), just dump it to post
    if current_block:
        post.extend(current_block)
        
    return '\n'.join(pre), '\n'.join(post)

text = """%frag
  fragment definitions...
end
%scf
  MaxIter 200
end
"""
pre, post = _split_post_blocks_new(text)
print("PRE:")
print(pre)
print("POST:")
print(post)
